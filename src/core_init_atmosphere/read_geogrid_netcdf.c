// Copyright (c) 2024,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
/* File: read_geogrid_netcdf.c

   Subroutine to read a sub-region (tile) from a NetCDF geographic data file.
   This allows init_atmosphere to use NetCDF-format source data in addition
   to the traditional WPS geogrid binary format.

   The NetCDF file is expected to contain a 2D or 3D variable on a regular
   lat-lon grid. The routine reads a rectangular sub-region specified by
   start indices and counts.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#define NC_ERR_CHECK(e, msg, status_ptr) \
    if ((e) != NC_NOERR) { \
        fprintf(stderr, "read_geogrid_netcdf: %s: %s\n", (msg), nc_strerror(e)); \
        *(status_ptr) = 1; \
        return 1; \
    }

/*  In Fortran, use the following as an interface for read_geogrid_netcdf:

 use iso_c_binding, only : c_char, c_int, c_float, c_ptr, c_loc

 interface
    subroutine read_geogrid_netcdf(fname, varname, rarray, nx, ny, nz, &
                                   start_x, start_y, start_z, status) bind(C)
       use iso_c_binding, only : c_char, c_int, c_float, c_ptr
       character (c_char), dimension(*), intent(in) :: fname
       character (c_char), dimension(*), intent(in) :: varname
       type (c_ptr), value :: rarray
       integer (c_int), intent(in), value :: nx
       integer (c_int), intent(in), value :: ny
       integer (c_int), intent(in), value :: nz
       integer (c_int), intent(in), value :: start_x
       integer (c_int), intent(in), value :: start_y
       integer (c_int), intent(in), value :: start_z
       integer (c_int), intent(inout) :: status
    end subroutine read_geogrid_netcdf
 end interface

*/

int read_geogrid_netcdf(
      char * fname,            /* Path to the NetCDF file */
      char * varname,          /* Name of the variable to read */
      float * rarray,          /* Output: array to be filled [nx * ny * nz] */
      int nx,                  /* Number of pixels to read in x (column) direction */
      int ny,                  /* Number of pixels to read in y (row) direction */
      int nz,                  /* Number of levels to read in z direction */
      int start_x,            /* Starting x index (0-based) */
      int start_y,            /* Starting y index (0-based) */
      int start_z,            /* Starting z index (0-based) */
      int * status)
{
   int ncid, varid, retval;
   int ndims;
   int dimids[NC_MAX_VAR_DIMS];
   size_t start[3], count[3];
   size_t dim_lens[3];
   int actual_nx, actual_ny, actual_nz;
   int i, j, k;
   float *tmpbuf = NULL;
   float fill_value;
   int has_fill;

   *status = 0;

   /* Open the NetCDF file */
   retval = nc_open(fname, NC_NOWRITE, &ncid);
   NC_ERR_CHECK(retval, "opening file", status);

   /* Get the variable ID */
   retval = nc_inq_varid(ncid, varname, &varid);
   NC_ERR_CHECK(retval, "finding variable", status);

   /* Get number of dimensions */
   retval = nc_inq_varndims(ncid, varid, &ndims);
   NC_ERR_CHECK(retval, "getting ndims", status);

   if (ndims < 2 || ndims > 3) {
      fprintf(stderr, "read_geogrid_netcdf: variable '%s' has %d dimensions, expected 2 or 3\n",
              varname, ndims);
      nc_close(ncid);
      *status = 1;
      return 1;
   }

   /* Get dimension IDs and lengths */
   retval = nc_inq_vardimid(ncid, varid, dimids);
   NC_ERR_CHECK(retval, "getting dimids", status);

   for (i = 0; i < ndims; i++) {
      retval = nc_inq_dimlen(ncid, dimids[i], &dim_lens[i]);
      NC_ERR_CHECK(retval, "getting dim length", status);
   }

   /*
    * Determine dimension ordering.
    * For 2D: assume (y, x) - i.e., (lat, lon) or (south_north, west_east)
    * For 3D: assume (z, y, x) or (y, x, z)
    *
    * We assume the most common convention: last dimension varies fastest,
    * and the ordering is (..., y, x) for 2D or (z, y, x) for 3D.
    *
    * Geogrid convention: data is stored row by row (y varies slowest for 2D).
    * NetCDF C convention: first dimension varies slowest.
    * So for a 2D array with dims (ny_total, nx_total), dim[0]=ny, dim[1]=nx.
    */

   if (ndims == 2) {
      /* (y, x) ordering */
      actual_ny = (int)dim_lens[0];
      actual_nx = (int)dim_lens[1];
      actual_nz = 1;

      /* Clamp the read region to the actual data extent */
      int read_nx = nx;
      int read_ny = ny;
      int sx = start_x;
      int sy = start_y;

      if (sx < 0) sx = 0;
      if (sy < 0) sy = 0;
      if (sx + read_nx > actual_nx) read_nx = actual_nx - sx;
      if (sy + read_ny > actual_ny) read_ny = actual_ny - sy;

      if (read_nx <= 0 || read_ny <= 0) {
         /* Entire tile is outside the data domain - fill with zeros */
         for (i = 0; i < nx * ny * nz; i++) rarray[i] = 0.0f;
         nc_close(ncid);
         return 0;
      }

      /* Read the sub-region */
      start[0] = (size_t)sy;
      start[1] = (size_t)sx;
      count[0] = (size_t)read_ny;
      count[1] = (size_t)read_nx;

      /* If we're reading less than the full tile (edge clamping), use a temp buffer */
      if (read_nx != nx || read_ny != ny) {
         tmpbuf = (float *)calloc((size_t)(nx * ny), sizeof(float));
         if (!tmpbuf) {
            fprintf(stderr, "read_geogrid_netcdf: memory allocation failed\n");
            nc_close(ncid);
            *status = 1;
            return 1;
         }

         /* Read into temporary buffer row by row, or as a block and copy */
         float *readbuf = (float *)malloc((size_t)(read_nx * read_ny) * sizeof(float));
         if (!readbuf) {
            free(tmpbuf);
            nc_close(ncid);
            *status = 1;
            return 1;
         }

         retval = nc_get_vara_float(ncid, varid, start, count, readbuf);
         if (retval != NC_NOERR) {
            fprintf(stderr, "read_geogrid_netcdf: reading data: %s\n", nc_strerror(retval));
            free(tmpbuf);
            free(readbuf);
            nc_close(ncid);
            *status = 1;
            return 1;
         }

         /* Copy into the full tile buffer, offset appropriately */
         int x_offset = (start_x < 0) ? -start_x : 0;
         int y_offset = (start_y < 0) ? -start_y : 0;
         for (j = 0; j < read_ny; j++) {
            for (i = 0; i < read_nx; i++) {
               tmpbuf[(j + y_offset) * nx + (i + x_offset)] = readbuf[j * read_nx + i];
            }
         }
         free(readbuf);

         /* Copy to output (geogrid convention: x varies fastest in memory) */
         for (i = 0; i < nx * ny; i++) {
            rarray[i] = tmpbuf[i];
         }
         free(tmpbuf);
      } else {
         /* Simple case: read directly into output array */
         retval = nc_get_vara_float(ncid, varid, start, count, rarray);
         NC_ERR_CHECK(retval, "reading data", status);
      }

   } else {
      /* 3D: assume (z, y, x) ordering */
      actual_nz = (int)dim_lens[0];
      actual_ny = (int)dim_lens[1];
      actual_nx = (int)dim_lens[2];

      int read_nx = nx;
      int read_ny = ny;
      int read_nz = nz;
      int sx = start_x;
      int sy = start_y;
      int sz = start_z;

      if (sx < 0) sx = 0;
      if (sy < 0) sy = 0;
      if (sz < 0) sz = 0;
      if (sx + read_nx > actual_nx) read_nx = actual_nx - sx;
      if (sy + read_ny > actual_ny) read_ny = actual_ny - sy;
      if (sz + read_nz > actual_nz) read_nz = actual_nz - sz;

      if (read_nx <= 0 || read_ny <= 0 || read_nz <= 0) {
         for (i = 0; i < nx * ny * nz; i++) rarray[i] = 0.0f;
         nc_close(ncid);
         return 0;
      }

      start[0] = (size_t)sz;
      start[1] = (size_t)sy;
      start[2] = (size_t)sx;
      count[0] = (size_t)read_nz;
      count[1] = (size_t)read_ny;
      count[2] = (size_t)read_nx;

      if (read_nx != nx || read_ny != ny || read_nz != nz) {
         int total = nx * ny * nz;
         int read_total = read_nx * read_ny * read_nz;
         tmpbuf = (float *)calloc((size_t)total, sizeof(float));
         float *readbuf = (float *)malloc((size_t)read_total * sizeof(float));
         if (!tmpbuf || !readbuf) {
            if (tmpbuf) free(tmpbuf);
            if (readbuf) free(readbuf);
            nc_close(ncid);
            *status = 1;
            return 1;
         }

         retval = nc_get_vara_float(ncid, varid, start, count, readbuf);
         if (retval != NC_NOERR) {
            fprintf(stderr, "read_geogrid_netcdf: reading 3D data: %s\n", nc_strerror(retval));
            free(tmpbuf);
            free(readbuf);
            nc_close(ncid);
            *status = 1;
            return 1;
         }

         int x_offset = (start_x < 0) ? -start_x : 0;
         int y_offset = (start_y < 0) ? -start_y : 0;
         int z_offset = (start_z < 0) ? -start_z : 0;
         for (k = 0; k < read_nz; k++) {
            for (j = 0; j < read_ny; j++) {
               for (i = 0; i < read_nx; i++) {
                  int dst_idx = (k + z_offset) * (nx * ny) + (j + y_offset) * nx + (i + x_offset);
                  int src_idx = k * (read_nx * read_ny) + j * read_nx + i;
                  tmpbuf[dst_idx] = readbuf[src_idx];
               }
            }
         }
         free(readbuf);
         for (i = 0; i < total; i++) rarray[i] = tmpbuf[i];
         free(tmpbuf);
      } else {
         retval = nc_get_vara_float(ncid, varid, start, count, rarray);
         NC_ERR_CHECK(retval, "reading 3D data", status);
      }
   }

   /* Check for fill values and replace with 0 */
   has_fill = 0;
   retval = nc_get_att_float(ncid, varid, "_FillValue", &fill_value);
   if (retval == NC_NOERR) has_fill = 1;
   if (!has_fill) {
      retval = nc_get_att_float(ncid, varid, "missing_value", &fill_value);
      if (retval == NC_NOERR) has_fill = 1;
   }

   if (has_fill) {
      int total = nx * ny * nz;
      for (i = 0; i < total; i++) {
         if (rarray[i] == fill_value) rarray[i] = 0.0f;
      }
   }

   /* Apply scale_factor and add_offset if present (CF convention) */
   {
      float scale = 1.0f, offset = 0.0f;
      int has_scale = 0, has_offset = 0;

      retval = nc_get_att_float(ncid, varid, "scale_factor", &scale);
      if (retval == NC_NOERR) has_scale = 1;

      retval = nc_get_att_float(ncid, varid, "add_offset", &offset);
      if (retval == NC_NOERR) has_offset = 1;

      if (has_scale || has_offset) {
         int total = nx * ny * nz;
         for (i = 0; i < total; i++) {
            rarray[i] = rarray[i] * scale + offset;
         }
      }
   }

   nc_close(ncid);
   return 0;
}
