#!/bin/csh

# Set  GRIB_API
setenv GRIB_API_INCLUDE -I/usr/local/ecmwf_tools/pgi-9.0-1/include
setenv GRIB_API_LIB     -L/usr/local/ecmwf_tools/pgi-9.0-1/lib

# Set netCDF
set ncdf_incs  = `nc-config --fflags`
set ncdf_libs  = `nc-config --flibs`

# Remove old files
\rm *.o
\rm *.a
\rm frgbtocdf_api

# libgm2em + libtra + libipo
pgf90 -c -O libgm2em.f ${ncdf_incs}
ar r libgm2em.a libgm2em.o
ranlib libgm2em.a

pgf90 -c -O libtra.f ${ncdf_incs}
ar r libtra.a libtra.o
ranlib libtra.a

pgf90 -c -O libipo.f ${ncdf_incs}
ar r libipo.a libipo.o
ranlib libipo.a

# cdfio + cdfplus
pgf90 -c -O -I${NETCDF}/include libcdfio.f ${ncdf_incs}
ar r libcdfio.a libcdfio.o
ranlib libcdfio.a

pgf90 -c -O -I${NETCDF}/include libcdfplus.f ${ncdf_incs}
ar r libcdfplus.a libcdfplus.o
ranlib libcdfplus.a

# fgrb2cdf_api
pgf90 -O -o fgrbtocdf_api ${GRIB_API_INCLUDE} fgrbtocdf_api.f ${ncdf_incs} -L. -lcdfio -lcdfplus -ltra -lipo -lgm2em ${ncdf_libs} ${GRIB_API_INCLUDE} ${GRIB_API_LIB}  -lgrib_api_f77 -lgrib_api -ljasper -lm

exit 0

