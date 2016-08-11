#!/bin/csh

# -------------------------------------------------
# Set some parameters
# -------------------------------------------------

# Input GRIB directory
set grbdir=/lhome/sprenger/lagranto.ecmwf/cdo

# Output netCDF directory
set cdfdir=/lhome/sprenger/lagranto.ecmwf/cdo

# Start and end date for conversion, and time step
set startdate = 20160201_00
set finaldate = 20160201_18
set timestep  = 6

# -------------------------------------------------
# Do the conversion
# -------------------------------------------------

# Incrrement finaldate by one timestep - to include finaldate
set finaldate=`newtime ${finaldate} ${timestep}` 

# Change to grib directory
cd ${cdfdir}

# Start loop over all dates
set date=${startdate}
loop:

# Convert an${date}_uvwt
\rm -f P${date}_tuvw
cdo -f nc -t ecmwf copy -invertlat -chname,W,OMEGA ${grbdir}/an${date}_tuvw P${date}_tuvw

# Convert an${date}_q
\rm -f P${date}_q
cdo -f nc -t ecmwf copy -invertlat ${grbdir}/an${date}_q P${date}_q

# Convert an${date}_ps
\rm -f P${date}_ps
cdo -f nc -t ecmwf copy -invertlat ${grbdir}/an${date}_ps P${date}_ps
ncap2   -O -s 'PS=0.01f*exp(LNSP)' P${date}_ps  P${date}_ps

# Merge all files
\rm -f P${date}
cdo -f nc merge P${date}_tuvw P${date}_q P${date}_ps P${date}
\rm -f P${date}_tuvw
\rm -f P${date}_q
\rm -f P${date}_ps

# Proceed to next date
set date=`newtime ${date} ${timestep}` 
if ( "${date}" != "${finaldate}" ) goto loop

exit 0
