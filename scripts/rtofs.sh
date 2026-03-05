#!/bin/bash 
#wcoss2:
#PBS -N newdrift
#PBS -o driftout
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev_transfer
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1

#wcoss2:
module load intel netcdf
module load prod_envir wgrib2
COMIN=$HOME/noscrub/model_intercompare/rtofs_cice/rtofs.$PDY/

#macos: COMIN=/Volumes/Data/rtofs/

#ursa:
#module load intel-oneapi-compilers
#module load hpc-x/2.18.1-icc
#module load netcdf-c/4.9.2
#module load netcdf-fortran/4.6.1
#export NETCDF=$NETCDF_FORTRAN_ROOT
#module list
#ursa: COMIN=$HOME/clim_data/rtofs/rtofs.$PDY/

#initialize
#drift_in -- file with full 6 values drifters, set to -99 for i,j,clat, clon
#  at readin/initialize, set i,j and clat/clon = ilat/ilon

set -xe

cp $HOME/rgdev/newdrift/fix/merged.nc drift_in.nc
#cp $HOME/rgdev/newdrift/fix/skiles_pts.nc drift_in.nc

#Loop:
#forecast hours 000 to 072 by 1
#forecast hours 072 to 192 by 3

EXDIR=${EXDIR:-$HOME/rgdev/devdrift/exec}
PDY=${PDY:-20260101}

hhh=000
# Pick up from partial run:
#cp drift_f010.nc drift_in.nc
#hhh=011
while [ $hhh -le 192 ] 
#while [ $hhh -le 024 ] 
#while [ $hhh -le 001 ] 
do
  fname=rtofs_glo_2ds_f${hhh}_ice.nc
  if [ ! -f ${COMIN}/$fname ] ; then
    echo could not find ${COMIN}/$fname
    exit 1
  #else
  #  ls -l ${COMIN}/$fname 
  fi

  echo \'${COMIN}/$fname\' > runin 
  echo \'drift_in.nc\'   >> runin
  echo \'out_${hhh}.nc\' >> runin
  if [ $hhh -lt 72 ] ; then
    export dt=1
  else
    export dt=3
  fi
  dtsec=`expr $dt \* 3600 `
  echo $dtsec >> runin 
  echo 1      >> runin
  echo 1      >> runin
  if [ $hhh -eq 000 ] ; then
    echo .FALSE. >> runin
  else
    echo .TRUE. >> runin
  fi
  cat $HOME/rgdev/devdrift/scripts/rtofs.vars >> runin
  echo runin | time $EXDIR/drifter

  cp out_${hhh}.nc drift_f${hhh}.nc
  mv out_${hhh}.nc drift_in.nc

  hhh=`expr $hhh + $dt`
  if [ $hhh -lt 10 ] ; then
    hhh=00$hhh
  elif [ $hhh -lt 100 ] ; then
    hhh=0$hhh
  fi

done
#endloop

exit
#mv outputs to $com
if [ -f drift_f192.nc ] ; then
  mkdir -p $COMOUT/$PDY
  mv *.nc ${PDY}.out $COMOUT/$PDY
fi
