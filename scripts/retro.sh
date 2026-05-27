#!/bin/bash 
#PBS -N devdrift
#PBS -o driftout
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev_transfer
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1

#Robert Grumbine
#27 May 2026

set -xe

PDY=${PDY:-20250601}

##wcoss2
#module load intel netcdf
#module load prod_envir wgrib2
#COMIN=$HOME/noscrub/model_intercompare/rtofs_cice/rtofs.$PDY/
#COMIN=$HOME/noscrub/retros/gfs.$PDY/00/model/ice/history/

#ursa:
 . /etc/profile.d/modules.sh
module load intel-oneapi-compilers
module load hpc-x/2.18.1-icc
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1
export NETCDF=$NETCDF_FORTRAN_ROOT
#COMIN=$HOME/clim_data/rtofs/rtofs.$PDY/
COMIN=$HOME/clim_data/stream4/gfs.${PDY}/00/model/ice/history/

#macos: COMIN=/Volumes/Data/rtofs/

#cd $HOME/rgdev/devdrift/sorc/

#initialize
#drift_in -- file with full 6 values drifters, set to -99 for i,j,clat, clon
#  at readin/initialize, set i,j and clat/clon = ilat/ilon
cp $HOME/rgdev/devdrift/fix/merged.nc drift_in.nc
#cp $HOME/rgdev/devdrift/fix/skiles_pts.nc drift_in.nc

#Loop:
#forecast hours 000 to 384 by 6

EXDIR=${EXDIR:-$HOME/rgdev/devdrift/exec}
PDY=${PDY:-20260501}

hhh=006
count=0
# Pick up from partial run:
#cp drift_f010.nc drift_in.nc
#hhh=011

while [ $hhh -le 384 ] 
#while [ $hhh -le 024 ] 
#while [ $hhh -le 000 ] 
do
  fname=gfs.t00z.6hr_avg.f${hhh}.nc
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
    export dt=6
  else
    export dt=6
  fi
  dtsec=`expr $dt \* 3600 `
  echo $dtsec >> runin 
  echo 1      >> runin
  echo 1      >> runin
  if [ $count -eq 0 ] ; then  # If count == 0, cold start
    echo .FALSE. >> runin
  else
    echo .TRUE. >> runin
  fi
  cat $HOME/rgdev/devdrift/scripts/ufs.vars >> runin
  echo runin | time $EXDIR/drifter

  cp out_${hhh}.nc drift_f${hhh}.nc
  mv out_${hhh}.nc drift_in.nc

  hhh=`expr $hhh + $dt`
  if [ $hhh -lt 10 ] ; then
    hhh=00$hhh
  elif [ $hhh -lt 100 ] ; then
    hhh=0$hhh
  fi
  count=`expr $count + 1`
done
#endloop


#mv outputs to $com
if [ -f drift_f024.nc ] ; then
  mkdir -p $COMOUT/$PDY
  mv *.nc ${PDY}.out $COMOUT/$PDY
fi
