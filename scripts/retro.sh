#!/bin/bash 
#PBS -N newdrift
#PBS -o driftout
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev_transfer
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1

#set -xe

#initialize
#drift_in -- file with full 6 values drifters, set to -99 for i,j,clat, clon
#  at readin/initialize, set i,j and clat/clon = ilat/ilon
module load intel netcdf
module load prod_envir wgrib2
#module list

#cd $HOME/rgdev/newdrift/sorc/

cp $HOME/rgdev/newdrift/fix/merged.nc drift_in.nc
#cp $HOME/rgdev/newdrift/fix/skiles_pts.nc drift_in.nc
#cp $HOME/rgdev/newdrift/fix/seaice_edge.nc drift_in.nc

#Loop:
#forecast hours 000 to 072 by 1
#forecast hours 072 to 192 by 3

EXDIR=${EXDIR:-$HOME/rgdev/newdrift/exec}
tag=${tag:-20250101}
#macos: base=/Volumes/Data/rtofs/
#hera: base=$HOME/clim_data/rtofs/rtofs.$tag/
#wcoss2:
#base=$HOME/noscrub/model_intercompare/rtofs_cice/rtofs.$tag/
base=$HOME/noscrub/retros/gfs.$tag/00/model/ice/history/

hhh=006
count=0
# Pick up from partial run:
#cp drift_f010.nc drift_in.nc
#hhh=011
while [ $hhh -le 240 ] 
#while [ $hhh -le 024 ] 
#while [ $hhh -le 000 ] 
do
  fname=gfs.ice.t00z.6hr_avg.f${hhh}.nc
  if [ ! -f ${base}/$fname ] ; then
    echo could not find ${base}/$fname
    exit 1
  #else
  #  ls -l ${base}/$fname 
  fi

  echo \'${base}/$fname\' > runin 
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
  echo runin | time $EXDIR/drifter_retro

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
if [ -f drift_f192.nc ] ; then
  mkdir -p $COMOUT/$tag
  mv *.nc ${tag}.out $COMOUT/$tag
fi
