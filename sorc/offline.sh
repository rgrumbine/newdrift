#!/bin/sh

set -xe

#initialize
#drift_in -- file with full 6 values drifters, set to -99 for i,j,clat, clon
#  at readin/initialize, set i,j and clat/clon = ilat/ilon

#cp ../fix/drift_ref.nc drift_in.nc
cp ../fix/merged.nc drift_in.nc

#Loop:
#forecast hours 000 to 072 by 1
#forecast hours 072 to 192 by 3

#base=$HOME/rtofs/
#macos: base=/Volumes/Data/rtofs/
#hera: base=$HOME/clim_data/rtofs/rtofs.20241101/
#wcoss2:
base=$HOME/noscrub/model_intercompare/rtofs_cice/rtofs.20250405/

#for hhh in 000
hhh=078
#while [ $hhh -le 192 ] 
while [ $hhh -le 078 ] 
do
  fname=rtofs_glo_2ds_f${hhh}_ice.nc

  #echo \'${base}/$fname\' > runin 
  echo \'cice.nc\' > runin
  echo \'drift_in.nc\' >> runin
  echo \'out_${hhh}.nc\' >> runin
  if [ $hhh -lt 72 ] ; then
    export dt=1
  else
    export dt=3
  fi
  echo $dt >> runin 
  echo 1   >> runin
  echo 1   >> runin
  echo runin | ./drifter 

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

#mv outputs to $com


#For inline:
#  subroutine that accepts initial (full) drifters, forcing, and dt
#    returns updated locations
#  accepts info on whether/where to write out

