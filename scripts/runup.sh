#!/bin/bash 
#PBS -N driftup2
#PBS -o driftup2
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1

#set -xe

cd $HOME/rgdev/newdrift/scripts

export tag=20260101
export COMOUT=$HOME/noscrub/newdrift
export end=`date +"%Y%m%d"`
#export end=20260101

while [ $tag -le $end ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time ./rtofs.sh > ${tag}.out
   # rm *.nc
    #time ./rtofs.sh 
  else
    echo zzz have $tag already
  fi
  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
