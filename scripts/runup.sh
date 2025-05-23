#!/bin/bash 
#PBS -N driftup
#PBS -o driftup
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1


cd $HOME/rgdev/newdrift/scripts

export tag=20250510
export COMOUT=$HOME/noscrub/newdrift
export end=`date +"%Y%m%d"`

while [ $tag -le $end ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time ./try.sh > ${tag}.out
  else
    echo zzz have $tag already
  fi
  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
