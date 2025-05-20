#!/bin/bash 
#PBS -N driftup
#PBS -o driftup
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1


cd $HOME/rgdev/newdrift/scripts

export tag=20250501
export COMOUT=$HOME/noscrub/newdrift

while [ $tag -le 20250518 ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time ./try.sh > ${tag}.out
  fi
  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
