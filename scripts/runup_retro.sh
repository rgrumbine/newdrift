#!/bin/bash 
#PBS -N driftup
#PBS -o driftup
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1

set -x

cd $HOME/rgdev/newdrift/scripts

export tag=20241211
export COMOUT=$HOME/noscrub/newdrift_retro
if [ ! -d $COMOUT ] ; then
  mkdir -p $COMOUT
fi

export end=`date +"%Y%m%d"`
export end=$tag
export end=20250117

while [ $tag -le $end ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time ./ufs.sh > ${tag}.out
  else
    echo zzz have $tag already
  fi

  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
