#!/bin/bash 
#PBS -N driftup
#PBS -o driftup
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1


cd $HOME/rgdev/newdrift/scripts

export tag=20250101
export COMOUT=$HOME/noscrub/newdrift_retro
if [ ! -d $COMOUT ] ; then
  mkdir -p $COMOUT
fi

export end=`date +"%Y%m%d"`
export end=$tag
export enc=20250110

while [ $tag -le $end ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time ./retro.sh > ${tag}.out
  else
    echo zzz have $tag already
  fi
  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
