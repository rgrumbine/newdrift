#!/bin/bash 
#ursa
#SBATCH -J newdrift
#SBATCH -e newdrift.err
#SBATCH -o newdrift.out
#SBATCH -t 7:55:00
#SBATCH -q batch
#SBATCH -A marine-cpu
#SBATCH -N 1
#SBATCH --mail-type FAIL
#SBATCH --mail-user USER@system

###WCOSS2
##PBS -N driftup
##PBS -o driftup
##PBS -j oe
##PBS -A ICE-DEV
##PBS -q dev
##PBS -l walltime=6:00:00
##PBS -l select=1:ncpus=1

set -xe
pid=$$
mkdir -p /scratch3/NCEPDEV/stmp/wx21rg/newdrift.$pid
cd /scratch3/NCEPDEV/stmp/wx21rg/newdrift.$pid

export tag=20251201
export COMOUT=$HOME/noscrub/newdrift
export end=`date +"%Y%m%d"`
export end=20251231

while [ $tag -le $end ]
do
  if [ ! -d $COMOUT/$tag ] ; then
    time $HOME/rgdev/newdrift/scripts/rtofs.sh > ${tag}.out
    #rm *.nc
    #time ./rtofs.sh 
  else
    echo zzz have $tag already
  fi
  tag=`expr $tag + 1`
  tag=`$HOME/bin/dtgfix3 $tag`
done
