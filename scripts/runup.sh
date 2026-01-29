#!/bin/bash 
#ursa
#SBATCH -J newdrift09
#SBATCH -e newdrift09.err
#SBATCH -o newdrift09.out
#SBATCH -t 7:55:00
#SBATCH -q batch
#SBATCH -A marine-cpu
#SBATCH -N 1
#SBATCH --mem=3g
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

export PDY=20250901
export COMOUT=$HOME/noscrub/newdrift
export end=`date +"%Y%m%d"`
export end=20250930

while [ $PDY -le $end ]
do
  if [ ! -d $COMOUT/$PDY ] ; then
    time $HOME/rgdev/newdrift/scripts/rtofs.sh > ${PDY}.out
    #rm *.nc
    #time ./rtofs.sh 
  else
    echo zzz have $PDY already
  fi
  PDY=`expr $PDY + 1`
  PDY=`$HOME/bin/dtgfix3 $PDY`
done
