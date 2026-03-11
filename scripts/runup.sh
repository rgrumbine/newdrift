#!/bin/bash 
###WCOSS2
#PBS -N driftup4
#PBS -o driftup4
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1

#ursa
#SBATCH -J devdrift09
#SBATCH -e devdrift09.err
#SBATCH -o devdrift09.out
#SBATCH -t 7:55:00
#SBATCH -q batch
#SBATCH -A marine-cpu
#SBATCH -N 1
#SBATCH --mem=3g
#SBATCH --mail-type FAIL
#SBATCH --mail-user USER@system

set -xe
pid=$$

##Wcoss2
#mkdir -p /lfs/h2/emc/ptmp/wx21rg/devdrift.$pid
#cd /lfs/h2/emc/ptmp/wx21rg/devdrift.$pid
##ursa
mkdir -p /scratch3/NCEPDEV/stmp/wx21rg/devdrift.$pid
cd /scratch3/NCEPDEV/stmp/wx21rg/devdrift.$pid

export PDY=20260226
export COMOUT=$HOME/noscrub/devdrift
export end=`date +"%Y%m%d"`
export end=20260226


while [ $PDY -le $end ]
do
  if [ ! -d $COMOUT/$PDY ] ; then
    time $HOME/rgdev/devdrift/scripts/rtofs.sh > ${PDY}.out
    #rm *.nc
  else
    echo zzz have $PDY already
  fi
  PDY=`expr $PDY + 1`
  PDY=`$HOME/bin/dtgfix3 $PDY`
done
