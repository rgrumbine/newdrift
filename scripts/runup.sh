#!/bin/bash 
###WCOSS2
#PBS -N driftup
#PBS -o driftup
#PBS -j oe
#PBS -A ICE-DEV
#PBS -q dev
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=1
##ursa
#SBATCH -J rtofsdrifta
#SBATCH -e rtofsdrifta.err
#SBATCH -o rtofsdrifta.out
#SBATCH -t 7:55:00
#SBATCH -q batch
#SBATCH -A marine-cpu
#SBATCH -N 1
#SBATCH --mem=3g
#SBATCH --mail-type FAIL
#SBATCH --mail-user robert.grumbine@noaa.gov

#Robert Grumbine
#27 May 2026

set -xe
pid=$$

##Wcoss2
mkdir -p /lfs/h2/emc/ptmp/wx21rg/devdrift.$pid
cd /lfs/h2/emc/ptmp/wx21rg/devdrift.$pid

##ursa
#mkdir -p /scratch3/NCEPDEV/stmp/wx21rg/devdrift.$pid
#cd /scratch3/NCEPDEV/stmp/wx21rg/devdrift.$pid

export PDY=20260601
export COMOUT=$HOME/noscrub/devdrift
export end=`date +"%Y%m%d"`
#export end=20260331

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
