#!/bin/bash 
##ursa
#SBATCH -J devdrift
#SBATCH -e devdrift.err
#SBATCH -o devdrift.out
#SBATCH -t 7:55:00
#SBATCH -q batch
#SBATCH -A marine-cpu
#SBATCH -N 1
#SBATCH --mem=3g
##Wcoss2
##PBS -N driftup
##PBS -o driftup
##PBS -j oe
##PBS -A ICE-DEV
##PBS -q dev
##PBS -l walltime=6:00:00
##PBS -l select=1:ncpus=1


set -x

cd $HOME/rgdev/devdrift/scripts

export PDY=20250601
export COMOUT=$HOME/noscrub/devdrift_retro
if [ ! -d $COMOUT ] ; then
  mkdir -p $COMOUT
fi

export end=`date +"%Y%m%d"`
export end=20260409
#export end=$PDY

while [ $PDY -le $end ]
do
  if [ ! -d $COMOUT/$PDY ] ; then
    time ./retro.sh > ${PDY}.out
  else
    echo zzz have $PDY already
  fi

  PDY=`expr $PDY + 1`
  PDY=`$HOME/bin/dtgfix3 $PDY`
done
