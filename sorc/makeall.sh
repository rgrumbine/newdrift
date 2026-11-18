#!/bin/bash

#Robert Grumbine
#28 October 2026

#modules :

## wcoss2:
module load PrgEnv-intel/8.3.3
module load netcdf/4.7.4
module load intel-classic/2022.2.0.262
module load intel/19.1.3.304
ln -sf mk.wcoss2 mk.this

## gaea:
#module load intel-classic/2023.2.0
#module load cray-hdf5/1.12.2.11
#module load cray-netcdf/4.9.0.9
#module load intel/2025.2
#module load cray-hdf5/1.14.3.7
#module load cray-netcdf/4.9.0.17
#export NETCDF=$NETCDF_DIR
#ln -sf mk.gaea mk.this

## ursa:
#module load intel-oneapi-compilers
#module load hpc-x/2.18.1-icc
#module load netcdf-c/4.9.2
#module load netcdf-fortran/4.6.2
#export NETCDF=$NETCDF_FORTRAN_ROOT
#ln -sf mk.ursa mk.this

## hera:
#module load hpc/1.2.0  hpc-intel/2022.1.2
#module load netcdf/4.7.0
#ln -sf mk.hera mk.this

## desk
#ln -sf mk.desk mk.this
#-------------------------------------------
set -x

echo zzz NETCDF 
echo `env | grep  NETCDF`

make
if [ ! -d ../exec ] ; then
  mkdir ../exec
fi
mv drifter in_drifter ../exec
