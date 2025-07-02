#!/bin/sh

#Wcoss2
#module load PrgEnv-intel/8.3.3
#module load intel-classic/2022.2.0.262
#module load netcdf/4.7.4

# Gaea
module load intel-classic/2023.2.0
module load cray-hdf5/1.12.2.11
module load cray-netcdf/4.9.0.9
export NETCDF=$NETCDF_DIR

#debug: ensure a new implementation is built
rm all
if [ ! -f all ] ; then
  ./makeall.sh
fi
if [ ! -f all ] ; then
  echo could not build executable
  exit 1
fi

export COMin=$HOME/noscrub/model_intercompare/rtofs_cice/
export PDY=20230718

if [ -f output.nc ] ; then
  rm output.nc
fi

#for fh in n00 f24 f48 f72 f96 f120 f144 f168 f192
for fh in f24
do
  fname=rtofs.${PDY}/rtofs_glo.t00z.${fh}.cice_inst
  ln -sf ${COMin}/$fname cice.nc
  cat input | ./all
done

