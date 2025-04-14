#!/bin/sh
#
f=00
#set -x
while [ $f -le 192 ]
do
  if [ -f drift_f0$f.nc ] ; then
    if [ ! -f f$f ] ; then
      python3 /u/robert.grumbine/rgdev/newdrift/sorc/scan_buoy.py drift_f0$f.nc >  f$f
      sort -nr -k 4 f$f > g$f
    fi
  fi

  f=`expr $f + 1`
  if [ $f -lt 10 ] ; then
    f=0$f
  fi
done 
