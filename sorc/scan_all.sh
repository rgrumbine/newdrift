#!/bin/sh
#
f=000
#set -x
while [ $f -le 192 ]
do
  if [ -f drift_f$f.nc ] ; then
    if [ ! -f f$f ] ; then
      python3 /u/robert.grumbine/rgdev/newdrift/sorc/scan_buoy.py drift_f$f.nc >  f$f
      sort -nr -k 4 f$f | grep -iv e > g$f
    fi
  fi

  f=`expr $f + 1`
  if [ $f -lt 100 ] ; then
    f=0$f
    if [ $f -lt 10 ] ; then
      f=0$f
    fi
  fi
done 
