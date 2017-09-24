#!/bin/bash

#This script prepares TURBOMOLE optimization with B97-D3/def2-SVP (+COSMO) (+GCP)
# (change defaults at the top of the script)
#usage: tmol_mkinp.sh <xyz file> <charge> [cosmo]

#############
title="B97-D3/def2-SV(P) OPT"
charge="0"
basis="def2-SV(P)"
method="b97-d"
#method_gcp="dft/svp"
#############

if [ ! -z $1 ]; then
    charge=$1
fi

cat > control << EOF
\$redund_inp
 metric 1
\$symmetry c1
\$coord    file=coord
\$user-defined bonds    file=coord
\$scfiterlimit  200
\$scfconv   5
\$scfdamp   start=1.500  step=0.050  min=0.050
\$statpt
  itrvec      0
  radmax      0.3
  radmin      1.0d-4
  tradius     0.3
  threchange  5.0d-4
  thrmaxdispl 3.0d-2
  thrmaxgrad  3.0d-2
  thrrmsdispl 1.0d-2
  thrrmsgrad  1.0d-2
\$end
EOF

x2t tm.xyz > coord
if [ $? -ne 0 ]; then
   echo "Abbormal termination from x2t" && exit 1
fi

define &> /dev/null << EOF
$title
y
ired
*
b
all $basis
c
$bsse_atoms
*
eht
y
$charge
y
n
dft
func
$method
on
*
dsp
on
*
ri
on
*
marij

*
q
EOF

if [ $? -ne 0 ]; then
   echo "Abbormal termination from define" && exit 1
fi

if [[ $2 == "cosmo" ]]; then
  cosmoprep  &> /dev/null<< EOF
$3










r all b
*

EOF
  if [ $? -ne 0 ]; then
     echo "Abbormal termination from cosmoprep" && exit 1
  fi
fi

# run opt
jobex -c 50 &> /dev/null 
