#!/bin/bash

#This script prepares TURBOMOLE optimization with B97-D3/def2-SVP (+COSMO) (+GCP)
# (change defaults at the top of the script)
#usage: tmol_mkinp.sh <xyz file> <charge> [cosmo]

#############
title="B97-D3/def2-SVP OPT"
charge="0"
basis="def2-SVP"
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
\$end
EOF

x2t tm.xyz > coord
if [ $? -ne 0 ]; then
   echo "Abbormal termination from x2t" && exit 1
fi

define << EOF
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
  cosmoprep << EOF
$3










r all b
*

EOF
  if [ $? -ne 0 ]; then
     echo "Abbormal termination from cosmoprep" && exit 1
  fi
fi

# run opt
jobex -c 100 -energy 6 -gcart 3 -gexp 3 
