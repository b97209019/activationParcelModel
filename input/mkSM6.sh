#! /bin/sh

function output(){
  cat << EOF > SM6R${R}N${N}
"SM6 R${R} N${N}"
  1                         'NCCN: number of CCN mode'
  `expr ${N} \* 100` 0.69314718055994529  .${R}e-6  'ZCCN[CM^-3],CCNSIGMA,CCNMODE[m]'
EOF
}

for R in 02 05 10 ; do
    for N in 01 03 10 ; do
        output
    done
done
