#! /bin/sh
function output(){
  cat << EOF > ../input/input.data
#	U3m/s,	aW,	DTs,	TSTART,	TSTOP,	printDts
	${W},	0,	1e-1,	0,	288000,	-60
#	P(hPa),	T(K),	RH(%),	Z(m),	ac,	aT
	${P}.	${T}.	100.	0.	${ac}	0.7
#	aerosoltype	aerosol_distribution_file_name	nbin
	lognormal	${ccn_type}			800
#	kappa_file
	kap${COM}
EOF
}
#	powerlaw	${ccn_type}
#      for ccn_type in ccn1D.data ccn2D.data ccn3D.data ccn4D.data; do
#      for ccn_type in SM1 SM2 SM3 SM4 SM5; do
#      for ccn_type in ck_A ck_B ck_C ck_ccn1Fit; do
#                for ccn_type in ccn1.data ccn2.data ccn3.data \
#                                ccn4.data ccn5.data ccn6.data \
#                                ccn7.data \
#                                ccn1D.data ccn2D.data ccn3D.data ccn4D.data ; do
#                                ccn1D.data ccn2D.data ccn3D.data ccn4D.data ; do

#acs='1.0 0.1 0.01 0.001'
#in foffice
#acs='1.0'
#in foffice2
acs='0.04'
#in fmine
#acs='0.005'
foruse='fit'
#foruse='test'
#foruse='gccn'

case $foruse in
fit)
    COMs='ASd NACLd HCLd Sd Nd'
    ambcase=1
    ccntycase=1
;;
test)
    COMs='t0 t1 t2 t3 t4 t5 t6 t7 t8'
    ambcase=2
    ccntycase=3
;;
gccn)
    COMs='ASd '
    ambcase=3
    ccntycase=2
;;
esac
case $ambcase in
1)
    Ws='0.1 0.21544347 0.46415888 1.0 2.15443469 4.64158883 10.'
    Ps='1000 900 800 700 600 500'
    Ts='-5 0 5 10 15 20 25 30'
;;
2)
    Ws='.1 .2 .5 1. 2. 5. 10 '
    Ps='925 850 750 650 550'
    Ts='-2 4 10 16'
;;
3)
    Ws='.1 .2 .5 1. 2. 5. 10 '
    Ps='1000 '
    Ts='16.85 '
;;
esac

case $ccntycase in
1)
    ccntys='SM1 SM2 SM3 SM4 SM5 '
    ccntys+='ccn1.data ccn2.data ccn3.data '
    ccntys+='ccn4.data ccn5.data '
    ccntys+='ccn6.data ccn7.data ccn8.data'
#    ccntys+='ccn8.data ccn9.data ccna.data '
;;
2)
    ccntys='DM1R0sxNx  DM1R7sxNx  DM2RxsxN4  DM3Rxs1Nx '
    ccntys+='DM1R1sxNx  DM1R8sxNx  DM2RxsxN5  DM3Rxs2Nx '
    ccntys+='DM1R2sxNx  DM1R9sxNx  DM2RxsxN6  DM3Rxs3Nx '
    ccntys+='DM1R3sxNx  DM2RxsxN0  DM2RxsxN7  DM3Rxs4Nx '
    ccntys+='DM1R4sxNx  DM2RxsxN1  DM2RxsxN8  DM3Rxs5Nx '
    ccntys+='DM1R5sxNx  DM2RxsxN2  DM2RxsxN9  DM3Rxs6Nx '
    ccntys+='DM1R6sxNx  DM2RxsxN3  DM3Rxs0Nx  DM3Rxs7Nx '
;;
3)
    ccntys='ccn1D.data ccn2D.data ccn3D.data ccn4D.data ccn5D.data '
    ccntys+='ccn1H.data ccn2H.data ccn3H.data ccn4H.data ccn5H.data '
    ccntys+='SM6R02N01  SM6R02N03  SM6R02N10 '
    ccntys+='SM6R05N01  SM6R05N03  SM6R05N10 '
    ccntys+='SM6R10N01  SM6R10N03  SM6R10N10 '
;;
esac

output="../output/exp_set_01_07_19_${foruse}ac${acs}_"

for COM in ${COMs}; do 
    for ac in ${acs}; do 
        for W in ${Ws}; do
            for P in ${Ps}; do
                for T in ${Ts}; do
                    for ccn_type in ${ccntys} ; do
                        output
                        ./mainAS
                    done
                done
            done
        done
    done
    mkdir ${output}${COM}
    mv ../output/*.csv ${output}${COM}/.
done
