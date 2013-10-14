#! /bin/ksh

rMin=$1
dr=$2
rMax=$3

laMin=$4
dla=$5
laMax=$6

nrun=5

r=$rMin
while [[ $r -le $rMax ]] ; do
	la=$laMin
	while [[ $la -le $laMax ]] ; do
		printf "r: %8f la:%8f\n" "$r" "$la"
	la=`echo "" | awk ' {print '$la'+'$dla'} '`
	done
r=`echo "" | awk ' {print '$r'+'$dr'} '`
done
