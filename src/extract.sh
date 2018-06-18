#!/bin/bash
rm -f final-ts/*

echo "calc label energy" > data-electronic-energies.txt
for i in $(find 03-ts-validated/ -maxdepth 1 -mindepth 1 ); 
do
	CALC=$(basename $i)
	LABEL=$(echo $CALC | awk -F '-' '{ print $NF }')
	LASTLOG="$(find 03-ts-validated/$CALC -path '*/ts-*' -name run.log | sort | tail -1)"
	test -f "$LASTLOG" || continue
	echo $CALC $LABEL $(grep 'FINAL SINGLE POINT ENERGY' $LASTLOG | tail -1 | sed 's/.* //')
done >> data-electronic-energies.txt

echo "calc label enthalpy" > data-gibbs-enthalpy.txt
for i in $(find 03-ts-validated/ -maxdepth 1 -mindepth 1 ); 
do
	CALC=$(basename $i)
	LABEL=$(echo $CALC | awk -F '-' '{ print $NF }')
	LASTLOG="$(find 03-ts-validated/$CALC -path '*/nf-*' -name run.log | sort | tail -1)"
	test -f "$LASTLOG" || continue
	echo $CALC $LABEL $(grep 'Final Gibbs free enthalpy' $LASTLOG | sed 's/ Eh//;s/.* //' | tail -1)
	cp $(echo $LASTLOG | sed 's/run.log$/inp.xyz/') final-ts/$CALC.xyz
done >> data-gibbs-enthalpy.txt

for i in 0*
do
	#ls $i > data-$i.txt
	find $i -maxdepth 1 -mindepth 1 -type d  -printf '%f\n' > data-$i.txt
done 

for i in OBSOLETE/*
do
        basename $(tar tzf $i | head -n1)
done > data-obsolete.txt

rm data-geo-ts.xyz
for i in final-ts/*; 
do 
	LABEL=$(echo $(basename $i) | sed 's/.xyz//')
	cat $i | sed "2s/.*/$LABEL/" >> data-geo-ts.xyz
done
