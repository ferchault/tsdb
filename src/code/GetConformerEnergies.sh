#!/bin/bash
find_success_conformers() {
        for i in $(find $1 -type d -name 'gogp*' | sort -n); do grep "HURRAY" $i/run.log &> /dev/null && echo $i; done
}

for dir in conformer/reactants/*/
do
	MOL=$(basename $dir)
	MOLDIR="RESULTS/conformer/reactants/$MOL"
	mkdir -p $MOLDIR
	for conformer in $(find_success_conformers $dir)
	do
		FIRST=$conformer
		CONFID=$(basename $conformer | sed 's/gogp//;s/-.*//')
		LAST=$(echo $conformer | sed 's/-[0-9][0-9]\//-00/')
		echo "reactant VERTICAL" $(grep -H "FINAL SINGLE POINT ENERGY" $FIRST/run.log | head -n1)
		echo "reactant RELAXED" $(grep -H "FINAL SINGLE POINT ENERGY" $LAST/run.log  | tail -n1)
		cp $FIRST/inp.xyz $MOLDIR/conf-${CONFID}-vertical.xyz
		cp $LAST/run.xyz $MOLDIR/conf-${CONFID}-relaxed.xyz
	done | sed 's/conformer.reactants.//;s/\/[^/]*-/ /;s/\/.* / /'
	grep -H conf $dir/cs-00/generate.log | sed 's/conformer.reactants./reactant UFF /;s/\/.*conf-/ /'
done | awk '{ $5 = ($2 == "UFF" ? $5/627.509: $5)} 1' | sed 's/UFF/UFF\/\/UFF/;s/VERTICAL/MP2\/6-311G\(d\)\/\/UFF/;s/RELAXED/MP2\/6-311G\(d\)\/\/MP2\/6-311G\(d\)/' > RESULTS/conformer-energies.txt

