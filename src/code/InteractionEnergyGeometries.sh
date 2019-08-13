#!/bin/bash
cat RESULTS/electronic-forward-barriers.txt | while read line; do 
	# python3 code/DoPlaceReactant.py conformer/reactants/C_A_C_E_A_0/gogp43-03/ D 3.2 sn2
 # sn2 A_B_A_B_C_D A_B_A_B_C_D MP2/6-311G(d)//MP2/6-311G(d) -5633.1349192384705 reactant 12 3.6351690000000003 MP2/6-311G(d)//MP2/6-311G(d) -5633.157872 0.02295276152926817
	TARGET=$(echo $line | cut -d " " -f 3)
	REACTANT=$(echo $line | cut -d " " -f 3 | sed 's/_.$/_0/')
	CONFORMER_INT=$(echo $line | cut -d " " -f 5)
	printf -v CONFORMER "%02d" $CONFORMER_INT
	Y=$(echo $line | cut -d " " -f 3 | sed 's/.*_//')
	DIST=$(echo $line | cut -d " " -f 6)
	RXN=$(echo $line | cut -d " " -f 1)
	
	mkdir -p RESULTS/conformer/reactants/$TARGET/
	python3 code/DoPlaceReactant.py $(find conformer/reactants/$REACTANT/gogp${CONFORMER}-*  -name "run.xyz"  | sort | tail -1 | sed 's/.run.xyz//') $Y $DIST $RXN > RESULTS/conformer/reactants/$TARGET/${RXN}-min.xyz
done
