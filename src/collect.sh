#!/bin/bash
BASEDIR="/mnt/ECHO/bookie/RESULTS"

rm -rf /mnt/ECHO/bookie/RESULTS/*


getlast() { find $1 -type d -name "$2-??" | sort| tail -1;} 

do_bucket() {
BUCKET=$1
SOLVENT=$2
RXN=$3

for CALC in $(cd $BUCKET; find . -mindepth 1 -maxdepth 1 -type d | sed 's/^..//')
do
	TARGETLABEL=$(echo $CALC | sed -r 's/.*[-_]*(sn2|e2ts)_//')
	CALCLABEL=$(echo $CALC | sed -r 's/[-_](sn2|e2ts)_/-/;s/(sn2|e2ts)_//')
	CALCDIR="$BASEDIR/$BUCKET/$CALCLABEL"

	mkdir -p $CALCDIR
	echo $CALC
	
	for i in $(python3 ~/workcopies/tsdb/src/workflow.py $RXN $SOLVENT /mnt/ECHO/bookie/$BUCKET/$CALC | grep 'COMPLETED' | sed 's/.*COMPLETED//')
	do
		# PROGRESS
		echo "$BUCKET $CALCLABEL $i" >> $BASEDIR/completed.txt
		
		# GEOMETRIES
		LASTDIR=$(getlast /mnt/ECHO/bookie/$BUCKET/$CALC $i)
		test -f $LASTDIR/inp.xyz && cp $LASTDIR/inp.xyz $CALCDIR/$i.xyz
		test -f $LASTDIR/run.xyz && cp $LASTDIR/run.xyz $CALCDIR/$i.xyz
		test -f $LASTDIR/reactant.xyz && cp $LASTDIR/reactant.xyz $CALCDIR/r.xyz
		test -f $LASTDIR/product.xyz && cp $LASTDIR/product.xyz $CALCDIR/p.xyz

		# ENERGIES
		LOGFILE=$LASTDIR/run.log
		test -f $LOGFILE || continue
		if [[ "$i" =~ ^nf.*$ ]]
		then
			# get free energy
			TS_kcal=$(grep 'sn= 3  qrot/sn=' $LOGFILE | tail -1 | awk '{ print $9 }')
			H_har=$(grep 'Total enthalpy                    ...' $LOGFILE  | tail -1 | awk '{ print $4}')
			G=$(echo "scale=5;($H_har)-($TS_kcal)/627.5095" | bc)
			echo $BUCKET $CALCLABEL $TARGETLABEL 'MP2/6-311G(d)//MP2/6-311G(d)' $i $G >> $BASEDIR/gibbs-free-energies.txt
		else
			# get total energy
			E_har=$(grep "FINAL SINGLE POINT ENERGY" $LOGFILE | tail -1| awk '{ print $5}')
			echo $BUCKET $CALCLABEL $TARGETLABEL 'MP2/6-311G(d)//MP2/6-311G(d)' $i $E_har >> $BASEDIR/total-electronic-energies.txt
		fi

		# HOMO/LUMO by Marco Bragato
		LOGFILE=$LASTDIR/run.log
		ne=$(grep " Number of Electrons    NEL             ...." $LOGFILE | head -1 | awk '{print $6}')
		lumo=$(( $ne / 2 | bc -l))
		homo=$(( $ne / 2 - 1 | bc -l))
		E_lumo=$(grep " $lumo   0.0000   " $LOGFILE | tail -1 | awk '{print $3}')
		E_homo=$(grep " $homo   2.0000   " $LOGFILE | tail -1 | awk '{print $3}')
		gap=$( echo $E_lumo- $E_homo | bc -l)
		echo $BUCKET $CALCLABEL $TARGETLABEL 'MP2/6-311G(d)//MP2/6-311G(d)' $i "$E_homo $E_lumo $gap" >> $BASEDIR/homo-lumo.txt


		#Dipole moment by Marco Bragato
		LOGFILE=$LASTDIR/run.log
		dipole_x_y_z=$(grep "Total Dipole" $LOGFILE | tail -1 | awk '{print $5 "     " $6 "     " $7}')
		dipole_tot=$(grep -A2 "Total Dipole" $LOGFILE | tail -1 | awk '{print $4}')
		echo $BUCKET $CALCLABEL $TARGETLABEL 'MP2/6-311G(d)//MP2/6-311G(d)' $i "$dipole_x_y_z $dipole_tot" >> $BASEDIR/dipole-moment.txt
		
	done
done

}

do_bucket e2 gasphase e2
#do_bucket e2i water e2
do_bucket sn2 gasphase sn2
#do_bucket sn2i water sn2

# Interaction energies
(cd /mnt/ECHO/bookie/conformer/reactants; . extract.sh)
cat /mnt/ECHO/bookie/conformer/reactants/*.scanpoints > /mnt/ECHO/bookie/conformer/reactants/tmp.summary
python3 code/InteractionEnergyScanReactants.py /mnt/ECHO/bookie/conformer/reactants/tmp.summary > $BASEDIR/interaction-energies.txt
./code/InteractionEnergyGeometries.sh
python3 code/GetBarriers.py

# Conformer energies and geometries
./code/GetConformerEnergies.sh


# Small atoms calculations
grep "MP2 TOTAL ENERGY" singleatoms/*/run.log | sed 's/singleatoms.//;s/_/ /;s/\/.*-/ MP2\/6-311G\(d\)\/\/MP2\/6-311G\(d\) -/;s/ Eh//' > $BASEDIR/atomisation-electronic.txt

#for fn in $(find interaction -name run.log );
#do
#        grep -H "FINAL SINGLE POINT ENERGY" $fn | tail -1;
#done | sed 's/.*interaction.//;s/.runs./ /;s/\// /g;s/run.* /MP2\/6-311G\(d\)\/\/MP2\/6-311G\(d\) /;s/reactants/reactant/' > RESULTS/interaction-energies.txt
# Interaction geometries
#for fn in $(find interaction -name run.log ); 
#do 
#	DIR=$(dirname $fn)
#	NEWFILE=$(echo $DIR | sed 's/^/RESULTS\//;s/runs.//;s/$/.xyz/')
#	NEWDIR=$(dirname $NEWFILE)
#	mkdir -p $NEWDIR
#	cp $DIR/inp.xyz $NEWFILE
#done

# Upload
#tar czf RESULTS.tgz RESULTS/
#scp RESULTS.tgz alchemy:share-nobackup/
#rm RESULTS.tgz
