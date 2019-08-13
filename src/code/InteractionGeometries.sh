#!/bin/bash
cat /mnt/ECHO/bookie/conformer/reactants/tmp.summary | while read line; do
	Y=$(echo $line | cut -d"/" -f3 | cut -d"-" -f2)
	TARGET=$(echo $line | cut -d"/" -f1 | sed "s/_0/_$Y/")
	CONFID=$(echo $line | cut -d"/" -f2 | sed "s/gogp//;s/-.*//")
	RXN=$(echo $line | cut -d"/" -f3 | sed 's/-.*//')
	DIST=$(echo $line | cut -d"/" -f3 | sed 's/.*-//')
	DIRNAME="RESULTS/conformer/reactants/$TARGET/conf-$CONFID"
	OUTNAME="$DIRNAME/$RXN-$DIST.xyz"
	SRCNAME=$(echo $line | sed 's/run.log.*/inp.xyz/')
	mkdir -p $DIRNAME
	cp "conformer/reactants/$SRCNAME" $OUTNAME
done
