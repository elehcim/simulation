#!/bin/bash
SIMNAME=sim60003
SNAP=snapshot_0065

periax=( 200.0 200.0  598.0 400.0 10.0  ) 
apoax=(  600.0 1000.0 600.0 600.0 600.0 )
radius=( 600.0 800.0  600.0 600.0 600.0 )

# Loop commad inspired from here: https://stackoverflow.com/a/28725562/1611927
# idx is simply looping in [0:len(periax)-1]
for idx in "${!periax[@]}"; do
	p=${periax[$idx]}
	a=${apoax[$idx]}
	r=${radius[$idx]}
	# TODO: in the future use logging or a way to use the gic name as computed by kick_ic.py 
	OUTFILE="${SIMNAME}_${SNAP}.kicked_p${p}_a${a}_r${r}_c8.15.gic"
	echo "Writing outfile $OUTFILE"
	LOGFILE="${OUTFILE%.gic}.log"
	echo $LOGFILE
	COMMAND="python kick_ic.py ~/sim/MoRIA/$SIMNAME/$SNAP --prefix=$SIMNAME -a $a -r $r -p $p"
	# save the command as part of the output
	echo $COMMAND > $LOGFILE
	# Run it
	eval $COMMAND >> $LOGFILE
done
