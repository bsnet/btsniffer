#!/bin/bash

LAP="219F87"      # your target LAP here
FOLDER="folder"   # your folder name here
TRACENAME="trace" # your trace name here

# You might need to tweak the following paths
DELAYFILE="${FOLDER}/synchtime${TRACENAME}.txt"
LOTRACEFILE="${FOLDER}/lo${TRACENAME}.dat"
HITRACEFILE="${FOLDER}/hi${TRACENAME}.dat"

EXECPATH="./sources/frame-processing"
TIMESTART=0
REPEAT=0

while [ true ]; do

  OUTPUTFILE="${FOLDER}/output${TRACENAME}_REPEAT${REPEAT}.txt"
  if [ -f $OUTPUTFILE ]; then
    echo "Output file ${OUTPUTFILE} already present, check please..."
    exit 1
  fi

  DELAYVALS=$(./getdelay.pl ${FOLDER}/synchtime${TRACENAME}.txt $TIMESTART)
  if [ "$DELAYVALS" = "" ]; then
    echo "Empty vals, quitting..."
    exit 1
  fi
  
  DELAY=$(echo $DELAYVALS | awk '{ print $2 }')
  
  echo "Starting at $TIMESTART with delay $DELAY"

  $EXECPATH/btdecoder $LOTRACEFILE $HITRACEFILE $DELAY $TIMESTART $LAP > $OUTPUTFILE
  
  if [ "$?" = "0" ]; then
    echo "Terminated successfully"
  else
    echo "Aborted"
    exit 1
  fi

  TIMEFOUND=0
  RES=$(cat $OUTPUTFILE | ./terminate.pl ${LAP})
  if [[ $RES =~ \(([0-9]+)\) ]]; then
    TIMEFOUND=${BASH_REMATCH[1]}
  fi
  
  if [ "$TIMEFOUND" = "0" ]; then
    echo "found at zero, quitting..."
    break
  fi
  
  TIMESTART=$(($TIMESTART+$TIMEFOUND))
  
  echo "Ended at $TIMEFOUND"

  REPEAT=$((REPEAT+1))

done
