#!/bin/bash

echo "Building synch processor..."
cd synch-processing ; ./compile.sh ; cd ..

echo "Building frame processor..."
cd frame-processing ; ./compile.sh ; cd ..

echo "Done"
