#!/bin/bash

DIRECTORY="../../../python27"

if [ ! -d "$DIRECTORY" ]; then
  mkdir $DIRECTORY
fi
 
virtualenv $DIRECTORY
# csh $DIRECTORY/bin/activate.csh
source $DIRECTORY/bin/activate
pip install -r dependencies.txt
deactivate