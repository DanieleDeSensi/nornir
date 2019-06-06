#!/bin/sh

# When a new parameter is added to the archdata.xml configuration file.
# 1. Change the CONFIGURATION_VERSION #define in parameters.cpp file.
# 2. Insert another condition in the if for checking if the new field
#    has been correctly created (e.g. 'grep "ticksPerNs" ... ')

if [ "$(id -u)" != "0" ]; then
	echo "============================ ATTENTION ============================" 
	echo "| The next steps need to be run with sudo (it is needed in order  |"
	echo "| to read some architecture specific information and to place the |"
	echo "| configuration files in system directories).                     |"
	echo "==================================================================="
fi

sudo $1/microbench/runmicrobenchs.sh $1

