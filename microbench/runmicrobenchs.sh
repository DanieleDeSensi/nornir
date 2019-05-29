#!/bin/bash

# When a new parameter is added to the archdata.xml configuration file.
# 1. Change the CONFIGURATION_VERSION #define in parameters.cpp file.
# 2. Insert another condition in the if for checking if the new field
#    has been correctly created (e.g. 'grep "ticksPerNs" ... ')

if [ "$(id -u)" != "0" ]; then
	echo "============================ ATTENTION ============================" 
	echo "| You need to run this command with sudo (this is needed in order |"
	echo "| to read some architecture specific information and to place the |"
	echo "| configuration files in system directories).                     |"
	echo "==================================================================="
	exit 1
fi

CURRENT_VERSION=$(grep "define CONFIGURATION_VERSION" ../src/parameters.cpp | cut -d ' ' -f 3 | tr -d '"')

CONFPATH_ROOT=$XDG_CONFIG_DIRS
if [ -z "$CONFPATH_ROOT" ];
then
	CONFPATH_ROOT="/etc/xdg"
else
	CONFPATH_ROOT=$(echo $XDG_CONFIG_DIRS | cut -d ':'  -f 1)
fi

echo "Creating "$CONFPATH_ROOT"/nornir folder."
mkdir -p $CONFPATH_ROOT"/nornir"

CONFPATH_FILE=$CONFPATH_ROOT"/nornir/archdata.xml"
CONFPATH_VOLTAGE=$CONFPATH_ROOT"/nornir/voltage.csv"
CONFPATH_VERFILE=$CONFPATH_ROOT"/nornir/version.csv"

if grep $CURRENT_VERSION "$CONFPATH_VERFILE" >/dev/null 2>&1 && \
   grep "ticksPerNs" "$CONFPATH_FILE" >/dev/null 2>&1 && \
   grep "idlePower" "$CONFPATH_FILE" >/dev/null 2>&1 && \
   [ -f $CONFPATH_VOLTAGE ] && \
   [ -f $CONFPATH_FILE ]; 
then
	echo "Configuration files for version "$CURRENT_VERSION " already created."
else
	echo "Creating configuration files for version "$CURRENT_VERSION
	echo $CURRENT_VERSION > $CONFPATH_VERFILE
	rm -rf $CONFPATH_FILE
	# Header
	echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" >> $CONFPATH_FILE
	echo "<archData>" >> $CONFPATH_FILE
	
	# TickPerNs
	echo -n "    " >> $CONFPATH_FILE
	./ticksPerNs >> $CONFPATH_FILE
	
	# Voltage table
	if [ ! -f $CONFPATH_VOLTAGE ]; then
		if [ ! -f voltageTable.txt ]; then
			./voltageTable
		fi
		cp voltageTable.txt $CONFPATH_VOLTAGE
	fi

        # Idle power
        echo -n "    " >> $CONFPATH_FILE
        ./idlePower >> $CONFPATH_FILE
	
	# Footer
	echo "</archData>" >> $CONFPATH_FILE
	chmod -R ugo+r $CONFPATH_ROOT"/nornir/"
fi

