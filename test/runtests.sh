#!/bin/bash

tar -xf validationdata.tar.gz
for t in *.cpp; do
# Unzip mammut simulation data
	cd ../src/external/mammut/test/archs && tar -xf repara.tar.gz
	cd ../../../../../test 
	./$(basename "$t" .cpp)
    exitvalue=$?
	cd ../src/external/mammut/test/archs && rm -rf repara
	cd ../../../../../test
    if [ $exitvalue -ne 0 ]; then
        break
    fi
done
rm -rf validationdata

exit $exitvalue
