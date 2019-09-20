#!/bin/bash

BINARY_DIR=$1
SOURCE_DIR=$2
TEST=$3

if [ ! -d ${BINARY_DIR}/test/validationdata ]; then
	tar -xf ${SOURCE_DIR}/test/validationdata.tar.gz -C ${BINARY_DIR}/test/
fi

if [ ! -d ${BINARY_DIR}/test/archconfig ]; then
	cp -r ${SOURCE_DIR}/test/archconfig ${BINARY_DIR}/test/
fi

mkdir -p ${BINARY_DIR}/test/mammut-test/archs/

rm -rf ${BINARY_DIR}/test/mammut-test/archs/* && tar -xf ${BINARY_DIR}/src/mammut_repo-prefix/src/mammut_repo/test/archs/repara.tar.gz -C ${BINARY_DIR}/test/mammut-test/archs/ && tar -xf ${BINARY_DIR}/src/mammut_repo-prefix/src/mammut_repo/test/archs/c8.tar.gz -C ${BINARY_DIR}/test/mammut-test/archs/
${TEST}



