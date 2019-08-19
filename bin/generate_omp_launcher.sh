#!/bin/sh
BUILD_DIR=$1
INSTALL_DIR=$2

echo "#!/bin/sh"                                                                              >  ${BUILD_DIR}/bin/nornir_openmp
echo "export LD_LIBRARY_PATH=${INSTALL_DIR}/lib/:\${LD_LIBRARY_PATH}"                         >> ${BUILD_DIR}/bin/nornir_openmp
echo "export LD_PRELOAD=${INSTALL_DIR}/lib/libomp_nornir.so:${INSTALL_DIR}/lib/libnornir.so"  >> ${BUILD_DIR}/bin/nornir_openmp
echo "COMMAND=\$1"                                                                            >> ${BUILD_DIR}/bin/nornir_openmp
echo "PARAMS=\$2"                                                                             >> ${BUILD_DIR}/bin/nornir_openmp
echo "export NORNIR_OMP_PARAMETERS=\${PARAMS}"                                                >> ${BUILD_DIR}/bin/nornir_openmp
echo "\${COMMAND}"                                                                            >> ${BUILD_DIR}/bin/nornir_openmp
