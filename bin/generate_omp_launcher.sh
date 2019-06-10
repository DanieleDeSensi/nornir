#!/bin/sh
BUILD_DIR=$1
INSTALL_DIR=$2

echo "#!/bin/sh"                                                                              >  ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "export LD_LIBRARY_PATH=${INSTALL_DIR}/lib/:\${LD_LIBRARY_PATH}"                         >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "export LD_PRELOAD=${INSTALL_DIR}/lib/libomp_nornir.so:${INSTALL_DIR}/lib/libnornir.so"  >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "COMMAND=\$1"                                                                            >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "PARAMS=\$2"                                                                             >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "export NORNIR_OMP_PARAMETERS=\${PARAMS}"                                                >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh
echo "\${COMMAND}"                                                                            >> ${BUILD_DIR}/bin/nornir_omp_launcher.sh 