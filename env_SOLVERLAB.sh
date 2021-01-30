#!/bin/bash

export CDMATH_INSTALL=@CMAKE_INSTALL_PREFIX@
export COREFLOWS_INSTALL=@CMAKE_INSTALL_PREFIX@
export PETSC_DIR=@PETSC_INSTALL@
export PETSC_ARCH=@PETSC_ARCH@
export PETSC_INCLUDES=@PETSC_INCLUDES_INSTALL@
export PETSC_LIBRARIES=@PETSC_DIR@/@PETSC_ARCH@/lib
export MEDFILE_ROOT_DIR=@MEDFILE_ROOT_DIR@
export MEDFILE_INCLUDE_DIRS=@MEDFILE_INCLUDE_DIRS@
export MEDFILE_LIBRARIES=@MEDFILE_LIBRARIES_INSTALL@
export MEDCOUPLING_ROOT_DIR=@MEDCOUPLING_ROOT_DIR@
export MEDCOUPLING_INCLUDE_DIR=@MEDCOUPLING_INCLUDE_DIR@
export MEDCOUPLING_LIBRARIES=@MEDCOUPLING_LIBRARIES@
export PV_LIB_DIR=@PV_LIB_DIR@
export PV_PYTHON_DIR=@PV_PYTHON_DIR@

if [ @SOLVERLAB_WITH_MPI@ = ON ]
then
  export mpirun=${PETSC_DIR}/${PETSC_ARCH}/bin/mpirun
  export mpiexec=${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec
fi

export COREFLOWS=$COREFLOWS_INSTALL/bin/Executable/COREFLOWSMainExe
export LD_LIBRARY_PATH=$COREFLOWS_INSTALL/lib:$CDMATH_INSTALL/lib:${PETSC_DIR}/lib:$PETSC_LIBRARIES:${MEDCOUPLING_LIBRARIES}:${MEDFILE_ROOT_DIR}/lib:$MEDFILE_LIBRARIES:${MEDFILE_C_LIBRARIES}:$PV_LIB_DIR:${LD_LIBRARY_PATH}
export PYTHONPATH=$COREFLOWS_INSTALL/lib/coreflows:$COREFLOWS_INSTALL/bin/coreflows:$COREFLOWS_INSTALL/lib/python:$CDMATH_INSTALL/lib/cdmath:$CDMATH_INSTALL/bin/cdmath:$CDMATH_INSTALL/bin/cdmath/postprocessing:$PV_PYTHON_DIR:${PETSC_DIR}/${PETSC_ARCH}/lib:${PETSC_DIR}/lib:$PETSC_LIBRARIES:${MEDCOUPLING_LIBRARIES}:${MEDFILE_ROOT_DIR}/lib:$MEDFILE_LIBRARIES:${MEDFILE_C_LIBRARIES}:${PYTHONPATH}

export SOLVERLAB_ROOT_DIR=$COREFLOWS_INSTALL
export SOLVERLABGUI=$COREFLOWS_INSTALL/bin/salome/COREFLOWS_Standalone.py
