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

export COREFLOWS=$COREFLOWS_INSTALL/bin/Executable/COREFLOWSMainExe
export LD_LIBRARY_PATH=$COREFLOWS_INSTALL/lib:$CDMATH_INSTALL/lib:${PETSC_DIR}/${PETSC_ARCH}/lib:$PETSC_LIBRARIES:${MEDCOUPLING_LIBRARIES}:${MEDFILE_ROOT_DIR}/lib:$MEDFILE_LIBRARIES:${MEDFILE_C_LIBRARIES}:$PV_LIB_DIR:${LD_LIBRARY_PATH}
export PYTHONPATH=$COREFLOWS_INSTALL/lib:$COREFLOWS_INSTALL/lib/coreflows:$COREFLOWS_INSTALL/bin/coreflows:$COREFLOWS_INSTALL/lib/python3.7/site-packages/salome:$CDMATH_INSTALL/lib/cdmath:$CDMATH_INSTALL/bin/cdmath:$CDMATH_INSTALL/bin/cdmath/postprocessing:$CDMATH_INSTALL:$PV_PYTHON_DIR:${PETSC_DIR}/${PETSC_ARCH}/lib:$PETSC_LIBRARIES:${MEDCOUPLING_LIBRARIES}:${MEDFILE_ROOT_DIR}/lib:$MEDFILE_LIBRARIES:${MEDFILE_C_LIBRARIES}:${PYTHONPATH}
export COREFLOWSGUI=$COREFLOWS_INSTALL/bin/salome/COREFLOWS_Standalone.py
