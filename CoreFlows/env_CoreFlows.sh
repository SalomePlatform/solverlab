#!/bin/bash

export CDMATH_DIR=@CDMATH_DIR@
source $CDMATH_DIR/env_CDMATH.sh

export CoreFlows_INSTALL=@CMAKE_INSTALL_PREFIX@
export PETSC_DIR=@PETSC_DIR@
export PETSC_ARCH=@PETSC_ARCH@
export PETSC_INCLUDES=@PETSC_INCLUDES_PATH@
export PETSC_LIBRARIES=@PETSC_LIBRARIES@

#------------------------------------------------------------------------------------------------------------------- 
export CoreFlows=$CoreFlows_INSTALL/bin/Executable/CoreFlowsMainExe
export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:${PETSC_DIR}/lib:/usr/lib64/:$CoreFlows_INSTALL/lib:$PETSC_LIBRARIES:${LD_LIBRARY_PATH}
export PYTHONPATH=$CoreFlows_INSTALL/lib:$CoreFlows_INSTALL/lib/CoreFlows_Python:$CoreFlows_INSTALL/bin/CoreFlows_Python:$CoreFlows_INSTALL/lib/python2.7/site-packages/salome:${PYTHONPATH}
export CoreFlowsGUI=$CoreFlows_INSTALL/bin/salome/CoreFlows_Standalone.py
