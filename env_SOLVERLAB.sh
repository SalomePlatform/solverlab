#!/bin/bash

# this bash file have to be sourced 

##### utilities #####

RED="\e[0;31m"
GREEN="\e[0;32m"
YELLOW="\e[0;33m"
BOLD="\e[1m"
NC="\e[0m" # No Color

function f_red {
  echo -e ${RED}${@}${NC}
}

function f_green {
  echo -e ${GREEN}${@}${NC}
}

function f_info {
  f_green "INFO : ${@}"
}

function f_error {
  f_red "ERROR : ${@}"
  # exit 1  # ko avoid if sourced file
}

function f_warning {
  f_red "WARNING : ${@}"
}

##### main #####

export SOLVERLAB_INSTALL=@CMAKE_INSTALL_PREFIX@
export PETSC_DIR=@PETSC_INSTALL@
export PETSC_ARCH=@PETSC_ARCH@
export PETSC_INCLUDES=@PETSC_INCLUDES_INSTALL@
export PETSC_LIBRARIES=@PETSC_DIR@/@PETSC_ARCH@/lib
export PETSC4PY_ROOT_DIR=@PETSC4PY_ROOT_DIR@
export SLEPC4PY_ROOT_DIR=@SLEPC4PY_ROOT_DIR@
export MEDFILE_ROOT_DIR=@MEDFILE_ROOT_DIR@
export MEDFILE_INCLUDE_DIRS=@MEDFILE_INCLUDE_DIRS@
export MEDFILE_LIBRARIES=@MEDFILE_LIBRARIES_INSTALL@
export MEDCOUPLING_ROOT_DIR=@MEDCOUPLING_ROOT_DIR@
export MEDCOUPLING_INCLUDE_DIR=@MEDCOUPLING_INCLUDE_DIR@
export MEDCOUPLING_LIBRARIES=@MEDCOUPLING_LIBRARIES@
export PV_LIB_DIR=@PV_LIB_DIR@
export PV_PYTHON_DIR=@PV_PYTHON_DIR@
export HDF5_ROOT=@HDF5_ROOT@

export SOLVERLAB=${SOLVERLAB_INSTALL}/bin/Executable/COREFLOWSMainExe

export LD_LIBRARY_PATH=\
${SOLVERLAB_INSTALL}/lib:\
${PETSC_DIR}/lib:\
${PETSC_DIR}/${PETSC_ARCH}/lib:\
${MEDCOUPLING_LIBRARIES}:\
${MEDFILE_ROOT_DIR}/lib:\
${HDF5_ROOT}/lib:\
${PV_LIB_DIR}:\
${LD_LIBRARY_PATH}

export PYTHONPATH=\
${SOLVERLAB_INSTALL}/bin:\
${SOLVERLAB_INSTALL}/bin/cdmath:\
${SOLVERLAB_INSTALL}/bin/cdmath/postprocessing:\
${SOLVERLAB_INSTALL}/bin/coreflows:\
${SOLVERLAB_INSTALL}/lib:\
${SOLVERLAB_INSTALL}/lib/cdmath:\
${SOLVERLAB_INSTALL}/lib/coreflows:\
${PETSC_DIR}/${PETSC_ARCH}/lib:\
${PETSC_DIR}/lib:\
${PETSC4PY_ROOT_DIR}:\
${PETSC4PY_ROOT_DIR}/lib/${PETSC_ARCH}:\
${SLEPC4PY_ROOT_DIR}:\
${SLEPC4PY_ROOT_DIR}/lib/${PETSC_ARCH}:\
${MEDCOUPLING_LIBRARIES}:\
${MEDCOUPLING_LIBRARIES}/python@Python_VERSION_MAJOR@.@Python_VERSION_MINOR@/site-packages/:\
${MEDFILE_ROOT_DIR}/lib:\
${HDF5_ROOT}/lib:\
${PV_PYTHON_DIR}:\
${PYTHONPATH}


if [ @SOLVERLAB_WITH_MPI@ = ON ] # test SOLVERLAB_WITH_MPI
then
  export mpirun=$PETSC_DIR/$PETSC_ARCH/bin/mpirun
  export mpiexec=$PETSC_DIR/$PETSC_ARCH/bin/mpiexec
  export MPI4PY_ROOT_DIR=@MPI4PY_ROOT_DIR@
  export PYTHONPATH=$MPI4PY_ROOT_DIR:$PYTHONPATH
fi

if [ @PRELOAD_NETCDF_ON_UBUNTU20@ = ON ]  # test PRELOAD_NETCDF_ON_UBUNTU20
then
    export  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libnetcdf.so
fi

export SOLVERLAB_ROOT_DIR=${SOLVERLAB_INSTALL}

f_info SOLVERLAB_ROOT_DIR=${SOLVERLAB_INSTALL}

# TODO 230308 Have to be tested because needs a preexisting python3 context with PyQt5/numpy/matplotlib/paraview etc. compatible!
# this is PROTOTYPING of dir PACKAGESPY to get solverlabGUI (standalone installation, not in salome context)

# PACKAGESPY in salome usage with PACKAGESPY_ROOT_DIR set
if [ @SOLVERLAB_WITH_GUI@ = ON ]; then
  if [ -z ${PACKAGESPY_ROOT_DIR} ]; then # test existing PACKAGESPY_ROOT_DIR to get solverlabGUI
    f_error "solverlabGUI installation needs env var PACKAGESPY_ROOT_DIR set"
  else
    if [ ! -d ${PACKAGESPY_ROOT_DIR} ]; then # test existing dir PACKAGESPY to get solverlabGUI standalone 
      f_error "solverlabGUI installation needs ${PACKAGESPY_ROOT_DIR} directory existing"
    fi
  fi
  if [ -z ${SOLVERLABGUI_ROOT_DIR} ]; then # test existing solverlabGUI
    f_error "solverlabGUI installation needs env var SOLVERLABGUI_ROOT_DIR set"
  else
    if [ ! -d ${SOLVERLABGUI_ROOT_DIR} ]; then # test existing SOLVERLABGUI_ROOT_DIR to get solverlabGUI
      f_error "solverlabGUI installation needs ${SOLVERLABGUI_ROOT_DIR} directory existing"
    fi 
  fi
  f_info SOLVERLABGUI_ROOT_DIR=${SOLVERLABGUI_ROOT_DIR}
fi


# solverlabGUI needs something like that
# export PYTHONPATH=${PACKAGESPY_ROOT_DIR}/packagespy:${PYTHONPATH}
# export PATH=${SOLVERLAB_ROOT_DIR}/solverlabGUI:${PATH}
# echo solverlabGUI installation Wambeke seem ok"
# main initial python script to get standalone GUI
# echo "To launch solverlabGUI type '\${SOLVERLAB_ROOT_DIR}/solverlabGUI --gui'"

