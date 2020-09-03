cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project (CDMATH)
set (CDMATH_VERSION_MAJOR 1)
set (CDMATH_VERSION_MINOR 0)

include(ExternalProject)#For PETSc, MED and MEDCoupling

# Project options
option (CDMATH_WITH_PETSC "Compile CDMATH with PETSc linking." OFF)
option (CDMATH_WITH_DOCUMENTATION "Generate documentation with doxygen." OFF)
option (CDMATH_WITH_PYTHON "Compile Python interface for CDMATH." OFF)
option (CDMATH_WITH_POSTPRO "Install postprocessing Python modules." OFF)
option (CDMATH_WITH_PACKAGE "Generate RPM, Debian and tarball packages." OFF)
option (CDMATH_WITH_TESTS "Compile unit testing." OFF)

#Path to installed libraries
set (PETSC_DIR            OFF CACHE STRING "PETSc library path" )
set (SLEPC_DIR            OFF CACHE STRING "SLEPc library path" )
set (F2CBLASLAPACK_DIR    OFF CACHE STRING "F2CBLASLAPACK library path" )
set (HDF5_ROOT            OFF CACHE STRING "HDF5 library path" )#For the compilation of MED with a specific version of HDF5
set (MEDFILE_ROOT_DIR     OFF CACHE STRING "MED library path" )
set (MEDCOUPLING_ROOT_DIR OFF CACHE STRING "MEDCoupling library path" )

#url of tarball librairies
set (DOWNLOAD_PETSC         http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.13.2.tar.gz
     CACHE STRING           "PETSc tarball path/url" )
set (DOWNLOAD_SLEPC         https://slepc.upv.es/download/distrib/slepc-3.13.2.tar.gz
     CACHE STRING           "SLEPC tarball path/url" )
set (DOWNLOAD_F2CBLASLAPACK ${CDMATH_SOURCE_DIR}/pre_requis/PETSc/f2cblaslapack-3.4.2.q4.tar.gz
     CACHE STRING           "F2CBLASLAPACK tarball path/url" )
set (DOWNLOAD_HDF5          ${CDMATH_SOURCE_DIR}/pre_requis/PETSc/hdf5-1.10.3.tar.gz  
     CACHE STRING           "HDF5 tarball path/url" )
set (DOWNLOAD_MED           http://files.salome-platform.org/Salome/other/med-4.0.0.tar.gz
     CACHE STRING           "MED tarball path/url")
set (DOWNLOAD_MEDCOUPLING   http://files.salome-platform.org/Salome/other/medCoupling-9.4.0.tar.gz
     CACHE STRING           "MEDCoupling tarball path/url" )


# Base directories
set (BASE_DIR ${CDMATH_SOURCE_DIR}/base)
set (MESH_DIR ${CDMATH_SOURCE_DIR}/mesh)
set (LINEARSOLVER_DIR ${CDMATH_SOURCE_DIR}/linearsolver)
set (CDMATH_SWIG_DIR ${CDMATH_SOURCE_DIR}/swig)
set (CDMATH_POSTPRO_DIR ${CDMATH_SOURCE_DIR}/postprocessing)

set (TESTS_DIR ${CDMATH_SOURCE_DIR}/tests)
list (APPEND CMAKE_MODULE_PATH "${CDMATH_SOURCE_DIR}/cmake_files")

# PETSc and HDF5
if (CDMATH_WITH_PETSC OR PETSC_DIR OR DEFINED ENV{PETSC_DIR} )

  if ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR)
    if   (NOT PETSC_DIR)
      set(PETSC_DIR $ENV{PETSC_DIR})
      message ( STATUS "Checking variable PETSC_DIR : " $ENV{PETSC_DIR} )
    else (NOT PETSC_DIR)
      message ( STATUS "Checking variable PETSC_DIR : " ${PETSC_DIR} )
    endif(NOT PETSC_DIR)

    find_package (PETSc 3.4 REQUIRED)
    petsc_get_version ()

    message ( STATUS "PETSc found. Version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" )
    set(PETSC_INSTALL ${PETSC_DIR})

    #Define and search slepc variables
    if   ( NOT SLEPC_DIR )
      if   ( DEFINED ENV{SLEPC_DIR} )
        set(SLEPC_DIR $ENV{SLEPC_DIR})
      else ( DEFINED ENV{SLEPC_DIR} )
        set(SLEPC_DIR ${PETSC_DIR}/${PETSC_ARCH})
      endif( DEFINED ENV{SLEPC_DIR} )
    endif( NOT SLEPC_DIR)

    message ( STATUS "Checking variable SLEPC_DIR" )
    if ( IS_DIRECTORY ${SLEPC_DIR}/include AND EXISTS ${SLEPC_DIR}/lib/libslepc.so)
      set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
      set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
      message( STATUS "SLEPc found at ${SLEPC_DIR}" )
    else()
      message( FATAL_ERROR "SLEPc not found at ${SLEPC_DIR}" )
    endif()

    #define hdf5 variables
    if   ( NOT HDF5_ROOT )#hdf5 is not defined in cmake arguments
      if  ( DEFINED ENV{HDF5_ROOT} )#hdf5 is defined in the environment
        set(HDF5_ROOT $ENV{HDF5_ROOT})
      else( DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment
        #HDF5 to be found in petsc external packages
        set(HDF5_ROOT ${PETSC_DIR}/${PETSC_ARCH})# define hint for hdf5/med installation in petsc external packages
      endif( DEFINED ENV{HDF5_ROOT} )
    endif( NOT HDF5_ROOT )
    set(HDF5_LIBRARY_DIR ${HDF5_ROOT}/lib)
    set(HDF5_INCLUDE_DIRS ${HDF5_ROOT}/include)
  else ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR )
    message ( STATUS "PETSC_DIR not set, searching PETSc and SLEPc in the system" )

    if ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )  #Case fedora/redhat system install
      message ( STATUS "PETSC includes found in /usr/include/petsc/" )
      message ( STATUS "PETSC library found in /usr/lib64" )
      set(PETSC_DIR /usr/)
      set(PETSC_INCLUDES  /usr/include/petsc /usr/include/petsc/petsc/mpiuni)
      set(PETSC_LIBRARIES /usr/lib64/libpetsc.so)

      set(PETSC_VERSION "3.8") #3.8 for fedora 26 , 3.9 for fedora 29 , 3.10 for fedora 30, , 3.12 for fedora 32
      
      set(CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, petsc-devel (>= 3.4)") # This is not fully working yet. 
      set(PETSC_INSTALL ${PETSC_DIR})

      #Define and search slepc variables
      if   ( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )
        message( STATUS "SLEPc includes found in /usr/include/slepc/" )
        message( STATUS "SLEPc library found in /usr/lib64/slepc/" )
        set(SLEPC_DIR /usr/)
        set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
        set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)
        set (CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, slepc-devel (>= 3.4)") # This is not fully working yet. 
      else ( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )
        message( FATAL_ERROR "SLEPc not found in the system" )
      endif( IS_DIRECTORY "/usr/include/slepc/" AND EXISTS "/usr/lib64/libslepc.so" )

      #HDF5 to be found in the system
      if   ( NOT HDF5_ROOT AND NOT DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment
        find_package(HDF5 REQUIRED)
        set(CPACK_DEBIAN_PACKAGE_DEPENDS "libhdf5-dev")
        set(CPACK_RPM_PACKAGE_REQUIRES "hdf5-devel")
      endif( NOT HDF5_ROOT AND NOT DEFINED ENV{HDF5_ROOT} )


    #elseif ( IS_DIRECTORY "/usr/lib/petsc/include" AND EXISTS "/usr/lib/petsc/lib/libpetsc_real.so" ) #Case ubuntu/debian system install
    #  message ( STATUS "PETSc includes found in /usr/lib/petsc/include")
    #  message ( STATUS "PETSc library  found in /usr/lib/petsc/lib/libpetsc_real.so")
    #  set(PETSC_DIR /usr/lib/petsc/)
    #  set(PETSC_INCLUDES  ${PETSC_DIR}/include ${PETSC_DIR}/include/petsc/mpiuni)
    #  set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc_real.so)

    #  find_package (PETSc 3.4 REQUIRED)
    #  petsc_get_version ()
    #  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev   (>= 3.4)") # This is not fully working yet. 

    #  Define and search slepc variables
    #  if   ( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )
    #    message( STATUS "SLEPc includes found in /usr/lib/slepc/include" )
    #    message( STATUS "SLEPc library  found in /usr/lib/slepc/lib/libslepc_real.so" )
    #    set(SLEPC_DIR /usr/lib/slepc/)
    #    set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
    #    set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc_real.so)
    #    set (CPACK_DEBIAN_PACKAGE_REQUIRES   "${CPACK_DEBIAN_PACKAGE_REQUIRES}, slepc-devel (>= 3.4)") # This is not fully working yet. 
    #  else ( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )
    #    message( FATAL_ERROR "SLEPc not found in the system" )
    #  endif( IS_DIRECTORY "/usr/lib/slepc/include" AND EXISTS "/usr/lib/slepc/lib/libslepc_real.so" )

    #elseif ( IS_DIRECTORY "/usr/local/lib/python2.7/dist-packages/petsc") #Case ubuntu/debian system pip install

    else ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )  # No petsc in system, do compile petsc along with slepc as an external package
      message (STATUS "PETSC not found in the system")
      message (STATUS "PETSC will be downloaded and compiled from ${DOWNLOAD_PETSC}" )
      message (STATUS "SLEPC will be downloaded and compiled from ${DOWNLOAD_SLEPC}" )
      message (STATUS "HDF5  will be downloaded and compiled from ${DOWNLOAD_HDF5}" )

      #extract tarball name
      string(LENGTH "${DOWNLOAD_PETSC}" tarball_url_length)#length of the tarball
      string(FIND   "${DOWNLOAD_PETSC}" "/" start_tarball_name REVERSE )# last occurence of "/"
      MATH(EXPR start_tarball_name "${start_tarball_name}+1")#start after the occurence of "/"
      MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
      string(SUBSTRING ${DOWNLOAD_PETSC} ${start_tarball_name} ${tarball_name_length} PETSC_TARBALL_NAME)

      set(PETSC_DIR ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${PETSC_TARBALL_NAME}) # Location of the final install 
      set(PETSC_INSTALL ${CMAKE_INSTALL_PREFIX}/share/petsc)#folder to copy petsc libraries and include files

      ExternalProject_Add (PETSc
      URL               ${DOWNLOAD_PETSC}   
      SOURCE_DIR        ${PETSC_DIR}
      BUILD_IN_SOURCE   TRUE
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${PETSC_INSTALL} --with-debugging=0 --with-mpi=0 --download-f2cblaslapack=${DOWNLOAD_F2CBLASLAPACK} --with-fc=0 --download-slepc=${DOWNLOAD_SLEPC} --download-hdf5=${DOWNLOAD_HDF5}
      BUILD_COMMAND     make
      TEST_COMMAND      make check
      INSTALL_COMMAND   make all
      INSTALL_DIR       ${PETSC_DIR}
      STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
      LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
      LOG_UPDATE        TRUE           # Wrap update in script to log output
      LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
      LOG_BUILD         TRUE           # Wrap build in script to log output
      LOG_TEST          TRUE           # Wrap test in script to log output
      LOG_INSTALL       TRUE           # Wrap install in script to log output
       )

      #file (STRINGS "${PETSC_DIR}/${PETSC_ARCH}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
      #message( STATUS "PETSc version is ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}" )
      message( STATUS "PETSc, SLEPc and HDF5 will be installed at ${PETSC_INSTALL}")

      string(SUBSTRING ${PETSC_TARBALL_NAME} 6 ${tarball_name_length} PETSC_VERSION)
      #define slepc variables
      set(SLEPC_DIR ${PETSC_INSTALL})
      set(SLEPC_INCLUDES ${SLEPC_DIR}/include)
      set(SLEPC_LIBRARIES ${SLEPC_DIR}/lib/libslepc.so)

      #define hdf5 variables
      if   ( NOT HDF5_ROOT )#hdf5 is not defined in cmake arguments
        if  ( DEFINED ENV{HDF5_ROOT} )#hdf5 is defined in the environment
          set(HDF5_ROOT $ENV{HDF5_ROOT})
        else( DEFINED ENV{HDF5_ROOT} )#hdf5 is not defined in the environment          
          set(HDF5_ROOT ${PETSC_INSTALL})#HDF5 to be found in petsc installation
        endif( DEFINED ENV{HDF5_ROOT} )
      endif( NOT HDF5_ROOT )
      set(HDF5_LIBRARY_DIR ${HDF5_ROOT}/lib)
      set(HDF5_INCLUDE_DIRS ${HDF5_ROOT}/include)

    endif ( IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so" )
  endif ( DEFINED ENV{PETSC_DIR} OR PETSC_DIR )

  if   ( ${PETSC_VERSION} VERSION_GREATER 3.5 )
    add_definitions(-DPETSC_VERSION_GREATER_3_5)
  endif( ${PETSC_VERSION} VERSION_GREATER 3.5 )
string(REPLACE ";" ":"  PETSC_INCLUDES_INSTALL "${PETSC_INCLUDES}")# use colon instead of semicolon in environment file env_CDMATH.sh
endif( CDMATH_WITH_PETSC OR PETSC_DIR  OR DEFINED ENV{PETSC_DIR} )

#MED
if( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})

  if   ( NOT MEDFILE_ROOT_DIR )
    set(MEDFILE_ROOT_DIR $ENV{MEDFILE_ROOT_DIR} )
  endif( NOT MEDFILE_ROOT_DIR )

  find_package (MEDFile REQUIRED)
  message (STATUS "MEDFile found in ${MEDFILE_ROOT_DIR}")

else( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})
  #string(FIND "${HDF5_LIBRARIES}" "libhdf5.so" pos)
  #string(SUBSTRING "${HDF5_LIBRARIES}" 0 ${pos} HDF5_LIBRARY_DIR)
  
  message(STATUS "MED will be downloaded and installed from ${DOWNLOAD_MED}")
  set(MACHINE PCLINUX)
  set(MEDFILE_DEFINITIONS "-D${MACHINE} -DMEDFILE_INSTALL_DOC=OFF")#Comment ne pas compiler la doc ?
  set(MEDFILE_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/share/med) # Location of the final install
  
  #extraction of the tarball archive name
  string(LENGTH "${DOWNLOAD_MED}" tarball_url_length)#length of the name of the tarball
  string(FIND "${DOWNLOAD_MED}" "/" start_tarball_name REVERSE )# last occurence of "/"
  MATH(EXPR start_tarball_name "${start_tarball_name}+1")#start after the occurence of "/"
  MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
  string(SUBSTRING ${DOWNLOAD_MED} ${start_tarball_name} ${tarball_name_length} MED_TARBALL_NAME)

  ExternalProject_Add (MED
        URL               ${DOWNLOAD_MED}  #location of med tarball
        SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${MED_TARBALL_NAME}
        BUILD_IN_SOURCE   TRUE
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${MEDFILE_ROOT_DIR}  --with-hdf5=${HDF5_ROOT} --with-hdf5-include=${HDF5_ROOT}/include --with-hdf5-lib=${HDF5_LIBRARY_DIR} --with-hdf5-bin=${HDF5_ROOT}/bin  --with-swig=yes 
        BUILD_COMMAND     make
        INSTALL_COMMAND   make install
        INSTALL_DIR       ${MEDFILE_ROOT_DIR}
        STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
        LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
        LOG_UPDATE        TRUE           # Wrap update in script to log output
        LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
        LOG_BUILD         TRUE           # Wrap build in script to log output
        LOG_TEST          TRUE           # Wrap test in script to log output
        LOG_INSTALL       TRUE           # Wrap install in script to log output
  )
  
  set(MEDFILE_INCLUDE_DIRS ${MEDFILE_ROOT_DIR}/include)# Nécessaire pour le medloader et les sous-dossiers mesh, et IJKMesh
  set(MEDFILE_C_LIBRARIES med medC)# Nécessaire pour le medloader

endif( MEDFILE_ROOT_DIR OR DEFINED ENV{MEDFILE_ROOT_DIR})

add_library(med   SHARED IMPORTED)
set_property(TARGET med  PROPERTY IMPORTED_LOCATION ${MEDFILE_ROOT_DIR}/lib/libmed.so)
add_library(medC   SHARED IMPORTED)
set_property(TARGET medC PROPERTY IMPORTED_LOCATION ${MEDFILE_ROOT_DIR}/lib/libmedC.so)

#MEDCoupling
if( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR})

  if( NOT MEDCOUPLING_ROOT_DIR )
    set(MEDCOUPLING_ROOT_DIR $ENV{MEDCOUPLING_ROOT_DIR} )
  endif( NOT MEDCOUPLING_ROOT_DIR )

  message (STATUS "Seeking MEDCoupling library in ${MEDCOUPLING_ROOT_DIR}")

  if( NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedloader.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedcoupling.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/include/MEDCoupling.hxx)
    message (FATAL_ERROR "MEDCoupling library not found in ${MEDCOUPLING_ROOT_DIR}")
  else( NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedloader.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedcoupling.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/include/MEDCoupling.hxx)
    message (STATUS "MEDCoupling library found in ${MEDCOUPLING_ROOT_DIR}")
  endif( NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedloader.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/lib/libmedcoupling.so OR NOT EXISTS ${MEDCOUPLING_ROOT_DIR}/include/MEDCoupling.hxx)

else( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR})
  message(STATUS "MEDCoupling will be downloaded and installed from ${DOWNLOAD_MEDCOUPLING}")
  set(MEDCOUPLING_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/share/medcoupling) # Location of the final install
  SET(ENV{HDF5_ROOT} ${HDF5_ROOT})#MEDCoupling install process seems to require an environment variable

  #extraction of the tarball archive name
  string(LENGTH "${DOWNLOAD_MEDCOUPLING}" tarball_url_length)#length of the tarball
  string(FIND   "${DOWNLOAD_MEDCOUPLING}" "/" start_tarball_name REVERSE )# last occurence of "/"
  MATH(EXPR start_tarball_name "${start_tarball_name}+1")#start after the occurence of "/"
  MATH(EXPR tarball_name_length "${tarball_url_length}-${start_tarball_name}-7")#name ends before .tar.gz that counts for 7 characters
  string(SUBSTRING ${DOWNLOAD_MEDCOUPLING} ${start_tarball_name} ${tarball_name_length} MEDCOUPLING_TARBALL_NAME)
  string(SUBSTRING ${MEDCOUPLING_TARBALL_NAME} 12 ${tarball_name_length} MEDCOUPLING_VERSION_NAME)

  message(STATUS "MEDCoupling version is ${MEDCOUPLING_VERSION_NAME}")

  ExternalProject_Add (MEDCoupling
        URL               ${DOWNLOAD_MEDCOUPLING}  #location of medcoupling tarball
        SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/pre_requis/${MEDCOUPLING_TARBALL_NAME}
#        BUILD_IN_SOURCE   FALSE
        CONFIGURE_COMMAND cmake <SOURCE_DIR>/${MEDCOUPLING_TARBALL_NAME} -DCMAKE_INSTALL_PREFIX=${MEDCOUPLING_ROOT_DIR} -DCONFIGURATION_ROOT_DIR=<SOURCE_DIR>/configuration-${MEDCOUPLING_VERSION_NAME} -DMEDFILE_ROOT_DIR=${MEDFILE_ROOT_DIR} -DMEDCOUPLING_ENABLE_PARTITIONER=OFF -DMEDCOUPLING_PARTITIONER_METIS=OFF -DMEDCOUPLING_PARTITIONER_SCOTCH=OFF -DMEDCOUPLING_ENABLE_RENUMBER=OFF -DMEDCOUPLING_BUILD_DOC=OFF -DHDF5_ROOT_DIR=${HDF5_ROOT}
        BUILD_COMMAND     make
        INSTALL_COMMAND   make install
        INSTALL_DIR       ${MEDCOUPLING_ROOT_DIR}
        STAMP_DIR         ${CMAKE_CURRENT_BINARY_DIR}/pre_requis
        LOG_DOWNLOAD      TRUE           # Wrap download in script to log output
        LOG_UPDATE        TRUE           # Wrap update in script to log output
        LOG_CONFIGURE     TRUE           # Wrap configure in script to log output
        LOG_BUILD         TRUE           # Wrap build in script to log output
        LOG_TEST          TRUE           # Wrap test in script to log output
        LOG_INSTALL       TRUE           # Wrap install in script to log output
  )

endif( MEDCOUPLING_ROOT_DIR OR DEFINED ENV{MEDCOUPLING_ROOT_DIR} )

set(MEDCOUPLING_INCLUDE_DIR ${MEDCOUPLING_ROOT_DIR}/include)
set(MEDCOUPLING_LIBRARIES   ${MEDCOUPLING_ROOT_DIR}/lib    )#for environment file env_CDMATH.sh 

add_library(medloader   SHARED IMPORTED )
set_property(TARGET medloader   PROPERTY IMPORTED_LOCATION ${MEDCOUPLING_ROOT_DIR}/lib/libmedloader.so)
add_library(medcoupling SHARED IMPORTED )
set_property(TARGET medcoupling PROPERTY IMPORTED_LOCATION ${MEDCOUPLING_ROOT_DIR}/lib/libmedcoupling.so)

if   (TARGET MED AND TARGET PETSc)
  ExternalProject_Add_StepDependencies( MED build PETSc )#PETSc doit être compilé avant MED car c'est PETSc qui compile HDF5
endif(TARGET MED AND TARGET PETSc)

if   (TARGET MED AND TARGET MEDCoupling)
  ExternalProject_Add_StepDependencies( MEDCoupling build MED )#MED doit être compilé avant MEDCoupling
endif(TARGET MED AND TARGET MEDCoupling)
 
if   (TARGET PETSc AND TARGET MEDCoupling)
  ExternalProject_Add_StepDependencies( MEDCoupling build PETSc )#PETSc doit être compilé avant MEDCoupling (car il contient hdf5)
endif(TARGET PETSc AND TARGET MEDCoupling)
 
# Paraview variables for env_CDMATH.sh
if   (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)
#find_package(VTK) 
#find_package(ParaView)  
  
  #extraction of the paraview version"
  IF   ( EXISTS "/usr/include/paraview/vtkPVConfig.h" )
    file(STRINGS /usr/include/paraview/vtkPVConfig.h vtkConfig)
  ELSE ( EXISTS "/usr/include/paraview/vtkPVConfig.h" )
    message(WARNING "Could not find ParaView configuration file vtkPVConfig.h. Postprocessing won't work")
  ENDIF( EXISTS "/usr/include/paraview/vtkPVConfig.h" )

  FOREACH(line ${vtkConfig})
    string(FIND "${line}" "#define PARAVIEW_VERSION_FULL " pos)
    IF(NOT ${pos} EQUAL -1)
      string(LENGTH ${line} line_length)#length of the tarball
      MATH(EXPR start_pv_version "${line_length}-6")#line ends with "x.y.z", that counts for 7 characters
      string(SUBSTRING ${line} ${start_pv_version} 5 PV_VERSION)
      break()
    ENDIF(NOT ${pos} EQUAL -1)
  ENDFOREACH(line vtkConfig)
  
  message(STATUS "ParaView version is ${PV_VERSION}" )
  IF   ( ${PV_VERSION} VERSION_GREATER 5.6.0 OR ${PV_VERSION} VERSION_EQUAL 5.6.0 )#Use python 3, use VERSION_GREATER_EQUAL if cmake >=3.7
    SET(PYTHON2OR3 "3")
  ELSE ( ${PV_VERSION} VERSION_GREATER 5.6.0 OR ${PV_VERSION} VERSION_EQUAL 5.6.0 )#Use python 2
    SET(PYTHON2OR3 "2")
  ENDIF( ${PV_VERSION} VERSION_GREATER 5.6.0 OR ${PV_VERSION} VERSION_EQUAL 5.6.0 )

    set (PV_LIB_DIR    /usr/lib/python${PYTHON2OR3}.7/dist-packages/paraview/:/usr/lib64/paraview/:/usr/lib/paraview/)
    set (PV_PYTHON_DIR /usr/lib/python${PYTHON2OR3}.7/dist-packages/paraview/:/usr/lib64/paraview/site-packages/:/usr/lib64/paraview/site-packages/paraview/:/usr/lib64/paraview/site-packages/vtk/:/usr/lib64/paraview/python${PYTHON2OR3}.7/site-packages/:/usr/lib64/paraview/python${PYTHON2OR3}.7/site-packages/paraview/:/usr/lib64/paraview/python${PYTHON2OR3}.7/site-packages/vtkmodules:/usr/lib/paraview/site-packages/paraview:/usr/lib/paraview/site-packages/paraview/vtk)

endif(CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)

# Swig interface
if   (CDMATH_WITH_PYTHON)
    IF   (${CMAKE_VERSION} VERSION_GREATER "3.12.0")
      find_package(Python ${PYTHON2OR3} REQUIRED COMPONENTS Interpreter Development )
      SET(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
    ELSE (${CMAKE_VERSION} VERSION_GREATER "3.12.0")
      find_package(PythonInterp ${PYTHON2OR3} REQUIRED )
      find_package(PythonLibs   ${PYTHON2OR3} REQUIRED )
      set(Python_LIBRARIES    ${PYTHON_LIBRARIES})
      set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
      set(Python_VERSION      ${PYTHON_VERSION_STRING})
    ENDIF(${CMAKE_VERSION} VERSION_GREATER "3.12.0")

  message(STATUS "Python version is ${Python_VERSION}")

  find_package (SWIG 3.0 REQUIRED)
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")
  set (CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES},   python-devel,      numpy, swig")
  add_subdirectory (${CDMATH_SWIG_DIR})
endif(CDMATH_WITH_PYTHON)

find_package (XDR REQUIRED)
add_subdirectory (${BASE_DIR})
add_subdirectory (${MESH_DIR})

if   (CDMATH_WITH_PETSC)
    add_subdirectory (${LINEARSOLVER_DIR})
endif(CDMATH_WITH_PETSC)

# Documentation
if (CDMATH_WITH_DOCUMENTATION)
  find_package (Doxygen)
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, doxygen, graphviz, mscgen")
  set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, doxygen, graphviz, mscgen")
  if   (NOT DOXYGEN_FOUND)
    message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
  endif(NOT DOXYGEN_FOUND)
  configure_file (Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
  							  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  if   ( IS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION doc/cdmath-dev-doc)
  endif( IS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html)
endif(CDMATH_WITH_DOCUMENTATION)

#Debug option : todo = sort flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -fprofile-arcs -ftest-coverage -lgcov")
if   (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  SET(CMAKE_EXE_LINKER_FLAGS "-pg -fprofile-arcs -ftest-coverage -lgcov")
  SET(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage")
endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")

# Tests
if   (CDMATH_WITH_TESTS)
  find_package (CPPUNIT REQUIRED)
  enable_testing ()
  add_subdirectory (${TESTS_DIR})#contains c++ and python tests
endif(CDMATH_WITH_TESTS)

# Packaging
if (CDMATH_WITH_PACKAGE)
  include (CPackLists.txt)
endif(CDMATH_WITH_PACKAGE)

# Postprocessing
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)
  add_subdirectory (postprocessing)
endif ()

# Configuration file
configure_file(
    ${PROJECT_SOURCE_DIR}/env_CDMATH.sh
    ${CMAKE_INSTALL_PREFIX}/env_CDMATH.sh
    @ONLY
)
													    

