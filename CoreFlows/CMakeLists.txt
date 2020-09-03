cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project (CoreFlows C CXX)
set (COREFLOWS_VERSION_MAJOR 0)
set (COREFLOWS_VERSION_MINOR 1)

option (CMAKE_BUILD_TYPE "Debug or Release compilation" Release)
option (COREFLOWS_WITH_DOCUMENTATION "Generate documentation for CoreFlows use." OFF)
option (COREFLOWS_WITH_PYTHON "Compile Python interface for COREFLOWS." OFF)
option (COREFLOWS_WITH_GUI "Compile Graphic user Interface for COREFLOWS." OFF)
option (COREFLOWS_WITH_PACKAGE "Generate RPM, Debian and tarball packages." OFF)

#Path to installed libraries
set (PETSC_DIR            OFF CACHE STRING "PETSc library path" )
set (SLEPC_DIR            OFF CACHE STRING "SLEPc library path" )

enable_testing()											   #
# Debug mode activates testing and profiling --------------------------------------------------------------#
if (CMAKE_BUILD_TYPE STREQUAL Debug)									   #
													   #
 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg --coverage")						   #
 SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg --coverage")					   #
 SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg --coverage")				   #
													   #
endif ()												   #
													   #
#----------------------------------------------------------------------------------------------------------#

list (APPEND CMAKE_MODULE_PATH "${CoreFlows_SOURCE_DIR}/cmake_files")                                       #  

# CDMATH ---------------------------------------------------------------------------------------------------#
find_package (CDMATH REQUIRED)                                                                              #
#-----------------------------------------------------------------------------------------------------------#


# PETSC ----------------------------------------------------------------------------------------------------#
message ( STATUS "Checking variable PETSC_DIR : " $ENV{PETSC_DIR} )

if   ( NOT PETSC_DIR AND NOT DEFINED ENV{PETSC_DIR} AND IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so")  #Case fedora/redhat system install
  message ( STATUS "PETSC includes found in /usr/include/petsc/" )
  message ( STATUS "PETSC library found in /usr/lib64" )
  set(PETSC_DIR /usr/)
  set(PETSC_ARCH lib64)
  set(PETSC_INCLUDES  /usr/include/petsc /usr/include/petsc/petsc/mpiuni)
  set(PETSC_INCLUDES_PATH  /usr/include/petsc:/usr/include/petsc/petsc/mpiuni)
  set(PETSC_LIBRARIES /usr/lib64/libpetsc.so)
  set(PETSC_VERSION "3.7") #3.7 for fedora 25/26 , 3.9 for fedora 29 , 3.10 for fedora 30, , 3.12 for fedora 32
  set(CPACK_RPM_PACKAGE_REQUIRES   "${CPACK_RPM_PACKAGE_REQUIRES}, petsc-devel (>= 3.4)") # This is not fully working yet. 

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

#elseif ( IS_DIRECTORY "/usr/lib/petsc/") #Case ubuntu/debian system install
#  message ( STATUS "PETSC found in /usr/lib/petsc/" )
#  set(PETSC_DIR /usr/lib/petsc/)
#  set(PETSC_ARCH )
#  find_package (PETSc 3.4 REQUIRED)
#  petsc_get_version ()
#  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev   (>= 3.4)") # This is not fully working yet. 

#elseif ( IS_DIRECTORY "/usr/local/lib/python2.7/dist-packages/petsc") #Case ubuntu/debian system pip install
#  message ( STATUS "PETSC found in /usr/local/lib/python2.7/dist-packages/petsc" )
#  set(PETSC_DIR /usr/local/lib/python2.7/dist-packages/petsc/)
#  set(PETSC_ARCH lib)
#  set(PETSC_INCLUDES       $PETSC_DIR/include /usr/include/openmpi)
#  set(PETSC_INCLUDES_PATH  $PETSC_DIR/include:/usr/include/openmpi)
#  set(PETSC_LIBRARIES      $PETSC_DIR/lib/libpetsc.so)
#  set(PETSC_VERSION "3.8") #value for Ubuntu 16.04 

else ( NOT PETSC_DIR AND NOT DEFINED ENV{PETSC_DIR} AND IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so")
  if(NOT PETSC_DIR)
    set(PETSC_DIR $ENV{PETSC_DIR})
  endif(NOT PETSC_DIR)

  find_package (PETSc 3.4 REQUIRED)
  petsc_get_version ()
  string(REPLACE ";" ":"  PETSC_INCLUDES_PATH "${PETSC_INCLUDES}")# use colon instead of semicolon in environment file env_CoreFlows.sh

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

endif( NOT PETSC_DIR AND NOT DEFINED ENV{PETSC_DIR} AND IS_DIRECTORY "/usr/include/petsc/" AND EXISTS "/usr/lib64/libpetsc.so")

if (${PETSC_VERSION} VERSION_GREATER 3.5)
  add_definitions(-DPETSC_VERSION_GREATER_3_5)
endif ()

#-----------------------------------------------------------------------------------------------------------#



# Swig interface -------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_PYTHON)                                                                                  #
   find_package (PYTHON REQUIRED)                                                                           #
   find_package (SWIG 3.0 REQUIRED)        

   set (PV_LIB_DIR /usr/lib/python2.7/dist-packages/:/usr/lib/python3.7/dist-packages/:/usr/lib64/paraview/)
   set (PV_PYTHON_DIR /usr/lib/python2.7/dist-packages:/usr/lib64/paraview/site-packages/:/usr/lib64/paraview/site-packages/paraview/:/usr/lib64/paraview/site-packages/vtk/:/usr/lib64/paraview/python2.7/site-packages/:/usr/lib64/paraview/python2.7/site-packages/paraview/:/usr/lib64/paraview/python2.7/site-packages/vtkmodules:/usr/lib64/paraview/python3.7/site-packages/:/usr/lib64/paraview/python3.7/site-packages/paraview/:/usr/lib64/paraview/python3.7/site-packages/vtkmodules)
   set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")     #
   set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python-devel, numpy, swig")              #
   set (CoreFlows_SWIG_DIR ${CoreFlows_SOURCE_DIR}/swig)                                                    #
endif ()                                                                                                    #
                                                                                                            #
#-----------------------------------------------------------------------------------------------------------#

# Base directories
set (CoreFlows_SRC ${CoreFlows_SOURCE_DIR}/Models ) 
set (CoreFlows_EXAMPLES ${CoreFlows_SOURCE_DIR}/examples)

set( CoreFlows_INCLUDES 
  ${CDMATH_INCLUDES}											    #
  ${MED_INCLUDES}											    #
  ${MEDCOUPLING_INCLUDES}										    #
  ${PETSC_INCLUDES}											    #
  ${CoreFlows_SRC}/inc											    #    
  )													    #

add_subdirectory (${CoreFlows_SRC})
add_subdirectory (${CoreFlows_EXAMPLES})
if (COREFLOWS_WITH_PYTHON)                                                                                  #
   add_subdirectory (${CoreFlows_SWIG_DIR})                                                                 #
endif ()                                                                                                    #

# Documentation --------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_DOCUMENTATION)                                                                           #
  find_package (Doxygen)                                                                                    #
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, doxygen, graphviz, mscgen")           #
  set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, doxygen, graphviz, mscgen")               #
  if (NOT DOXYGEN_FOUND)                                                                                    #
    message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")      #
  endif (NOT DOXYGEN_FOUND)                                                                                 #
  configure_file (Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)                              #
  FILE(GLOB_RECURSE CoreFlows_Documentation ${CoreFlows_SOURCE_DIR}/Documentation *.dox)                    #
  SET(Documentation_SRC                                                                                     #
    ${PROJECT_BINARY_DIR}/Doxyfile                                                                          #
    ${CoreFlows_Documentation}                                                                              #
      )                                                                                                     #
  add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE}                                                  #
  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)			                                            #
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/CoreFlows OPTIONAL)                  #
  INSTALL(FILES ${CoreFlows_SOURCE_DIR}/Documentation/CoreFlows.pdf DESTINATION share/doc/CoreFlows OPTIONAL)	    #
endif ()                                                                                                    #
                                                                                                            #
# ------------- --------------------------------------------------------------------------------------------#

# Graphic interface-----------------------------------------------------------------------------------------#
     													    #
if (COREFLOWS_WITH_GUI)											    					  														    #
 set (CoreFlows_GUI_DIR ${CoreFlows_SOURCE_DIR}/gui)                                                        #
 add_subdirectory(${CoreFlows_GUI_DIR})			

endif ()                                                                                                    #
													    #
# ------------- --------------------------------------------------------------------------------------------#

#--------------------- COMPILATION MAIN --------------------------------------------------------------------#
													    #
INCLUDE_DIRECTORIES(						   					    #
  ${CoreFlows_INCLUDES}											    #
  )													    #
													    #
SET(CoreFlowsMain_SOURCES										    #
    ${CoreFlows_SRC}/src/Fluide.cxx									    #
    ${CoreFlows_SRC}/src/DiffusionEquation.cxx								    #
    ${CoreFlows_SRC}/src/StationaryDiffusionEquation.cx							    #
    ${CoreFlows_SRC}/src/ProblemFluid.cxx								    #
    ${CoreFlows_SRC}/src/IsothermalTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/utilitaire_algebre.cxx								    #
    ${CoreFlows_SRC}/src/DriftModel.cxx									    #
    ${CoreFlows_SRC}/src/SinglePhase.cxx								    #
    ${CoreFlows_SRC}/src/ProblemCoreFlows.cxx								    #
    ${CoreFlows_SRC}/src/FiveEqsTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/TransportEquation.cxx								    #
  )													    #
													    #
SET(MAIN_SOURCES											    #
  main.cxx												    #
  )													    #
SET(_extra_lib_CoreFlows ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} )					    #
ADD_EXECUTABLE(CoreFlowsMainExe ${MAIN_SOURCES})							    #
SET_TARGET_PROPERTIES(CoreFlowsMainExe PROPERTIES COMPILE_FLAGS "")					    #
TARGET_LINK_LIBRARIES(CoreFlowsMainExe CoreFlows  ${_extra_lib_CoreFlows})				    #
INSTALL(TARGETS CoreFlowsMainExe DESTINATION bin/Executable)						    #
													    #
configure_file(
    ${PROJECT_SOURCE_DIR}/env_CoreFlows.sh
    ${CMAKE_INSTALL_PREFIX}/env_CoreFlows.sh
    @ONLY
)
													    
													    
# ------------- --------------------------------------------------------------------------------------------#


# Packaging
if (COREFLOWS_WITH_PACKAGE)
  include (CPackLists.txt)
endif ()
