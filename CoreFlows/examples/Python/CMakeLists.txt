# Graphic interface---------------------------------------------------------------------------------------------#
														#
if (COREFLOWS_WITH_GUI)													    			 
						    								#
	# Find KERNEL 
	# ==============
	SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
	IF( EXISTS ${KERNEL_ROOT_DIR} )
	  LIST(APPEND CMAKE_MODULE_PATH "${KERNEL_ROOT_DIR}/salome_adm/cmake_files")
	  LIST(APPEND CMAKE_PREFIX_PATH "${KERNEL_ROOT_DIR}/salome_adm/cmake_files")
	  INCLUDE(SalomeMacros)
	  FIND_PACKAGE(SalomeKERNEL REQUIRED)
	ELSE( EXISTS ${KERNEL_ROOT_DIR} )
	  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR!")
	ENDIF( EXISTS ${KERNEL_ROOT_DIR} )
	
	# Find SALOME GUI 
	# ==============
	SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR} CACHE PATH "Path to the Salome GUI")
	IF(EXISTS ${GUI_ROOT_DIR})
	  LIST(APPEND CMAKE_MODULE_PATH "${GUI_ROOT_DIR}/adm_local/cmake_files")
	  LIST(APPEND CMAKE_PREFIX_PATH "${GUI_ROOT_DIR}/adm_local/cmake_files")
	  FIND_PACKAGE(SalomeGUI)
	ELSE(EXISTS ${GUI_ROOT_DIR})
	  MESSAGE(FATAL_ERROR "We absolutely need a Salome GUI, please define GUI_ROOT_DIR!")
	ENDIF(EXISTS ${GUI_ROOT_DIR})
	
	# Platform setup
	# ==============
	INCLUDE(SalomeSetupPlatform)   # From KERNEL
	    
	# Prerequisites
	# =============
	
	# Mandatory products
	IF(  $ENV{PARAVIEW_VERSION} STRLESS "5")
	         FIND_PACKAGE(SalomePython REQUIRED)
	ELSE()
	         FIND_PACKAGE(SalomePythonInterp REQUIRED)
	         FIND_PACKAGE(SalomePythonLibs REQUIRED)
	ENDIF()
	
	SET(_all_SCRIPTS
	    main_tests.py
	    )
	
	SALOME_INSTALL_SCRIPTS("${_all_SCRIPTS}" ${SALOME_INSTALL_SCRIPT_PYTHON}/tests)
	
	SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PROJECT_BINARY_DIR}/src) 
	SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
	
	ADD_TEST(PythonTests ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/main_tests.py)
	SET_TESTS_PROPERTIES(PythonTests PROPERTIES ENVIRONMENT "${tests_env}")
						    								#
endif (COREFLOWS_WITH_GUI)                                                                                                        #
#---------------------------------------------------------------------------------------------------------------#



############# python example tests

function(CreatePythonTest pythonFile)
    add_test(${pythonFile} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${pythonFile})
    SET_PROPERTY(TEST ${pythonFile} PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling:${PETSC_DIR}/${PETSC_ARCH}/lib)
    SET_PROPERTY(TEST ${pythonFile} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/CoreFlows_Python:${CMAKE_INSTALL_PREFIX}/bin/CoreFlows_Python:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/cdmath:${CDMATH_DIR}/bin/cdmath:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling)
     message("setting python example test ${pythonFile}")
endfunction(CreatePythonTest)

# copy all *.py tests files before build
file(GLOB pythonTestFiles "*.py" )
file(COPY ${pythonTestFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
)


CreatePythonTest(DiffusionEquation_1DHeatedRod.py)
CreatePythonTest(DriftModel_1DBoilingAssembly.py)
CreatePythonTest(DriftModel_1DBoilingChannel.py)
CreatePythonTest(DriftModel_1DDepressurisation.py)
CreatePythonTest(DriftModel_1DPorosityJump.py)
CreatePythonTest(DriftModel_1DPressureLoss.py)
CreatePythonTest(DriftModel_1DRiemannProblem.py)
CreatePythonTest(DriftModel_1DVidangeReservoir.py)
CreatePythonTest(DriftModel_2BranchesBoilingChannels.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravity.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravityTriangles.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravityBarriers.py)
CreatePythonTest(DriftModel_2DBoilingChannelBarrier.py)
CreatePythonTest(DriftModel_2DInclinedBoilingChannelBarrier.py)
CreatePythonTest(DriftModel_2DInclinedBoilingChannel.py)
CreatePythonTest(DriftModel_2DPorosityJump.py)
CreatePythonTest(DriftModel_2DPressureLoss.py)
CreatePythonTest(DriftModel_2DVidangeReservoir.py)
CreatePythonTest(DriftModel_2DVidangeReservoirUnstructured.py)
CreatePythonTest(DriftModel_3DBoilingChannelBarrier.py)
CreatePythonTest(FiveEqsTwoFluid_1DBoilingAssembly.py)
CreatePythonTest(FiveEqsTwoFluid_1DBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(FiveEqsTwoFluid_2DInclinedBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid_2DInclinedSedimentation.py)
CreatePythonTest(FiveEqsTwoFluid_2DVidangeReservoir.py)
CreatePythonTest(IsothermalTwoFluid_1DSedimentation.py)
CreatePythonTest(IsothermalTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(IsothermalTwoFluid_2DVidangeReservoir.py)
CreatePythonTest(SinglePhase_1DDepressurisation.py)
CreatePythonTest(SinglePhase_1DHeatedAssembly.py)
CreatePythonTest(SinglePhase_1DHeatedChannel.py)
CreatePythonTest(SinglePhase_1DRiemannProblem.py)
CreatePythonTest(SinglePhase_1DWaterHammer.py)
CreatePythonTest(SinglePhase_2BranchesHeatedChannels.py)
CreatePythonTest(SinglePhase_2DHeatedChannelInclined.py)
CreatePythonTest(SinglePhase_2DLidDrivenCavity.py)
CreatePythonTest(SinglePhase_2DLidDrivenCavity_unstructured.py)
CreatePythonTest(SinglePhase_2DSphericalExplosion_unstructured.py)
CreatePythonTest(SinglePhase_3DSphericalExplosion_unstructured.py)
CreatePythonTest(SinglePhase_2DThermalDiffusion.py)
CreatePythonTest(SinglePhase_2DVidangeReservoir.py)
CreatePythonTest(SinglePhase_2DWallHeatedChannel_ChangeSect.py)
CreatePythonTest(SinglePhase_3DHeatDrivenCavity.py)
CreatePythonTest(SinglePhase_3DVortexTube_NoCone_NoViscosity.py)
CreatePythonTest(SinglePhase_3DVortexTube_WithCone_NoViscosity.py)
CreatePythonTest(TransportEquation_1DHeatedChannel.py)
CreatePythonTest(StationaryDiffusionEquation_2DEF.py)
CreatePythonTest(StationaryDiffusionEquation_2DEF_Neumann.py)
CreatePythonTest(StationaryDiffusionEquation_2DFV_StructuredTriangles.py)
CreatePythonTest(StationaryDiffusionEquation_2DFV_EquilateralTriangles.py)
CreatePythonTest(StationaryDiffusionEquation_2DFV_StructuredSquares.py)
CreatePythonTest(StationaryDiffusionEquation_2DFV_StructuredSquares_Neumann.py)
CreatePythonTest(StationaryDiffusionEquation_3DEF.py)
CreatePythonTest(StationaryDiffusionEquation_3DFV_StructuredCubes.py)
CreatePythonTest(StationaryDiffusionEquation_3DFV_StructuredTetrahedra.py)
CreatePythonTest(StationaryDiffusionEquation_3DEF_RoomCooling.py)
CreatePythonTest(StationaryDiffusionEquation_3DVF_RoomCooling_StructuredCubes.py)
CreatePythonTest(StationaryDiffusionEquation_3DVF_RoomCooling_UnstructuredTetras.py)

############# python convergence tests

function(CreatePythonConvergenceTest model convergenceTest pythonSolver MESHES_AND_PICTURES_PATH)
    FILE(GLOB MY_MESHES_AND_PICTURES  MESHES_AND_PICTURES_PATH )
    file(COPY convergence_${model}_${convergenceTest}.py ../${pythonSolver} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE )

    add_test(convergence_${model}_${convergenceTest} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/convergence_${model}_${convergenceTest}.py)
    SET_PROPERTY(TEST convergence_${model}_${convergenceTest} PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling:${PETSC_DIR}/${PETSC_ARCH}/lib:${PV_LIB_DIR})
    SET_PROPERTY(TEST convergence_${model}_${convergenceTest} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/CoreFlows_Python:${CMAKE_INSTALL_PREFIX}/bin/CoreFlows_Python:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/cdmath:${CDMATH_DIR}/bin/cdmath:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling:${CDMATH_DIR}/bin/cdmath/postprocessing:${PV_PYTHON_DIR})
     message("setting python convergence test convergence_${model}_${convergenceTest}")
endfunction(CreatePythonConvergenceTest)

add_subdirectory (Convergence/StationaryDiffusion)

