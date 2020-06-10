# nom du projet tests
PROJECT(tests)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (CDMATH_WITH_TESTS)
  ADD_SUBDIRECTORY(cdmath)

  INCLUDE_DIRECTORIES(
    ${MEDCOUPLING_INCLUDE_DIR}
    ${BASE_DIR}/inc
    ${MESH_DIR}/inc
    cdmath
    )

  SET(_extra_lib_TESTS)
  if (CDMATH_WITH_PETSC)
      INCLUDE_DIRECTORIES(
        ${LINEARSOLVER_DIR}/inc
        ${PETSC_INCLUDES}
        ${SLEPC_INCLUDES}
        )
      SET(_extra_lib_TESTS linearsolver ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
  endif (CDMATH_WITH_PETSC)

  SET(MAIN_SOURCES
    MainTests.cxx
    )

  ADD_EXECUTABLE(TestCDMATH ${MAIN_SOURCES})

  if (CDMATH_WITH_PETSC)
     SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS -DCDMATH_WITH_PETSC)
  else(CDMATH_WITH_PETSC)
     SET_TARGET_PROPERTIES(TestCDMATH PROPERTIES COMPILE_FLAGS "")
  endif (CDMATH_WITH_PETSC)

  TARGET_LINK_LIBRARIES(TestCDMATH TestMesh base mesh ${CPPUNIT_LIBS} ${_extra_lib_TESTS})

  ADD_TEST(UnitTestCDMATH TestCDMATH)
  SET_PROPERTY(TEST UnitTestCDMATH PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib")

  INSTALL(TARGETS TestCDMATH DESTINATION bin/tests)
ENDIF (CDMATH_WITH_TESTS)

ADD_SUBDIRECTORY(ressources)

IF (CDMATH_WITH_PYTHON)
  ADD_SUBDIRECTORY(swig)
ENDIF (CDMATH_WITH_PYTHON)

FILE(COPY doc/ DESTINATION convergence)
INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/convergence DESTINATION share/doc/)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)
  # Examples tests
  set (examples_dir ${CDMATH_SOURCE_DIR}/tests/examples)
  add_subdirectory (${examples_dir})
  add_subdirectory (validation)

  add_custom_target (tests COMMAND ctest -O tests.log)
  add_custom_target (check COMMAND ctest -E 'validation|Example')# may be replace ctest -E with ctest -LE
  add_custom_target (examples COMMAND ctest -R Example)# may be replace ctest -R with ctest -L
  add_custom_target (validation COMMAND ctest -R validation)# may be replace ctest -R with ctest -L

  add_custom_target (poisson COMMAND ctest -R Poisson)# may be replace ctest -R with ctest -L
  add_custom_target (diffusion COMMAND ctest -R Diffusion)# may be replace ctest -R with ctest -L
  add_custom_target (transport COMMAND ctest -R Transport)# may be replace ctest -R with ctest -L
  add_custom_target (wave COMMAND ctest -R WaveSystem)# may be replace ctest -R with ctest -L
  add_custom_target (upwind COMMAND ctest -R Upwind)# may be replace ctest -R with ctest -L
  add_custom_target (pstag COMMAND ctest -R PStag)# may be replace ctest -R with ctest -L
  add_custom_target (staggered COMMAND ctest -R Staggered)# may be replace ctest -R with ctest -L
  add_custom_target (centered COMMAND ctest -R Centered)# may be replace ctest -R with ctest -L
  add_custom_target (triangles COMMAND ctest -R triangles)# may be replace ctest -R with ctest -L
  add_custom_target (checkerboard COMMAND ctest -R checkerboard)# may be replace ctest -R with ctest -L
  add_custom_target (squares COMMAND ctest -R squares)# may be replace ctest -R with ctest -L
  add_custom_target (scaling COMMAND ctest -R '_scaling' -E no_scaling)# may be replace ctest -R with ctest -L
  add_custom_target (noscaling COMMAND ctest -R no_scaling)# may be replace ctest -R with ctest -L
  add_custom_target (2D COMMAND ctest -R 2D)# may be replace ctest -R with ctest -L
  add_custom_target (3D COMMAND ctest -R 3D)# may be replace ctest -R with ctest -L
  add_custom_target (shock COMMAND ctest -R Shock)# may be replace ctest -R with ctest -L
  add_custom_target (vortex COMMAND ctest -R Vortex)# may be replace ctest -R with ctest -L
  add_custom_target (riemann COMMAND ctest -R RiemannProblem)# may be replace ctest -R with ctest -L
endif ()

