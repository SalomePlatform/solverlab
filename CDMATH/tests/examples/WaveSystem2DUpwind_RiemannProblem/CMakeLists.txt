
SET(MESH_MED
  ../../ressources/squareWithTriangles.med
  ../../ressources/squareWithSquares.med
  ../../ressources/squareWithBrickWall.med
  ../../ressources/squareWithCheckerboardSquares.med
  ../../ressources/squareWithDeformedQuadrangles.med
  ../../ressources/squareWithHexagons.med
  )

file(COPY ${MESH_MED} DESTINATION ${CMAKE_BINARY_DIR}/tests/examples/WaveSystemUpwind)
install(FILES ${MESH_MED} DESTINATION share/examples/WaveSystemUpwind)

if (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    SET(IMPLICIT_SCHEME  0 )

    SET(MESH_FILE  ../../ressources/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_triangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_squares PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithBrickWall.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_brickwall PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_brickwall APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithCheckerboardSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_checkerboard PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_checkerboard APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithDeformedQuadrangles.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_deformedQuadrangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_deformedQuadrangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_deformedQuadrangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_hexagons PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindExplicit_SQUARE_hexagons APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(IMPLICIT_SCHEME  1 )

    SET(MESH_FILE  ../../ressources/meshSquare.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_triangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_triangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_triangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_squares ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_squares PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_squares APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithBrickWall.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_brickwall ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_brickwall PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_brickwall APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithCheckerboardSquares.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_checkerboard ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_checkerboard PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_checkerboard APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithDeformedQuadrangles.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_deformedQuadrangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_deformedQuadrangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_deformedQuadrangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

    SET(MESH_FILE  ../../ressources/squareWithHexagons.med  )

    ADD_TEST(ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_hexagons ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/WaveSystemUpwind.py  ${MESH_FILE} ${IMPLICIT_SCHEME})
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_hexagons PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")
    SET_PROPERTY(TEST ExampleWaveSystem_2DRiemannProblem_UpwindImplicit_SQUARE_hexagons APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CDMATH_BINARY_DIR}/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath")

endif (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


