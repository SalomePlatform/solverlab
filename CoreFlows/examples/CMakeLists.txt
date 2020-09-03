project(test)

file(GLOB MESHES_TO_INSTALL resources )
install(DIRECTORY ${MESHES_TO_INSTALL} DESTINATION share/examples)

add_subdirectory (C)

if (COREFLOWS_WITH_PYTHON )

  add_subdirectory (Python)

endif (COREFLOWS_WITH_PYTHON )

add_custom_target (tests COMMAND ctest -O testsCoreFlows.log)

add_custom_target (SinglePhase COMMAND ctest -R SinglePhase)# may be replace ctest -R with ctest -L
add_custom_target (DriftModel  COMMAND ctest -R DriftModel) # may be replace ctest -R with ctest -L
add_custom_target (IsothermalTwoFluid COMMAND ctest -R IsothermalTwoFluid)# may be replace ctest -R with ctest -L
add_custom_target (   FiveEqsTwoFluid COMMAND ctest -R    FiveEqsTwoFluid)# may be replace ctest -R with ctest -L

add_custom_target (DiffusionEquation COMMAND ctest -R DiffusionEquation -E StationaryDiffusionEquation)# may be replace ctest -R with ctest -L
add_custom_target (diffusion         COMMAND ctest -R DiffusionEquation)# may be replace ctest -R with ctest -L
add_custom_target (TransportEquation COMMAND ctest -R TransportEquation)# may be replace ctest -R with ctest -L
add_custom_target (transport         COMMAND ctest -R TransportEquation)# may be replace ctest -R with ctest -L
add_custom_target (StationaryDiffusionEquation COMMAND ctest -R StationaryDiffusionEquation)# may be replace ctest -R with ctest -L

add_custom_target (convergence COMMAND ctest -R convergence)# may be replace ctest -R with ctest -L

add_custom_target (fv         COMMAND ctest -R FV)# may be replace ctest -R with ctest -L
add_custom_target (FV         COMMAND ctest -R FV)# may be replace ctest -R with ctest -L
add_custom_target (fe         COMMAND ctest -R FE)# may be replace ctest -R with ctest -L
add_custom_target (FE         COMMAND ctest -R FE)# may be replace ctest -R with ctest -L
add_custom_target (1D         COMMAND ctest -R 1D)# may be replace ctest -R with ctest -L
add_custom_target (2D         COMMAND ctest -R 2D)# may be replace ctest -R with ctest -L
add_custom_target (3D         COMMAND ctest -R 3D)# may be replace ctest -R with ctest -L
add_custom_target (Dirichlet  COMMAND ctest -R Dirichlet)# may be replace ctest -R with ctest -L
add_custom_target (Neumann    COMMAND ctest -R Neumann)# may be replace ctest -R with ctest -L

