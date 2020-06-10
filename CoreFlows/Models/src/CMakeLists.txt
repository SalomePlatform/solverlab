
INCLUDE_DIRECTORIES(
  ${PETSC_INCLUDES} 
  ${CDMATH_INCLUDES} 
  ${CDMATH_INCLUDES}/med											    #
  ${CDMATH_INCLUDES}/medcoupling											    #
  ${CoreFlows_SRC}/inc      
  )

SET(src_models_CXX
    DriftModel.cxx
    IsothermalTwoFluid.cxx
    SinglePhase.cxx
    FiveEqsTwoFluid.cxx
    ProblemCoreFlows.cxx
    TransportEquation.cxx
    DiffusionEquation.cxx
    StationaryDiffusionEquation.cxx
    Fluide.cxx
    ProblemFluid.cxx
    utilitaire_algebre.cxx
  )

ADD_LIBRARY(CoreFlows SHARED ${src_models_CXX})
target_link_libraries(CoreFlows ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES})

INSTALL(TARGETS CoreFlows DESTINATION lib)
