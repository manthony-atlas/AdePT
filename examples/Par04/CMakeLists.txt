# SPDX-FileCopyrightText: 2021 CERN
# SPDX-License-Identifier: Apache-2.0

if(NOT TARGET G4HepEm::g4HepEm)
  message(STATUS "Disabling TestEm3 (needs G4HepEm)")
  return()
endif()

file(GLOB sources "${PROJECT_SOURCE_DIR}/examples/Par04/src/*.cc")

# Par04
add_executable(examplePar04 examplePar04.cpp AdeptIntegration.cpp AdeptIntegration.cu electrons.cu gammas.cu relocation.cu ${sources})

message(STATUS "srcs: ${sources}")

target_include_directories(examplePar04 PRIVATE ${PROJECT_SOURCE_DIR}/examples/Par04/include ${PROJECT_SOURCE_DIR}/examples/Par04)
target_link_libraries(examplePar04 PRIVATE AdePT CopCore::CopCore VecGeom::vecgeom VecGeom::vecgeomcuda_static VecGeom::vgdml ${Geant4_LIBRARIES} G4HepEm::g4HepEmData G4HepEm::g4HepEmInit G4HepEm::g4HepEmRun)
target_compile_options(examplePar04 PRIVATE "$<$<AND:$<OR:$<CONFIG:Debug>>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>")

set_target_properties(examplePar04 PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
