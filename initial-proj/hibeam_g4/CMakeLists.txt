cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(hibeam_g4)

include(CheckCXXCompilerFlag)
set(HAVE_ISYSTEM true)
set (CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS Geom)
if(ROOT_FOUND)
  list(APPEND EXTRA_LIBRARIES ROOT::Geom)
  message(STATUS "ROOT found.")
endif()

SET(VGM_INSTALL $ENV{VGM_INSTALL})
find_package(VGM REQUIRED)
if (VGM_FOUND)
  message(STATUS "VGM :" ${VGM_INCLUDE_DIRS})
  include_directories(${VGM_INCLUDE_DIRS})
  list(APPEND EXTRA_LIBRARIES ${VGM_LIBRARIES} )
endif()
  
find_package(Geant4 11.0.0 REQUIRED)
if(Geant4_FOUND)
  message(STATUS "Geant4 :" ${Geant4_INCLUDE_DIRS})
endif()

message(STATIS "libraries :" ${EXTRA_LIBRARIES})
# Locate sources and headers for this project
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

# Add the executable, and link it to the Geant4 libraries
set(library_name _${PROJECT_NAME})

add_executable(hibeam_g4 hibeam_g4.cc ${sources} ${headers})
target_link_libraries(hibeam_g4 ${library_name} ${Geant4_LIBRARIES} ${EXTRA_LIBRARIES})

set(SourceLib ${PROJECT_SOURCE_DIR}/src/PrimaryEvent.cc ${PROJECT_SOURCE_DIR}/src/MultiHits.cc)
add_library(${library_name} SHARED ${SourceLib} ${headers})
target_link_libraries(${library_name})

# Copy config files and macros to the build directory.
set(HG4_SCRIPTS
  example.config
  icons.mac
  gui.mac
  init.mac 
  init_vis.mac 
  run_gps.mac 
  run_mcpl.mac 
  run_scattering.mac 
  set_colors.mac 
  vis.mac
  )

foreach(_script ${HG4_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macros/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
