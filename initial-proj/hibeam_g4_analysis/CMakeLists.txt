cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(hibeam_g4_analysis)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS Core Tree)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Create a shared library.
add_library(hibeam_g4 SHARED src/hibeam_g4.cxx src/hibeam_g4_source.cxx)
target_include_directories(hibeam_g4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hibeam_g4 PUBLIC ROOT::Core ROOT::Tree)

add_library(hibeam_dict SHARED src/SECHit.cxx src/TPCHit.cxx src/CVHit.cxx src/ScintHit.cxx src/ScintAna.cxx src/SimpleSciHit.cxx src/ProtoTPCHit.cxx)
target_include_directories(hibeam_dict PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(hibeam_dict PUBLIC ROOT::Core ROOT::Tree)

# Generate the dictionary and add to an existing target:
root_generate_dictionary(hibeam_d include/TPCHit.hh include/SECHit.hh include/CVHit.hh include/ScintHit.hh include/ScintAna.hh include/SimpleSciHit.hh include/ProtoTPCHit.hh MODULE hibeam_dict LINKDEF include/Linkdef.h)

# Create the main program using the library.
add_executable(hibeam_ana hibeam_ana.cxx)
target_link_libraries(hibeam_ana PUBLIC hibeam_g4 hibeam_dict ROOT::Core ROOT::Tree)


