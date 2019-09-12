# $Id: CMakeLists.txt 14280 2009-10-16 13:48:46Z loose $

cmake_minimum_required(VERSION 3.5)

project(LofarStMan)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

if(PORTABLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native") 
endif(PORTABLE)

enable_testing()

find_package(Casacore COMPONENTS casa measures tables REQUIRED)
include_directories(${CASACORE_INCLUDE_DIR})

add_library(lofarstman SHARED
  src/LofarStMan.cc
  src/LofarColumn.cc
  src/Register.cc
)

target_include_directories(lofarstman
	PUBLIC
	$<INSTALL_INTERFACE:include>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(lofarstman PRIVATE cxx_std_11)

target_link_libraries(lofarstman PUBLIC ${CASACORE_LIBRARIES})

# bin_scripts(makeFLAGwritable)
add_subdirectory(test)
include(GNUInstallDirs)

install(TARGETS lofarstman
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
