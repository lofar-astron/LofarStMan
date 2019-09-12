# $Id: CMakeLists.txt 14273 2009-10-16 10:08:29Z loose $

include_directories($<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

add_executable(tfix tfix.cc)
target_compile_features(tfix PRIVATE cxx_std_11)
target_link_libraries(tfix lofarstman)
add_executable(tLofarStMan tLofarStMan.cc)
target_compile_features(tLofarStMan PRIVATE cxx_std_11)
target_link_libraries(tLofarStMan lofarstman)
add_executable(tIOPerf tIOPerf.cc)
target_compile_features(tIOPerf PRIVATE cxx_std_11)
target_link_libraries(tIOPerf lofarstman)

add_test(tfix tfix)
add_test(tLofarStMan tLofarStMan)
add_test(tIOPerf tIOPerf)
