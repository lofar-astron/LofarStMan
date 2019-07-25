# $Id: CMakeLists.txt 14280 2009-10-16 13:48:46Z loose $

lofar_package(LofarStMan 1.0 DEPENDS Common)

include(LofarFindPackage)
lofar_find_package(Casacore COMPONENTS casa measures tables REQUIRED)

add_subdirectory(include/LofarStMan)
add_subdirectory(src)
add_subdirectory(test)
