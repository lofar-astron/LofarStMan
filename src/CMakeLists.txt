# $Id: CMakeLists.txt 14273 2009-10-16 10:08:29Z loose $

include(LofarPackageVersion)

lofar_add_library(lofarstman
  Package__Version.cc
  LofarStMan.cc
  LofarColumn.cc
  Register.cc
  )

lofar_add_bin_program(versionlofarstman versionlofarstman.cc)
lofar_add_bin_program(fixlsmeta fixlsmeta.cc)

lofar_add_bin_scripts(makeFLAGwritable)
