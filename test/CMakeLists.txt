# $Id: CMakeLists.txt 14273 2009-10-16 10:08:29Z loose $

include(LofarCTest)

lofar_add_test(tLofarStMan tLofarStMan.cc)
lofar_add_test(tfix tfix.cc)
lofar_add_test(tIOPerf tIOPerf.cc)
