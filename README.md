# LofarStMan

This package provides a casacore storage manager to read raw correlator output of the LOFAR telescope.

To use it, install it and make sure that the shared library `liblofarstman.so` is put in the linker path, e.g. by setting
```
export LD_LIBRARY_PATH=/path/to/liblofarstman.so:$LD_LIBRARY_PATH
```
If the library is installed in a system location, this is not necessary.
When the shared library is found, programs like `taql` and [DP3](https://git.astron.nl/RD/DP3/) can read the MeasurementSet for further processing.
