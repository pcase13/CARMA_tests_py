# CARMA_tests_py
Python plotting scripts for the CARMA tests cases, netcdf4 version

## Use:
This project includes example output files for each of the implemented tests. To use on one's own version of the output, copy the appropriate read_[test].py into your /run/carma directory and run:

`python read_[test].py`

This will show the figure as well as save it to a png under the name read_[test].png. 

Note: [test] should be replaced with the name of the test in all above examples.

## Required python libraries:
[netCDF4](http://unidata.github.io/netcdf4-python/)

[hdf5](https://www.h5py.org/)

[numpy](http://www.numpy.org/)

[matplotlib](https://matplotlib.org/)
