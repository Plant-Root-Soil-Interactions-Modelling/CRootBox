# CRootBox

The fastest way to try CRootBox is to read the examples. Just uncomment the desired example in the `examples/main.cpp` file and compile and run it.

```bash
cmake .
make
cd examples && ./test_crootbox
```

`cmake . ` runs CMake which configures the CRootBox libraries. `make` builds the libraries and the C++ example. `cd examples && ./test_crootbox` runs the example.

# Folder sructure

/src			CRootBox C++ codes
/examples 		Some examples how to use the CRootBox
/modelparameter		Some root parameter, and a plant parameter files
/scripts 		Pyhthon scripts for visualization with Paraview, and Matlab scripts for parameter export
/results 		Nice result images


# Documentation

Create the documentation by running doxygen in the folder
$ doxygen doxy_config

The documentation should now be located in the folder /doc

# Python bindings

A Python library is automatically built by CMake if Python 3 and boost-python is installed on your system.

# References

Please cite one or all of the following papers if you make use of CRootBox in your publication.

Andrea Schnepf, Daniel Leitner, Magdalena Landl, Guillaume Lobet, Trung Hieu Mai, Shehan Morandage, Cheng Sheng, Mirjam Zörner, Jan Vanderborght, Harry Vereecken; CRootBox: a structural–functional modelling framework for root systems, Annals of Botany, Volume 121, Issue 5, 18 April 2018, Pages 1033–1053, https://doi.org/10.1093/aob/mcx221

STATISTICAL CHARACTERIZATION OF THE ROOT SYSTEM ARCHITECTURE MODEL CROOTBOX
Andrea Schnepf; Katrin Huber; Magdalena Landl; Félicien Meunier; Lukas Petrich; Volker Schmidt
doi: 10.2136/vzj2017.12.0212; Date posted: May 04, 2018
