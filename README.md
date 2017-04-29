# CRootBox

The fastest way to try CRootBox is to read the examples. Just uncomment the example in the main.cpp file to try, compile and run it. 

The code should compile with any c++11 compiler, e.g. for g++:
$ g++ *.cpp -std=c++11
$ ./a.out


# Folder sructure

/			CRootBox C++ codes
/examples 		Some examples how to use the CRootBox
/modelparameter		Some root parameter, and a plant parameter files
/scripts 		Pyhthon scripts for visualization with Paraview, and Matlab scripts for parameter export
/results 		Nice result images


# Documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation should now be located in the folder /doc

To build the shared library py_rootbox for coupling with Python pleaser refer to 'python building guide.txt'

