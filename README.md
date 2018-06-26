# CRootBox

The fastest way to try CRootBox is to read the examples. Just uncomment the example in the main.cpp file to try, compile and run it. 

The code should compile with any c++11 compiler, e.g. for g++:

    g++ *.cpp -std=c++11
    ./a.out


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

# References
Andrea Schnepf, Daniel Leitner, Magdalena Landl, Guillaume Lobet, Trung Hieu Mai, Shehan Morandage, Cheng Sheng, Mirjam Zörner, Jan Vanderborght, Harry Vereecken; CRootBox: a structural–functional modelling framework for root systems, Annals of Botany, Volume 121, Issue 5, 18 April 2018, Pages 1033–1053, https://doi.org/10.1093/aob/mcx221
