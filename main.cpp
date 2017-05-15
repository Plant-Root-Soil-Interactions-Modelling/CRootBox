// Copyright (C) 2016 Daniel Leitner and Andrea Schnepf. See //license.txt for details.

#include "RootSystem.h"
#include "analysis.h"

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "examples/example1.h"
#include "examples/example1_wb_dgf.h"
#include "examples/example2.h"
#include "examples/example3.h"
#include "examples/example4.h"
#include "examples/example5.h"

#include "examples/benchmarks.h"

#include "examples/shehan_SoilCore.h"
#include "examples/shehan_RhizoTubes.h"
#include "examples/shehan_Trenches.h"

#include "examples/Exudation/example_exudation.h"

//#include "DumuxRootSystem.h" // suggested coupling
//#include "examples/example_dumux.h"

/**
 * Starts an examples (from the examples folder)
 */
int main(int argc, char* argv[])
{
    string name="";

    if (argc>1) {
        name= argv[1];
    }

    // example1(); // open parameter file, and output VTP
    // example1_wb_dgf(); // root growth inside a big box to simulate soil surface, open parameter file, and output VTP
    // example2(); // like example 1, but with put geometry
    // example3(); // more than 1 plant
    // example4(); // rhizotubes (an example for a more complex geomety)
    // example5(); // hydrotropism

    // benchmarks();

    if (argc>1) {
        cout<<"starting simulation: "<< name <<"\n";
        shehan_SoilCore(name, false); // put true here to export geometry
    } else {
        shehan_SoilCore(); // with default values
    }

    // shehan_SoilCore("wheat",true);

    // shehan_RhizoTubes("wheat",true);

    // shehan_Trenches("wheat",true);

    // example_dumux(); // tests the suggested dumux coupling

    // example_exudation();

    return(0);

}



