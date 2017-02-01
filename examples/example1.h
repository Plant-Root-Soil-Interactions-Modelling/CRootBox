/**
 * Example 1
 *
 * 1) Opens plant and root parameters from a file
 * 2) Simulates root growth
 * 3) Outputs a VTP (for vizualisation in ParaView)
 *    In Paraview: use tubePLot.py script for fancy visualisation (Macro/Add new macro...), apply after opening file
 *
 *  Additionally, exports the line segments as .txt file to import into Matlab for postprocessing
 */
using namespace std;

void example1()
{
    RootSystem rootsystem;

    string name = "Maize_Pheno1_Leitner_et_al_2014";

    /*
     * Open plant and root parameter from a file
     */
    rootsystem.openFile(name);
    rootsystem.writeParameters(std::cout);

    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    double simtime = 60; // 20, 40, 60 days
    double dt = 60; // try other values here
    int N = round(simtime/dt);

    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp",RootSystem::ot_segments); // use ot_polylines for nicer visualization, ot_segments for animations

    /*
     * Export segments in RSML format
     */
    rootsystem.write(name+".rsml");

    /*
     * Export segments for Matlab analysis
     */
    SegmentAnalyser analysis(rootsystem);
    analysis.write(name+".txt");

    /*
      Total length and surface
     */
    double l = analysis.getSummed(RootSystem::st_length);
    std::cout << "Root system length " << l << " cm \n";

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";

}
