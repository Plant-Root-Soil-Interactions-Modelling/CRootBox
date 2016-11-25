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

void example1_wb()
{
    RootSystem rootsystem;

    string name = "Cauliflower_Vansteenkiste_et_al_2014";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);
    rootsystem.writeParameters(std::cout);

    // THETA by hand
    rootsystem.getRootParameter(2)->theta = M_PI/4.;
    rootsystem.getRootParameter(2)->thetas = 0; // no std

    /*
     * Set geometry
     */
    //creates a square 27*27 cm containter with height 1.5 cm (used in parametrisation experiment)
    SDF_PlantBox box(900,900,900);
    rootsystem.setGeometry(&box);
    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    double simtime = 10; // 20, 40, 60 days
    double dt = 1; // try other values here
    int N = round(simtime/dt);

    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp",RootSystem::ot_segments); // use ot_polylines for nicer visualization, ot_segments for animations

    /*
     * Export segments in DGF format TODO finish and test this method
     */
    // rootsystem.write(name+".dgf");

    /*
     * Export segments for Matlab analysis
     */
    AnalysisSDF analysis(rootsystem);
    analysis.write(name+".txt");

    /*
      Total length and surface
     */
    double l = analysis.getSummed(RootSystem::st_length);
    std::cout << "Visible Length " << l << " cm \n";

    /*
     * Export visible rhizotrons segments
     */
    cout << "writing VTP... \n";
    analysis.write(name+"_croped.vtp");
    cout << "writing RootBox segments... \n";
    analysis.write(name+"_all.txt");

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";

}
