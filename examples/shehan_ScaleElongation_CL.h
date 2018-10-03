#include <sstream>
#include <fstream>

/**
 * Example for Shehan
 *
 * Scale elongation rate according to some 1d tabular data
 * including carbon limitation
 *
 * see shehan_ScaleElongation.h for definition of class ScaleElongation and function readCSV
 *
 * proof of concept, please adjust
 */

namespace CRootBox {

/**
 * The maximal possible total root system length increment per day for time t, and duration dt
 * (based on carbon, or leaf area)
 */
double length_increment(double t, double dt) {
    return 5. / dt; // 1 cm per day
}

/**
 * Example how to implement varying elongation rates
 */
void shehan_ScaleElongation_CL()
{
    using namespace std;

    RootSystem rootsystem;

    string name = "Zea_mays_5_Leitner_2014";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);

    /*
     * Create the tables for water content and temperature and the scale elongation function
     */
    // grid coordinates
    int n = 7;
    std::vector<double> z_{0., 15., 25., 50., 70., 100., 140.};

    // temperature data
    std::cout << "reading temperature data ...";
    auto field_temp = readCSV("/home/daniel/workspace/CRootBox/examples/wheat.csv",';',2,1);
    std::vector<double> temp = field_temp.at(0);

    // water content data (located between the grid coordinates, i.e.l n-1 data points), where do the data come from? read from file? or copy paste into the code?
    std::vector<double> wc(n-1);
    for (size_t i = 0; i<wc.size(); i++) {
        wc[i] = 0.5;
    }

    // create scale elongation function
    Grid1D water_content = Grid1D(n, z_, wc);
    Grid1D temperature = Grid1D(n, z_, temp );
    ScaleElongation se = ScaleElongation(&water_content, &temperature);

    // create carbon scaling (on top)
    ProportionalElongation pe = ProportionalElongation();
    pe.setBaseLookUp(&se);

    // "manually" set the scale elongation function
    for (int i = 1; i < 7; i++) {
        rootsystem.getRootTypeParameter(i)->se = &pe;
    }

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 7*30; // days
    double dt = 0.5 * 1./24.;  // dt is half an hour
    size_t N = round(simtime/dt);

    assert(N<=field_temp.size()); // check if enough data are available

    for (size_t i=0; i<N; i++) {

        rootsystem.simulate(dt, length_increment(i * dt, dt), &pe, false);

        // update field data:
        temperature.data =  field_temp.at(i);
        // water_content.data = ... (type is vector<double>)

    }

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    auto rl = rootsystem.getScalar(RootSystem::st_length);
    double tl = std::accumulate(rl.begin(), rl.end(), 0);
    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes, " << tl << " cm total length \n";
}

} // end namespace CRootBox
