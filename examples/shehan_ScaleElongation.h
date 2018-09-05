/**
 * Example for Shehan
 *
 * Scale elongation rate according to some 1d tabular data
 *
 * proof of concept, please adjust
 */

namespace CRootBox {

/**
 * The elongation rate is calculated in dependence of water content and temperature,
 * put the empirical relationship into the private function f in line 43.
 *
 */
class ScaleElongation :public SoilLookUp
{

public:

    ScaleElongation(Grid1D* wc, Grid1D* temp)
    {
        water_content = wc;
        temperature = temp;
    }

    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const ///< Returns a scalar property of the soil, 1. per default
    {
        double wc = water_content->getValue(pos,root);
        double temp = temperature->getValue(pos,root);
        return f(wc,temp);
    }

    virtual std::string toString() const { return "ScaleElongation based on water content and temperature"; } ///< Quick info about the object for debugging

    Grid1D* water_content;
    Grid1D* temperature;

private:

    // this function describes the dependency the elongation scale to
    double f(double wc, double temp) const
    {
        return 0.5*(temp-minT)/(maxT-minT)+0.5; // <- TODO put the magic function here
    }

    double minT = -10;
    double maxT = 40;

};

/**
 * Example how to implement varying elongation rates
 */
void shehan_ScaleElongation()
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
    int n = 100;
    double topZ = 0;
    double botZ = -100; // cm
    std::vector<double> z_(n);
    for (size_t i = 0; i<z_.size(); i++) {
        z_[i] = botZ + i*(topZ-botZ)/(n-1);
    }
    // temperature data (located between the grid coordinates, i.e.l n-1 data points), where do the data come from? read from file? or copy paste into the code?
    std::vector<double> temperature(n-1);
    for (size_t i = 0; i<temperature.size(); i++) {
        temperature[i] = 0.5;
    }
    // water content data (located between the grid coordinates, i.e.l n-1 data points), where do the data come from? read from file? or copy paste into the code?
    std::vector<double> water_content(n-1);
    for (size_t i = 0; i<water_content.size(); i++) {
        water_content[i] = 0.5;
    }
    // create scale elongation function
    Grid1D wc = Grid1D(n, z_, water_content);
    Grid1D temp = Grid1D(n, z_, temperature);
    ScaleElongation se = ScaleElongation(&wc, &temp);
    // "manually" set the scale elongation function
    for (int i=1; i<7; i++) {
        rootsystem.getRootTypeParameter(i)->se = &se;
    }

    /**
     * Initialize
     */
    rootsystem.initialize(); //it is important to call initialize() after setGeometry()

    /*
     * Simulate
     */
    double simtime = 30; // e.g. 30 or 60 days
    double dt = simtime;
    int N = round(simtime/dt);
    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";
}

} // end namespace CRootBox
