#include <sstream>
#include <fstream>

/**
 * Example for Shehan
 *
 * Scale elongation rate according to some 1d tabular data
 *
 * proof of concept, please adjust
 */

namespace CRootBox {

/*
 *  reads a csv into a vector^2
 */
std::vector<std::vector<double>> readCSV(std::string name, char delimeter, int ignore_rows, int ignore_cols)
{
    std::ifstream file(name);
    std::vector<std::vector<double>> data(0);
    int rc = 0;
    std::string line = "sfdfsdffs";
    std::getline(file, line);
    std::cout << line << rc << "!!!!\n";

    while (getline(file, line)) {
        std::cout << line << rc << "\n";
        rc++;
        if (rc>ignore_rows) {
            std::stringstream rss(line);
            std::vector<double> row(0);
            int cc = 0;
            std::string val;
            while (std::getline(rss, val, delimeter)) {
                cc++;
                if (cc>ignore_cols) {
                    std::stringstream convert(val);
                    double d;
                    convert >> d;
                    row.push_back(d);
                }
            }
            if (rss.bad()) {
                std::cout << "bad file";
            } else if (!rss.eof()) {
                std::cout << "format error2";
            } else {
                std::cout << "eof";
            }
            data.push_back(row);
        }
    }
    if (file.bad()) {
        std::cout << "bad file";
    } else if (!file.eof()) {
        std::cout << "format error1";
    } else {
        std::cout << "eof";
    }
    std::cout << std::flush;

    return data;
}

/**
 * The elongation rate is calculated in dependence of water content and temperature,
 * The empirical relationship is located in the private function f.
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
    int n = 7;
    std::vector<double> z_{0., 15., 25., 50., 70., 100., 140.};

    // temperature data
    std::cout << "reading temperature data ...";
    auto field_temp = readCSV("wheat.csv",';',2,1);
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
    double simtime = 7*30; // days
    double dt = 0.5 * 1./24.;  // days
    size_t N = round(simtime/dt);

    assert(N<=field_temp.size()); // check if enough data are available

    for (size_t i=0; i<N; i++) {

        rootsystem.simulate(dt);

        // update field data by
        // water_content.data = ... (type is vector<double>)
        temperature.data =  field_temp.at(i);

    }

    /**
     * Export results (as vtp)
     */
    rootsystem.write(name+".vtp");

    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";
}

} // end namespace CRootBox
