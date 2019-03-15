#include "gauss_legendre.h"

#include <functional>

namespace CRootBox {

/**
 * Example Exudation
 *
 * 1) Opens plant and root parameters from a file
 * 2) Simulates root growth
 * 3) Outputs a VTP (for vizualisation in ParaView)
 *    In Paraview: use tubePLot.py script for fancy visualisation (Macro/Add new macro...), apply after opening file
 *
 *  Computes analytical solution of moving point/line sources based on Carslaw and Jaeger
 */

/**
 * Parameters of the exudation
 */
class ExudationParameters {
public:
    // static
    double M=1e-5;
    double Dt=1e-5; // cm2/s
    double Dl=Dt;
    double theta=0.3;
    double R=1;
    double lambda_=1e-6;
    double l = 0.1; // cm

    // update for each root
    double age_r;
    Vector3d tip;
    Vector3d v;
    Root* r;

    // update for each position
    Vector3d pos;

    // sample pointsp per cm
    int N = 5;
};

/**
 * returns an interpolated point of root r with age a
 */
Vector3d pointAtAge(Root* r, double a) {
    a = std::max(0.,a);
    double et = r->getNodeETime(0)+a; // age -> emergence time
    size_t i=0;
    while (i<r->getNumberOfNodes()) {
        if (r->getNodeETime(i)>et) { // first index bigger than emergence time, interpolate i-1, i
            break;
        }
        i++;
    }
    if (i == r->getNumberOfNodes()) { // this happens if a root has stopped growing
        return r->getNode(i-1);
    }
    Vector3d n1 = r->getNode(i-1);
    Vector3d n2 = r->getNode(i);
    double t = (et - r->getNodeETime(i - 1)) / (r->getNodeETime(i) - r->getNodeETime(i - 1)); // t in (0,1]
    //    std::cout << "root " << r->id << " age between (" << r->getNodeETime(i - 1) << " and " << r->getNodeETime(i) << "], at " << et << "\n";
    //    std::cout << " pos " << (n1.times(1. - t)).plus(n2.times(t)).toString() << "\n";
    return (n1.times(1. - t)).plus(n2.times(t));
}

/**
 * segment moving line source, root is represented by a straight segments
 */
double integrandSMLS(double t, double t2,  void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    if (p->age_r-t-t2 < 0) { // makes sense, does it?
        return 0.;
    }

    Vector3d nl = pointAtAge(p->r, p->age_r-t-t2);
    Vector3d xtip = pointAtAge(p->r, p->age_r-t);
    Vector3d x = p->pos.minus(xtip).plus(xtip.minus(nl));

    double c = -p->R / ( 4*p->Dl*t );
    double d = 8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dl*p->Dl*p->Dl*t*t*t);

    return (p->M*sqrt(p->R))/d *exp(c*x.times(x) - p->lambda_/p->R * t); // Eqn (11)
}


/**
 * segment moving point source, root is represented by a straight segments
 */
double integrandSMPS(double t, void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    Vector3d xtip = pointAtAge(p->r, p->age_r-t);
    Vector3d x = p->pos.minus(xtip);

    double c = -p->R / ( 4*p->Dl*t );
    double d = 8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dl*p->Dl*p->Dl*t*t*t);

    return (p->M*sqrt(p->R))/d *exp(c*x.times(x) - p->lambda_/p->R * t); // Eqn (11)
}

/**
 * moving point source, root is represented by a single straight line (substituted)
 */
double integrandMPS2(double t, void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    Vector3d xtip = p->tip.plus(p->v.times(t)); // for t=0 at tip, at t=age_r at base, as above
    Vector3d x = p->pos.minus(xtip);

    double c = -p->R / ( 4*p->Dl*t );
    double d = 8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dl*p->Dl*p->Dl*t*t*t);

    return (p->M*sqrt(p->R))/d *exp(c*x.times(x) - p->lambda_/p->R * t); // Eqn (11)
}

/**
 * moving point source, root is represented by a single straight line
 */
double integrandMPS(double t, void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    Vector3d xtip = p->tip.plus(p->v.times(p->age_r-t));
    Vector3d x = p->pos.minus(xtip);

    double c = -p->R / (4*p->Dl*(p->age_r-t));
    double d = 8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dt*p->Dt*p->Dl*(p->age_r-t)*(p->age_r-t)*(p->age_r-t));

    return (p->M*sqrt(p->R)) / d * exp(c*x.times(x)  - p->lambda_/p->R*(p->age_r - t)); // Eqn (11)
}

/**
 * TODO doc me!
 *
 * type: 0 moving point source (root is a single straight line)
 * type: 1 moving point source (root is represented by segments)
 * type: 2 moving point source (root is a single straight line), same substitution as type 0 (for debugging)
 * type: 3 moving line source (root is represented by segments)
 *
 */
std::vector<double> getExudateConcentration(RootSystem& rootsystem, ExudationParameters& params, int X, int Y, int Z, double width, double depth,
    int type = 0) {

    const auto& roots = rootsystem.getRoots();
    const auto& nodes = rootsystem.getNodes();

    double (*integrand)(double, void*);

    if (type == 0) {
        std::cout << "getExudateConcentration() with linear moving points source\n";
        integrand = integrandMPS;
    }
    if (type == 1) {
        std::cout << "getExudateConcentration() with segment wise moving point source\n";
        integrand = integrandSMPS;
    }
    if (type == 2) {
        std::cout << "getExudateConcentration() with segment wise moving point source\n";
        integrand = integrandMPS2;
    }

    auto tipsI = rootsystem.getRootTips();
    std::vector<Vector3d> tips;
    for (auto i : tipsI) {
        tips.push_back(nodes[i]);
    }
    auto basesI = rootsystem.getRootBases();
    std::vector<Vector3d> bases;
    for (auto i : basesI) {
        bases.push_back(nodes[i]);
    }

    double simtime = rootsystem.getSimTime();

    std::vector<double> ages; // age per root
    for (const auto& r : roots) {
        if (r->getNumberOfNodes()>1) {
            ages.push_back(simtime - r->getNodeETime(0));
        }
    }
    assert(roots.size()==ages.size());

    std::vector<double> allc = std::vector<double> (X*Y*Z); // many double
    for (size_t i=0; i<allc.size(); i++) {
        allc[i] = 0;
    }

    for (size_t i = 0; i< roots.size(); i++) {

        std::cout << i << "/" << roots.size()-1 << "\n";

        if (ages.at(i)>0) {

            double vx = - (tips[i].x-bases[i].x)/ages[i];
            double vy = - (tips[i].y-bases[i].y)/ages[i];
            double vz = - (tips[i].z-bases[i].z)/ages[i];

            params.v = Vector3d(vx,vy,vz);
            params.tip = tips[i];
            params.age_r = ages[i];
            params.r = roots[i];

            double r = roots[i]->param.r;
            double age = params.l/r; // days
            int n = std::ceil(params.N*roots[i]->length);

            std::cout << "root " << i << " N = "<< n << "\n";

            for (int x=0; x<X; x++) {
                for(int y=0; y<Y; y++) {
                    for (int z=0; z<Z; z++ ) {

                        params.pos = Vector3d(((double(x)-double(X)/2.)/X)*width,
                            ((double(y) - double(Y) / 2.) / Y) * width,
                            ((-double(z)) / Z) * depth);

                        int ind = x*Y*Z+y*Z+z;
                        // int ind = z*Y*X+y*X+x; // one is c, one is Fortran ordering

                        if (type<3) {
                            allc[ind] += gauss_legendre(n, integrand, &params, 0, params.age_r);
                        } else {
                            allc[ind] += gauss_legendre_2D_cube(n, integrandSMLS, &params, 0, params.age_r, 0, age);
                        }

                    }
                }
            }

        } // if

    }

    return allc;
}


/**
 *
 */
void example_exudation()
{
    RootSystem rootsystem;
    std::string name = "Zea_mays_1_Leitner_2010";

    /*
     * Plant and root parameter from a file
     */
    rootsystem.openFile(name);
    rootsystem.writeParameters(std::cout);

    // THETA by hand
    rootsystem.getRootTypeParameter(2)->theta = M_PI/4.;
    rootsystem.getRootTypeParameter(2)->thetas = 0; // no std

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
    double simtime = 10.; // 20, 40, 60 days
    double dt = 1.; // try other values here
    int N = round(simtime/dt);
    for (int i=0; i<N; i++) {
        rootsystem.simulate(dt);
    }

    std::cout << "Number of roots " << rootsystem.getRoots().size() << "\n";

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp");

    /*
     * Exudation
     */
    ExudationParameters params;

    std::vector<double> allC = getExudateConcentration(rootsystem, params, 10, 10, 100, 30, 50);

    //    /*
    //     * write as txt file
    //     */
    //    std::ofstream fos;
    //    fos.open(name+".txt");
    //    for (size_t i=0; i<allc.size(); i++) {
    //    	fos << allc[i] << "\n";
    //    }
    //    fos.close();
    //
    //    cout << "fin \n";

}

} // end namespace CRootBox
