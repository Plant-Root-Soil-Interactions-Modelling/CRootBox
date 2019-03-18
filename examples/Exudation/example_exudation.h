#include "gauss_legendre.h"

#include <functional>

namespace CRootBox {

double to32(double x) { return sqrt(x*x*x); }

/**
 * Parameters of the exudation model
 */
class ExudationParameters {
public:

    // Model parameters (same for all roots)
    double M=1e-5;
    double Dt=1e-5; // cm2/s
    double Dl=Dt;
    double theta=0.3;
    double R=1;
    double lambda_=1e-6;
    double l = 0.1; // cm
    double simtime = 10; // days

    // Resolution
    int N = 5; // sample points per cm along the root (for moving source)
    int N3 = 5; // sample points for 3d integration when root stopped growing
    double intRange = 1.; // integration range for stopped roots

    // Update for each root
    double age; // days
    double stopTime; // days
    Vector3d tip;
    Vector3d v;
    Root* r;
    std::vector<double> G = std::vector<double>(N3*N3*N3);
    bool initG = false;

    // Update integral position
    Vector3d pos;

    double g=0; // internal

};

/**
 * Returns an interpolated point of root r at age a
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
        std::cout << "pointAtAge(): warning... ";
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
 * simplistic integration in 3d
 */
double integrate3d(Vector3d p, int N, double range, double (*f)(Vector3d, double, void*), void* param) {
    ExudationParameters* p_ = (ExudationParameters*) param;
    double w = range/2;
    double dx3 = (range/N)*(range/N)*(range/N);
    double c = 0;
    for (int i=0; i<N; i++) {
        for (int j=0; i<N; i++) {
            for (int k=0; i<N; i++) {
                p_->g = p_->G[i*N*N+j*N+k];
                c += f(Vector3d(p.x-w+i*2*w, p.y-w+j*2*w, p.z-w+k*2*w), p_->simtime, param )*dx3;
            }
        }
    }
    return c;
}

/**
 * Eqn 13
 */
double integrand13(Vector3d x, double t, void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    double c =to32(p->R)*p->g / (4*p->Dl*M_PI*to32(p->simtime-p->stopTime));

    return c;
}


/**
 * segment moving line source, root is represented by a straight segments
 */
double integrandSMLS(double t, double l, void* param) {
    ExudationParameters* p = (ExudationParameters*) param;

    Vector3d xtip = pointAtAge(p->r, p->age-t);

    double tl = p->r->getLength( p->age-t ); // tip
    if (tl-l<0) { // if root smaller l
        return 0.;
    }
    double agel = p->r->getAge(tl-l);
    Vector3d nl = pointAtAge(p->r, agel);

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

    Vector3d xtip = pointAtAge(p->r, p->age-t);
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

    Vector3d xtip = p->tip.plus(p->v.times(t)); // for t=0 at tip, at t=age at base, as above
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

    Vector3d xtip = p->tip.plus(p->v.times(p->age-t));
    Vector3d x = p->pos.minus(xtip);

    double c = -p->R / (4*p->Dl*(p->age));
    double d = 8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dt*p->Dt*p->Dl*(p->age-t)*(p->age-t)*(p->age-t));

    return (p->M*sqrt(p->R)) / d * exp(c*x.times(x)  - p->lambda_/p->R*(p->age - t)); // Eqn (11)
}

/**
 * TODO doc me!
 *
 * type: 0 (integrandMPS) moving point source (root is a single straight line)
 * type: 1 (integrandSMPS) moving point source (root is represented by segments)
 * type: 2 (integrandMPS2) moving point source (root is a single straight line), same substitution as type 1 (for debugging)
 * type: 3 (integrandSMLS) moving line source (root is represented by segments)
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

    // tips
    auto tipsI = rootsystem.getRootTips();
    std::vector<Vector3d> tips;
    for (auto i : tipsI) {
        tips.push_back(nodes[i]);
    }

    // bases
    auto basesI = rootsystem.getRootBases();
    std::vector<Vector3d> bases;
    for (auto i : basesI) {
        bases.push_back(nodes[i]);
    }

    // ages
    std::vector<double> ages; // age per root
    for (const auto& r : roots) {
        if (r->getNumberOfNodes()>1) {
            ages.push_back(r->getNodeETime(r->getNumberOfNodes()-1) - r->getNodeETime(0));
        }
    }

    // stopTime
    params.simtime = rootsystem.getSimTime();
    std::vector<double> stopTime; // age per root
    for (const auto& r : roots) {
        if (r->getNumberOfNodes()>1) {
            double sTime = r->getNodeETime(r->getNumberOfNodes()-1);
            if (r->active) {
                stopTime.push_back(0);
            } else {
                stopTime.push_back(sTime);
            }
        }
    }

    std::vector<double> allc = std::vector<double> (X*Y*Z); // many double
    std::fill(allc.begin(), allc.end(), 0); // all zero

    for (size_t i = 0; i< roots.size(); i++) {

        std::cout << "Root #" << i+1 << "/" << roots.size() << "\n";

        if (ages.at(i)>0) {

            params.v = Vector3d(tips[i].minus(bases[i]).times(-1./ages[i]));
            params.tip = tips[i];
            params.age = ages[i];
            params.r = roots[i];
            params.stopTime = stopTime[i];

            int n = std::ceil(params.N*roots[i]->length);
            std::cout << "root " << i << " N = "<< n << "\n";

            for (int x=0; x<X; x++) {
                for(int y=0; y<Y; y++) {
                    for (int z=0; z<Z; z++ ) {

                        params.pos = Vector3d(((double(x)-double(X)/2.)/X)*width,
                            ((double(y) - double(Y) / 2.) / Y) * width,
                            ((-double(z)) / Z) * depth);

                        // if (sdfRS.getDist(pos)<intRange) { // todo

                        int ind = x*Y*Z+y*Z+z;

                        // different flavors of Eqn (11)
                        double g = 0;
                        if (type<3) {
                            g = gauss_legendre(n, integrand, &params, 0, params.age);
                        } else {
                            g = gauss_legendre_2D_cube(n, integrandSMLS, &params, 0, params.age, 0, params.l);
                        }
                        allc[ind] += g;

                        if (params.stopTime>0) { // root stopped growing

                            // calculate G
                            if (!params.initG) { // called once per root
                                int N3 = params.N3;
                                double w = params.intRange/2;
                                auto p = params.pos;
                                for (int i=0; i<N3; i++) {
                                    for (int j=0; i<N3; i++) {
                                        for (int k=0; i<N3; i++) {
                                            params.pos = Vector3d(p.x-w+i*2*w,p.y-w+j*2*w,p.z-w+k*2*w);
                                            if (type<3) {
                                                params.G[i*N3*N3+j*N3+k] = gauss_legendre(n, integrand, &params, 0, params.age);
                                            } else {
                                                params.G[i*N3*N3+j*N3+k] = gauss_legendre_2D_cube(n, integrandSMLS, &params, 0, params.age, 0, params.l);
                                            }
                                        }
                                    }
                                }
                                params.pos = p;
                                params.initG = true;
                            }

                            // integrate Eqn (13)
                            allc[ind] += integrate3d(params.pos, params.N3, params.intRange, integrand13, &params);

                        }

                        // }

                    }
                }
            }

            params.initG = false; // next root

        } // if ages.at(i)>0

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




/* Unused code snipets
 *
 // lengths
    std::vector<std::vector<double>> rootLength = std::vector<std::vector<double>>(0);
    std::vector<std::vector<double>> cumLength = std::vector<std::vector<double>>(0);
    for (const auto& r : roots) {
        std::vector<double> rl = std::vector<double>(0); // root length
        rl.push_back(0.);
        for (size_t i=1; i<r->getNumberOfNodes(); i++) {
            rl.push_back( (r->getNode(i).minus(r->getNode(i-1))).length() );
        }
        rootLength.push_back(rl);
        std::vector<double> cl = rl; // cumulative length
        std::partial_sum(rl.begin(), rl.end(), cl.begin(), std::plus<double>());
        cumLength.push_back(cl);
//        std::cout << "Length: \n";
//        for (size_t j=0; j<rl.size(); j++) { std::cout << rl[j] << ", "; }
//        std::cout << "\n";
//        for (size_t j=0; j<cl.size(); j++) { std::cout << cl[j] << ", "; }
//        std::cout << "\n";
    }
 *
 *
 */


