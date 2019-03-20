#include <functional>

#include "gauss_legendre.h"

#include "soil.h"
#include "sdf_rs.h"



namespace CRootBox {

class ExudationModel {
public:

    enum IntegrationType { mps_straight = 0, mps = 1, mls = 2 };

    /*
     * Model parameters (same for all roots)
     */
    double Q = 1e-5;
    double Dl = 1e-5; // cm2 / day
    double theta = 0.3;
    double R = 1;
    double k = 1e-6;
    double l = 0.1; // cm

    /*
     *  Numerical parameters
     */
    EquidistantGrid3D grid;
    int type = mps;
    int n0 = 5; // integration points per cm
    double thresh13 = 1e-15; // threshold for Eqn 13
    bool calc13 = true;
    double observationRadius = 5;

    ExudationModel(double width, double depth, int n, RootSystem& rs) :ExudationModel(width, width, depth, n, n, n, rs) {
    }

    ExudationModel(double length, double width, double depth, int nx, int ny, int nz, RootSystem& rs) :grid(EquidistantGrid3D(length, width, depth, nx, ny, nz)) {

        dx3 = (length/nx)*(width/ny)*(depth/nz); // for integration of eqn 13
        simtime = rs.getSimTime();
        roots = rs.getRoots();

        for (const auto& r : roots) {
            if (r->getNumberOfNodes()>1) {
                // root age (until root stopped growing)
                double a = r->getNodeETime(r->getNumberOfNodes()-1) - r->getNodeETime(0);
                age.push_back(a);
                // time when the root stopped growing
                double sTime = r->getNodeETime(r->getNumberOfNodes()-1);
                if (r->active) {
                    stopTime.push_back(0);
                } else {
                    stopTime.push_back(sTime);
                }
                // root tip
                auto t = r->getNode(r->getNumberOfNodes()-1);
                tip.push_back(t);
                // direction towards root base
                auto base = r->getNode(0);
                v.push_back(base.minus(t).times(1./a));
            }

            sdfs.push_back(SDF_RootSystem(*r, observationRadius));

        }

    }

    /**
     * For each root for each grid point
     */
    std::vector<double> calculate() {

        std::fill(grid.data.begin(), grid.data.end(), 0); // set data to zero
        g_.resize(grid.data.size()); // saves last root contribution

        for (size_t ri = 0; ri< roots.size(); ri++) {

            std::cout << "Root #" << ri << "/" << roots.size() << ", age "<< age[ri] << "\n";

            if (age[ri]>0) {

                // per root (passed to integrands)
                r_ = roots[ri]; // eq 11
                n_ = int(n0*r_->length); // number of integration points eq 11
                age_ = age[ri]; // eq 11
                v_ = v[ri]; // for mps_straight, eq 11
                tip_ = tip[ri]; // for mps_straight, eq 11

                st_ = stopTime[ri]; // eq 13
                st_ *= calc13;

                // EQN 11
                for (size_t i = 0; i<grid.nx; i++) {
                    for(size_t j = 0; j<grid.ny; j++) {
                        for (size_t k = 0; k<grid.nz; k++) {

                            x_ = grid.getGridPoint(i,j,k); // integration point

                            if (-sdfs[ri].getDist(x_)<observationRadius) {

                                size_t lind = i*(grid.ny*grid.nz)+j*grid.nz+k;

                                // different flavors of Eqn (11)
                                double c = eqn11(0, age_, 0, l);
                                grid.data[lind] += c;
                                g_[lind] = c;

                            }

                        }
                    }
                }

                // EQN 13
                if (st_>0) { // has stopped growing
                    std::cout << "13!";
                    for (size_t i = 0; i<grid.nx; i++) {
                        std::cout << "*";
                        for(size_t j = 0; j<grid.ny; j++) {
                            for (size_t k = 0; k<grid.nz; k++) {

                                size_t lind = i*(grid.ny*grid.nz)+j*grid.nz+k;
                                if (g_[lind] > thresh13) {

                                    x_ = grid.getGridPoint(i,j,k);
                                    if (-sdfs[ri].getDist(x_)<observationRadius) {

                                        // Eqn (13)
                                        grid.data[lind] += integrate13();

                                    }

                                }
                            }
                        }
                    }
                    std::cout << "13\n";
                }


            } // if ages.at(i)>0

        }

        return grid.data;
    }

    double eqn11(double x0, double xend, double y0, double yend) {
        switch (type) {
        case mps_straight: {
            return gauss_legendre(n_, integrandMPS_straight, this, x0, xend);
        }
        case mps: {
            return gauss_legendre(n_, integrandMPS, this, x0, xend);
        }
        case mls: {
            return gauss_legendre_2D_cube(n_, integrandMLS, this, x0, xend, y0, yend);
        }
        }
        std::cout << "Unknown integration type \n";
        return 0.;
    }

    // simplistic integration in 3d
    double integrate13() {
        double c = 0;
        for (size_t i = 0; i<grid.nx; i++) {
            for(size_t j = 0; j<grid.ny; j++) {
                for (size_t k = 0; k<grid.nz; k++) {
                    Vector3d y = grid.getGridPoint(i,j,k);
                    size_t lind = i*(grid.ny*grid.nz)+j*grid.nz+k;
                    c += integrand13(y,lind, simtime)*dx3;
                }
            }
        }
        return c;
    }

    // integrand Eqn 13
    double integrand13(Vector3d& y, size_t lind, double t) {
        double dt = t-st_;
        double c = to32(R)*g_[lind] / to32(4*Dl*M_PI*dt);
        Vector3d z = x_.minus(y);
        return c*exp(-R/(4*Dl*dt) * z.times(z) - k*dt/R);
    }

    // Returns the linearly interpolated position along the root r at age a
    static Vector3d pointAtAge(Root* r, double a) {
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
            std::cout << "pointAtAge(): warning age is older than the root \n";
            return r->getNode(i-1);
        }
        Vector3d n1 = r->getNode(i-1);
        Vector3d n2 = r->getNode(i);
        double t = (et - r->getNodeETime(i - 1)) / (r->getNodeETime(i) - r->getNodeETime(i - 1)); // t in (0,1]
        return (n1.times(1. - t)).plus(n2.times(t));
    }

    static double to32(double x) { return sqrt(x*x*x); }

    static double to3(double x) { return x*x*x; }

    // point source, root is represented by a single straight line (substituted)
    static double integrandMPS_straight(double t, void* param) {
        ExudationModel* p = (ExudationModel*) param;
        double c = -p->R / ( 4*p->Dl*t );
        double d = 8*(p->theta)*ExudationModel::to32(M_PI*p->Dl*t);

        Vector3d xtip = p->tip_.plus(p->v_.times(t)); // for t=0 at tip, at t=age at base, as above
        Vector3d z = p->x_.minus(xtip);

        return ((p->Q)*sqrt(p->R))/d *exp(c*z.times(z) - p->k/p->R * t); // Eqn (11)
    }

    // moving line source, root is represented by a straight segments
    static double integrandMLS(double t, double l, void* param) {
        ExudationModel* p = (ExudationModel*) param;
        double c = -(p->R) / ( 4*(p->Dl)*t );
        double d = 8*(p->theta)*ExudationModel::to32(M_PI*p->Dl*t);

        Vector3d xtip = ExudationModel::pointAtAge(p->r_, p->age_-t);
        double tl = p->r_->getLength( p->age_-t ); // tip
        if (tl-l<0) { // if root smaller l
            return 0.;
        }
        double agel = p->r_->getAge(tl-l);
        Vector3d nl = p->ExudationModel::pointAtAge(p->r_, agel); // todo ???
        Vector3d z = p->x_.minus(xtip).plus(xtip.minus(nl));

        return ((p->Q)*sqrt(p->R))/d *exp(c*z.times(z) - p->k/p->R * t); // Eqn (11)
    }

    // moving point source, root is represented by a straight segments
    static double integrandMPS(double t, void* param) {
        ExudationModel* p = (ExudationModel*) param;
        double c = -p->R / ( 4*p->Dl*t );
        double d = 8*(p->theta)*ExudationModel::to32(M_PI*p->Dl*t);

        Vector3d xtip = ExudationModel::pointAtAge(p->r_, p->age_-t);
        Vector3d z = p->x_.minus(xtip);

        return ((p->Q)*sqrt(p->R))/d *exp(c*z.times(z) - p->k/p->R * t); // Eqn (11)
    }

    // Root system
    std::vector<Root*> roots;
    std::vector<double> age; // age of root, until it stopped growing
    std::vector<double> stopTime; // time when root stopped growing, 0 if it has not
    std::vector<Vector3d> tip;
    std::vector<Vector3d> v; // direction from tip towards root base
    double simtime;
    double dx3 = 1;
    std::vector<SDF_RootSystem> sdfs; // direction from tip towards root base


    // Set before integrating
    Vector3d x_ = Vector3d(); // integration point
    int n_ = 0;
    Root* r_ = nullptr; // current root
    double age_ = 0;
    Vector3d tip_ = Vector3d();
    Vector3d v_ = Vector3d();
    double st_ = 0; // stop time (eqn 13)
    std::vector<double> g_;

};


}

