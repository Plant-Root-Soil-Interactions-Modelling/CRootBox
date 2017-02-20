#include "gauss_legendre.h"

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

class Exudation_params {
public:
	// static
	double M=1e-5;
	double Dt=1e-5;  //cm2/s
	double Dl=Dt;
	double theta=0.3;
	double R=1;
	double lambda_=1e-6;

	// update for each root
	double age_r;
	Vector3d tip;
	Vector3d v;

	// update for each position
	Vector3d pos;
};

//
// Define integrand
//
double integrand(double t,void* param) {
	Exudation_params* p = (Exudation_params*) param;

	double dn = 4*p->R*p->Dl*(p->age_r-t);

	double x1 = (p->pos.x - p->tip.x)*p->R -p->v.x*(p->age_r-t);
    double exp_x = -x1*x1/dn;

    double y1= (p->pos.y-p->tip.y)*p->R -p->v.y*(p->age_r-t);
	double exp_y = -y1*y1/dn;

    double z1 = (p->pos.z-p->tip.z)*p->R -p->v.z*(p->age_r-t);
    double exp_z = -z1*z1/dn;

    return p->M/(8*p->theta*sqrt(M_PI*M_PI*M_PI*p->Dt*p->Dt*p->Dl*(p->age_r-t)*(p->age_r-t)*(p->age_r-t)))*
    		exp( exp_x + exp_y + exp_z  - p->lambda_*(p->age_r-t)/p->R );
}


using namespace std;

void example_exudation()
{
    RootSystem rootsystem;

    string name = "Zea_mays_1_Leitner_2010";

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

    cout << "Number of roots " << rootsystem.getRoots().size() << "\n";

    /*
     * Export final result (as vtp)
     */
    rootsystem.write(name+".vtp",RootSystem::ot_segments); // use ot_polylines for nicer visualization, ot_segments for animations

    /*
     * Exudation
     */
    Exudation_params params;
    double L = 1;

    std::vector<Root*> roots = rootsystem.getRoots();
    auto tips = rootsystem.getRootTips(roots); // nodes
    auto bases = rootsystem.getRootBases(roots);

    vector<double> ages;
    for (const auto& r : roots) {
    	if (r->getNumberOfNodes()>1) {
    		ages.push_back(simtime - r->getNodeETime(0));
    	}
    }
    assert(roots.size()==ages.size());


    int X = 10;
    int Y = 10;
    int Z = 10;
    double width = 50;
    double depth = 100;

    vector<double> allc = vector<double> (X*Y*Z); // many double
    for (size_t i=0; i<allc.size(); i++) {
    	allc[i] = 0;
    }

    for (size_t i = 0; i< roots.size(); i++) {

    	cout << i << "\n";

    	if (ages.at(i)>0) {

            double vx = - (tips[i].x-bases[i].x)/ages[i];
            double vy = - (tips[i].y-bases[i].y)/ages[i];
            double vz = - (tips[i].z-bases[i].z)/ages[i];
            Vector3d v = Vector3d(vx,vy,vz);
			v.normalize();
			double vl = v.length();
            double l = min(L,vl);   // if root is not yet as long as L

            params.v = v;
            params.tip = tips[i];
			params.age_r = ages[i];

            for (int x=0; x<X; x++) {
            	for(int y=0; y<Y; y++) {
            		for (int z=0; z<Z; z++ ) {

            			params.pos = Vector3d(((double(x)-double(X)/2.)/X)*width,
            					((double(y)-double(Y)/2.)/Y)*width,
								((-double(z))/Z)*depth);

            			int ind = x*Y*Z+y*Z+z;

            			allc[ind] += gauss_legendre(5,integrand,&params,0,params.age_r);
            			// gauss_legendre_2D_cube()
            		}
            	}
            }


    	}

    }

    /*
     * write as txt file
     */
    std::ofstream fos;
    fos.open(name+".txt");
    for (size_t i=0; i<allc.size(); i++) {
    	fos << allc[i] << "\n";
    }
    fos.close();

    cout << "fin \n";


}
