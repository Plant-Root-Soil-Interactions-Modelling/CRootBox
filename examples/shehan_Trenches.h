/**
 * Shehans experimental set-up
 *
 * (C) Virtual analysis of trenches
 *
 * simulates the growing root systems (6*37) in the field
 * uses trenches to analyse the number of roots in a coarse grid
 *
 */
using namespace std;


/**
 * Creates the geometry of the trenches
 *
 * \return a vector of 9 infinitely large trenches
 */
vector<SDF_HalfPlane*> fieldTrenches() {
	int N = 9; // number of trenches
	vector<SDF_HalfPlane*> trenches;
	double dist1 = 18; // cm horizontal distance between plants
	double dist2 = 6; 	// cm vertical distance between the trenches
	for (int i=0; i<N; i++) {
		Vector3d o(1.5*dist1, 5.*dist2+i*dist2,  -160);  // actually, only y component is important
		Vector3d p1(1.5*dist1, 5.*dist2+i*dist2,  0);
		Vector3d p2(3.5*dist1, 5.*dist2+i*dist2,  -160);
		trenches.push_back(new SDF_HalfPlane(o,p1,p2));
	}
	return trenches;
}

/**
 * Creates a field of multiply root system relevant for the trenches (4*21)
 *
 * @param name     filename of the parameter file
 * @param geom     confining geometry
 */
vector<RootSystem*> initializeRootSystems2(string name, SignedDistanceFunction* geom = new SignedDistanceFunction())
{
	auto gen = mt19937(chrono::system_clock::now().time_since_epoch().count());
	auto UD = uniform_real_distribution<double>(0,1); // random stuff, does it work now?
	int M=6;
	int N=37;
	double dist1 = 18; // [cm] (M-1)*18 = 90
	double dist2 = 3; // [cm] (N-1)*3 = 108
	vector<RootSystem*> allRS;
	for (int i=1; i<(M-1); i++) {
		for (int j=8;j<(N-8); j++) {
			RootSystem* rs = new RootSystem();
			allRS.push_back(rs);
			rs->openFile(name);
			rs->setGeometry(geom);
			rs->getRootTypeParameter(4)->theta = 80./180.*M_PI; // fix insertion angle of the basal roots
			rs->getRootSystemParameter()->seedPos = Vector3d(dist1*i,dist2*j,-3); // set position of seed [cm]
			double s = UD(gen);
			rs->setSeed(s); // randomly select a seed
			rs->initialize(4,5);
		}
	}
	return allRS;
}



/**
 * Trench analysis in the virtual field experiment
 */
void shehan_Trenches(string name = "wheat", bool exportVTP = false)
{

	/*
	 * Initialize
	 */
	auto allRS = initializeRootSystems2(name); // no geometry, trenches are destructive and created after roots are grown

	/*
	 * Simulate
	 */
	vector<double> times = {0, 30, 60, 90, 120, 150, 180, 210, 240};
	simulateRS(times, allRS);

	/*
	 * Analysis: cut the roots along the trenches
	 */
	double dist1 = 18; // cm
	int m = 4;
	int n = 32;
	SDF_PlantBox box_(36.,1.e9,160.);
	SDF_RotateTranslate box(&box_, Vector3d(18.+ 27., 0., 0.) );
	vector<SDF_HalfPlane*> trenches = fieldTrenches();
	vector<vector<vector<double>>> finalMT;

	for (size_t t=0; t<times.size()-1; t++) {

		std::cout<<"create empty matrix";
		// create empty n*m matrix
		vector<vector<double>> finalmatrix(n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				finalmatrix[i].push_back(0.);
			}
		}

		// loop over trenches
		int ii=0;
		for (auto& tr : trenches) {
			std::cout << "\nANALYSE TRENCH " << ++ii <<"\n";

			// cut along the trench
			std::cout << "cut \n";
			SegmentAnalyser analyser = getResult(allRS,times.at(t+1));
			SegmentAnalyser cut = analyser.cut(*tr);
			cut.crop(&box); // cut with bounding box

			// split into grid
			std::cout <<"grid \n";
			vector<vector<SegmentAnalyser>> anamatrix1 = cut.distribution2(0,160, tr->o.x-dist1,tr->o.x,n,m);
			vector<vector<SegmentAnalyser>> anamatrix2 = cut.distribution2(0,160, tr->o.x,tr->o.x+dist1,n,m);

			// save root count into matrix
			std::cout<<"count \n";
			for (int i=0; i<n; i++) {
				for (int j=0; j<m; j++) {
					finalmatrix[i][j] += anamatrix1[i][j].segments.size();
					finalmatrix[i][j] += anamatrix2[i][j].segments.size();
				}
			}

//			/*
//			 * Export trenches for last time step
//			 */
//			if (t==times.size()-2) {
//				string vtpname = name + "_trench_"+std::to_string(ii)+".vtp";
//				cut.write(vtpname);
//				if (ii==1) { // only once
//					//analyser.at(t)->write(name+std::to_string(int(t))+".vtp");
//				}
//			}

		}

		finalMT.push_back(finalmatrix); // a lot of finals we have here (one matrix per time step)
	}

	/*
	 * Export the final matrix (one file per matrix)
	 */
	double N = double(trenches.size())*2.*(4.5*5.);
	// to calculate the mean number of roots per cm^2 over the trenches
	for (int j=0; j<m; j++) {
		string matname = name+"_trench_matrix"+std::to_string(j+1)+".txt"; // grid number
		std::ofstream fos;
		fos.open(matname.c_str());
		for (int i=0; i<n ; i++) {
			for (size_t t=0; t<times.size()-1; t++) {
				std::cout << std::fixed << std::setprecision(4)<< finalMT.at(t)[i][j]/N << "\t";
				fos << std::fixed << std::setprecision(4)<< finalMT.at(t)[i][j]/N << "\t";
			}
			std::cout << "\n";
			fos << "\n";
		}
		fos.close();
		std::cout << "\n\n";
	}

	//	/*
	//	 * Export rootsystems (around 4.5GB)
	//	 */
	//	for (size_t i=0; i<times.size()-1; i++) {
	//		string vtpname = name + std::to_string(i+1)+".vtp";
	//		analyser.at(i)->write(vtpname);
	//	}

	if (exportVTP) {
		vector<SignedDistanceFunction*> tr_;
		for (const auto& tr : trenches) { // copy vector
			tr_.push_back(tr);
		}
		string gname = name + "_trench.py";
		SDF_Union trenchgeometry = SDF_Union(tr_);
		allRS[0]->setGeometry(&trenchgeometry); // just for writing
		allRS[0]->write(gname);
	}

}


