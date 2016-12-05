/**
 * Example 3b (for Shehan)
 *
 * Creates N*N root systems at different locations with rhizotube geometry
 * Uses the AnalysisSDF class to write the VTP, and RootBox segments into one file.
 *
 */
using namespace std;

vector<RootSystem*> example3b()
{
  auto gen = default_random_engine(chrono::system_clock::now().time_since_epoch().count());
  auto UD = uniform_real_distribution<double>(0,1); // random stuff

  string name = "wheat";

  /*
   * Rhizotubes as obstacles (from Example_Rhizotubes.m)
   */
  // Box
  double boxX = 96.;
  double boxY = 126.;
  double boxZ = 130.;
  SDF_PlantBox box(boxX,boxY,boxZ);
  // Single Rhizotube
  double r = 2.*3.2; // cm
  SDF_PlantContainer rhizotube(r,r,96.,false);
  SDF_RotateTranslate rhizoX(&rhizotube, 90., SDF_RotateTranslate::yaxis, Vector3d(boxX/2.,0.,0));
  // Experimental setting
  vector<SignedDistanceFunction*> rhizotubes_;
  const int tubeN=6;
  double y_[tubeN] = { -30, -18, -6, 6, 18, 30 };
  double z_[tubeN]= { -10, -20, -40, -60, -80, -120 };
  for (int i=0; i<tubeN; i++) {
      rhizotubes_.push_back(new SDF_RotateTranslate(&rhizoX, 0, SDF_RotateTranslate::xaxis, Vector3d(0,y_[i],z_[i])));
  }
  SDF_Union rhizotubes(rhizotubes_);
  // Final geometry
  SDF_Difference geometry(&box, &rhizotubes);
  SDF_RotateTranslate geometry2(&geometry, 0, SDF_RotateTranslate::yaxis, Vector3d(boxX/2.,boxY/2.,0));

  /*
   * Creates and initializes N*N root systems
   */
  int N=2;
  double dist = 20; // [cm]
  vector<RootSystem*> allRS;
  for (int i=0; i<N; i++) {
      for (int j=0;j<N; j++) {
          RootSystem* rs = new RootSystem();
          allRS.push_back(rs);
          rs->openFile(name);
          rs->getRootTypeParameter(4)->theta = 80./180.*M_PI; // fix insertion angle of the basal roots
          rs->rsparam.seedPos = Vector3d(dist*i+dist/2.,dist*j+dist/2,-3); // set position of seed [cm]
          rs->setGeometry(&geometry2);
          rs->setSeed(double(UD(gen))); // randomly select a seed
          rs->initialize(4,5);
      }
  }

  /*
   * Simulate
   */
  double simtime = 120; // days
  for (auto rs : allRS) { // simulate all
      rs->simulate(simtime);
      cout << "Finished with a total of " << rs->getNumberOfNodes()<< " nodes\n";
  }

  /*
   * Export results into a single file containing all root systems
   */
  AnalysisSDF allresults;
  for (const auto& rs : allRS) { // copy all to the segment analyser
      allresults.addSegments(*rs);
  }

  // Export VTP
  allresults.write(name+"_all.vtp");

  // Export RootBox segments
  allresults.write(name+"_all.txt");

  /*
   * Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)
   */
  allRS.at(0)->write(name+".py");

  return allRS;
}
