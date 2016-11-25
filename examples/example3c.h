/**
 * Example 3c (for Shehan)
 *
 * runs Example3b and then only shows the segments that are visible friom within
 *
 */
using namespace std;

void example3c()
{
  auto allRS = example3b(); // run example3b and use roots systems for the following analysis

  string name = "wheat";

  /*
   * Rhizotubes a little larger than in example3b
   */
  // Box
  double boxX = 96.;
  double boxY = 126.;
  double boxZ = 130.;
  SDF_PlantBox box(boxX,boxY,boxZ);
  // A single Rhizotube
  double r = 2.5*3.2; // cm, larger than in example 3b
  SDF_PlantContainer rhizotube(r,r,96.,false);
  SDF_RotateTranslate rhizoX(&rhizotube, 90., SDF_RotateTranslate::yaxis, Vector3d(boxX/2.,0.,0));
  // The experimental setting
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
  SDF_Complement visible(&geometry2);

  /*
   * Create the analyser classes
   */
  AnalysisSDF allresults;
  for (const auto& rs : allRS) {
      allresults.addSegments(*rs);
  }

  /*
   * Crop to visible segments
   */
  allresults.crop(&visible);

  /*
   * Total length and surface
   */
  double l = allresults.getSummed(RootSystem::st_length);
  double s = allresults.getSummed(RootSystem::st_surface);
  int nor = allresults.getNumberOfRoots();
  std::cout << "Visible Length " << l << " cm \n";
  std::cout << "Visible Surface " << s << " cm2 \n";
  std::cout << "Number of roots " << nor << "\n";

  /*
   * Export visible rhizotrons segments
   */
  cout << "writing VTP... \n";
  allresults.write(name+"_croped.vtp");
  cout << "writing RootBox segments... \n";
  allresults.write(name+"_all.txt");


}
