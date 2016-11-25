/**
 * Example testing the suggested coupling to Dumux
 * (suggested coupling)
 * TODO
 */
using namespace std;

void example_dumux()
{

    DumuxRootSystem rootsystem;

    /*
     * Set root parameters
     */
    char name[] = "modelparameter/sorghum.rparam";
    ifstream fis;
    fis.open(name);
    rootsystem.readParameters(fis);
    fis.close();
    //rootsystem.writeParameters(cout); // show parameters

    /*
     * Set plant parameters
     */
    RootSystemParameter plant;
    //plant.read(); // TODO
    rootsystem.plantparam = plant; // not implemented, just default

    /*
     * Set geometry (optionally)
     */
    rootsystem.setGeometry(new SDF_PlantContainer());

    /*
     * Initialize
     */
    rootsystem.initialize();

    /*
     * Simulate
     */
    //    rootsystem.baseRoots.at(0)->param.write(cout); // write base root

    double simtime = 30;
    double dt = 1;
    int N = round(simtime/dt);

    int dumuxID = 100; // assume someting for debugging
    rootsystem.getInitialNodes(); //
    vector<Vector3d> nodes = rootsystem.getNewNodes();
    cout << "Initially, "<< nodes.size() <<" node: "<<nodes.at(0).toString() <<"\n";
    vector<int> dids;
    dids.push_back(dumuxID);
    rootsystem.setDumuxIds(dids);

    for (int i=0; i<N; i++) {

        // simulate
        rootsystem.simulate(dt);

        // nodes and ids
        rootsystem.findNewNodes(double(i)*dt,double(i+1)*dt);
        nodes = rootsystem.getNewNodes();
        vector<int> dids = rootsystem.getDumuxIds();

        // debug
        cout << nodes.size() <<" new nodes: connected to dumux_ids ";
        for (auto const& i :dids) {
            cout << i <<", ";
        }
        cout <<"\n";

        // now add to foamgrid and obtain dumux ids
        vector<int> newids;
        for (size_t i=0; i<dids.size(); i++) {
            dumuxID++;
            newids.push_back(dumuxID);
        }
        //

        rootsystem.setDumuxIds(newids);

    }

    /**
     * Export results (as vtp or rsml)
     */
    //char rsmlname[] = "results.rsml";    //   filename
    char vtpname[] = "sorghum.vtp";    //   filename
    ofstream fos;
    //fos.open(rsmlname);
    fos.open(vtpname);
    rootsystem.writeVTP(fos);
    //rootsystem.writeRSML(fos);
    fos.close();

    cout << "   nodes " << rootsystem.getNumberOfNodes();

}
