// usage:
// ./DT_Reader [command] ["input filename"] ["output filename"]
// commands:
//      qa -        build Quality Assurance histograms
//      qvector -   build Q-vector histograms
//      flow -      build Flow histograms

#include "Reader.h"
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
    if(argc<4)
    {
        cout << "Error: incorrect number of arguments" << endl;
        return 1;
    }
    auto t = new Reader(argv[2]);
    std::string command = argv[1];
    if( command == "qa" )
    {
        t->BuildQAHistograms(argv[3]);
        return 0;
    }
    if( command == "qvector" )
    {
        t->BuildQvector3SeHistograms(argv[3]);
        return 0;
    }
    if( command == "flow" )
    {
        t->BuildFlow3SeHistograms(argv[3]);
        return 0;
    }
    
    cout << "Error: unknow command" << endl;
    return 2;
}