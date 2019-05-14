// usage:
// ./DT_Reader ["input filename"] [command] [configuration] ["output filename"]
// commands:
//      eventqa -   build Quality Assurance histograms of Event varriables
//      trackqa -   build Quality Assurance histograms of Track varriables
//      qvector -   build Q-vector histograms
//      flow -      build Flow histograms
// configiration (optional for flow & qvector)
//      [1/0][1/0][1/0]
//      first -     0/1 - channel selection off/on
//      second -    0/1 - Using Fw-Adc, Fw-Z
//      third -     0/1 - Using All spectators or only protons 

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
    auto reader = new Reader(argv[1]);
    std::string command = argv[2];
    std::string config = argv[3];
    int channelSelection = (int)config.at(0)-'0';
    int fwZ = (int)config.at(1)-'0';
    int protons = (int)config.at(2)-'0';
    switch (argc)
    {
        case 4:
            if( command == "eventqa" )
            {
                reader->BuildEventQaHistograms(argv[3]);
                return 0;
            }
            if( command == "trackqa" )
            {
                reader->BuildTrackQaHistograms(argv[3]);
                return 0;
            }
            if( command == "qvector" )
            {
                reader->BuildQvector3SeHistograms(argv[3]);
                return 0;
            }
            if( command == "flow" )
            {
                reader->BuildFlow3SeHistograms(argv[3]);
                return 0;
            }
            break;
        case 5:
            if( command == "qvector" )
            {
                reader->BuildQvector3SeHistograms(argv[4], channelSelection, fwZ, protons);
                return 0;
            }
            if( command == "flow" )
            {
                reader->BuildFlow3SeHistograms(argv[4], channelSelection, fwZ, protons);
                return 0;
            }
            break;    
        default:
            cout << "Error: incorrect number of arguments" << endl;
            return 1;
            break;
    }
    
    cout << "Error: unknow command" << endl;
    return 2;
}