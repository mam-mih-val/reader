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
        std::cerr << "Error 1: incorrect number of arguments" << endl;
        std::cerr << argc << " arguments were given, minimum 4 is required" << endl;
        return 1;
    }
    std::string command = argv[1];
    bool channelSelection;
    std::string signal="adc";
    float minSignal=0;
    float maxSignal=999;
    auto reader = new Reader(argv[ argc-2 ]);
    if(argc>4)
    {
        for(int i=2; i<argc-2; i+=2)
        {
            std::string flag=argv[i];
            if(flag=="--signal")
            {
                signal=argv[i+1];
                continue;
            }
            if(flag=="--perchannel")
            {
                channelSelection=std::atoi(argv[i+1]);
                continue;
            }
            if(flag=="--min")
            {
                minSignal=std::atof(argv[i+1]);
                continue;
            }
            if(flag=="--max")
            {
                maxSignal=std::atof(argv[i+1]);
                continue;
            }
            std::cerr << "Error 3: Unknown flag" << endl;
            std::cerr << "Only --signal, --perchannel, --min, --max flags are supported" << endl;
            return 3;
        }
    }
    if( command=="eventqa" )
    {
        reader->BuildEventQaHistograms(argv[argc-1]);
        return 0;
    }
    if( command=="trackqa" )
    {
        reader->BuildTrackQaHistograms(argv[argc-1]);
        return 0;
    }
    if( command=="qvector" )
    {
        cout << "Configuration:" << endl;
        cout << "Channel Selection: " << channelSelection << endl;
        cout << "Signal type: " << signal << endl;
        cout << "Minimal signal: " << minSignal << endl;
        cout << "Maximal signal: " << maxSignal << endl;
        reader->BuildQvector3SeHistograms(argv[argc-1], channelSelection, signal, minSignal, maxSignal);
        return 0;
    }
    if( command=="flow" )
    {
        cout << "Configuration:" << endl;
        cout << "Channel Selection: " << channelSelection << endl;
        cout << "Signal type: " << signal << endl;
        cout << "Minimal signal: " << minSignal << endl;
        cout << "Maximal signal: " << maxSignal << endl;
        reader->BuildFlow3SeHistograms(argv[argc-1], channelSelection, signal, minSignal, maxSignal);
        return 0;
    }
    std::cerr << "Error 2: unknow command" << endl;
    std::cerr << command << " is not known. Only eventqa, trackqa, qvector, flow commands are supported" << endl;
    return 2;
}