#include "Reader.h"

int main(int argc, char** argv)
{
    auto t = new Reader(argv[1]);
    switch (argc)
    {
        case 2:
            t->BuildQvectorHistograms("Qvector");
            break;
        case 3:
            t->BuildQvectorHistograms(argv[2]);
            break;
        default:
            cout << "Error: DT_Reader get 1 or 2 arguments" << endl;
            return 1;
    }
    return 0;
}