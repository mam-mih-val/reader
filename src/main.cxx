#include "Reader.h"

int main()
{
    Reader* t = new Reader("~/output.root");
    t->BuildQvectorHistograms("Qvector");
	//t->BuildQAHistograms("QA");cd
    return 0;
}