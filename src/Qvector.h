#pragma once
#include "HADES_constants.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

class Qvector
{
    private:
    Float_t fQ[2][2];
    DataTreeEvent* fEvent;
    public:
    Qvector();
    ~Qvector() {};
    void        Estimate(DataTreeEvent* fEvent, bool bSubEvent=0);
    void        Recenter(Float_t* fCorrection);
    Float_t     GetComponent(int i, int j=0);
    Float_t     GetPsiEP(int j=0);
};
