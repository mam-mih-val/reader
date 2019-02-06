#include "Qvector.h"

Qvector::Qvector()
{
    fEvent = new DataTreeEvent;
    for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            fQ[i][j]=0;
        }
    }
}

void Qvector::Estimate(DataTreeEvent* fEvent, bool bSubEvent)
{
    for (int i=0;i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            fQ[i][j]=0;
        }
    }
    Int_t iNPSDModules = fEvent->GetNPSDModules();
    DataTreePSDModule* fModule;
    Float_t fPhi, fChargeModule;
    if (!bSubEvent) 
    {
        Float_t fChargeSum;
        for(unsigned int i=0; i<iNPSDModules;i++)
        {
            fModule = fEvent->GetPSDModule(i);
            if( fModule->GetId()<0 )
                continue;
            fPhi = fModule->GetPhi();
            fChargeModule = fModule->GetEnergy();
            fQ[0][0]+=fChargeModule*cos(fPhi);
            fQ[1][0]+=fChargeModule*sin(fPhi);
            fChargeSum+=fChargeModule;
        }
        fQ[0][0]/=fChargeSum;
        fQ[1][0]/=fChargeSum;
    }
    else
    {
        Float_t fChargeSum[2];
        fChargeSum[0]=0; fChargeSum[1]=0;
        vector<DataTreePSDModule*> vModules;
        for(unsigned int i=0; i<iNPSDModules;i++)
        {
            fModule = fEvent->GetPSDModule(i);
            if( fModule->GetId()<0 )
                continue;
            vModules.push_back(fModule);
        }
        random_shuffle( vModules.begin(),vModules.end() );
        for(unsigned int i=0;i<vModules.size();i++)
        {
            fChargeModule = vModules[i]->GetEnergy();
            fPhi = vModules[i]->GetPhi();
            int p = i%2;
            fQ[0][p]+=fChargeModule*cos(fPhi);
            fQ[1][p]+=fChargeModule*sin(fPhi);
            fChargeSum[p]+=fChargeModule;
        }
        for (int i=0;i<2;i++)
        {
            for (int j=0;j<2;j++)
            {
                fQ[i][j]/=fChargeSum[j];
            }
        }
    }
}

void Qvector::Recenter(Float_t* fCorrection)
{
    for(int j=0;j<2;j++)
    {
        fQ[0][j]-=fCorrection[0];
        fQ[1][j]-=fCorrection[1];
    }
}

Float_t Qvector::GetComponent(int i, int j)
{
    return fQ[i][j];
}

Float_t Qvector::GetPsiEP(int j)
{
    return atan2(fQ[1][j],fQ[0][j]);
}
