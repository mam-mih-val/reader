#pragma once

#include "TH1F.h"
#include "TFile.h"

#include "HADES_constants.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

using std::vector;
using std::cout;
using std::endl;

class Centrality
{
	private:
	TH1F* hCentralityPercentile;
	public:
	Centrality() {};
	Centrality(TString FileName);
	~Centrality() {};
	float	GetCentralityClass(DataTreeEvent* fEvent);
	int		GetNumClasses() { return hCentralityPercentile->GetNbinsX(); }
	void	LoadCentralityPercentile(TString FileName);
};