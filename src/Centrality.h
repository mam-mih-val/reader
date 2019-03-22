#pragma once

#include "TH1F.h"
#include "TFile.h"

#include "HADES_constants.h"

#define  DATATREE_SHINE
#include "DataTreeEvent.h"

class Centrality
{
	private:
	TH1F* hCentralityPercentile;
	public:
	Centrality() {};
	Centrality(TString FileName);
	~Centrality() {};
	float	GetCentralityClass(DataTreeEvent* fEvent);
	void	LoadCentralityPercentile(TString FileName);
};