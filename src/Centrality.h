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
	DataTreeEvent* fEvent;
	Centrality() {};
	public:
	Centrality(DataTreeEvent* _fEvent, TString FileName);
	~Centrality() {};
	float	GetCentralityClass();
	int		GetNumClasses();
	void	LoadCentralityPercentile(TString FileName);
};