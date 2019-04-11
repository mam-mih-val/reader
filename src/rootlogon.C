{
    gSystem->Load("/home/maikhail/DataTree/build/libDataTree.so");
	gROOT->ProcessLine(".L Centrality.cxx+");
    gROOT->ProcessLine(".L Selector.cxx+");
    gROOT->ProcessLine(".L Qvector.cxx+");
    gROOT->ProcessLine(".L Qvector3SE.cxx++");
    gROOT->ProcessLine(".L Flow.cxx+");
    gROOT->ProcessLine(".L Flow3SE.cxx++");
    // gROOT->ProcessLine(".L TrackQA.cxx+");
    // gROOT->ProcessLine(".L EventQA.cxx+");
    gROOT->ProcessLine(".L Reader.cxx++");
}
