{
    gSystem->Load("/home/DataTree/build/libDataTree.so");
    gROOT->ProcessLine(".L Selector.cxx");
    gROOT->ProcessLine(".L Reader.cxx");
}
