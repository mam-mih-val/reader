void script()
{
    Reader* t = new Reader("~/output.root");
    t->BuildFlowHistograms("Flow");
	//t->BuildQAHistograms("QA");
}