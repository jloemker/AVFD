#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Event.h"
#include "Particle.h"
#include "CalculateFlowCME.h"
using namespace std;

//Int_t centID;

void runSingleCentrality(Int_t centID){
//Int_t centID;
//for(Int_t centID = 0; centID < 8; centID++){
	//Loop over all centrality bins
	std::cout << "Centrality ID: " << centID <<std::endl;
	        //Int_t centID = {0,1,2,3,4,5,6,7}; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70

	//should be according to the number of files per centrality
	const Int_t nSplit = 10;
	//TFile *f[nCentBin][nSplit];

	CalculateFlowCME *fQC = new CalculateFlowCME("CalculateFlowCME");
	fQC->UserCreateOutputObjects();
	//fQC->SetmaxPtCut(3);
	//fQC->SetminPtCut(0.2);
	//fQC->SetminNtrackCut(500);
	//fQC->SetmaxEtaCut(0.8);
	fQC->SetmaxPtCut(5);
	fQC->SetminPtCut(0.2);
	fQC->SetminNtrackCut(500);
	fQC->SetmaxEtaCut(5);
	fQC->SetdoQA(kTRUE);

	//fQC->SetEtaGapNeg(-0.1);
	//fQC->SetEtaGapPos(0.1);
	fQC->SetEtaGapNeg(0);
	fQC->SetEtaGapPos(0);
	
	Int_t nTotalEvent = 0;
	for (Int_t k = 0; k < nSplit; k++) {
		string directory;
		//string directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.44TeV/Centrality0-5/Baseline/job-%d/particle_distribution_final/%d.dat",ithJob,ithFile);
		if (centID == 0 || centID == 1) {
			directory = Form("/user/jlomker/project/AVFD/5.44TeV/a-0.1/tree_Baseline_5.44TeV_Cent%d_%d_%d.root", (centID)*5, (centID+1)*5, k);
		} else {
			directory = Form("/user/jlomker/project/AVFD/5.44TeV/a-0.1/tree_Baseline_5.44TeV_Cent%d_%d_%d.root", (centID-1)*10, (centID)*10, k);
		}
		TFile *f = TFile::Open(directory.c_str(), "READ");
		
		if (!f) {
			std::cout<<"Input file: "<<directory<<" does not exist"<<"\n";
			continue;
		} else {
			std::cout<<"Input file: "<<directory<<" exists"<<"\n";
		}
        
		TTree* tree = (TTree*)f->Get("events");

		// create a pointer to an event object. This will be used
		// to read the branch values.
		Event *event = new Event();

		// get the branch and set the branch address
		TBranch *bnevent = tree->GetBranch("event");
		bnevent->SetAddress(&event);

		Long64_t nevent = tree->GetEntries();
		cout<<"total event no. is "<<nevent<<"\n";

		for (Long64_t j = 0; j < nevent; j++) {
		//for (Long64_t j = 0; j < 1000; j++) {
			if (j%100 == 0) std::cout<<"Event: "<<j<<"\n";

			bnevent->GetEntry(j);
			
			if (event) {
				fQC->SetEvent(event);
				
				if (centID == 0 || centID == 1) { // fill 0-5% and 5-10% to the same bin as the histograms in CalculateFlowCME does not distinguish 0-5 and 5-10%
					fQC->SetCentralityEBE(gRandom->Uniform(0,10));
				} else {
					fQC->SetCentralityEBE(gRandom->Uniform((centID-1)*10,(centID)*10));
				}
				
				fQC->UserExec();
				nTotalEvent++;
			}
		}
		
		f->Close();
	}
  	
	cout<<"nTotalEvent ======== "<<nTotalEvent<<endl;
	fQC->Terminate();

	// Save list holding histogram with weights:
	TFile *ResultsFile;
	if (centID == 0 || centID == 1) {
		ResultsFile = new TFile(Form("result/AnalysisResults_a-0.1_5.44TeV_Cent%d_%d.root", (centID)*5, (centID+1)*5), "RECREATE");
	} else {
		ResultsFile = new TFile(Form("result/AnalysisResults_a-0.1_5.44TeV_Cent%d_%d.root", (centID-1)*10, (centID)*10), "RECREATE");
	}  

	ResultsFile->WriteObject(fQC->GetQAList(),"QAList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowQCList(),"FlowQCList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowGFList(),"FlowGFList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowFromBWList(),"FlowFromBWList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMEList(),"CMEList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMWList(),"CMWList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMWQAList(),"CMWQAList","SingleKey");
//	}	
}
