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

void runSingleCentrality(Int_t centID, Int_t dirID, Int_t pT){
//for(Int_t centID = 0; centID < 8; centID++){ 
	std::cout << "Centrality ID: " << centID << "Directory ID: " << dirID <<std::endl;
        //Int_t centID = {0,1,2,3,4,5,6,7}; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70
        //DirID = 0: baseline, 1: a-0.1
        int val1 = (centID-1)*10;
        int val2 = (centID)*10;
        if (centID == 0 || centID == 1) {
        val1 = (centID)*5;
        val2 = (centID+1)*5;
        }
	Double_t pTmin = 0.2;//to pick the full pT range with random integer
	Double_t pTmax = 5.0;
	if(pT == 1){//for small and low pT range
	pTmin = 0.2;
	pTmax = 1;
	}
	if(pT == 2){//for small high pT range
	pTmin = 3.0;
	pTmax = 5.0;
	}
	//should be according to the number of files per centrality
	const Int_t nSplit = 10;
	//TFile *f[nCentBin][nSplit];
	CalculateFlowCME *fQC = new CalculateFlowCME("CalculateFlowCME");
	fQC->UserCreateOutputObjects();
	//fQC->SetmaxPtCut(3);
	//fQC->SetminPtCut(0.2);
	//fQC->SetminNtrackCut(500);
	//fQC->SetmaxEtaCut(0.8);
	fQC->SetmaxPtCut(pTmax);
	fQC->SetminPtCut(pTmin);
	fQC->SetminNtrackCut(500);
	fQC->SetmaxEtaCut(0.8);
	fQC->SetdoQA(kTRUE);

	//fQC->SetEtaGapNeg(-0.1);
	//fQC->SetEtaGapPos(0.1);
	fQC->SetEtaGapNeg(0);
	fQC->SetEtaGapPos(0);
	
	Int_t nTotalEvent = 0;
	for (Int_t k = 0; k < nSplit; k++) {
		string directory;
	//here comes the input directory	
	//string directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.44TeV/Centrality0-5/Baseline/job-%d/particle_distribution_final/%d.dat",ithJob,ithFile);
		//directory = Form("tree_5.44TeV_Cent%d_%d_%d.root", val1, val2, k);//submit script puts me into dirID/...
		directory = Form("/data/alice/jlomker/AVFD/Centrality-%d/dirID-0/tree_5.44TeV_Cent%d_%d_%d.root",centID, val1,val2,k);
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
	// Here comes the output directory
	ResultsFile = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-%d/full/Results_5.44TeV_pTrange_%d_Cent%d_%d.root",dirID, pT, val1, val2), "RECREATE");	  

	ResultsFile->WriteObject(fQC->GetQAList(),"QAList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowQCList(),"FlowQCList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowGFList(),"FlowGFList","SingleKey");
	ResultsFile->WriteObject(fQC->GetFlowFromBWList(),"FlowFromBWList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMEList(),"CMEList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMWList(),"CMWList","SingleKey");
	ResultsFile->WriteObject(fQC->GetCMWQAList(),"CMWQAList","SingleKey");
//	}	
}
