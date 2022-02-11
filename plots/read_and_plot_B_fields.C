#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TProfile2D.h"
using namespace std;

void read_and_plot_B_fields() {

	const Int_t nCentBin = 8;
	Double_t CentBinEdge[nCentBin+1] = {0, 5, 10, 20, 30, 40, 50, 60, 70};

	TProfile2D *feBxfield2D[nCentBin];
	TProfile2D *feByfield2D[nCentBin];
	TProfile2D *feExfield2D[nCentBin];
	TProfile2D *feEyfield2D[nCentBin];
	
	for (int centID = 0; centID < 8; centID++) {
		if (centID == 0 || centID == 1) {
			feBxfield2D[centID] = new TProfile2D(Form("feBxfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), Form("feBxfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), 260, -13, 13, 260, -13, 13);
			feByfield2D[centID] = new TProfile2D(Form("feByfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), Form("feByfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), 260, -13, 13, 260, -13, 13);
			feExfield2D[centID] = new TProfile2D(Form("feExfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), Form("feExfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), 260, -13, 13, 260, -13, 13);
			feEyfield2D[centID] = new TProfile2D(Form("feEyfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), Form("feEyfield2D_cent%d_%d",  (centID)*5, (centID+1)*5), 260, -13, 13, 260, -13, 13);
		} else {
			feBxfield2D[centID] = new TProfile2D(Form("feBxfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), Form("feBxfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), 260, -13, 13, 260, -13, 13);
			feByfield2D[centID] = new TProfile2D(Form("feByfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), Form("feByfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), 260, -13, 13, 260, -13, 13);
			feExfield2D[centID] = new TProfile2D(Form("feExfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), Form("feExfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), 260, -13, 13, 260, -13, 13);
			feEyfield2D[centID] = new TProfile2D(Form("feEyfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), Form("feEyfield2D_cent%d_%d",  (centID-1)*10, (centID)*10), 260, -13, 13, 260, -13, 13);
		}
	}
		

	TProfile *favgeBxfield = new TProfile("favgeBxfield", "favgeBxfield", nCentBin, CentBinEdge);
	TProfile *favgeByfield = new TProfile("favgeByfield", "favgeByfield", nCentBin, CentBinEdge);
	TProfile *favgeExfield = new TProfile("favgeExfield", "favgeExfield", nCentBin, CentBinEdge);
	TProfile *favgeEyfield = new TProfile("favgeEyfield", "favgeEyfield", nCentBin, CentBinEdge);

	Int_t centID = 0; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70%

	int nJob = 2000;
	for (int centID = 0; centID < 8; centID++) {
		for (int ithJob = 0; ithJob <= nJob; ithJob++) {
			string directory;
			if (centID == 0 || centID == 1) {
				directory = Form("/dcache/alice/jlomker/sim/NoBField/5.02TeV/Centrality%d_%d/job-%d/Result/event-1/EMField.dat", (centID)*5, (centID+1)*5, ithJob);
				//directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.02TeV/Centrality%d-%d/BaselineForBField/job-%d/Result/event-1/EMField.dat", (centID)*5, (centID+1)*5, ithJob);
				//directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.44TeV/Centrality%d-%d/BaselineForBField/job-%d/Result/event-1/EMField.dat", (centID)*5, (centID+1)*5, ithJob);
			} else {
				directory = Form("/dcache/alice/jlomker/sim/NoBField/5.02TeV/Centrality%d_%d/job-%d/Result/event-1/EMField.dat", (centID-1)*10, (centID)*10, ithJob);
				//directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.02TeV/Centrality%d-%d/BaselineForBField/job-%d/Result/event-1/EMField.dat", (centID-1)*10, (centID)*10, ithJob);
				//directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.44TeV/Centrality%d-%d/BaselineForBField/job-%d/Result/event-1/EMField.dat", (centID-1)*10, (centID)*10, ithJob);
			}
			std::ifstream file_dat;

			file_dat.open(directory.c_str());
			//file_dat.open("/home/alidock/simulation/2020_AVFD_2.76TeV/EMField.dat");

			if (!file_dat.is_open()) {
			  std::cout<<"Cannot find "<<directory<<std::endl;
			  continue;
			} else {
			  if (ithJob%100==0) std::cout<<"Processing job #"<<ithJob<<std::endl;
			}


			Double_t xPos = -999, yPos = -999, eBx = -999, eBy = -999, eEx = -999, eEy = -999;
			Double_t emptyloc = -1, emptyloc2 = -1;

			while (file_dat >> xPos >> yPos >> eBx >> eBy >> emptyloc >> eEx >> eEy >> emptyloc2) {
				//cout<< xPos << "\t" << yPos << "\t" << eBx<< "\t" << eBy <<"\t"<< emptyloc <<"\t"<< eEx <<"\t"<< eEy  <<"\t"<< emptyloc2 << endl;
				
				feBxfield2D[centID]->Fill(xPos, yPos, eBx);
				feByfield2D[centID]->Fill(xPos, yPos, eBy);
				feExfield2D[centID]->Fill(xPos, yPos, eEx);
				feEyfield2D[centID]->Fill(xPos, yPos, eEy);
				if (centID == 0 || centID == 1) {
					favgeBxfield->Fill(centID*5+2.5,eBx);
					favgeByfield->Fill(centID*5+2.5,eBy);
					favgeExfield->Fill(centID*5+2.5,eEx);
					favgeEyfield->Fill(centID*5+2.5,eEy);
				} else {
					favgeBxfield->Fill((centID-1)*10+5,eBx);
					favgeByfield->Fill((centID-1)*10+5,eBy);
					favgeExfield->Fill((centID-1)*10+5,eEx);
					favgeEyfield->Fill((centID-1)*10+5,eEy);
				}
			}

			file_dat.close();
		}
	}

	TFile *file = new TFile("BField/B_fields_5.02TeV_lifetime_0.root", "RECREATE");
	for (int centID = 0; centID < 8; centID++) {
		feBxfield2D[centID]->Write();
		feByfield2D[centID]->Write();
		feExfield2D[centID]->Write();
		feEyfield2D[centID]->Write();
	}
	favgeBxfield->Write();
	favgeByfield->Write();
	favgeExfield->Write();
	favgeEyfield->Write();
	file->Close();

}
