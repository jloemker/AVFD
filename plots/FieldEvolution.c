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

//********************************************************
//********************************************************
//  Trying to create snapshots of the EM - Field evolution
//  - Assumes Global tau dependence
//*******************************************************
//*******************************************************

//general initializations in main and arguments of the funtions

//reading input file
void read_initial_fields() {//maybe also the derivative.. or calculate derivative for a proper evolution 
 int nJob = 1000;
        for (int centID = 2; centID < 4; centID++) {//here I have to specify my input !
                for (int ithJob = 0; ithJob <= nJob; ithJob++) {
                        string directory;
                        if (centID == 0 || centID == 1) {
                                directory = Form("/dcache/alice/jlomker/sim/NoBField/5.02TeV/Centrality%d_%d/job-%d/Result/event-1/EMField.dat", (centID)*5, (centID+1)*5, ithJob);
 			} else {
                                directory = Form("/dcache/alice/jlomker/sim/NoBField/5.02TeV/Centrality%d_%d/job-%d/Result/event-1/EMField.dat", (centID-1)*10, (centID)*10, ithJob);
			}
 			std::ifstream file_dat;
                        file_dat.open(directory.c_str());
			if (!file_dat.is_open()) {
				std::cout<<"Cannot find "<<directory<<std::endl;
			  	continue;
			} else {
			  if (ithJob%100==0) std::cout<<"Processing job #"<<ithJob<<std::endl;
			}
			Double_t xPos = -999, yPos = -999, eBx = -999, eBy = -999, eEx = -999, eEy = -999;
			Double_t emptyloc = -1, emptyloc2 = -1;

			while (file_dat >> xPos >> yPos >> eBx >> eBy >> emptyloc >> eEx >> eEy >> emptyloc2) {
				/*feBxfield2D[centID]->Fill(xPos, yPos, eBx);//we need this for extra dim [tau] -> don't think this is the way to go...bc we only reduce the filed stepwise, without updaing the positions... we indeed hae to add the derivative ... yes that will take a bit ... 
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
				}*/
			}//end of while
			file_dat.close();
		}//end of ithJob
	}//end of centrality loop
}//end read initial

//writing output file that contains the evolution b(tau)
void write_evolution(double_t tau){
 TFile *file = new TFile("BField/Field_evolution_.root", "RECREATE");
 for (int centID = 0; centID < 8; centID++) {
	feBxfield2D[centID][tauID]->Write();
	feByfield2D[centID][tauID]->Write();
	feExfield2D[centID][tauID]->Write();
	feEyfield2D[centID][tauID]->Write();
 }
	favgeBxfield->Write();
	favgeByfield->Write();
	favgeExfield->Write();
	favgeEyfield->Write();
	file->Close();
}//end write


//calculation tools
double B_Tau(double tau)
{//tauB is lifetime of the field (?) -> old plots 0.2fm/c, new plots 0 fm/c
 return TauB*TauB/(tau*tau+TauB*TauB);
}

double dB_dTau(double tau)
{
 return -2*tau*TauB*TauB/(tau*tau+TauB*TauB)/(tau*tau+TauB*TauB);
}


//******************************************************
//
//			MAIN
//
//******************************************************		

void FieldEvolution(){
	//EM_B.. from file input read em files

	//from evolution code with idx = ieta*nx*ny + ix*ny + iy;
	//every point was stored in the .dat file ->use Shi's script for every tau and calc B_at_tau ... that easy ?!
	double B_at_Tau = B_Tau(tau);
	double dBdTau_at_Tau = dB_dTau(tau);
	double map_Bx = EM_Bx[idx] * B_at_Tau;
	double map_By = EM_By[idx] * B_at_Tau;
	double map_Bz = EM_Bz[idx] * B_at_Tau;
	double map_Ex = EM_Ex[idx] * B_at_Tau;
	double map_Ey = EM_Ey[idx] * B_at_Tau;
	double map_Ez = EM_Ez[idx] * B_at_Tau;

}
