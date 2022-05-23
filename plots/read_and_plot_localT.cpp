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

void read_and_plot_localT() {
	const Int_t tauInit = 1;// bc ta0 = =0.tauInit
	const Int_t nCentBin = 8;
	const Int_t nTau = 300;//from 0.6 to 1 in dtau = 0.01
	const Double_t dTau = 0.01;
	Double_t CentBinEdge[nCentBin+1] = {0, 5, 10, 20, 30, 40, 50, 60, 70};
	Double_t Tau[nTau+1] = {0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1}; 
	TProfile2D *Tlocal2D[nCentBin][nTau];//make them all 2D and set centID = 3 to one value
	
	for (int tauID = 0; tauID <= nTau; tauID++){
		for (int centID = 0; centID < 8; centID++) {
			if (centID == 0 || centID == 1) {
				Tlocal2D[centID][tauID] = new TProfile2D(Form("Tlocal2D_tau%d_cent%d_%d", tauID,  (centID)*5, (centID+1)*5), Form("Tlocal2D__tau%d_cent%d_%d", tauID,  (centID)*5, (centID+1)*5), 46, -13, 13, 46, -13, 13);
			} else {
				Tlocal2D[centID][tauID] = new TProfile2D(Form("Tlocal2D_tau%d_cent%d_%d", tauID,  (centID-1)*10, (centID)*10), Form("Tlocal2D_tau%d_cent%d_%d", tauID,  (centID-1)*10, (centID)*10), 46, -13, 13, 46, -13, 13);
			}
		}
	}	

	TProfile *Tavge = new TProfile("Tavge", "Tavge", nCentBin, CentBinEdge);
	Int_t centID = 0; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70%
	int nJob = 1000;
	for (int centID = 0; centID < 8; centID++) {
		for (int ithJob = 0; ithJob <= nJob; ithJob++) {
			string directory;
			if (centID == 0 || centID == 1) {
				directory = Form("/dcache/alice/jlomker/sim/TestEM/tau_init_0.%d/BField0.2/Centrality%d_%d/job-%d/Result/event-1/VISHNUmovie_evolution.dat",tauInit, (centID)*5, (centID+1)*5, ithJob);
			} else {
				directory = Form("/dcache/alice/jlomker/sim/TestEM/tau_init_0.%d/BField0.2/Centrality%d_%d/job-%d/Result/event-1/VISHNUmovie_evolution.dat",tauInit, (centID-1)*10, (centID)*10, ithJob);
			}
			std::ifstream file_dat;

			file_dat.open(directory.c_str());

			if (!file_dat.is_open()) {
			  std::cout<<"Cannot find "<<directory<<std::endl;
			  continue;
			} else {
			  if (ithJob%100==0) std::cout<<"Processing job #"<<ithJob<<std::endl;
			}


			Double_t tau = -999, x = -999, y = -999, T = -999;

			while (file_dat >> tau>> x >> y >> T) {
		//cout<< tau << "\t" << x << "\t" << y << "\t" << T << endl;
				int tauID = (tau - 0.6)*100;
		//cout<<tauID<<centID<<endl;
				if(tauID <= nTau){
					Tlocal2D[centID][tauID]->Fill(x, y, T);
					if (centID == 0 || centID == 1) {
						Tavge->Fill(centID*5+2.5,T);
					} else {
						Tavge->Fill((centID-1)*10+5,T);
					}
				}
			}

			file_dat.close();
		}
	}

	TFile *file = new TFile(Form("BField/Tevolution_5.02TeV_tauInitial0.%d_Blifetime_0.2.root", tauInit), "RECREATE");
	for (int centID = 0; centID < 8; centID++) {
		for (int tauID = 0; tauID <= nTau; tauID++){
			Tlocal2D[centID][tauID]->Write();
		}
	}
	Tavge->Write();
	file->Close();
	}
