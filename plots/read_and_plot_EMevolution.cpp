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


void read_and_plot_EMevolution() {
	const Int_t nTau = 300;
	const Int_t nCentBin = 8;
	const Int_t nEta = 3; //0: n_s < 0; 1: n_s = 0; 2: n_s > 0
	Double_t CentBinEdge[nCentBin+1] = {0, 5, 10, 20, 30, 40, 50, 60, 70};

	TProfile2D *feBxfield2D[nCentBin][nEta][nTau];//add calculation nTau
	TProfile2D *feByfield2D[nCentBin][nEta][nTau];//add definition nEta
	TProfile2D *feExfield2D[nCentBin][nEta][nTau];
	TProfile2D *feEyfield2D[nCentBin][nEta][nTau];

	TProfile *lrf_Bx[nCentBin][nEta];//add TProfiles for B_xyz/E_xyz(t) for cent and eta
	TProfile *lrf_By[nCentBin][nEta];
	TProfile *lrf_Bz[nCentBin][nEta];
	 
	TProfile *lrf_Ex[nCentBin][nEta];
	TProfile *lrf_Ey[nCentBin][nEta];
	TProfile *lrf_Ez[nCentBin][nEta];

	for (int centID = 0; centID < 8; centID++) {
		for (int etaID = 0; etaID < 4; etaID++) {
				for (int tauID = 0; tauID < nTau; tauID++) {
					feBxfield2D[centID][etaID][tauID] = new TProfile2D(Form("feBxfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), Form("feBxfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), 100, -13, 13, 100, -13, 13);
					feByfield2D[centID][etaID][tauID] = new TProfile2D(Form("feByfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), Form("feByfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), 100, -13, 13, 100, -13, 13);
					feExfield2D[centID][etaID][tauID] = new TProfile2D(Form("feExfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), Form("feExfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), 100, -13, 13, 100, -13, 13);
					feEyfield2D[centID][etaID][tauID] = new TProfile2D(Form("feEyfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), Form("feEyfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, centID, centID + 1), 100, -13, 13, 100, -13, 13);
				}
				TProfile* lrf_Bx = new TProfile(Form("lrf_Bxfield_eta%d_cent%d0_%d0", etaID, centID, centID +1 ), Form("lrf_Bxfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);
				TProfile* lrf_By = new TProfile(Form("lrf_Byfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), Form("lrf_Byfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);
				TProfile* lrf_Bz = new TProfile(Form("lrf_Bzfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), Form("lrf_Bzfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);

				TProfile* lrf_Ex = new TProfile(Form("lrf_Exfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), Form("lrf_Exfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);
				TProfile* lrf_Ey = new TProfile(Form("lrf_Eyfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), Form("lrf_Eyfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);
				TProfile* lrf_Ez = new TProfile(Form("lrf_Ezfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), Form("lrf_Ezfield_eta%d_cent%d0_%d0", etaID, centID, centID + 1), nTau, 0, 3);
		}
	}
	

	//Int_t centID = 0; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70%
	int nJob = 1000;
	double tauInitial = 0.6; 
	for (int centID = 0; centID < 8; centID++) {
		for (int ithJob = 0; ithJob <= nJob; ithJob++) {
			string directory;
			directory = Form("/dcache/alice/jlomker/sim/TestEM/tau_init_0.6/BField0.2/Centrality%d_%d/job-%d/Result/event-1/check_lrf_EMfields.dat", (centID-1)*10, (centID)*10, ithJob);
			std::ifstream file_dat;

			file_dat.open(directory.c_str());
			if (!file_dat.is_open()) {
			  std::cout<<"Cannot find "<<directory<<std::endl;
			  continue;
			} else {
			  if (ithJob%100==0) std::cout<<"Processing job #"<<ithJob<<std::endl;
			}
			TFile *skip0;
			Double_t tau = -999, xPos = -999, yPos = -999, eta = -999, eEx = -999, eEy = -999, eEz = -999, eBx = -999, eBy = -999, eBz = -999;//E and B in [1/fm^2]
			//cout << file_dat <<endl;
			//skip0 = file_dat.seekg(1,ios::beg);
			//fgets(line, 50, input) != NULL
			file_dat.ignore(1);
			while (file_dat >> tau >> xPos >> yPos >> eta >> eEx >> eEy >> eEz >> eBx >> eBy >> eBz){
		//	if(file_dat.eof()){break;}
			cout<< xPos << "\t" << yPos << "\t" << eBx<< "\t" << eBy <<"\t" <<"\t"<< eEx <<"\t"<< eEy  <<"\t"<< endl;
				/*tau  x  y  eta  Ex[1/fm^2]  Ey[1/fm^2]  Ez[1/fm^2]  Bx[1/fm^2]  By[1/fm^2]  Bz[1/fm^2]
 				   0   0  0   0 =-ns		0	   0		0	   0		0	   0
 				   0   0  0   1 = 0		0	   0		0 	   0		0	   0
    				   0   0  0   2 = ns		0	   0 		0          0		0	   0
 				dynamically calculate the tauID based on data in dTau = 0.1 steps ...allocate ID's up to tau = 2 
				*/
				int etaID = 1;//set to eta_ns = 0
				int tauID = (tau - tauInitial) * 100; //dTau is then 0.01 .. not sure if this is ideal though...maybe 0.02 steps are better
				if(eta < 0){etaID = 0;}
				if(eta > 0){etaID = 2;}
				
				feBxfield2D[centID][etaID][tauID]->Fill(xPos, yPos, eBx);
				feByfield2D[centID][etaID][tauID]->Fill(xPos, yPos, eBy);
				feExfield2D[centID][etaID][tauID]->Fill(xPos, yPos, eEx);
				feEyfield2D[centID][etaID][tauID]->Fill(xPos, yPos, eEy);

				lrf_Bx[centID][etaID]->Fill(tau, eBx);//add TProfiles for B_xyz/E_xyz(t) for cent and eta
				lrf_By[centID][etaID]->Fill(tau, eBy);
				lrf_Bz[centID][etaID]->Fill(tau, eBz);
				lrf_Ex[centID][etaID]->Fill(tau, eEx);
				lrf_Ey[centID][etaID]->Fill(tau, eEy);
				lrf_Ez[centID][etaID]->Fill(tau, eEz);			
			}
			file_dat.close();
		}
	}

	TFile *file = new TFile("BField/EM_evolution_5.02TeV_initial0.6_lifetime_0.2.root", "RECREATE");
	for (int centID = 0; centID < 8; centID++) {
		for (int etaID = 0; etaID < 4; etaID++) {
			for (int tauID = 0; tauID < nTau; tauID++) {
				feBxfield2D[centID][etaID][tauID]->Write();
				feByfield2D[centID][etaID][tauID]->Write();
				feExfield2D[centID][etaID][tauID]->Write();
				feEyfield2D[centID][etaID][tauID]->Write();
			}
			lrf_Bx[centID][etaID]->Write();
			lrf_By[centID][etaID]->Write();
			lrf_Bz[centID][etaID]->Write();

			lrf_Ex[centID][etaID]->Write();
			lrf_Ey[centID][etaID]->Write();
			lrf_Ez[centID][etaID]->Write();

		}
	}
	file->Close();
	
}
