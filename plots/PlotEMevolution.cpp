#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotEMevolution(Int_t start, Int_t end) {//argument for time stamps ?
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	const int tauInitial = 6; 
	TFile* f = new TFile(Form("/project/alice/users/jlomker/AVFD/plots/BField/EM_evolution_5.02TeV_initial0.%d_lifetime_0.2.root",tauInitial));
	TProfile2D* Ex2D;
	TProfile2D* Ey2D;
	TProfile2D* Bx2D;
	TProfile2D* By2D;

	TProfile* Ex, Ey, Ez;
	TProfile* Bx, By, Bz;

	for (Int_t c = 2; c < 3; c++) {//one canvas vor every eta // canvas split in Bx,y upper and Ex Ey lower  
		for (Int_t etaID = 0; etaID < 3; etaID++){
			for (Int_t tauID = start; tauID <= end; tauID++) {
				Float_t tau = tauInitial/10 + tauID * 0.010000;//snap0 = 0.6 -> snapi = 0.6 + i*0.01
				Float_t eta = 0.0;
				if (etaID == 0) { eta = -0.85;}
				if (etaID == 2) { eta = 0.85;};
				Ex2D = (TProfile2D*)f->Get(Form("feExfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
				Ey2D = (TProfile2D*)f->Get(Form("feEyfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
				Bx2D = (TProfile2D*)f->Get(Form("feBxfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
				By2D = (TProfile2D*)f->Get(Form("feByfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));

				TCanvas* c1 = new TCanvas("c1", "Fields");//,200,10,700,500);
				c1->SetTitle(Form("For #tau_{0} = 0.%d with #tau(B) = 0.2 at #eta_{s} = %f and #tau = %f ",tauInitial, eta, tau));
				c1->Divide(2, 2);
				c1->cd(1);//gStyle->SetPalette(1);//to set colors ... very intense with 1
				Bx2D->SetTitle(Form("B_{x}"));// (#tau = % f) for #eta_{ s } = % f", tau, eta));
				Bx2D->Draw("COLZ");//->SetContour(10);//to group intensities more roughly
				c1->cd(2);				
				By2D->SetTitle(Form("B_{y}"));// (#tau = % f) for #eta_{ s } = % f", tau, eta));
				By2D->Draw("COLZ");
				c1->cd(3);
				Ex2D->SetTitle(Form("E_{x}"));// (#tau = % f) for #eta_{ s } = % f", tau, eta));
				Ex2D->Draw("COLZ");
				c1->cd(4);
				Ey2D->SetTitle(Form("E_{y}"));// (#tau = % f) for #eta_{ s } = % f", tau, eta));
				Ey2D->Draw("COLZ");
				c1->SaveAs(Form("BField/BField0.2/EMevolution_tauInit%d_tauID%d_eta%d_cent%d0_%d0.pdf", tauInitial, tauID, etaID, c, c + 1));//double tau..
				Bx2D->Reset();
				By2D->Reset();
				Ex2D->Reset();
				Ey2D->Reset();
				delete c1;
			}//here for the 1D plots ...
		}
	}
}
