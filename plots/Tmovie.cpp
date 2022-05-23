#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void Tmovie(Int_t start, Int_t end) {//argument for time stamps ?
    const Int_t tauInit = 6;
    gROOT->SetBatch();
    gStyle->SetOptStat(0);
    TFile* f = new TFile(Form("/project/alice/users/jlomker/AVFD/plots/BField/Tevolution_5.02TeV_tauInitial0.%d_Blifetime_0.2.root",tauInit));
    TProfile2D* T1local2D;

    for (Int_t c = 2; c < 3; c++) {
        for (Int_t snap1 = start; snap1 < end; snap1++) {
	    double tau1 = 0.600000 + snap1 * 0.010000;
            T1local2D = (TProfile2D*)f->Get(Form("Tlocal2D_tau%d_cent%d0_%d0;1", snap1, c, c + 1));
            TCanvas* c1 = new TCanvas("c1", "Fields");//,200,10,700,500);
            T1local2D->SetTitle(Form("local temperature(#tau = %f)", tau1));
            T1local2D->Draw("COLZ");//->SetContour(10);//to group intensities more roughly
            c1->SaveAs(Form("BField/BField0.2/movieT/tau_initial0.%d/cent%d0_%d0/Tevolution_tauID%d.pdf",tauInit, c, c + 1, snap1));//double tau..
            T1local2D->Reset();
            delete c1;
        }
    }
}

