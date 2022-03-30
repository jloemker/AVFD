#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotEMevolution(Int_t start, Int_t end) {//add legends and make shape better 
    gROOT->SetBatch();
    gStyle->SetOptStat(0);
    const int tauInitial = 1;//becomes 0.tauInitial
    TFile* f = new TFile(Form("/project/alice/users/jlomker/AVFD/plots/BField/EM_evolution_5.02TeV_initial0.%d_lifetime_0.2.root", tauInitial));
    TProfile2D* Ex2D;
    TProfile2D* Ey2D;
    TProfile2D* Bx2D;
    TProfile2D* By2D;

    TProfile* Ex;
    TProfile* Ey;
    TProfile* Ez;
    TProfile* Bx;
    TProfile* By;
    TProfile* Bz;

    TProfile* cEx;
    TProfile* cEy;
    TProfile* cEz;
    TProfile* cBx;
    TProfile* cBy;
    TProfile* cBz;

    for (Int_t c = 2; c < 3; c++) {//one canvas vor every eta // canvas split in Bx,y upper and Ex Ey lower
        for (Int_t etaID = 0; etaID < 3; etaID++) {
            for (Int_t tauID = start; tauID <= end; tauID++) {
                Float_t tau = tauInitial / 10 + tauID * 0.010000;//snap0 = 0.6 -> snapi = 0.6 + i*0.01
                Float_t eta = 0.0;
                if (etaID == 0) { eta = -0.85; }
                if (etaID == 2) { eta = 0.85; };
                Ex2D = (TProfile2D*)f->Get(Form("feExfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
                Ey2D = (TProfile2D*)f->Get(Form("feEyfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
                Bx2D = (TProfile2D*)f->Get(Form("feBxfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));
                By2D = (TProfile2D*)f->Get(Form("feByfield2D_eta%d_tau%d_cent%d0_%d0", etaID, tauID, c, c + 1));

                TCanvas* c1 = new TCanvas("c1", "Fields");//,200,10,700,500);
                c1->Divide(2, 2);
                c1->SetTitle(Form("For #tau_{0} = 0.%d with #tau(B) = 0.2 at #eta_{s} = %f and #tau = %f ", tauInitial, eta, tau));
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
                c1->SaveAs(Form("BField/BField0.2/movieEM/EMevolution_tauInit%d_tauID%d_eta%d_cent%d0_%d0.pdf", tauInitial, tauID, etaID, c, c + 1));//double tau..
                Bx2D->Reset();
                By2D->Reset();
                Ex2D->Reset();
                Ey2D->Reset();
                delete c1;
            }//here for the 1D plots ... Ex,y,z for every eta in one plot same for Bx,y,z || Form("lrf_Bxfield_eta%d_cent%d0_%d0", etaID, centID, centID +1 ),
            
            Ex = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            Ey = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            Ez = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
	    cEx = (TProfile*)Ex->Clone("cEx");
	    cEy = (TProfile*)Ey->Clone("cEy");
            cEz = (TProfile*)Ez->Clone("cEz");

            auto E1 = new TLegend(0.5, 0.6, 0.9, 0.9);
            E1->SetHeader(Form("Pb-Pb, 5.02TeV, Cent. %d0-%d0", c, c + 1), "C");
            E1->AddEntry(Ex, "E_{x}", "l");
            E1->AddEntry(Ey, "E_{y}", "l");
            E1->AddEntry(Ez, "E_{z}", "l");

            Bx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            By = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            Bz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            cBx = (TProfile*)Bx->Clone("BBx");
            cBy = (TProfile*)By->Clone("cBy");
            cBz = (TProfile*)Bz->Clone("cBz");

            auto B1 = new TLegend(0.5, 0.1, 0.9, 0.4);
            B1->SetHeader(Form("Pb-Pb, 5.02TeV, Cent. %d0-%d0", c, c + 1), "C");
            B1->AddEntry(Bx, "B_{x}", "l");
            B1->AddEntry(By, "B_{y}", "l");
            B1->AddEntry(Bz, "B_{z}", "l");

            Double_t eta = 0.0;
            if (etaID == 0) { eta = -0.85; }
            if (etaID == 2) { eta = 0.85; };

	    auto Eta = new TLegend();//0, 0.8,0.4,1);
	    Eta->SetHeader(Form("#eta = %f",eta),"C");
            Eta->AddEntry(cBx,"Initialized at 0.6 fm/c"," ");
	    Eta->AddEntry(cBy,"B field lifetime 0.2 fm/c", " ");
	    TCanvas* e1 = new TCanvas("e1", "EFields");//,200,10,700,500); //legend spacing pipapo...
            e1->Divide(2,2);
	    e1->SetTopMargin(0.2);
	    e1->SetTitle(Form("For #tau_{0} = 0.%d with #tau(B) = 0.2 in centrality %d0 - %d0", tauInitial,c,c+1));
	    e1->cd(1);
	    e1->SetLeftMargin(0.2);
	    Ex->SetTitle("");
            Ex->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            Ex->GetXaxis()->SetTitle("#tau [fm/c]");
	    Ex->GetYaxis()->SetTitleOffset(1.5);
            Ex->SetLineColor(1);
            Ex->Draw();
            Ey->SetLineColor(2);
            Ey->Draw("SAME");
            Ez->SetLineColor(3);
            Ez->Draw("SAME");
            E1->Draw();
            e1->cd(2);
            e1->SetLeftMargin(0.2);
	    cEx->SetTitle("");
            cEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            cEx->GetXaxis()->SetTitle("#tau [fm/c]");
	    cEx->GetYaxis()->SetTitleOffset(1.5);
            cEx->SetLineColor(1);
            cEx->GetXaxis()->SetRangeUser(0.5,1);
            cEx->Draw();
            cEy->SetLineColor(2);
            cEy->Draw("SAME");
            cEz->SetLineColor(3);
            cEz->Draw("SAME");
	    Eta->Draw();
            e1->cd(3);
            e1->SetLeftMargin(0.2);
	    By->SetTitle("");
            By->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            By->GetXaxis()->SetTitle("#tau [fm/c]");
	    By->GetYaxis()->SetTitleOffset(1.5);
            By->SetLineColor(2);
            By->Draw();
            Bx->SetLineColor(1);
            Bx->Draw("SAME");
            Bz->SetLineColor(3);
            Bz->Draw("SAME");
            B1->Draw();
            e1->cd(4);
            e1->SetLeftMargin(0.2);
	    cBx->SetTitle("");
            cBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
	    cBx->GetYaxis()->SetTitleOffset(1.5);
            cBx->GetXaxis()->SetTitle("#tau [fm/c]");
            cBx->SetLineColor(1);
            cBx->GetXaxis()->SetRangeUser(0.5,1);
            cBx->Draw();
            cBy->SetLineColor(2);
            cBy->Draw("SAME");
            cBz->SetLineColor(3);
            cBz->Draw("SAME");
	    Eta->Draw();
            e1->SaveAs(Form("BField/BField0.2/EMFields_evolution_tauInit%d_eta%d_cent%d0_%d0.pdf", tauInitial,etaID, c, c + 1));//double tau..
            Ex->Reset();
            Ey->Reset();
            Ez->Reset();
            Bx->Reset();
            By->Reset();
            Bz->Reset();

            cEx->Reset();
            cEy->Reset();
            cEz->Reset();
            cBx->Reset();
            cBy->Reset();
            cBz->Reset();
            delete e1;
        }
    }
}

