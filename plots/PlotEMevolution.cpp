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

    TProfile* posEx;
    TProfile* posEy;
    TProfile* posEz;
    TProfile* posBx;
    TProfile* posBy;
    TProfile* posBz;

    TProfile* negEx;
    TProfile* negEy;
    TProfile* negEz;
    TProfile* negBx;
    TProfile* negBy;
    TProfile* negBz;

    for (Int_t c = 1; c < 7; c++) {//one canvas vor every eta // canvas split in Bx,y upper and Ex Ey lower
        //for (Int_t etaID = 0; etaID < 3; etaID++) {
        Int_t etaID = 1;
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

	    negEx = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            negEy = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            negEz = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            posEx = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent%d0_%d0", 2, c, c + 1));
            posEy = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent%d0_%d0", 2, c, c + 1));
            posEz = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent%d0_%d0", 2, c, c + 1));
	    //For full plot
            auto C = new TLegend(0.62,0.13,0.89,0.29);
            C->AddEntry(Ex, "x", "l");
            C->AddEntry(Ey, "y", "l");
            C->AddEntry(Ez, "z", "l");
	    //For split plots
            auto cc = new TLegend(0.33,0.47,0.45,0.55);
            cc->AddEntry(Ex, "x", "l");
            cc->AddEntry(Ey, "y", "l");
            cc->AddEntry(Ez, "z", "l");

            Bx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            By = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            Bz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0", etaID, c, c + 1));
            cBx = (TProfile*)Bx->Clone("cBx");
            cBy = (TProfile*)By->Clone("cBy");
            cBz = (TProfile*)Bz->Clone("cBz");

            negBx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0",0, c, c + 1));
            negBy = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0",0, c, c + 1));
            negBz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0",0, c, c + 1));
            posBx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0",2, c, c + 1));
            posBy = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0",2, c, c + 1));
            posBz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0",2, c, c + 1));

            auto I = new TLegend(0.46, 0.3, 0.85, 0.55);
	    I->SetTextSize(0.038);
            I->SetHeader("AVFD Simulation","C");
	    I->SetEntrySeparation(0.006);
	    I->SetTextSize(0.034);
	    I->AddEntry((TObject*)0, "Pb-Pb @ #sqrt{s} = 5.02TeV","");
	    I->SetEntrySeparation(0.005);
	    I->AddEntry((TObject*)0, Form("Centrality %d0%%-%d0%%", c, c + 1), " ");
            I->SetEntrySeparation(0.005);
	    I->AddEntry((TObject*)0, "Initialized at #tau_{0} = 0.1 fm/c"," ");
            I->SetEntrySeparation(0.005);
	    I->AddEntry((TObject*)0, "With #tau_{B} = 0.2 fm/c", " ");
	    I->SetTextAlign(13);//I->SetEntrySeparation(0.0);//I->SetColumnSeparation(0);
	    I->SetBorderSize(0);
            Double_t eta = 0.0;

            auto i = new TLegend(0.58,0.46,0.79,0.59);
            i->SetTextSize(0.015);
            i->SetHeader("AVFD Simulation","C");
            //i->SetEntrySeparation(0.006);
            //i->SetTextSize(0.034);
            i->AddEntry((TObject*)0, "Pb-Pb @ #sqrt{s} = 5.02TeV","");
            //i->SetEntrySeparation(0.005);
            i->AddEntry((TObject*)0, Form("Centrality %d0%%-%d0%%", c, c + 1), " ");
            //i->SetEntrySeparation(0.005);
            i->AddEntry((TObject*)0, "Initialized at #tau_{0} = 0.1 fm/c"," ");
            //i->SetEntrySeparation(0.005);
            i->AddEntry((TObject*)0, "With #tau_{B} = 0.2 fm/c", " ");
            i->SetTextAlign(13);//I->SetEntrySeparation(0.0);//I->SetColumnSeparation(0);


            if (etaID == 0) { eta = -0.85; }
            if (etaID == 2) { eta = 0.85; };

	    //make E and B full separate for eta 0 and then cut for the 3 etas
	    auto Eta = new TLegend();//0, 0.8,0.4,1);
	    Eta->SetHeader(Form("#eta = %f",eta),"C");
            Eta->AddEntry(cBx,"Initialized at #tau_{0} = 0.1 fm/c"," ");
	    Eta->AddEntry(cBy,"#tau_{B} = 0.2 fm/c", " ");
	    
	    //canvas for full distribution
	    TCanvas* e1 = new TCanvas("e1", "EFields",600,250); //legend spacing pipapo...
	    e1->Divide(2,1);
	    e1->SetLeftMargin(4);
	    e1->SetRightMargin(0);
	    e1->cd(1);
	    e1->SetLeftMargin(4);
	    e1->SetRightMargin(0.0);
	    Ey->SetTitle("");
            Ey->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
	    //Ex->GetYaxis()->SetTitleSize(0.15);
	    //Ex->GetYaxis()->SetLabelSize(0.1);
	    Ey->GetYaxis()->SetTitleOffset(1.5);
            Ey->GetXaxis()->SetTitle("#tau [fm/c]");
	    Ey->GetXaxis()->SetRangeUser(0.0,3);
            Ey->SetLineColor(2);
            Ey->Draw("HIST L");//or Hist C
            Ex->SetLineColor(1);
	    Ex->Draw("HIST L SAME");
            Ez->SetLineColor(3);
            Ez->Draw("HIST L SAME");
            C->Draw();
	  
            e1->cd(2);
            e1->SetLeftMargin(4);
	    By->SetTitle("");
            By->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            By->GetXaxis()->SetTitle("#tau [fm/c]");
	    By->GetYaxis()->SetTitleOffset(1.5);
	    //By->GetYaxis()->SetLabelOffset(1.0);
	    By->GetXaxis()->SetRangeUser(0.0,3);
            By->SetLineColor(2);
            By->Draw("HIST L");
            Bx->SetLineColor(1);
            Bx->Draw("HIST L SAME");
            Bz->SetLineColor(3);
            Bz->Draw("HIST L SAME");
	    I->Draw();
            e1->SaveAs(Form("Evolution/EMFields_fullevolution_tauInit%d_cent%d0_%d0.png", tauInitial, c, c + 1));//double tau..
            Ex->Reset();
            Ey->Reset();
            Ez->Reset();
            Bx->Reset();
            By->Reset();
            Bz->Reset();
	    delete e1;

	    Double_t ymin = -0.033;
	    Double_t ymax = 0.033;
	    Double_t bmin = -0.0075;//the full -B peak should be in the first plot!
	    Double_t bmax = 0.0025;
	    //Canvas for cuts in different eta
            TCanvas* e2 = new TCanvas("e2", "EFields");//,200,10,700,500); //legend spacing pipapo...
            e2->SetLeftMargin(0.2);
	    e2->Divide(3,2, 1,0);
            //e2->SetLeftMargin(0.3);
	    e2->cd(1);
	    e2->SetBottomMargin(0);
            e2->SetLeftMargin(0.2);
	    e2->SetTopMargin(1);
            negEx->SetTitleOffset(4);
	    negEx->SetTitle("#eta_{s} = -0.85");
            negEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            negEx->GetXaxis()->SetTitle("#tau [fm/c]");
            negEx->GetYaxis()->SetTitleOffset(1.9);
	    negEx->GetYaxis()->SetRangeUser(ymin,ymax);
            negEx->SetLineColor(kBlue+1);
            negEx->GetXaxis()->SetRangeUser(0.0,0.8);
            negEx->Draw("HIST L");
            negEy->SetLineColor(kRed+1);
            negEy->Draw("HIST L SAME");
            negEz->SetLineColor(kGreen+1);
            negEz->Draw("HIST L SAME");
            
	    e2->cd();
	    e2->cd(2);
	    e2->SetBottomMargin(0);
            e2->SetLeftMargin(0.2);
	    e2->SetTopMargin(1);
            cEx->SetTitle("#eta_{s} = 0");
	    cEx->SetTitleOffset(4);
            cEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            cEx->GetXaxis()->SetTitle("#tau [fm/c]");
            cEx->GetYaxis()->SetTitleOffset(1.9);
	    cEx->GetYaxis()->SetRangeUser(ymin,ymax);
            cEx->SetLineColor(1);
            cEx->GetXaxis()->SetRangeUser(0.0,0.8);
            cEx->Draw("Hist L");
            cEy->SetLineColor(2);
            cEy->Draw("HIST L SAME");
            cEz->SetLineColor(3);
            cEz->Draw("HIST L SAME");
	    e2->cd();
            e2->cd(3);
	    e2->SetBottomMargin(0);
            e2->SetLeftMargin(0.2);
	    e2->SetTopMargin(1);
            posEx->SetTitle("#eta_{s} = 0.85");
	    posEx->SetTitleOffset(4);
            posEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            posEx->GetXaxis()->SetTitle("#tau [fm/c]");
            posEx->GetYaxis()->SetTitleOffset(1.9);
	    posEx->GetYaxis()->SetRangeUser(ymin,ymax);
            posEx->SetLineColor(1);
            posEx->GetXaxis()->SetRangeUser(0.0,0.8);
            posEx->Draw("HIST L");
            posEy->SetLineColor(2);
            posEy->Draw("HIST L SAME");
            posEz->SetLineColor(3);
            posEz->Draw("HIST L SAME");
	    
	    e2->cd();
            e2->cd(4);
	    e2->SetTopMargin(0);
            e2->SetLeftMargin(0.2);
            negBx->SetTitle("");
            negBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            negBx->GetYaxis()->SetTitleOffset(1.9);
            negBx->GetXaxis()->SetTitle("#tau [fm/c]");
	    negBx->GetYaxis()->SetRangeUser(bmin,bmax);
            negBx->SetLineColor(1);
            negBx->GetXaxis()->SetRangeUser(0.0,0.8);
            negBx->Draw("HIST L");
            negBy->SetLineColor(2);
            negBy->Draw("HIST L SAME");
            negBz->SetLineColor(3);
            negBz->Draw("HIST L SAME");
            e2->cd();
	    cc->Draw();
	    e2->cd();
	    e2->cd(5);
            e2->SetLeftMargin(0.2);
	    e2->SetTopMargin(0);
            cBx->SetTitle("");
            cBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            cBx->GetYaxis()->SetTitleOffset(1.9);
            cBx->GetXaxis()->SetTitle("#tau [fm/c]");
	    cBx->GetYaxis()->SetRangeUser(bmin,bmax);
            cBx->SetLineColor(1);
            cBx->GetXaxis()->SetRangeUser(0.0,0.8);
            cBx->Draw("HIST L");
            cBy->SetLineColor(2);
            cBy->Draw("HIST L SAME");
            cBz->SetLineColor(3);
            cBz->Draw("HIST L SAME");
            e2->cd();
	    e2->cd(6);
            gPad->SetTickx(1);
	    e2->SetLeftMargin(0.2);
	    e2->SetTopMargin(0);
            posBx->SetTitle(" ");
            posBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            posBx->GetYaxis()->SetTitleOffset(1.9);
            posBx->GetXaxis()->SetTitle("#tau [fm/c]");
	    posBx->GetYaxis()->SetRangeUser(bmin,bmax);
            posBx->SetLineColor(1);
            posBx->GetXaxis()->SetRangeUser(0.0,0.8);
            posBx->Draw("HIST L");
            posBy->SetLineColor(2);
            posBy->Draw("HIST L SAME");
            posBz->SetLineColor(3);
            posBz->Draw("HIST L SAME");
	    e2->cd();
	    i->Draw();
	    e2->SaveAs(Form("Evolution/EMFields_Splitsevolution_tauInit%d_cent%d0_%d0.pdf", tauInitial, c, c + 1));//double tau..
            cEx->Reset();
            cEy->Reset();
            cEz->Reset();
            posEx->Reset();
            posEy->Reset();
            posEz->Reset();
            negEx->Reset();
            negEy->Reset();
            negEz->Reset();
	    cBx->Reset();
            cBy->Reset();
            cBz->Reset();
            posBx->Reset();
            posBy->Reset();
            posBz->Reset();
            negBx->Reset();
            negBy->Reset();
            negBz->Reset();
 	    delete e2;
    //    }
    }
}
