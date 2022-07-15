#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void ThesisEMevolution() {//add legends and make shape better 
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
    TProfile* Ex2;
    TProfile* Ey2;
    TProfile* Ez2;
    TProfile* Bx2;
    TProfile* By2;
    TProfile* Bz2;
    TProfile* cEx2;
    TProfile* cEy2;
    TProfile* cEz2;
    TProfile* cBx2;
    TProfile* cBy2;
    TProfile* cBz2;

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

    for (Int_t c = 3; c < 4; c++) {//one canvas vor every eta // canvas split in Bx,y upper and Ex Ey lower
        //for (Int_t etaID = 0; etaID < 3; etaID++) {
        Int_t etaID = 1;
            
            Ex = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent30_40", etaID));
            Ey = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent30_40", etaID));
            Ez = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent30_40", etaID));
	    cEx = (TProfile*)Ex->Clone("cEx");
	    cEy = (TProfile*)Ey->Clone("cEy");
            cEz = (TProfile*)Ez->Clone("cEz");

            Ex2 = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent60_70", etaID));
            Ey2 = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent60_70", etaID));
            Ez2 = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent60_70", etaID));
            cEx2 = (TProfile*)Ex2->Clone("cEx2");
            cEy2 = (TProfile*)Ey2->Clone("cEy2");
            cEz2 = (TProfile*)Ez2->Clone("cEz2");


	    negEx = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            negEy = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            negEz = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent%d0_%d0", 0, c, c + 1));
            posEx = (TProfile*)f->Get(Form("lrf_Exfield_eta%d_cent%d0_%d0", 2, c, c + 1));
            posEy = (TProfile*)f->Get(Form("lrf_Eyfield_eta%d_cent%d0_%d0", 2, c, c + 1));
            posEz = (TProfile*)f->Get(Form("lrf_Ezfield_eta%d_cent%d0_%d0", 2, c, c + 1));
	    //For full plot
            auto C = new TLegend(0.2,0.4,0.46,0.5);
            C->SetHeader("#bf{Centrality 30%-40%}", "C");
	    C->SetNColumns(3);
            C->AddEntry(Ex, "x ", "l");
            C->AddEntry(Ey, "y ", "l");
            C->AddEntry(Ez, "z ", "l");
            C->SetBorderSize(0);
            C->SetFillStyle(0);

	    auto C1 = new TLegend(0.2,0.8,0.46,0.9);
	    C1->SetHeader("#bf{Centrality 60%-70%}", "C");
 	    C1->SetNColumns(3);
            C1->AddEntry(Ex2, "x ", "l");
            C1->AddEntry(Ey2, "y ", "l");
            C1->AddEntry(Ez2, "z ", "l");
	    C1->SetBorderSize(0);
	    C1->SetFillStyle(0);
	    //For split plots
            auto cc = new TLegend(0.1,0.85,0.9,0.95);
	    cc->SetTextSize(0.035);
	    cc->SetNColumns(4);
	    cc->AddEntry((TObject*)0, "#bf{Centrality: 30%-40%}","" );
            cc->AddEntry(Ex, "x", "l");
            cc->AddEntry(Ey, "y", "l");
            cc->AddEntry(Ez, "z", "l");
	    cc->SetBorderSize(0);
	    cc->SetFillStyle(0);

	    auto Ce = new TLegend(0.05,0.4,0.35,0.5);
	    Ce->SetTextSize(0.033);
	    Ce->SetHeader("#bf{#eta_{s} = - 0.85}", "C");
	    Ce->SetBorderSize(0);
	    Ce->SetFillStyle(0);
            auto Ce2 = new TLegend(0.35,0.4,0.65,0.5);
            Ce2->SetTextSize(0.033);
            Ce2->SetHeader("#bf{#eta_{s} = 0.0}", "C");
            Ce2->SetBorderSize(0);
            Ce2->SetFillStyle(0);
            auto Ce3 = new TLegend(0.65,0.4,0.95,0.5);
            Ce3->SetTextSize(0.033);
            Ce3->SetHeader("#bf{#eta_{s} = 0.85} ", "C");
            Ce3->SetBorderSize(0);
            Ce3->SetFillStyle(0);


            Bx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent30_40", etaID));
            By = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent30_40", etaID));
            Bz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent30_40", etaID));
            cBx = (TProfile*)Bx->Clone("cBx");
            cBy = (TProfile*)By->Clone("cBy");
            cBz = (TProfile*)Bz->Clone("cBz");
            Bx2 = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent60_70", etaID));
            By2 = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent60_70", etaID));
            Bz2 = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent60_70", etaID));
            cBx2 = (TProfile*)Bx2->Clone("cBx2");
            cBy2 = (TProfile*)By2->Clone("cBy2");
            cBz2 = (TProfile*)Bz2->Clone("cBz2");


            negBx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0",0, c, c + 1));
            negBy = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0",0, c, c + 1));
            negBz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0",0, c, c + 1));
            posBx = (TProfile*)f->Get(Form("lrf_Bxfield_eta%d_cent%d0_%d0",2, c, c + 1));
            posBy = (TProfile*)f->Get(Form("lrf_Byfield_eta%d_cent%d0_%d0",2, c, c + 1));
            posBz = (TProfile*)f->Get(Form("lrf_Bzfield_eta%d_cent%d0_%d0",2, c, c + 1));

            auto I = new TLegend(0.1, 0.9, 0.9, 1);
	    //I->SetTextSize(0.08);
            I->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV initialized at #tau_{0} = 0.1 fm/c with #tau_{B} = 0.2 fm/c","C");
	    //I->SetTextSize(0.034);
            //I->SetEntrySeparation(0.005);
	    //I->AddEntry((TObject*)0, "Initialized at #tau_{0} = 0.1 fm/c"," ");
            //I->SetEntrySeparation(0.005);
	    //I->AddEntry((TObject*)0, "With #tau_{B} = 0.2 fm/c", " ");
	    //I->SetTextAlign(13);//I->SetEntrySeparation(0.0);//I->SetColumnSeparation(0);
	    I->SetBorderSize(0);
	    I->SetFillStyle(0);
            Double_t eta = 0.0;

            auto i = new TLegend(0.58,0.46,0.79,0.59);
            i->SetTextSize(0.015);
            i->SetHeader("AVFD Simulation","C");
            i->AddEntry((TObject*)0, "Pb-Pb @ #sqrt{s} = 5.02TeV","");
            i->AddEntry((TObject*)0, Form("Centrality %d0%%-%d0%%", c, c + 1), " ");
            i->AddEntry((TObject*)0, "Initialized at #tau_{0} = 0.1 fm/c"," ");
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
	    TCanvas* c1 = new TCanvas("c1", "EFields",700,400); //legend spacing pipapo...
	    TPad *e1 = new TPad("e1","e1",0,0.5,0.5,1);
	    e1->SetLeftMargin(0.2);
	    e1->SetRightMargin(0.8);
	    e1->SetBottomMargin(0.);
	    e1->SetTopMargin(0.2);
	    e1->Draw();
	    e1->cd();
	    Ey2->SetTitle(" ");
            Ey2->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
	    Ey2->GetYaxis()->SetTitleOffset(0.99);
            Ey2->GetYaxis()->SetTitleSize(0.08);
            Ey2->GetYaxis()->SetLabelSize(0.07);
	    Ey2->GetYaxis()->SetRangeUser(-0.0051,0.0059);
	    Ey2->GetXaxis()->SetRangeUser(0.0,3);
            Ey2->SetLineColor(kRed+3);
            Ey2->DrawCopy("HIST L");//or Hist C
            Ex2->SetLineColor(kBlue+3);
	    Ex2->DrawCopy("HIST L SAME");
            Ez2->SetLineColor(kGreen+3);
            Ez2->DrawCopy("HIST L SAME");
		
            c1->cd();
            TPad *e11 = new TPad("e1","e1",0,0.,0.5,0.5);
            e11->SetLeftMargin(0.2);
	    e11->SetTopMargin(0.);
            e11->SetRightMargin(0.8);
            e11->SetBottomMargin(0.2);
            e11->Draw();
            e11->cd();
            Ey->SetTitle(" ");
            Ey->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            Ey->GetYaxis()->SetTitleSize(0.08);
	    Ey->GetYaxis()->SetLabelSize(0.07);
            Ey->GetXaxis()->SetTitleSize(0.08);
            Ey->GetXaxis()->SetLabelSize(0.07);
            Ey->GetYaxis()->SetTitleOffset(0.99);
            Ey->GetYaxis()->SetRangeUser(-0.0051,0.0059);
            Ey->GetXaxis()->SetTitle("#tau [fm/c]");
            Ey->GetXaxis()->SetRangeUser(0.0,3);
	    Ey->SetLineColor(kRed+1);
            Ey->DrawCopy("HIST L");//or Hist C
            Ex->SetLineColor(kBlue+1);
            Ex->DrawCopy("HIST L SAME");
            Ez->SetLineColor(kGreen+1);
            Ez->DrawCopy("HIST L SAME");
            
            c1->cd();
	    C->Draw();
	    TPad *e12 = new TPad("e12","e12",0.5,0,1,1);
            e12->SetLeftMargin(0.1);
            e12->SetRightMargin(0.1);
            e12->SetBottomMargin(0.1);
	    e12->SetTopMargin(0.1);
	    e12->Draw();
            e12->cd();
	    By2->SetTitle(" ");
            By2->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            By2->GetXaxis()->SetTitle("#tau [fm/c]");
	    By2->GetYaxis()->SetTitleOffset(1.1);//By->GetYaxis()->SetLabelOffset(1.0);	    //By2->GetYaxis()->SetRangeUser(-0.8,0.1);
            By2->GetYaxis()->SetTitleSize(0.04);
            By2->GetXaxis()->SetTitleSize(0.04);
            By2->GetYaxis()->SetLabelSize(0.035);
            By2->GetXaxis()->SetLabelSize(0.035);
	    By2->GetXaxis()->SetRangeUser(0.0,3);
            By2->SetLineColor(kRed+3);
            //By2->SetLineStyle(2);
            By2->Draw("HIST L");//or Hist C
            Bx2->SetLineColor(kBlue+3);
            Bx2->Draw("HIST L SAME");
            Bz2->SetLineColor(kGreen+3);
            Bz2->Draw("HIST L SAME");
            By->SetLineColor(kRed+1);
            By->Draw("HIST L Same");//or Hist C
            Bx->SetLineColor(kBlue+1);
            Bx->Draw("HIST L SAME");
            Bz->SetLineColor(kGreen+1);
            Bz->Draw("HIST L SAME");
	    c1->cd();
	    C1->Draw();
	    I->Draw(); 
	    TPad *e13 = new TPad("e13","e13", 0.65,0.12,0.93,0.7);
	    e13->Draw();
	    e13->cd();
	    Bx2->SetTitle("");
            Bx2->GetYaxis()->SetRangeUser(-0.01,0.02);
	    Bx2->GetXaxis()->SetRangeUser(0,1);
            Bx2->SetLineColor(kBlue+3);
	    Bx2->SetLineWidth(1);
            Bx2->Draw("HIST L");//or Hist C
            Bz2->SetLineColor(kGreen+3);
	    Bz2->SetLineWidth(1);
            Bz2->Draw("HIST L SAME");
            Bx->SetLineColor(kBlue+1);
	    Bx->SetLineWidth(1);
            Bx->Draw("HIST L SAME");
	    Bz->SetLineWidth(1);
            Bz->SetLineColor(kGreen+1);
	    Bz->Draw("HIST L SAME");

            c1->SaveAs("Evolution/ThesisEvol.png");//double tau..
            Ex->Reset();
            Ey->Reset();
            Ez->Reset();
            Bx->Reset();
            By->Reset();
            Bz->Reset();
            Ex2->Reset();
            Ey2->Reset();
            Ez2->Reset();
            Bx2->Reset();
            By2->Reset();
            Bz2->Reset();

	    delete c1;

	    Double_t ymin = -0.033;
	    Double_t ymax = 0.033;
	    Double_t bmin = -0.0075;//the full -B peak should be in the first plot!
	    Double_t bmax = 0.0025;
	    //Canvas for cuts in different eta
            TCanvas* c2 = new TCanvas("e2", "EFields");//,200,10,700,500); //legend spacing pipapo...
            TPad* e20 = new TPad("e20","e20",0.0,0.5,0.35,0.9);
	    e20->SetLeftMargin(0.23);
            e20->SetRightMargin(0.);
	    e20->SetBottomMargin(0.);
	    e20->SetTopMargin(0.15);
	    e20->Draw();
	    e20->cd();
            //negEx->SetTitleSize(0.4);
	    negEx->SetTitle("");
            negEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            negEx->GetXaxis()->SetTitle("#tau [fm/c]");
            negEx->GetYaxis()->SetTitleOffset(1.4);
	    negEx->GetYaxis()->SetTitleSize(0.08);
	    negEx->GetYaxis()->SetLabelSize(0.07);
	    negEx->GetXaxis()->SetTitleOffset(0.4);
	    negEx->GetYaxis()->SetRangeUser(ymin,ymax);
            negEx->SetLineColor(kBlue+1);
            negEx->GetXaxis()->SetRangeUser(0.0,0.82);
            negEx->Draw("HIST L");
            negEy->SetLineColor(kRed+1);
            negEy->Draw("HIST L SAME");
            negEz->SetLineColor(kGreen+1);
            negEz->Draw("HIST L SAME");
            
	    c2->cd();
	    I->Draw();
            TPad* e21 = new TPad("e21","e21",0.35,0.5,0.65,0.9);
            e21->SetLeftMargin(0.);
            e21->SetRightMargin(0.);
            e21->SetBottomMargin(0);
	    e21->SetTopMargin(0.15);
            e21->Draw();
            e21->cd();
            cEx->SetTitle(" ");
	    cEx->SetTitleOffset(1);
            //cEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            cEx->GetXaxis()->SetTitle("#tau [fm/c]");
            cEx->GetXaxis()->SetTitleOffset(0.9);
	    cEx->GetYaxis()->SetRangeUser(ymin,ymax);
            cEx->SetLineColor(kBlue+1);
            cEx->GetXaxis()->SetRangeUser(0.0,0.82);
            cEx->Draw("Hist L");
            cEy->SetLineColor(kRed+2);
            cEy->Draw("HIST L SAME");
            cEz->SetLineColor(kGreen+1);
            cEz->Draw("HIST L SAME");

	    c2->cd();
            TPad* e22 = new TPad("e22","e22",0.65,0.5,0.95,0.9);
            e22->SetLeftMargin(0.);
            e22->SetRightMargin(0.1);
            e22->SetBottomMargin(0);
	    e22->SetTopMargin(0.15);
            e22->Draw();
            e22->cd();
            posEx->SetTitle(" ");
	    posEx->SetTitleOffset(4);
           // posEx->GetYaxis()->SetTitle("E^{lrf} [1/fm^{2}]");
            posEx->GetXaxis()->SetTitle("#tau [fm/c]");
            posEx->GetXaxis()->SetTitleOffset(0.9);
	    posEx->GetYaxis()->SetRangeUser(ymin,ymax);
            posEx->SetLineColor(kBlue+1);
            posEx->GetXaxis()->SetRangeUser(0.0,0.82);
            posEx->Draw("HIST L");
            posEy->SetLineColor(kRed+1);
            posEy->Draw("HIST L SAME");
            posEz->SetLineColor(kGreen+1);
            posEz->Draw("HIST L SAME");
	    
	    c2->cd();
	    TPad* e23 = new TPad("e23","e23",0.0,0.,0.35,0.45);
            e23->SetLeftMargin(0.23);
            e23->SetRightMargin(0.);
            e23->SetBottomMargin(0.2);
	    e23->SetTopMargin(0.1);
            e23->Draw();
            e23->cd();	 
            negBx->SetTitle("");
            negBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            negBx->GetYaxis()->SetTitleOffset(1.5);
            negBx->GetYaxis()->SetTitleSize(0.07);
	    negBx->GetYaxis()->SetLabelSize(0.06);
	    negBx->GetYaxis()->SetNdivisions(6);
	    negBx->GetXaxis()->SetTitleSize(0.07);//negBx->GetXaxis()->SetNdivisions(2);
	    negBx->GetXaxis()->SetLabelSize(0.06);
            negBx->GetXaxis()->SetTitleOffset(0.8);
            negBx->GetXaxis()->SetTitle("#tau [fm/c]");
	    negBx->GetYaxis()->SetRangeUser(bmin,bmax);
            negBx->SetLineColor(kBlue+1);
            negBx->GetXaxis()->SetRangeUser(0.0,0.82);
            negBx->Draw("HIST L");
            negBy->SetLineColor(kRed+1);
            negBy->Draw("HIST L SAME");
            negBz->SetLineColor(kGreen+1);
            negBz->Draw("HIST L SAME");
	    	    
	    c2->cd();
	    cc->Draw();
            TPad* e24 = new TPad("e24","e24",0.35,0.0,0.65,0.45);
            e24->SetLeftMargin(0.);
            e24->SetRightMargin(0.);
            e24->SetBottomMargin(0.2);
            e24->SetTopMargin(0.1);
            e24->Draw();
            e24->cd();
            cBx->SetTitle("");
            //cBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            cBx->GetXaxis()->SetTitle("#tau [fm/c]");
            cBx->GetXaxis()->SetTitleSize(0.07);
	    cBx->GetYaxis()->SetNdivisions(6);
            cBx->GetXaxis()->SetLabelSize(0.06);
            cBx->GetXaxis()->SetTitleOffset(0.8);
	    cBx->GetYaxis()->SetRangeUser(bmin,bmax);
            cBx->SetLineColor(kBlue+1);
            cBx->GetXaxis()->SetRangeUser(0.0,0.82);
            cBx->Draw("HIST L");
            cBy->SetLineColor(kRed+1);
            cBy->Draw("HIST L SAME");
            cBz->SetLineColor(kGreen+1);
            cBz->Draw("HIST L SAME");
            
	    c2->cd();
	    TPad* e25 = new TPad("e25","e25",0.65,0.,0.95,0.45);
            e25->SetLeftMargin(0.);
            e25->SetRightMargin(0.1);
            e25->SetBottomMargin(0.2);
            e25->SetTopMargin(0.1);
            e25->Draw();
            e25->cd();
            posBx->SetTitle(" ");
            //posBx->GetYaxis()->SetTitle("B^{lrf} [1/fm^{2}]");
            posBx->GetXaxis()->SetTitle("#tau [fm/c]");
            posBx->GetXaxis()->SetTitleSize(0.07);
	    posBx->GetYaxis()->SetNdivisions(6);
            posBx->GetXaxis()->SetLabelSize(0.06);
            posBx->GetXaxis()->SetTitleOffset(0.8);
	    posBx->GetYaxis()->SetRangeUser(bmin,bmax);
            posBx->SetLineColor(kBlue+1);
            posBx->GetXaxis()->SetRangeUser(0.0,0.82);
            posBx->Draw("HIST L");
            posBy->SetLineColor(kRed+1);
            posBy->Draw("HIST L SAME");
            posBz->SetLineColor(kGreen+1);
            posBz->Draw("HIST L SAME");
	    
	    c2->cd();
	    Ce->Draw();
	    Ce2->Draw();
	    Ce3->Draw();
	    c2->SaveAs("Evolution/EMFields_Thesis_Splitsevolution.pdf");//double tau..
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
 	    delete c2;
    //    }
    }
}
