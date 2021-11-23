#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

//easy for ratio plot: https://root.cern/doc/v610/ratioplotOld_8C_source.html

void PlotMacro(){
gStyle -> SetOptStat(0);
//gROOT->SetBatch();
//

//**************Quick and Dirty *
//				*
//*******************************
//small sets in project test dir with all or pos: /project/alice/users/jlomker/AVFD/test/dirID-0/all/AnalysisResults_Cent20_30.root
//for large sets all/pos: /data/alice/jlomker/AVFD/result/dirID-0/all/AnalysisResults_Cent20_30.root
//data/alice/jlomker/AVFD/result/dirID-0/pos
//AnalysisResults_5.44TeV_Cent20_30.root

//plot small sample 2 for all and pos on top of eachother:
TFile* file_all = new TFile("/project/alice/users/jlomker/AVFD/test/dirID-0/all/v1/AnalysisResults_Cent20_30.root");
TFile* file_pos = new TFile("/project/alice/users/jlomker/AVFD/test/dirID-0/pos/eta_pt/Analysis_pTrange_0_Cent20_30.root");
//Getting the list(s) of files for all
TList *qa = (TList*) file_all->Get("QAList;1");//Kinematics from POI`s
TList *flowQC = (TList*) file_all->Get("FlowQCList;1");//Matrices for calculations ?
TList *cme = (TList*) file_all->Get("CMEList;1");//(SS/OS/Delta) delta, gamma 
//for pos particles 
TList *posFlowQC = (TList*) file_pos->Get("FlowQCList;1");//Matrices for calculations ?

//TH1F *v2_all = (TH1F*) flowQC->FindObject("fFlowQCFinalEtaDifHist[2][0][0]");
TH1F *v2_pos = (TH1F*) posFlowQC->FindObject("fFlowQCFinalEtaDifDeltaHist[2][0]");

//TH1F *v3_all = (TH1F*) flowQC->FindObject("fFlowQCFinalEtaDifHist[2][1][0]");
TH1F *v3_pos = (TH1F*) posFlowQC->FindObject("fFlowQCFinalEtaDifDeltaHist[2][1]");

TCanvas *v2 = new TCanvas("v2","pT",400,400);
v2->SetLeftMargin(0.2);
v2_pos->GetYaxis() -> SetTitleOffset(2.0);
/*v2_all ->GetYaxis()->SetTitle("differential flow v_{1}");
v2_all ->GetXaxis()->SetTitle("p_{T} [GeV]");
v2_all ->SetLineColor(kBlack);
v2_all->GetXaxis()->SetRangeUser(0.,5.);
v2_all->Draw();*/
v2_pos ->GetYaxis()->SetTitle("differential flow v_{1}");
//v2_pos ->GetXaxis()->SetTitle("p_{T} [GeV]");
v2_pos->GetXaxis()->SetTitle("eta");
v2_pos ->GetXaxis()->SetRangeUser(-1.,1.);
v2_pos->SetLineColor(kRed);
v2_pos->Draw("same");
auto legend = new TLegend();
legend->SetHeader("Small sample Xe-Xe: 5.44TeV, 20-30%","C");
//legend->AddEntry(v2_all,"Charged particles", "l");
legend->AddEntry(v2_pos,"Pos particles", "l");
//legend->Draw();
//v2 ->SaveAs("test_eta_delta_v1.pdf");

TCanvas *v3 = new TCanvas("v3","pT",400,400);
v3->SetLeftMargin(0.2);
v3_pos->GetYaxis() -> SetTitleOffset(2.0);
/*
v3_all ->GetYaxis()->SetTitle("differential flow v_{2}");
v3_all ->GetXaxis()->SetTitle("p_{T} [GeV]");
v3_all ->SetLineColor(kBlack);
v3_all->GetXaxis()->SetRangeUser(0.,5.);
v3_all->Draw();*/
v3_pos ->GetYaxis()->SetTitle("differential flow v_{2}");
//v3_pos ->GetXaxis()->SetTitle("p_{T} [GeV]");
v3_pos ->GetXaxis()->SetTitle("eta");
v3_pos ->GetXaxis()->SetRangeUser(-1.,1.);
v3_pos->SetLineColor(kRed);
v3_pos->Draw("same");
auto legend2 = new TLegend();
legend2->SetHeader("Small sample Xe-Xe: 5.44TeV, 20-30%","C");
//legend2->AddEntry(v3_all,"Charged particles", "l");
legend2->AddEntry(v3_pos,"Pos particles", "l");
//legend2->Draw();
//v3 ->SaveAs("test_small_eta_delta_v2.pdf");






//Plots for large sample
TFile* FpT = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Results_5.44TeV_pTrange_0_Cent20_30.root");//Full pT range
TFile* SpT = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Results_5.44TeV_pTrange_1_Cent20_30.root");//Small pT range
TFile* LpT = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Results_5.44TeV_pTrange_2_Cent20_30.root");//Large pT range

//TList *QA = (TList*) File_all->Get("QAList;1");//Kinematics from POI`s
TList *full = (TList*) FpT->Get("FlowQCList;1");
TList *small = (TList*) SpT->Get("FlowQCList;1");
TList *large = (TList*) LpT->Get("FlowQCList;1");

//Differential plots*******************************************************************
//										      *
//        TH1D *fFlowQCFinal**DifHist[Ncharge][fCRCMaxnCen][fFlowNHarm][fFlowQCNCov]; *
//        TH1D *fFlowQCFinal**DifDeltaHist[fCRCMaxnCen][fFlowNHarm];		      *
//*************************************************************************************

//************
//Plots vs eta 
//************        

///For vn vs eta differential flow for positive 0 and negative 1 particles in 3 pTranges
TH1F *v1_full = (TH1F*) full->FindObject("fFlowQCFinalEtaDifHist[0][2][0][0]");
TH1F *v1_small = (TH1F*) small->FindObject("fFlowQCFinalEtaDifHist[0][2][0][0]");
TH1F *v1_large = (TH1F*) large->FindObject("fFlowQCFinalEtaDifHist[0][2][0][0]");

TCanvas *V1 = new TCanvas("V1","Eta",400,400);
V1->SetLeftMargin(0.2);
V1->SetBottomMargin(0.2);
v1_full->GetYaxis()->SetTitle("differential flow v_{1}");//change n
v1_full->GetXaxis()->SetTitle("#eta");
v1_full->SetLineColor(kGreen);
v1_full->GetXaxis()->SetRangeUser(-0.9,0.9);
v1_full->GetXaxis()->SetTitleOffset(2.0);
v1_full->Draw();
v1_small->SetLineColor(kRed);
v1_small->Draw("SAME");
v1_large->SetLineColor(kBlue);
v1_large->Draw("SAME");
auto Legend1 = new TLegend();
Legend1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend1->AddEntry(v1_full,"pT: 0.2-5 GeV", "l");
Legend1->AddEntry(v1_small,"pT: 0.2-1 GeV", "l");
Legend1->AddEntry(v1_large,"pT: 3.0-5 GeV");
Legend1->Draw();
V1->SaveAs("Pos_diffv1_eta.pdf");//change Pos/Neg and vn

///For differental delta vn vs eta in 3 pTranges
TH1F *delv1_full = (TH1F*) full->FindObject("fFlowQCFinalEtaDifDeltaHist[2][0]");
TH1F *delv1_small = (TH1F*) small->FindObject("fFlowQCFinalEtaDifDeltaHist[2][0]");
TH1F *delv1_large = (TH1F*) large->FindObject("fFlowQCFinalEtaDifDeltaHist[2][0]");

TCanvas *delV1 = new TCanvas("DeltaV1","Eta",400,400);
delV1->SetLeftMargin(0.2);
delV1->SetBottomMargin(0.2);
delv1_full->GetYaxis()->SetTitle("differential #Delta v_{1}");//change n
delv1_full->GetXaxis()->SetTitle("#eta");
delv1_full->SetLineColor(kGreen);
delv1_full->GetXaxis()->SetRangeUser(-0.9,0.9);
delv1_full->GetXaxis()->SetTitleOffset(2.0);
delv1_full->Draw();
delv1_small->SetLineColor(kRed);
delv1_small->Draw("SAME");
delv1_large->SetLineColor(kBlue);
delv1_large->Draw("SAME");
auto Legend2 = new TLegend();
Legend2->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend2->AddEntry(delv1_full,"pT: 0.2-5 GeV", "l");
Legend2->AddEntry(delv1_small,"pT: 0.2-1 GeV", "l");
Legend2->AddEntry(delv1_large,"pT: 3.0-5 GeV");
Legend2->Draw();
delV1->SaveAs("Del_diffv1_eta.pdf");//change vn


//************
//Plots vs pT
//************

//For vn vs pT differential flow for pos and neg 
//in each pt range the vn from pos & neg particle
//
TH1F *pos_full = (TH1F*) full->FindObject("fFlowQCFinalPtDifHist[0][2][0][0]");
TH1F *neg_full = (TH1F*) full->FindObject("fFlowQCFinalPtDifHist[1][2][0][0]");
TH1F *pos_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[0][2][0][0]");
TH1F *neg_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[1][2][0][0]");
TH1F *pos_large = (TH1F*) large->FindObject("fFlowQCFinalPtDifHist[0][2][0][0]");
TH1F *neg_large = (TH1F*) large->FindObject("fFlowQCFinalPtDifHist[1][2][0][0]");

//TH1F *pos_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[0][2][0][0]");
TCanvas *Vn = new TCanvas("Vn","Pt",400,400);
Vn->SetLeftMargin(0.2);
Vn->SetBottomMargin(0.2);
pos_full->GetYaxis()->SetTitle("differential flow v_{1}");//change n
pos_full->GetXaxis()->SetTitle("p_{T} [GeV]");
pos_full->SetLineColor(kGreen);
pos_full->GetXaxis()->SetRangeUser(0.2,5.0);
pos_full->GetXaxis()->SetTitleOffset(1.0);
pos_full->Draw();
neg_full->SetLineColor(kRed);
neg_full->Draw("SAME");
auto Legend3 = new TLegend();
Legend3->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend3->AddEntry(pos_full,"Positive particles", "l");
Legend3->AddEntry(neg_full,"Negative particles", "l");
Legend3->Draw();
Vn->SaveAs("diffv1_pT_0.pdf");//change Pos/Neg and vn

TCanvas *Vn1 = new TCanvas("Vn1","Pt1",400,400);
Vn1->SetLeftMargin(0.2);
Vn1->SetBottomMargin(0.2);
pos_small->GetYaxis()->SetTitle("differential flow v_{1}");//change n
pos_small->GetXaxis()->SetTitle("p_{T}");
pos_small->SetLineColor(kGreen);
pos_small->GetXaxis()->SetRangeUser(0.2,1.0);
pos_small->GetXaxis()->SetTitleOffset(1.0);
pos_small->Draw();
neg_small->SetLineColor(kRed);
neg_small->Draw("SAME");
auto Legend4 = new TLegend();
Legend4->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend4->AddEntry(pos_full,"Positive particles", "l");
Legend4->AddEntry(neg_full,"Negative particles", "l");
Legend4->Draw();
Vn1->SaveAs("diffv1_pT_1.pdf");//change Pos/Neg and vn

TCanvas *Vn2 = new TCanvas("Vn2","Pt2",400,400);
Vn2->SetLeftMargin(0.2);
Vn2->SetBottomMargin(0.2);
pos_large->GetYaxis()->SetTitle("differential flow v_{1}");//change n
pos_large->GetXaxis()->SetTitle("p_{T} [GeV]");
pos_large->SetLineColor(kGreen);
pos_large->GetXaxis()->SetRangeUser(3.0,5.0);
pos_large->GetXaxis()->SetTitleOffset(1.0);
pos_large->Draw();
neg_large->SetLineColor(kRed);
neg_large->Draw("SAME");
auto Legend5 = new TLegend();
Legend5->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend5->AddEntry(pos_large,"Positive particles", "l");
Legend5->AddEntry(neg_large,"Negative particles", "l");
Legend5->Draw();
Vn2->SaveAs("diffv1_pT_2.pdf");//change Pos/Neg and vn


//for deltavn vs pT plots
//in each pt range the vn from pos & neg particle: full small large to find objects in pT ranges 0.2-5, 0.2-1 and 3.0 to 5
TH1F *delvn = (TH1F*) full->FindObject("fFlowQCFinalPtDifDeltaHist[2][0]");

TCanvas *delVn = new TCanvas("DeltaVn","Pt",400,400);
delVn->SetLeftMargin(0.2);
delVn->SetBottomMargin(0.2);
delvn->GetYaxis()->SetTitle("differential #Delta v_{1}");//change n
delvn->GetXaxis()->SetTitle("p_{T} [GeV]");
delvn->SetLineColor(kGreen);
delvn->GetXaxis()->SetRangeUser(0.2,5.0);
delvn->GetXaxis()->SetTitleOffset(2.0);
delvn->GetXaxis()->SetTitleOffset(1.0);
delvn->Draw();
//delV1_small->SetLineColor(kRed); 0.2 to 1 Color code for the 2 pT ranges
//delv1_large->SetLineColor(kBlue); 3 to 5
auto Legend6 = new TLegend();
Legend6->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend6->AddEntry(delvn,"pT: 0.2-5 GeV", "l");
//Legend->AddEntry(delV1_small,"pT: 0.2-1 GeV", "l");
//Legend->AddEntry(delv1_large,"pT: 3.0-5 GeV")
Legend6->Draw();
delVn->SaveAs("Del_diffv1_pT_0.pdf");//change Pos/Neg and vn

/*

TFile *Published = new TFile("/project/alice/users/jlomker/AVFD/plots/published.root");
TList *p = (TList*) Published->Get("Table 5");
 
TH1F *h1 = (TH1F*) p->FindObject("");
TH1F *h2 = (TH1F*) FlowQC->FindObject("fFlowQCFinalPtDifHist[2][0][0]");


//define canvas
TCanvas *c = new TCanvas("c", "canvas", 800, 800);
//upper plot in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetGridx();         // Vertical grid
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pad
h1->SetStats(0);          // No statistics on upper plot
h1->Draw();               // Draw h1
h2->Draw("same");         // Draw h2 on top of h1

h1->GetYaxis()->SetLabelSize(0.);
TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
axis->SetLabelSize(15);
axis->Draw();
//lower plot in pad
c->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();       // pad2 becomes the current pad
//define raito plot
TH1F *h3 = (TH1F*)h1->Clone("h3");
h3->SetLineColor(kBlack);
h3->SetMinimum(0.8);  // Define Y ..
h3->SetMaximum(1.35); // .. range
h3->Sumw2();
h3->SetStats(0);      // No statistics on lower plot
h3->Divide(h2);
h3->SetMarkerStyle(21);
h3->Draw("ep");       // Draw the ratio plot
//h1 settings
h1->SetLineColor(kBlue+1);
h1->SetLineWidth(2);
// Y axis h1 plot settings
h1->GetYaxis()->SetTitleSize(20);
h1->GetYaxis()->SetTitleFont(43);
h1->GetYaxis()->SetTitleOffset(1.55);
//h2 setting
h2->SetLineColor(kRed);
h2->SetLineWidth(2);

// Ratio plot (h3) settings
h3->SetTitle(""); // Remove the ratio title 
// Y axis ratio plot settings
h3->GetYaxis()->SetTitle("ratio h1/h2 ");
h3->GetYaxis()->SetNdivisions(505);
h3->GetYaxis()->SetTitleSize(20);
h3->GetYaxis()->SetTitleFont(43);
h3->GetYaxis()->SetTitleOffset(1.55);
h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetYaxis()->SetLabelSize(15);
// X axis ratio plot settings
h3->GetXaxis()->SetTitleSize(20);
h3->GetXaxis()->SetTitleFont(43);
h3->GetXaxis()->SetTitleOffset(4.);
h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetXaxis()->SetLabelSize(15);

c->SaveAs("ratio.pdf");*/
}

