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
TH1F *v2_pos = (TH1F*) posFlowQC->FindObject("fFlowQCFinalEtaDifHist[0][2][0][0]");

//TH1F *v3_all = (TH1F*) flowQC->FindObject("fFlowQCFinalEtaDifHist[2][1][0]");
TH1F *v3_pos = (TH1F*) posFlowQC->FindObject("fFlowQCFinalEtaDifHist[0][2][1][0]");

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
legend->Draw();
v2 ->SaveAs("test_eta_v1.pdf");

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
legend2->Draw();
v3 ->SaveAs("test_small_eta_v2.pdf");

//Plots for large sample
TFile* File_all = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/all/v1/AnalysisResults_5.44TeV_Cent20_30.root");
TFile* File_pos = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/pos/v1/AnalysisResults_5.44TeV_Cent20_30.root");
TList *QA = (TList*) File_all->Get("QAList;1");//Kinematics from POI`s
TList *FlowQC = (TList*) File_all->Get("FlowQCList;1");
TList *PosFlowQC = (TList*) File_pos->Get("FlowQCList;1");
//pt differential plots
TH1F *V2_all = (TH1F*) FlowQC->FindObject("fFlowQCFinalPtDifHist[2][0][0]");
TH1F *V2_pos = (TH1F*) PosFlowQC->FindObject("fFlowQCFinalPtDifHist[2][0][0]");

TH1F *V3_all = (TH1F*) FlowQC->FindObject("fFlowQCFinalPtDifHist[2][1][0]");
TH1F *V3_pos = (TH1F*) PosFlowQC->FindObject("fFlowQCFinalPtDifHist[2][1][0]");

TCanvas *V2 = new TCanvas("V2","pT",400,400);
V2_all ->GetYaxis()->SetTitle("differential flow v_{1}");
V2_all ->GetXaxis()->SetTitle("p_{T} [GeV]");
V2_all ->SetLineColor(kBlack);
V2_all->GetXaxis()->SetRangeUser(0.,5.);
V2_all->Draw();
V2_pos ->GetYaxis()->SetTitle("differential flow v_{1}");
V2_pos ->GetXaxis()->SetTitle("p_{T} [GeV]");
V2_pos ->GetXaxis()->SetRangeUser(0.,5.);
V2_pos->SetLineColor(kRed);
V2_pos->Draw("same");
auto Legend = new TLegend();
Legend->SetHeader("Large sample Xe-Xe: 5.44TeV, 20-30%","C");
Legend->AddEntry(V2_all,"Charged particles", "l");
Legend->AddEntry(V2_pos,"Pos particles", "l");
Legend->Draw();
V2 ->SaveAs("charged_vs_pos_pt_v1.pdf");

TCanvas *V3 = new TCanvas("V3","pT",400,400);
V3_all ->GetYaxis()->SetTitle("differential flow v_{2}");
V3_all ->GetXaxis()->SetTitle("p_{T} [GeV]");
V3_all ->SetLineColor(kBlack);
V3_all->GetXaxis()->SetRangeUser(0.,5.);
V3_all->Draw();
V3_pos ->GetYaxis()->SetTitle("differential flow v_{2}");
V3_pos ->GetXaxis()->SetTitle("p_{T} [GeV]");
V3_pos ->GetXaxis()->SetRangeUser(0.,5.);
V3_pos->SetLineColor(kRed);
V3_pos->Draw("same");
auto Legend2 = new TLegend(0.1,0.7,0.48,0.9);
Legend2->SetHeader("Large sample Xe-Xe: 5.44TeV, 20-30%","C");
Legend2->AddEntry(V3_all,"Charged particles", "l");
Legend2->AddEntry(V3_pos,"Pos particles", "l");
Legend2->Draw();
V3 ->SaveAs("charged_vs_pos_pt_v2.pdf");

//Comparison plots with real data
//define TH1F here  How the hack can I get the values from the file ?!
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

