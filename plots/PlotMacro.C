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
TFile* File = new TFile("/project/alice/users/jlomker/AnalysisResult_a-0.1_5.44TeV.root");

//Getting the list(s) of files
TList *QA = (TList*) File->Get("QAList;1");//Kinematics from POI`s
TList *FlowQC = (TList*) File->Get("FlowQCList;1");//Matrices for calculations ?
TList *FlowGF = (TList*) File->Get("FLowGFList;1");//Matrices for calculations ?
TList *FlowFromBW = (TList*) File->Get("FlowFromBWList;1");//Flow harmonics v2 -> v4
TList *CME = (TList*) File->Get("CMEList;1");//(SS/OS/Delta) delta, gamma 
TList *CMW = (TList*) File->Get("CMWList;1");//Different methods for flow harmonics, eta, etc..
TList *CMWQA = (TList*) File->Get("CMWQAList;1");//2 more (empty) directories

//Transverse momentum example
//Extracting the histogram(s) from the list             
TH1F *pT = (TH1F*) QA->FindObject("fPtChargedParticlesDistribution");
//Define the Canvas for the plot
TCanvas *c1 = new TCanvas("c1","pT",400,400);
// Some things to make things more professional 
pT ->GetYaxis()->SetTitle("Transverse momentum p_{T}");
pT ->GetXaxis()->SetTitle("GeV");
pT ->SetTitle("Charged particles");
//pT->Draw();
//c1->SaveAs("pT.pdf");


//double canvas:
TCanvas *duo = new TCanvas("duo","duo",200,300);
duo->Divide(1,2,0,0);

duo->cd(1);
duo->cd(1)->SetTopMargin(0.1);
duo->cd(1)->SetLeftMargin(0.2);
duo->cd(1)->SetRightMargin(0.1);
duo->cd(1)->Draw();

TH1F *delDelta = (TH1F*) CME->FindObject("DeltaD11");
delDelta->SetTitle("");
//some beauty issues left - marker style and consistency with the legend ... 
delDelta->SetLineColor(2);
delDelta->SetLineStyle(1);
delDelta->SetLineWidth(3);
delDelta->GetYaxis()->SetRangeUser(-0.002,0.004);
delDelta->GetYaxis()->SetNdivisions(3);
delDelta->GetYaxis()->SetTitle("#Delta #delta_{1}");
delDelta->GetYaxis()->SetTitleSize(0.05);
delDelta->GetYaxis()->CenterTitle();
delDelta->GetXaxis()->SetRangeUser(0,70);
delDelta->GetXaxis()->SetTickLength(0);
delDelta->GetXaxis()->SetLabelSize(0);
delDelta->Draw();
//Adding the legend -> some beauty issues left
auto legend = new TLegend(0.3,0.35);
legend->SetHeader("AVFD Xe-Xe, #sqrt{s_{NN}} = 4.55TeV","C");
legend->AddEntry("delDelta", " n_{5}/s = 0.1 - LCC = 0");

duo->cd(2);
duo->cd(2)->SetLeftMargin(0.2);
duo->cd(2)->SetBottomMargin(0.2);
duo->cd(2)->SetRightMargin(0.1);
duo->cd(2)->Draw();
TH1F *delGamma = (TH1F*) CME->FindObject("DeltaG112");
delGamma->SetTitle("");
delGamma->SetLineColor(2);
delGamma->SetLineStyle(1);
delGamma->SetLineWidth(3);
delGamma->GetYaxis()->SetRangeUser(0,0.0022);
delGamma->GetYaxis()->SetNdivisions(3);
delGamma->GetYaxis()->SetTitle("#Delta #gamma_{1,1}");
delGamma->GetYaxis()->SetTitleSize(0.05);
delGamma->GetYaxis()->CenterTitle();
delGamma->GetXaxis()->SetRangeUser(0, 70);
delGamma ->GetXaxis()->SetTitle("centrality, %");
delGamma ->GetXaxis()->CenterTitle();
delGamma ->GetXaxis()->SetTitleSize(0.05);
delGamma ->Draw();

legend->AddEntry("delGamma", "n_{5}/s = 0.1 - LCC = 0");
legend->SetBorderSize(0);
legend->Draw();
duo->SaveAs("delD_delG.pdf");
}

