#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotMacro(){
//gStyle -> SetOptStat(0);
//gROOT->SetBatch();
TFile* File = new TFile("/project/alice/users/jlomker/AnalysisResult_a-0.1_5.44TeV.root");

//v2 = File.Get('FlowFromBWList/V2IntPro');//or V2IntProQC
//gpp= File.Get('CMEList/')

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
pT->Draw();
c1->SaveAs("pT.pdf");

//Delta delta
TCanvas *c2 = new TCanvas("c2","Delta delta",400,400);
TH1F *delDelta = (TH1F*) CME->FindObject("DeltaD11");
//delDelta ->GetYaxis()->SetTitle("#Delta #delta"); not reallz sure if it is what i think
delDelta ->GetXaxis()->SetTitle("Centrality");
delDelta->Draw();
c2->SaveAs("delD11.pdf");

//Delta gamma
TCanvas *c3 = new TCanvas("c3","Delta gamma", 400, 400);
TH1F *delGamma = (TH1F*) CME->FindObject("DeltaG112");
//delgamma ->GetYaxis()->SetTitle("#Delta #Gamma"); not reallz sure if it is what i think
delGamma ->GetXaxis()->SetTitle("Centrality");
delGamma ->Draw();
c3->SaveAs("delG112.pdf");
}

