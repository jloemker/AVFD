#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotMacro(){
//gStyle -> SetOptStat(0);
//gROOT->SetBatch();
TFile* File = new TFile("/user/jlomker/project/AVFD/result/AnalysisResult_a-0.1_5.44TeV.root");

//v2 = File.Get('FlowFromBWList/V2IntPro');//or V2IntProQC
//gpp= File.Get('CMEList/')
//                 
//Getting the list(s) of files
TList *QA = (TList*) File->Get("QAList;1");
//Extracting the histogram(s) from the list             
TH1F *pT = (TH1F*) QA->FindObject("fPtChargedParticlesDistribution");
//Define the Canvas for the plot
TCanvas *c1 = new TCanvas("c1","pT",400,400);
//
// Some things to make things more professional 
//pT -> GetYaxis()->SetTitle("#eta AK8_{2}");
//pT ->GetXaxis()->SetTitle("#eta AK8_{1}");
//pT -> SetTitle("#eta AK8_{1} vs #eta AK8_{2}");
pT->Draw();
c1->SaveAs("pT.pdf");
/*
TCanvas *c2 = new TCanvas("c2","Delta delta",400,400);
TH1F *delDelta = (TH1F*) File->Get("CMEList;1/DeltaD11");
delDelta->Draw();
c2->SaveAs("delDelta.pdf");

TCanvas *c3 = new TCanvas("c3","Delta gamma", 400, 400);
TH1F *delGamma = (TH1F*) File->Get("CMEList;1/DeltaG112");
delGamma ->Draw();
c3->SaveAs("delGamma.pdf");
*/
}

