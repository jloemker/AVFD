#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotMacro(){
//gStyle -> SetOptStat(0);
//gROOT->SetBatch();
TFile * File = new TFile("/user/jlomker/project/AVFD/result/AnalysisResult_a-0.1_5.44TeV.root");

//v2 = File.Get('FlowFromBWList/V2IntPro');//or V2IntProQC
//gpp= File.Get('CMEList/')
              
TCanvas *c1 = new TCanvas("c1","pT",400,400);                  
TDirectory *dir1 = (TDirectory*) File->Get("QAList;1");
TH1F *pT = (TH1F*) dir1->Get("fPtChargedParticlesDistribution");
//h_eta -> GetYaxis()->SetTitle("#eta AK8_{2}");
//h_eta ->GetXaxis()->SetTitle("#eta AK8_{1}");
//h_eta -> SetTitle("#eta AK8_{1} vs #eta AK8_{2}");
pT->Draw();
c1->SaveAs("pT.pdf");

TCanvas *c2 = new TCanvas("c2","Delta delta",400,400);
TH1F *delDelta = (TH1F*) File->Get("CMEList;1/DeltaD11");
delDelta->Draw();
c2->SaveAs("delDelta.pdf");

TCanvas *c3 = new TCanvas("c3","Delta gamma", 400, 400);
TH1F *delGamma = (TH1F*) File->Get("CMEList;1/DeltaG112");
delGamma ->Draw();
c3->SaveAs("delGamma.pdf");
}
