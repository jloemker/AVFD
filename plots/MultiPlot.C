#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void MultiPlot(){
gStyle -> SetOptStat(0);
for(Int_t cen = 1; cen<7; cen++){

TFile* f = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-0/new/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",cen,cen+1));
TList *l = (TList*) f->Get("FlowQCList;1");

for(Int_t harm = 0; harm<1;harm++){
//histograms
//eta and pt pos/neg

TH1F *PvpT = (TH1F*) l->FindObject(Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",cen,harm));
TH1F *NvpT = (TH1F*) l->FindObject(Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",cen,harm));
TH1F *DvpT = (TH1F*) l->FindObject(Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",cen,harm));

TH1F *PvEta = (TH1F*) l->FindObject(Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",cen,harm));
TH1F *NvEta = (TH1F*) l->FindObject(Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",cen,harm));
TH1F *DvEta = (TH1F*) l->FindObject(Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",cen,harm));
/*
if(harm == 0){//figure out how to get delta v1 ...{
TProfile *PvpT = (TProfile*) f->FindObject(Form("fFlowQCQv1[0][%d][0]",cen));
TProfile *NvpT = (TProfile*) f->FindObject(Form("fFlowQCQv1[1][%d][0]",cen));
TProfile *PvEta = (TProfile*) f->FindObject(Form("fFlowQCQv1[0][%d][1]",cen));
TProfile *NvEta = (TProfile*) f->FindObject(Form("fFlowQCQv1[1][%d][1]",cen));
}*/

//pos and negative and delta in eta and pt for pos and neg

//pos and neg in pT
TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1., 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetLeftMargin(0.2);
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pad//and filled with pos
PvpT->GetXaxis()->SetRangeUser(0.,5);
PvpT->SetLineColor(kRed-1);
PvpT->SetLineWidth(2);
PvpT->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
PvpT->SetStats(0);
PvpT->Draw();
NvpT->SetLineColor(kBlue-1);
NvpT->SetLineWidth(2);
NvpT->Draw("same");  
auto L1 = new TLegend();
L1->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
L1->AddEntry(PvpT,Form("AVFD, pos, Cent %d0-%d0",cen,cen+1), "l");
L1->AddEntry(NvpT,Form("AVFD, neg, Cent %d0-%d0",cen,cen+1), "l");
L1->Draw();
c->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.5);
pad2->SetTopMargin(0);
pad2->SetLeftMargin(0.2);
pad2->SetRightMargin(4.);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd(); 
DvpT->GetXaxis()->SetRangeUser(0.,5.);
DvpT->SetLineColor(kGreen-1);
DvpT->SetLineWidth(2);
DvpT->GetXaxis()->SetTitle("p_{T} [GeV]");
DvpT->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
DvpT->Draw();
c->SaveAs(Form("v%d_pT_cen%d.pdf",harm+1,cen));

TCanvas *cEta = new TCanvas("ceta", "canvas", 800, 800);//upper plot in pad1
TPad *padeta = new TPad("padeta", "padeta", 0, 0.5, 1., 1.0);
padeta->SetBottomMargin(0); // Upper and lower plot are joined
padeta->SetLeftMargin(0.2);
padeta->Draw();             // Draw the upper pad: pad1
padeta->cd();               // pad1 becomes the current pad//and filled with pos
PvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
PvEta->SetLineColor(kRed+1);
PvEta->SetLineWidth(2);
PvEta->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
PvEta->SetStats(0);
PvEta->Draw();
NvEta->SetLineColor(kBlue+1);
NvEta->SetLineWidth(2);
NvEta->Draw("same");
L1->Draw();
cEta->cd();          // Go back to the main canvas before defining pad2
TPad *padeta2 = new TPad("padeta2", "padeta2", 0, 0., 1, 0.5);
padeta2->SetTopMargin(0);
padeta2->SetLeftMargin(0.2);
padeta2->SetRightMargin(4.);
padeta2->SetBottomMargin(0.2);
padeta2->SetGridx(); // vertical grid
padeta2->Draw();
padeta2->cd();
DvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
DvEta->GetXaxis()->SetTitle("#eta");
DvEta->SetLineColor(kGreen+1);
DvEta->SetLineWidth(2);
DvEta->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
DvEta->Draw();
cEta->SaveAs(Form("v%d_Eta_cen%d.pdf",harm+1,cen));
}//Harm loop
}//Cent loop
//delta hists for {1,2,3} and {4,5,6} in eta and pT for different centralities
/*
for(Int_t hr = 1; hr<4;hr ++){

for(Int_t ce = 1; ce<4; ce++){

TFile* fi = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-0/new/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",ce,ce+1));
TList *li = (TList*) fi->Get("FlowQCList;1");

//single canvas plot delta vpT
TCanvas *DPT = new TCanvas("DPT","dpT",400,400);
DPT->SetLeftMargin(0.2);
v1->GetYaxis() -> SetTitleOffset(2.0);
v1->GetYaxis()->SetTitle("differential flow v_{1}");
v1->GetXaxis()->SetTitle("p_{T} [GeV]");
v1->GetXaxis()->SetRangeUser(0.5,3);
v1->SetLineColor(kRed);
v1->Draw("same");
auto L2 = new TLegend(0.2,0.7,0.48,0.9);
L2->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
L2->AddEntry(PvpT,Form("AVFD, pos, Cent %d0-%d0%",cen,cen+1), "l");
L2->AddEntry(NvpT,From("AVFD, neg, Cent %d0-%d0%",cen,cen+1), "l");
L2->Draw();

DPT->SaveAs(Form("sq_coverrptv1_%d.pdf",ce));
*/


//}//Cent LOOP
//}//Harm LOOP



}//end void
