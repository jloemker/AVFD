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
//Cent30_40_5.
TFile* high = new TFile("/project/alice/users/jlomker/AVFD/test/dirID-0/pos/eta_pt/Analysis_pTrange_0_eta_0_Cent20_30.root");
//Getting the list(s) of files for all
TList *qa = (TList*) file_all->Get("QAList;1");//Kinematics from POI`s
TList *flowQC = (TList*) file_all->Get("FlowQCList;1");//Matrices for calculations ?
TList *cme = (TList*) file_all->Get("CMEList;1");//(SS/OS/Delta) delta, gamma 
//for pos particles 
TList *posFlowQC = (TList*) file_pos->Get("FlowQCList;1");//Matrices for calculations ?
TList *higher = (TList*) high->Get("FlowQCList;1");

//TH1F *v2_all = (TH1F*) flowQC->FindObject("fFlowQCFinalEtaDifHist[2][0][0]");
TH1F *v1 = (TH1F*) higher->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *v1e = (TH1F*) higher->FindObject("fFlowQCFinalEtaDifHist[0][2][1][0]");
//TH1F *v3_all = (TH1F*) flowQC->FindObject("fFlowQCFinalEtaDifHist[2][1][0]");
TH1F *v3_pos = (TH1F*) posFlowQC->FindObject("fFlowQCFinalEtaDifDeltaHist[2][1]");

TCanvas *V2 = new TCanvas("v2","pT",400,400);
V2->SetLeftMargin(0.2);
v1->GetYaxis() -> SetTitleOffset(2.0);
v1->GetYaxis()->SetRangeUser(-0.3,0.3);
v1->GetXaxis()->SetRangeUser(0,5);
v1->SetLineColor(kRed);
v1->Draw("same");
auto legend = new TLegend();
legend->SetHeader("Small sample Xe-Xe: 5.44TeV, 20-30%","C");
V2 ->SaveAs("testptv2.pdf");


TCanvas *V1 = new TCanvas("v2","pT",400,400);
V1->SetLeftMargin(0.2);
v1e->GetYaxis() -> SetTitleOffset(2.0);
v1e ->GetYaxis()->SetTitle("differential flow v_{n}");
v1e->GetXaxis()->SetTitle("#eta");
v1e->GetYaxis()->SetRangeUser(-0.3,0.3);
v1e->GetXaxis()->SetRangeUser(-1,1);
v1e->SetLineColor(kRed);
v1e->Draw("same");
V1 ->SaveAs("testetav2.pdf");

////**********************************************************/////
////*********************************************************/////
//							    /////
//		This is useful 				   ////
//********************************************************///
//*******************************************************//		

//Plots for large sample in small eta window [-0.8,0.8]
TFile* pT0eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_0_eta_0_Cent20_30.root");//
TFile* pT1eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_1_eta_0_Cent20_30.root");//
TFile* pT2eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_2_eta_0_Cent20_30.root");//
TFile* pT3eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_3_eta_0_Cent20_30.root");//
TFile* pT4eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_4_eta_0_Cent20_30.root");//
TFile* pT5eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_5_eta_0_Cent20_30.root");//
//Plots for large sample in large eta window [-3,3]
TFile* pT0eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_0_eta_1_Cent20_30.root");//
TFile* pT1eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_1_eta_1_Cent20_30.root");//
TFile* pT2eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_2_eta_1_Cent20_30.root");//
TFile* pT3eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_3_eta_1_Cent20_30.root");//
TFile* pT4eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_4_eta_1_Cent20_30.root");//
TFile* pT5eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_5_eta_1_Cent20_30.root");//
//TList *QA = (TList*) File_all->Get("QAList;1");//Kinematics from POI`s
TList *eta0pT0 = (TList*) pT0eta0->Get("FlowQCList;1");
TList *eta0pT1 = (TList*) pT1eta0->Get("FlowQCList;1");
TList *eta0pT2 = (TList*) pT2eta0->Get("FlowQCList;1");
TList *eta0pT3 = (TList*) pT3eta0->Get("FlowQCList;1");
TList *eta0pT4 = (TList*) pT4eta0->Get("FlowQCList;1");
TList *eta0pT5 = (TList*) pT5eta0->Get("FlowQCList;1");

TList *eta1pT0 = (TList*) pT0eta1->Get("FlowQCList;1");
TList *eta1pT1 = (TList*) pT1eta1->Get("FlowQCList;1");
TList *eta1pT2 = (TList*) pT2eta1->Get("FlowQCList;1");
TList *eta1pT3 = (TList*) pT3eta1->Get("FlowQCList;1");
TList *eta1pT4 = (TList*) pT4eta1->Get("FlowQCList;1");
TList *eta1pT5 = (TList*) pT5eta1->Get("FlowQCList;1");

//Differential plots*******************************************************************
//										      *
//        TH1D *fFlowQCFinal**DifHist[Ncharge][fCRCMaxnCen][fFlowNHarm][fFlowQCNCov]; *
//        TH1D *fFlowQCFinal**DifDeltaHist[fCRCMaxnCen][fFlowNHarm];		      *
//*************************************************************************************

//************
//Plots vs eta 
//************        

///For vn vs eta differential flow for positive 0 particles in 6 pTranges
for(Int_t harm = 0; harm<1; harm++){
Int_t vn = harm+1;
TString Y = Form("differential flow v_{%d}",vn);
for(Int_t charge = 0; charge<1; charge++){
TString DifEtaHist = Form("fFlowQCFinalEtaDifHist[%d][2][%d][0]",charge,harm); 
//vn_etapT
TH1F *vn_00 = (TH1F*) eta0pT0->FindObject(DifEtaHist);
TH1F *vn_01 = (TH1F*) eta0pT1->FindObject(DifEtaHist);
TH1F *vn_02 = (TH1F*) eta0pT2->FindObject(DifEtaHist);
TH1F *vn_03 = (TH1F*) eta0pT3->FindObject(DifEtaHist);
TH1F *vn_04 = (TH1F*) eta0pT4->FindObject(DifEtaHist);
TH1F *vn_05 = (TH1F*) eta0pT5->FindObject(DifEtaHist);

TH1F *vn_10 = (TH1F*) eta1pT0->FindObject(DifEtaHist);
TH1F *vn_11 = (TH1F*) eta1pT1->FindObject(DifEtaHist);
TH1F *vn_12 = (TH1F*) eta1pT2->FindObject(DifEtaHist);
TH1F *vn_13 = (TH1F*) eta1pT3->FindObject(DifEtaHist);
TH1F *vn_14 = (TH1F*) eta1pT4->FindObject(DifEtaHist);
TH1F *vn_15 = (TH1F*) eta1pT5->FindObject(DifEtaHist);

//Histogram for eta in [-0.8,0.8]
TCanvas *Eta0 = new TCanvas("Eta0","Eta0",400,400);
Eta0->SetLeftMargin(0.2);
Eta0->SetBottomMargin(0.2);
vn_02->GetYaxis()->SetTitle(Y);//change n
vn_02->GetXaxis()->SetTitle("#eta");
vn_02->SetLineColor(kGreen+1);
vn_02->GetXaxis()->SetRangeUser(-0.9,0.9);
vn_02->GetXaxis()->SetTitleOffset(2.0);
//vn_00->GetYaxis()->SetRangeUser(-2,2);
vn_02->Draw();
vn_01->SetLineColor(kRed+1);
vn_01->Draw("SAME");
vn_02->SetLineColor(kBlue+1);
vn_02->Draw("SAME");
vn_03->SetLineColor(kViolet+1);
vn_03->Draw("SAME");
vn_04->SetLineColor(kOrange+1);
vn_04->Draw("SAME");
vn_05->SetLineColor(kBlack+1);
vn_05->Draw("SAME");
auto P_eta0 = new TLegend();
P_eta0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_eta0->AddEntry(vn_00,"Positive, pT: 0.2-5 GeV", "l");
P_eta0->AddEntry(vn_01,"Positive, pT: 0.2-0.7 GeV", "l");
P_eta0->AddEntry(vn_02,"Positive, pT: 1-2 GeV","l");
P_eta0->AddEntry(vn_03,"Positive, pT: 1-3 GeV", "l");
P_eta0->AddEntry(vn_04,"Positive, pT: 1.8-2.5 GeV","l");
P_eta0->AddEntry(vn_05,"Positive, pT: 0.2-3 GeV", "l");
auto N_eta0 = new TLegend();
N_eta0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_eta0->AddEntry(vn_00,"Negative, pT: 0.2-5 GeV", "l");
N_eta0->AddEntry(vn_01,"Negative, pT: 0.2-0.7 GeV", "l");
N_eta0->AddEntry(vn_02,"Negative, pT: 1-2 GeV","l");
N_eta0->AddEntry(vn_03,"Negative, pT: 1-3 GeV", "l");
N_eta0->AddEntry(vn_04,"Negative, pT: 1.8-2.5 GeV","l");
N_eta0->AddEntry(vn_05,"Negative, pT: 0.2-3 GeV", "l");
if(charge == 0){
P_eta0->Draw();
TString a = Form("Pos_diffv%d_Eta0_a2_pT.pdf",vn);
Eta0->SaveAs(a);
}
if(charge == 1){
N_eta0->Draw();
TString b = Form("Neg_diffv%d_Eta0_a2_pT.pdf",vn);
Eta0->SaveAs(b);
}

//For eta [-3,3]
TCanvas *Eta1 = new TCanvas("Eta1","Eta1",400,400);
Eta1->SetLeftMargin(0.2);
Eta1->SetBottomMargin(0.2);
vn_12->GetYaxis()->SetTitle(Y);
vn_12->GetXaxis()->SetTitle("#eta");
vn_12->SetLineColor(kGreen+1);
vn_12->GetXaxis()->SetRangeUser(-3.1,3.1);
vn_12->GetXaxis()->SetTitleOffset(2.0);
//vn_10->GetYaxis()->SetRangeUser(-2,2);
vn_12->Draw();
vn_11->SetLineColor(kRed+1);
vn_11->Draw("SAME");
vn_12->SetLineColor(kBlue+1);
vn_12->Draw("SAME");
vn_13->SetLineColor(kViolet+1);
vn_13->Draw("SAME");
vn_14->SetLineColor(kOrange+1);
vn_14->Draw("SAME");
vn_15->SetLineColor(kBlack+1);
vn_15->Draw("SAME");
auto P_eta1 = new TLegend();
P_eta1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_eta1->AddEntry(vn_10,"Positive, pT: 0.2-5 GeV", "l");
P_eta1->AddEntry(vn_11,"Positive, pT: 0.2-0.7 GeV", "l");
P_eta1->AddEntry(vn_12,"Positive, pT: 1-2 GeV","l");
P_eta1->AddEntry(vn_13,"Positive, pT: 1-3 GeV", "l");
P_eta1->AddEntry(vn_14,"Positive, pT: 1.8-2.5 GeV","l");
P_eta1->AddEntry(vn_15,"Positive, pT: 0.2-3 GeV", "l");
P_eta1->Draw();
auto N_eta1 = new TLegend();
N_eta1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_eta1->AddEntry(vn_10,"Negative, pT: 0.2-5 GeV", "l");
N_eta1->AddEntry(vn_11,"Negative, pT: 0.2-0.7 GeV", "l");
N_eta1->AddEntry(vn_12,"Negative, pT: 1-2 GeV","l");
N_eta1->AddEntry(vn_13,"Negative, pT: 1-3 GeV", "l");
N_eta1->AddEntry(vn_14,"Negative, pT: 1.8-2.5 GeV","l");
N_eta1->AddEntry(vn_15,"Negative, pT: 0.2-3 GeV", "l");
if(charge == 0){
P_eta1->Draw();
TString c = Form("Pos_diffv%d_Eta1_a1_pT.pdf",vn);
Eta1->SaveAs(c);
}
if(charge == 1){
N_eta1->Draw();
TString d = Form("Neg_diffv%d_Eta1_a1_pT.pdf",vn);
Eta1->SaveAs(d);
}


//******************
//////Plots vs pT  *
//////************ *
////For vn vs pT differential flow for pos and neg
////in each pt range the vn from pos & neg particle for eta 0 and 1
//
TString DifPtHist = Form("fFlowQCFinalPtDifHist[%d][2][%d][0]",charge,harm);
//    v_etapt
TH1F *v_00 = (TH1F*) eta0pT0->FindObject(DifPtHist);
TH1F *v_01 = (TH1F*) eta0pT1->FindObject(DifPtHist);
TH1F *v_02 = (TH1F*) eta0pT2->FindObject(DifPtHist);
TH1F *v_03 = (TH1F*) eta0pT3->FindObject(DifPtHist);
TH1F *v_04 = (TH1F*) eta0pT4->FindObject(DifPtHist);
TH1F *v_05 = (TH1F*) eta0pT5->FindObject(DifPtHist);

TH1F *v_10 = (TH1F*) eta1pT0->FindObject(DifPtHist);
TH1F *v_11 = (TH1F*) eta1pT1->FindObject(DifPtHist);
TH1F *v_12 = (TH1F*) eta1pT2->FindObject(DifPtHist);
TH1F *v_13 = (TH1F*) eta1pT3->FindObject(DifPtHist);
TH1F *v_14 = (TH1F*) eta1pT4->FindObject(DifPtHist);
TH1F *v_15 = (TH1F*) eta1pT5->FindObject(DifPtHist);
//Pt0
TCanvas *Pt0 = new TCanvas("Pt0","Pt0",400,400);
Pt0->SetLeftMargin(0.2);
Pt0->SetBottomMargin(0.2);
v_00->GetYaxis()->SetTitle(Y);
v_00->GetXaxis()->SetTitle("p_{T} [GeV]");
v_00->SetLineColor(kGreen-1);
v_00->GetXaxis()->SetRangeUser(0.2,5);
v_00->GetXaxis()->SetTitleOffset(2.0);//v_10->GetYaxis()->SetRangeUser(-2,2);
v_00->Draw();
v_10->SetLineColor(kGreen+2);
v_10->Draw("SAME");
auto P_pt0 = new TLegend();
P_pt0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt0->AddEntry(v_00,"Positive, pT: 0.2-5 GeV, #eta: 0.8", "l");
P_pt0->AddEntry(v_10,"Positive, pT: 0.2-5 GeV, #eta: 3", "l");
auto N_pt0 = new TLegend();
N_pt0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt0->AddEntry(vn_00,"Negative, pT: 0.2-5 GeV, #eta: 0.8", "l");
N_pt0->AddEntry(vn_10,"Negative, pT: 0.2-5 GeV, #eta: 3", "l");
if(charge == 0){
P_pt0->Draw();
TString e = Form("Pos_diffv%d_pT_0_Eta01.pdf",vn);
Pt0->SaveAs(e);
}
if(charge == 1){
N_pt0->Draw();
TString f = Form("Neg_diffv%d_pT_0_Eta01.pdf",vn);
Pt0->SaveAs(f);
}
//Pt1
TCanvas *Pt1 = new TCanvas("Pt1","Pt1",400,400);
Pt1->SetLeftMargin(0.2);
Pt1->SetBottomMargin(0.2);
v_01->GetYaxis()->SetTitle(Y);
v_01->GetXaxis()->SetTitle("p_{T} [GeV]");
v_01->SetLineColor(kRed-1);
v_01->GetXaxis()->SetRangeUser(0.2,0.7);
v_01->GetXaxis()->SetTitleOffset(2.0);//v_10->GetYaxis()->SetRangeUser(-2,2);
v_01->Draw();
v_11->SetLineColor(kRed+2);
v_11->Draw("SAME");
auto P_pt1 = new TLegend();
P_pt1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt1->AddEntry(v_01,"Positive, pT: 0.2-0.7 GeV, #eta: 0.8", "l");
P_pt1->AddEntry(v_11,"Positive, pT: 0.2-9.7 GeV, #eta: 3", "l");
auto N_pt1 = new TLegend();
N_pt1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt1->AddEntry(vn_01,"Negative, pT: 0.2-0.7 GeV, #eta: 0.8", "l");
N_pt1->AddEntry(vn_11,"Negative, pT: 0.2-0.7 GeV, #eta: 3", "l");
if(charge == 0){
TString g = Form("Pos_diffv%d_pT_1_Eta01.pdf",vn);
P_pt1->Draw();
Pt1->SaveAs(g);
}
if(charge == 1){
TString h = Form("Neg_diffv%d_pT_1_Eta01.pdf",vn);
N_pt1->Draw();
Pt1->SaveAs(h);
}
//Pt2
TCanvas *Pt2 = new TCanvas("Pt2","Pt2",400,400);
Pt2->SetLeftMargin(0.2);
Pt2->SetBottomMargin(0.2);
v_02->GetYaxis()->SetTitle(Y);
v_02->GetXaxis()->SetTitle("p_{T} [GeV]");
v_02->SetLineColor(kBlue-1);
v_02->GetXaxis()->SetRangeUser(1,2);
v_02->GetXaxis()->SetTitleOffset(2.0);//v_10->GetYaxis()->SetRangeUser(-2,2);
v_02->Draw();
v_12->SetLineColor(kBlue+2);
v_12->Draw("SAME");
auto P_pt2 = new TLegend();
P_pt2->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt2->AddEntry(v_02,"Positive, pT: 1-2 GeV, #eta: 0.8", "l");
P_pt2->AddEntry(v_12,"Positive, pT: 1-2 GeV, #eta: 3", "l");
auto N_pt2 = new TLegend();
N_pt2->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt2->AddEntry(vn_02,"Negative, pT: 1-2 GeV, #eta: 0.8", "l");
N_pt2->AddEntry(vn_12,"Negative, pT: 1-2 GeV, #eta: 3", "l");
if(charge == 0){
P_pt2->Draw();
TString i = Form("Pos_diffv%d_pT_2_Eta01.pdf",vn);
Pt2->SaveAs(i);
}
if(charge == 1){
TString j = Form("Neg_diffv%d_pT_2_Eta01.pdf",vn);
N_pt2->Draw();
Pt2->SaveAs(j);
}

TCanvas *Pt3 = new TCanvas("Pt3","Pt3",400,400);
Pt3->SetLeftMargin(0.2);
Pt3->SetBottomMargin(0.2);
v_03->GetYaxis()->SetTitle(Y);
v_03->GetXaxis()->SetTitle("p_{T} [GeV]");
v_03->SetLineColor(kViolet-1);
v_03->GetXaxis()->SetRangeUser(1,3.1);
v_03->GetXaxis()->SetTitleOffset(2.0);//v_10->GetYaxis()->SetRangeUser(-2,2);
v_03->Draw();
v_13->SetLineColor(kViolet+2);
v_13->Draw("SAME");
auto P_pt3 = new TLegend();
P_pt3->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt3->AddEntry(v_03,"Positive, pT: 1-3 GeV, #eta: 0.8", "l");
P_pt3->AddEntry(v_13,"Positive, pT: 1-3 GeV, #eta: 3", "l");
auto N_pt3 = new TLegend();
N_pt3->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt3->AddEntry(vn_03,"Negative, pT: 1-3 GeV, #eta: 0.8", "l");
N_pt3->AddEntry(vn_13,"Negative, pT: 1-3 GeV, #eta: 3", "l");
if(charge == 0){
TString k = Form("Pos_diffv%d_pT_3_Eta01.pdf",vn);
P_pt3->Draw();
Pt3->SaveAs(k);
}
if(charge == 1){
TString l = Form("Neg_diffv%d_pT_3_Eta01.pdf",vn);
N_pt3->Draw();
Pt3->SaveAs(l);
}

TCanvas *Pt4 = new TCanvas("Pt4","Pt4",400,400);
Pt4->SetLeftMargin(0.2);
Pt4->SetBottomMargin(0.2);
v_04->GetYaxis()->SetTitle(Y);
v_04->GetXaxis()->SetTitle("p_{T} [GeV]");
v_04->SetLineColor(kOrange-1);
v_04->GetXaxis()->SetRangeUser(1.8,2.5);
v_04->GetXaxis()->SetTitleOffset(2.0);//v_10->GetYaxis()->SetRangeUser(-2,2);
v_04->Draw();
v_14->SetLineColor(kOrange+2);
v_14->Draw("SAME");
auto P_pt4 = new TLegend();
P_pt4->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt4->AddEntry(v_04,"Positive, pT: 1.8-2.5 GeV, #eta: 0.8", "l");
P_pt4->AddEntry(v_14,"Positive, pT: 1.8-2.5 GeV, #eta: 3", "l");
auto N_pt4 = new TLegend();
N_pt4->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt4->AddEntry(vn_04,"Negative, pT: 1.8-2.5 GeV, #eta: 0.8", "l");
N_pt4->AddEntry(vn_14,"Negative, pT: 1.8-2.5 GeV, #eta: 3", "l");
if(charge == 0){
P_pt4->Draw();
TString m = Form("Pos_diffv%d_pT_4_Eta01.pdf",vn);
Pt4->SaveAs(m);
}
if(charge == 1){
TString n = Form("Neg_diffv%d_pT_4_Eta01.pdf",vn);
N_pt4->Draw();
Pt4->SaveAs(n);
}

TCanvas *Pt5 = new TCanvas("Pt5","Pt5",400,400);
Pt5->SetLeftMargin(0.2);
Pt5->SetBottomMargin(0.2);
v_05->GetYaxis()->SetTitle(Y);
v_05->GetXaxis()->SetTitle("p_{T} [GeV]");
v_05->SetLineColor(kPink+2);
v_05->GetXaxis()->SetRangeUser(0.2,3);
v_05->GetXaxis()->SetTitleOffset(2.0);//v_05->GetYaxis()->SetRangeUser(0.2,3);
v_05->Draw();
v_15->SetLineColor(kBlack);
v_15->Draw("SAME");
auto P_pt5 = new TLegend();
P_pt5->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
P_pt5->AddEntry(v_03,"Positive, pT: 0.2-3 GeV, #eta: 0.8", "l");
P_pt5->AddEntry(v_13,"Positive, pT: 0.2-3 GeV, #eta: 3", "l");
auto N_pt5 = new TLegend();
N_pt5->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
N_pt5->AddEntry(vn_05,"Negative, pT: 0.2-3 GeV, #eta: 0.8", "l");
N_pt5->AddEntry(vn_15,"Negative, pT: 0.2-3 GeV, #eta: 3", "l");
if(charge == 0){
TString o = Form("Pos_diffv%d_pT_5_Eta01.pdf",vn);
P_pt5->Draw();
Pt5->SaveAs(o);
}
if(charge == 1){
TString p = Form("Neg_diffv%d_pT5_Eta01_all_pT.pdf",vn);
N_pt5->Draw();
Pt5->SaveAs(p);
}


}//charge

//For differential \Delta vn vs Eta
TString DelY = Form("differential flow #Delta v_{%d} = v_{%d}(+q)-v_{%d}(-q)",vn,vn,vn);
TString DifDelEtaHist = Form("fFlowQCFinalEtaDifDeltaHist[2][%d]",harm);
//        vn_etapT
TH1F *Del_vn_00 = (TH1F*) eta0pT0->FindObject(DifDelEtaHist);
TH1F *Del_vn_01 = (TH1F*) eta0pT1->FindObject(DifDelEtaHist);
TH1F *Del_vn_02 = (TH1F*) eta0pT2->FindObject(DifDelEtaHist);
TH1F *Del_vn_03 = (TH1F*) eta0pT3->FindObject(DifDelEtaHist);
TH1F *Del_vn_04 = (TH1F*) eta0pT4->FindObject(DifDelEtaHist);
TH1F *Del_vn_05 = (TH1F*) eta0pT5->FindObject(DifDelEtaHist);

TH1F *Del_vn_10 = (TH1F*) eta1pT0->FindObject(DifDelEtaHist);
TH1F *Del_vn_11 = (TH1F*) eta1pT1->FindObject(DifDelEtaHist);
TH1F *Del_vn_12 = (TH1F*) eta1pT2->FindObject(DifDelEtaHist);
TH1F *Del_vn_13 = (TH1F*) eta1pT3->FindObject(DifDelEtaHist);
TH1F *Del_vn_14 = (TH1F*) eta1pT4->FindObject(DifDelEtaHist);
TH1F *Del_vn_15 = (TH1F*) eta1pT5->FindObject(DifDelEtaHist);
//small eta range
TCanvas *DelEta0 = new TCanvas("DelEta0","DelEta0",400,400);
DelEta0->SetLeftMargin(0.2);
DelEta0->SetBottomMargin(0.2);
Del_vn_00->GetYaxis()->SetTitle(DelY);
Del_vn_00->GetXaxis()->SetTitle("#eta");
Del_vn_00->SetLineColor(kGreen+3);
Del_vn_00->GetXaxis()->SetRangeUser(-0.9,0.9);
Del_vn_00->GetXaxis()->SetTitleOffset(2.0);
Del_vn_00->Draw();

Del_vn_01->SetLineColor(kRed);
Del_vn_01->Draw("SAME");
Del_vn_02->SetLineColor(kBlue+1);
Del_vn_02->Draw("SAME");

Del_vn_03->SetLineColor(kViolet+1);
Del_vn_03->Draw("SAME");
Del_vn_04->SetLineColor(kOrange-1);
Del_vn_04->Draw("SAME");
Del_vn_05->SetLineColor(kRed+3);
Del_vn_05->Draw("SAME");

auto Del_E0 = new TLegend();
Del_E0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_E0->AddEntry(Del_vn_00,"pT: 0.2-5 GeV", "l");
Del_E0->AddEntry(Del_vn_01,"pT: 0.2-0.7 GeV", "l");
Del_E0->AddEntry(Del_vn_02,"pT: 1-2 GeV","l");
Del_E0->AddEntry(Del_vn_03,"pT: 1-3 GeV", "l");
Del_E0->AddEntry(Del_vn_04,"pT: 1.8-2.5 GeV","l");
Del_E0->AddEntry(Del_vn_05,"pT: 0.2-3 GeV", "l");
Del_E0->Draw();
TString q = Form("Delta_diffv%d_Eta0_a1_pT.pdf",vn);
DelEta0->SaveAs(q);
//larger eta range
TCanvas *DelEta1 = new TCanvas("DelEta1","DelEta1",400,400);
DelEta1->SetLeftMargin(0.2);
DelEta1->SetBottomMargin(0.2);
Del_vn_10->GetYaxis()->SetTitle(DelY);
Del_vn_10->GetXaxis()->SetTitle("#eta");
Del_vn_10->SetLineColor(kGreen+3);
Del_vn_10->GetXaxis()->SetRangeUser(-3.1,3.1);
Del_vn_10->GetXaxis()->SetTitleOffset(2.0);
Del_vn_10->Draw();

Del_vn_11->SetLineColor(kRed);
Del_vn_11->Draw("SAME");
Del_vn_12->SetLineColor(kBlue+1);
Del_vn_12->Draw("SAME");

Del_vn_13->SetLineColor(kViolet+1);
Del_vn_13->Draw("SAME");
Del_vn_14->SetLineColor(kOrange-1);
Del_vn_14->Draw("SAME");
Del_vn_15->SetLineColor(kRed+3);
Del_vn_15->Draw("SAME");
auto Del_E1 = new TLegend();
Del_E1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_E1->AddEntry(Del_vn_10,"pT: 0.2-5 GeV", "l");
Del_E1->AddEntry(Del_vn_11,"pT: 0.2-0.7 GeV", "l");
Del_E1->AddEntry(Del_vn_12,"pT: 1-2 GeV","l");
Del_E1->AddEntry(Del_vn_13,"pT: 1-3 GeV", "l");
Del_E1->AddEntry(Del_vn_14,"pT: 1.8-2.5 GeV","l");
Del_E1->AddEntry(Del_vn_15,"pT: 0.2-3 GeV", "l");
Del_E1->Draw();
TString r = Form("Delta_diffv%d_Eta1_a1_pT.pdf",vn);
DelEta1->SaveAs(r);

//************
////Plots vs pT
////************
//
///*
// *
// //For vn vs pT differential flow for pos and neg
// //in each pt range the vn from pos & neg particle for eta 0 and 1
//
//for \Delta vn vs pt for the 2 eta ranges

TString DifDelPtHist = Form("fFlowQCFinalPtDifDeltaHist[2][%d]",harm);
//        vn_etapT
TH1F *Del_v_00 = (TH1F*) eta0pT0->FindObject(DifDelPtHist);
TH1F *Del_v_01 = (TH1F*) eta0pT1->FindObject(DifDelPtHist);
TH1F *Del_v_02 = (TH1F*) eta0pT2->FindObject(DifDelPtHist);
TH1F *Del_v_03 = (TH1F*) eta0pT3->FindObject(DifDelPtHist);
TH1F *Del_v_04 = (TH1F*) eta0pT4->FindObject(DifDelPtHist);
TH1F *Del_v_05 = (TH1F*) eta0pT5->FindObject(DifDelPtHist);

TH1F *Del_v_10 = (TH1F*) eta1pT0->FindObject(DifDelPtHist);
TH1F *Del_v_11 = (TH1F*) eta1pT1->FindObject(DifDelPtHist);
TH1F *Del_v_12 = (TH1F*) eta1pT2->FindObject(DifDelPtHist);
TH1F *Del_v_13 = (TH1F*) eta1pT3->FindObject(DifDelPtHist);
TH1F *Del_v_14 = (TH1F*) eta1pT4->FindObject(DifDelPtHist);
TH1F *Del_v_15 = (TH1F*) eta1pT5->FindObject(DifDelPtHist);

//Pt0
TCanvas *DelPt0 = new TCanvas("DelPt0","DelPt0",400,400);
DelPt0->SetLeftMargin(0.2);
DelPt0->SetBottomMargin(0.2);
Del_v_00->GetYaxis()->SetTitle(DelY);
Del_v_00->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_00->SetLineColor(kGreen+3);
Del_v_00->GetXaxis()->SetRangeUser(0.2,5);
Del_v_00->GetXaxis()->SetTitleOffset(2.0);
Del_v_00->Draw();
Del_v_10->SetLineColor(kGreen-2);
Del_v_10->Draw("SAME");
auto Del_Pt0 = new TLegend();
Del_Pt0->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt0->AddEntry(Del_v_00,"pT: 0.2-5 GeV, eta: 0.8", "l");
Del_Pt0->AddEntry(Del_v_10,"pT: 0.2-5 GeV, eta: 3", "l");
Del_Pt0->Draw();
TString s = Form("Delta_diffv%d_pT_0_E01.pdf",vn);
DelPt0->SaveAs(s);
//Pt1
TCanvas *DelPt1 = new TCanvas("DelPt1","DelPt1",400,400);
DelPt1->SetLeftMargin(0.2);
DelPt1->SetBottomMargin(0.2);
Del_v_01->GetYaxis()->SetTitle(DelY);
Del_v_01->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_01->SetLineColor(kRed+3);
Del_v_01->GetXaxis()->SetRangeUser(-0.2,0.7);
Del_v_01->GetXaxis()->SetTitleOffset(2.0);
Del_v_01->Draw();
Del_v_11->SetLineColor(kRed-2);
Del_v_11->Draw("SAME");
auto Del_Pt1 = new TLegend();
Del_Pt1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt1->AddEntry(Del_v_01,"pT: 0.2-0.7 GeV, eta: 0.8", "l");
Del_Pt1->AddEntry(Del_v_11,"pT: 0.2-0.7 GeV, eta: 3", "l");
Del_Pt1->Draw();
TString t = Form("Delta_diffv%d_pT_1_E01.pdf",vn);
DelPt0->SaveAs(t);
//Pt2
TCanvas *DelPt2 = new TCanvas("DelPt2","DelPt2",400,400);
DelPt0->SetLeftMargin(0.2);
DelPt0->SetBottomMargin(0.2);
Del_v_02->GetYaxis()->SetTitle(DelY);
Del_v_02->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_02->SetLineColor(kBlue+3);
Del_v_02->GetXaxis()->SetRangeUser(1.,2.1);
Del_v_02->GetXaxis()->SetTitleOffset(2.0);
Del_v_02->Draw();
Del_v_12->SetLineColor(kBlue-2);
Del_v_12->Draw("SAME");
auto Del_Pt2 = new TLegend();
Del_Pt2->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt2->AddEntry(Del_v_02,"pT: 1-2 GeV, eta: 0.8", "l");
Del_Pt2->AddEntry(Del_v_12,"pT: 1-2 GeV, eta: 3", "l");
Del_Pt2->Draw();
TString u = Form("Delta_diffv%d_pT_2_E01.pdf",vn);
DelPt2->SaveAs(u);
//Pt3
TCanvas *DelPt3 = new TCanvas("DelPt3","DelPt3",400,400);
DelPt3->SetLeftMargin(0.2);
DelPt3->SetBottomMargin(0.2);
Del_v_03->GetYaxis()->SetTitle(DelY);
Del_v_03->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_03->SetLineColor(kViolet+3);
Del_v_03->GetXaxis()->SetRangeUser(1,3.1);
Del_v_03->GetXaxis()->SetTitleOffset(2.0);
Del_v_03->Draw();
Del_v_13->SetLineColor(kViolet-2);
Del_v_13->Draw("SAME");
auto Del_Pt3 = new TLegend();
Del_Pt3->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt3->AddEntry(Del_v_03,"pT: 1-3 GeV, eta: 0.8", "l");
Del_Pt3->AddEntry(Del_v_13,"pT: 1-3 GeV, eta: 3", "l");
Del_Pt3->Draw();
TString v = Form("Delta_diffv%d_pT_3_E01.pdf",vn);
DelPt3->SaveAs(v);
//Pt4
TCanvas *DelPt4 = new TCanvas("DelPt4","DelPt4",400,400);
DelPt4->SetLeftMargin(0.2);
DelPt4->SetBottomMargin(0.2);
Del_v_04->GetYaxis()->SetTitle(DelY);
Del_v_04->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_04->SetLineColor(kOrange+3);
Del_v_04->GetXaxis()->SetRangeUser(1.8,2.5);
Del_v_04->GetXaxis()->SetTitleOffset(2.0);
Del_v_04->Draw();
Del_v_14->SetLineColor(kOrange-2);
Del_v_14->Draw("SAME");
auto Del_Pt4 = new TLegend();
Del_Pt4->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt4->AddEntry(Del_v_04,"pT: 1.8-2.5 GeV, eta: 0.8", "l");
Del_Pt4->AddEntry(Del_v_14,"pT: 1.8-2.5 GeV, eta: 3", "l");
Del_Pt4->Draw();
TString w = Form("Delta_diffv%d_pT_4_E01.pdf",vn);
DelPt4->SaveAs(w);
//Pt5
TCanvas *DelPt5 = new TCanvas("DelPt5","DelPt5",400,400);
DelPt5->SetLeftMargin(0.2);
DelPt5->SetBottomMargin(0.2);
Del_v_05->GetYaxis()->SetTitle(DelY);
Del_v_05->GetXaxis()->SetTitle("p_{T} [GeV]");
Del_v_05->SetLineColor(kBlack+1);
Del_v_05->GetXaxis()->SetRangeUser(0.2,3);
Del_v_05->GetXaxis()->SetTitleOffset(2.0);
Del_v_05->Draw();
Del_v_15->SetLineColor(kPink-1);
Del_v_15->Draw("SAME");
auto Del_Pt5 = new TLegend();
Del_Pt5->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Del_Pt5->AddEntry(Del_v_05,"pT: 0.2-3 GeV, eta: 0.8", "l");
Del_Pt5->AddEntry(Del_v_15,"pT: 0.2-3 GeV, eta: 3", "l");
Del_Pt5->Draw();
TString x = Form("Delta_diffv%d_pT_5_E01.pdf",vn);
DelPt5->SaveAs(x);

}//harm

//Next: clear the trash below and make finish the difference plot

/*
TCanvas *V1_1 = new TCanvas("V1_1","Eta1_1",400,400);
V1_1->SetLeftMargin(0.2);
V1_1->SetBottomMargin(0.2);
v1_small->GetYaxis()->SetTitle("differential flow v_{1}");//change n
v1_small->GetXaxis()->SetTitle("#eta");
v1_small->GetXaxis()->SetRangeUser(-0.9,0.9);
v1_small->GetXaxis()->SetTitleOffset(2.0);
v1_small->SetLineColor(kGreen+1);
v1_small->Draw();
auto Legend1_1 = new TLegend();
Legend1_1->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");TLegend1_1->AddEntry(v1_small,"Positive, pT: 0.2-3 GeV");
Legend1_1->Draw();
//V1_1->SaveAs("Pos_diffv1_eta_pT1.pdf");//change and vn

TCanvas *V12 = new TCanvas("V12","Eta2",400,400);
V12->SetLeftMargin(0.2);
V12->SetBottomMargin(0.2);
v1_large->GetYaxis()->SetTitle("differential flow v_{2}");//change n
v1_large->GetXaxis()->SetTitle("#eta");
v1_large->GetXaxis()->SetRangeUser(-0.9,0.9);
v1_large->GetXaxis()->SetTitleOffset(2.0);
v1_large->SetLineColor(kBlue);
v1_large->Draw();
auto Legend12 = new TLegend();
Legend12->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend12->AddEntry(v1_large,"Positive, pT: 3.0-5 GeV");
Legend12->Draw();
//V12->SaveAs("Pos_diffv2_eta_pT2.pdf");//change vn


///For vn vs eta differential flow for negative 1 particles in 3 pTranges
TH1F *v11_full = (TH1F*) full->FindObject("fFlowQCFinalEtaDifHist[1][2][1][0]");
TH1F *v11_small = (TH1F*) small->FindObject("fFlowQCFinalEtaDifHist[1][2][1][0]");
TH1F *v11_large = (TH1F*) large->FindObject("fFlowQCFinalEtaDifHist[1][2][1][0]");

TCanvas *V11 = new TCanvas("V11","Eta1",400,400);
V11->SetLeftMargin(0.2);
V11->SetBottomMargin(0.2);
v11_full->GetYaxis()->SetTitle("differential flow v_{2}");//change 11_full->GetXaxis()->SetTitle("#eta");
v11_full->SetLineColor(kGreen+3);
v11_full->GetXaxis()->SetRangeUser(-0.9,0.9);
v11_full->GetXaxis()->SetTitleOffset(2.0);
//v11_full->GetYaxis()->SetRangeUser(0.04,0.09);
v11_full->Draw();
v11_small->SetLineColor(kRed);
v11_small->Draw("SAME");
v11_large->SetLineColor(kBlue);
v11_large->Draw("SAME");
auto Legend11 = new TLegend();
Legend11->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend11->AddEntry(v11_full,"Negative, pT: 0.2-5 GeV", "l");
Legend11->AddEntry(v11_small,"Negative, pT: 0.2-1 GeV", "l");
Legend11->AddEntry(v11_large,"Negative, pT: 3.0-5 GeV");
Legend11->Draw();
//V11->SaveAs("Neg_diffv2_eta_pT012.pdf");//change Pos/Neg and vn

TCanvas *V112 = new TCanvas("V112","Eta12",400,400);
V112->SetLeftMargin(0.2);
V112->SetBottomMargin(0.2);
v11_large->GetYaxis()->SetTitle("differential flow v_{2}");//change n
v11_large->GetXaxis()->SetTitle("#eta");
v11_large->SetLineColor(kGreen+3);
v11_large->GetXaxis()->SetRangeUser(-0.9,0.9);
v11_large->GetXaxis()->SetTitleOffset(2.0);
v11_large->SetLineColor(kBlue);
v11_large->Draw();
auto Legend112 = new TLegend();
Legend112->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend112->AddEntry(v11_large,"Negative, pT: 3.0-5 GeV");
Legend112->Draw();
//V112->SaveAs("Neg_diffv2_eta_pT2.pdf");//change Pos/Neg and vn

//direct comparison pos to negative particles in eta
TCanvas *Vpn = new TCanvas("Vpn","Etapn",400,400);
Vpn->SetLeftMargin(0.2);
Vpn->SetBottomMargin(0.2);
v11_full->GetYaxis()->SetTitle("differential flow v_{2}");//change n
v11_full->GetXaxis()->SetTitle("#eta");
v11_full->SetLineColor(kViolet);
v11_full->GetXaxis()->SetRangeUser(-0.9,0.9);
v11_full->GetXaxis()->SetTitleOffset(1.0);
v11_full->GetYaxis()->SetRangeUser(0.04,0.08);
v11_full->Draw();
v1_full->SetLineColor(kOrange+2);
v1_full->GetYaxis()->SetRangeUser(0.01,0.04);
v1_full->Draw("SAME");
auto Legendnp = new TLegend();
Legendnp->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legendnp->AddEntry(v11_full,"Negative, pT: 0.2-5 GeV", "l");
Legendnp->AddEntry(v1_full,"Positive, pT: 0.2-5 GeV", "l");
Legendnp->Draw();
//Vpn->SaveAs("Np_diffv2_eta_pT0.pdf");//change Pos/Neg and vn




///For differental delta vn vs eta in 3 pTranges
TH1F *delv1_full = (TH1F*) full->FindObject("fFlowQCFinalEtaDifDeltaHist[2][1]");
TH1F *delv1_small = (TH1F*) small->FindObject("fFlowQCFinalEtaDifDeltaHist[2][1]");
TH1F *delv1_large = (TH1F*) large->FindObject("fFlowQCFinalEtaDifDeltaHist[2][1]");

TCanvas *delV1 = new TCanvas("DeltaV1","Eta",400,400);
delV1->SetLeftMargin(0.2);
delV1->SetBottomMargin(0.2);
delv1_full->GetYaxis()->SetTitle("differential #Delta v_{2}");//change n
delv1_full->GetXaxis()->SetTitle("#eta");
delv1_full->SetLineColor(kGreen+3);
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
//delV1->SaveAs("Delta_diffv2_eta_pT012.pdf");//change vn

TCanvas *delV12 = new TCanvas("DeltaV12","Eta2",400,400);
delV1->SetLeftMargin(0.2);
delV1->SetBottomMargin(0.3);
delv1_large->GetYaxis()->SetTitle("differential #Delta v_{2}");//change n
delv1_large->GetXaxis()->SetTitle("#eta");
delv1_large->SetLineColor(kGreen+2);
delv1_large->GetXaxis()->SetRangeUser(-0.9,0.9);
delv1_large->GetXaxis()->SetTitleOffset(1);
delv1_large->SetLineColor(kBlue);
delv1_large->Draw();
auto Legend22 = new TLegend();
Legend22->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend22->AddEntry(delv1_large,"pT: 3.0-5 GeV");
Legend22->Draw();
//delV12->SaveAs("Delta_diffv2_eta_pT2.pdf");//change vn


*/



//************
//Plots vs pT
//************

/*
 *
//For vn vs pT differential flow for pos and neg 
//in each pt range the vn from pos & neg particle
//
TH1F *pos_full = (TH1F*) full->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *neg_full = (TH1F*) full->FindObject("fFlowQCFinalPtDifHist[1][2][1][0]");
TH1F *pos_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *neg_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[1][2][1][0]");
TH1F *pos_large = (TH1F*) large->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *neg_large = (TH1F*) large->FindObject("fFlowQCFinalPtDifHist[1][2][1][0]");

//TH1F *pos_small = (TH1F*) small->FindObject("fFlowQCFinalPtDifHist[0][2][0][0]");
TCanvas *Vn = new TCanvas("Vn","Pt",400,400);
Vn->SetLeftMargin(0.2);
Vn->SetBottomMargin(0.2);
pos_full->GetYaxis()->SetTitle("differential flow v_{2}");//change n
pos_full->GetXaxis()->SetTitle("p_{T} [GeV]");
pos_full->SetLineColor(kOrange);
pos_full->GetXaxis()->SetRangeUser(0.2,5.0);
pos_full->GetXaxis()->SetTitleOffset(1.0);
pos_full->Draw();
neg_full->SetLineColor(kViolet);
neg_full->Draw("SAME");
auto Legend3 = new TLegend();
Legend3->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend3->AddEntry(pos_full,"Positive particles", "l");
Legend3->AddEntry(neg_full,"Negative particles", "l");
Legend3->Draw();
//Vn->SaveAs("diffv2_pT_0.pdf");//change Pos/Neg and vn

TCanvas *Vn1 = new TCanvas("Vn1","Pt1",400,400);
Vn1->SetLeftMargin(0.2);
Vn1->SetBottomMargin(0.2);
pos_small->GetYaxis()->SetTitle("differential flow v_{2}");//change n
pos_small->GetXaxis()->SetTitle("p_{T}");
pos_small->SetLineColor(kOrange);
pos_small->GetXaxis()->SetRangeUser(0.2,1.0);
pos_small->GetXaxis()->SetTitleOffset(1.0);
pos_small->Draw();
neg_small->SetLineColor(kViolet);
neg_small->Draw("SAME");
auto Legend4 = new TLegend();
Legend4->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend4->AddEntry(pos_full,"Positive particles", "l");
Legend4->AddEntry(neg_full,"Negative particles", "l");
Legend4->Draw();
//Vn1->SaveAs("diffv2_pT_1.pdf");//change Pos/Neg and vn

TCanvas *Vn2 = new TCanvas("Vn2","Pt2",400,400);
Vn2->SetLeftMargin(0.2);
Vn2->SetBottomMargin(0.2);
pos_large->GetYaxis()->SetTitle("differential flow v_{2}");//change n
pos_large->GetXaxis()->SetTitle("p_{T} [GeV]");
pos_large->SetLineColor(kOrange);
pos_large->GetXaxis()->SetRangeUser(3.0,5.0);
pos_large->GetXaxis()->SetTitleOffset(1.0);
pos_large->Draw();
neg_large->SetLineColor(kViolet);
neg_large->Draw("SAME");
auto Legend5 = new TLegend();
Legend5->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend5->AddEntry(pos_large,"Positive particles", "l");
Legend5->AddEntry(neg_large,"Negative particles", "l");
Legend5->Draw();
//Vn2->SaveAs("diffv2_pT_2.pdf");//change Pos/Neg and vn


//for deltavn vs pT plots
//in each pt range the vn from pos & neg particle: full small large to find objects in pT ranges 0.2-5, 0.2-1 and 3.0 to 5
TH1F *delvn = (TH1F*) full->FindObject("fFlowQCFinalPtDifDeltaHist[2][1]");
TH1F *delvn1 = (TH1F*) small->FindObject("fFlowQCFinalPtDifDeltaHist[2][1]");
TH1F *delvn2 = (TH1F*) large->FindObject("fFlowQCFinalPtDifDeltaHist[2][1]");

TCanvas *delVn = new TCanvas("DeltaVn","Pt",400,400);
delVn->SetLeftMargin(0.2);
delVn->SetBottomMargin(0.2);
delvn->GetYaxis()->SetTitle("differential #Delta v_{2}");//change n
delvn->GetXaxis()->SetTitle("p_{T} [GeV]");
delvn->SetLineColor(kGreen+1);
delvn->GetXaxis()->SetRangeUser(0.2,5.0);
delvn->GetXaxis()->SetTitleOffset(2.0);
delvn->GetXaxis()->SetTitleOffset(1.0);
delvn->Draw();
auto Legend6 = new TLegend();
Legend6->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend6->AddEntry(delvn,"pT: 0.2-5 GeV", "l");
Legend6->Draw();
//delVn->SaveAs("Delta_diffv2_pT_0.pdf");//change vn

TCanvas *delVn1 = new TCanvas("DeltaVn1","Pt1",400,400);
delVn1->SetLeftMargin(0.2);
delVn1->SetBottomMargin(0.2);
delvn1->GetYaxis()->SetTitle("differential #Delta v_{2}");//change n
delvn1->GetXaxis()->SetTitle("p_{T} [GeV]");
delvn1->GetXaxis()->SetRangeUser(0.2,1.0);
delvn1->GetXaxis()->SetTitleOffset(2.0);
delvn1->GetXaxis()->SetTitleOffset(1.0);
delvn1->SetLineColor(kRed+2);
delvn1->Draw();
auto Legend61 = new TLegend();
Legend61->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend61->AddEntry(delvn1,"pT: 0.2-1 GeV", "l");
Legend61->Draw();
//delVn1->SaveAs("Delta_diffv2_pT_1.pdf");//change vn

TCanvas *delVn2 = new TCanvas("DeltaVn2","Pt2",400,400);
delVn2->SetLeftMargin(0.2);
delVn2->SetBottomMargin(0.2);
delvn2->GetYaxis()->SetTitle("differential #Delta v_{2}");//change n
delvn2->GetXaxis()->SetTitle("p_{T} [GeV]");
delvn2->GetXaxis()->SetRangeUser(3,5);
delvn2->GetXaxis()->SetTitleOffset(2.0);
delvn2->GetXaxis()->SetTitleOffset(1.0);
delvn2->SetLineColor(kBlue);
delvn2->Draw();
auto Legend62 = new TLegend();
Legend62->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
Legend62->AddEntry(delvn2,"pT: 3-5 GeV", "l");
Legend62->Draw();
//delVn2->SaveAs("Delta_diffv2_pT_2.pdf");//change vn
*/

/*

TFile *Published = new TFile("/project/alice/users/jlomker/AVFD/plots/published.root");
TList *p = (TList*) Published->Get("Table 5");
 
TH1F *h1 = (TH1F*) p->FindObject("");
TH1F *h2 = (TH1F*) FlowQC->FindObject("fFlowQCFinalPtDifHist[2][0][0]");//define canvas
TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
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
axis->Draw();//lower plot in pad
c->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();       // pad2 becomes the current pad//define raito plot
TH1F *h3 = (TH1F*)h1->Clone("h3");
h3->SetLineColor(kBlack);
h3->SetMinimum(0.8);  // Define Y ..
h3->SetMaximum(1.35); // .. range
h3->Sumw2();
h3->SetStats(0);      // No statistics on lower plot
h3->Divide(h2);
h3->SetMarkerStyle(21);
h3->Draw("ep");       // Draw the ratio plot//h1 settings
h1->SetLineColor(kBlue+1);
h1->SetLineWidth(2);// Y axis h1 plot settings
h1->GetYaxis()->SetTitleSize(20);
h1->GetYaxis()->SetTitleFont(43);
h1->GetYaxis()->SetTitleOffset(1.55);//h2 setting
h2->SetLineColor(kRed);
h2->SetLineWidth(2);// Ratio plot (h3) settings
h3->SetTitle(""); // Remove the ratio title // Y axis ratio plot settings
h3->GetYaxis()->SetTitle("ratio h1/h2 ");
h3->GetYaxis()->SetNdivisions(505);
h3->GetYaxis()->SetTitleSize(20);
h3->GetYaxis()->SetTitleFont(43);
h3->GetYaxis()->SetTitleOffset(1.55);
h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetYaxis()->SetLabelSize(15);// X axis ratio plot settings
h3->GetXaxis()->SetTitleSize(20);
h3->GetXaxis()->SetTitleFont(43);
h3->GetXaxis()->SetTitleOffset(4.);
h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetXaxis()->SetLabelSize(15);

c->SaveAs("ratio.pdf");*/
}

