//add ratio & difference plot (and the subfile check error check here)

void RatioPlot(){
TFile *Published = new TFile("/project/alice/users/jlomker/AVFD/plots/published.root");
//the vocab. is a bit unhandy ...
TFile* pT0eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_0_eta_0_Cent20_30.root");//
TFile* pT1eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_1_eta_0_Cent20_30.root");//
TFile* pT2eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_2_eta_0_Cent20_30.root");//
TFile* pT3eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_3_eta_0_Cent20_30.root");//
TFile* pT4eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_4_eta_0_Cent20_30.root");//
TFile* pT5eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_5_eta_0_Cent20_30.root");//

TFile* pT0eta1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_0_eta_1_Cent20_30.root");//

TKey *p = (TKey*) Published->Get("Table 5;1");
TList *eta0pT0 = (TList*) pT0eta0->Get("FlowQCList;1");
TList *eta0pT1 = (TList*) pT1eta0->Get("FlowQCList;1");
TList *eta0pT2 = (TList*) pT2eta0->Get("FlowQCList;1");
TList *eta0pT3 = (TList*) pT3eta0->Get("FlowQCList;1");
TList *eta0pT4 = (TList*) pT4eta0->Get("FlowQCList;1");
//TList *eta0pT5 = (TList*) pT5eta0->Get("FlowQCList;1");

TList *eta0pT5 = (TList*) pT5eta0->Get("FlowQCList;1");
//20 to 30 % [positive = 0][centrality = 2][harmonic = 1 = (v2)][0]
TH1F *h1 = (TH1F*) eta0pT0->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");//check 1 in full range->has 36 Entries 

Double_t h2_x[18] = {0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,4.5};
Double_t h2_y[18] = {0.0336282, 0.0442939, 0.0563404,0.067766, 0.0780708, 0.0879414, 0.0966931, 0.105392, 0.120248, 0.138635,0.154528, 0.166843, 0.179147,0.187231,0.194451 ,0.192949,0.18926,0.170559};
Double_t h2_err[18] = {0.00149402, 0.001888284, 0.002114885,0.002380251, 0.002362924, 0.001851419, 0.002107053, 0.002307215, 0.002392648, 0.00281213, 0.003222799, 0.003611608,0.00404789, 0.00447, 0.00454807,0.00527164, 0.00623852, 0.00644067};
TH1F *h2 = new TH1F("ALICE Fit: Pol 2","ALICE vs AVFD",18,h2_x);
TH1F *hErr = new TH1F("ErrFit","ErrFit",18,h2_x);
for(int i=1; i < 19; ++i) {
   h2->SetBinContent(i,h2_y[i-1]);
   h2->SetBinError(i,h2_err[i-1]);
   hErr->SetBinContent(i,h2_y[i-1]+h2_err[i-1]);
}
//fit hist 2 with poln and then difference from fit (poln) to my data 
TF1 *func = new TF1("func","pol 2",0.,5);
//TF1 *funcErr = new TF1("funcErr","pol 3",0.,5);

h2->Fit("func","R");
TF1 *func_res2 = h2->GetFunction("func");
Double_t chi2_h2 = func_res2->GetChisquare();
Double_t NDF = func_res2 ->GetNDF();
cout<<"h2 ALICE: chi2/ndf "<<chi2_h2/NDF<<endl;

gStyle->SetOptFit(100);
gStyle->SetStatX(0.6);
gStyle->SetStatW(0.2);
gStyle->SetLabelOffset(1.);
gStyle->SetLabelFont(32);
TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1., 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetLeftMargin(4.);
pad1->SetRightMargin(4.);
pad1->SetGridx();         // Vertical grid
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pad//h1->SetStats(0);
h1->GetXaxis()->SetRangeUser(0.,5);
h1->GetYaxis()->SetRangeUser(0.,0.4);
h1->SetLineColor(kGreen-1);
h1->SetLineWidth(2);
h1->GetYaxis()->SetTitle("differential flow v_{2}");
h1->SetStats(0);
h1->Draw();
h2->SetLineColor(kBlue-1);
h2->SetLineWidth(2);
h2->Draw("same");         // Draw h2 on top of h1
auto L_Ratio = new TLegend(0.1,0.7,0.48,0.9);
L_Ratio->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
L_Ratio->AddEntry(h1,"AVFD,pos, pT: 0.2 - 5 GeV, #eta: 0.8", "l");
L_Ratio->AddEntry(h2,"ALICE,charged, pT: 0.25 - 4.5 GeV, #eta: 0.8", "l");
L_Ratio->Draw();
/* 0 was -5
TGaxis *axis = new TGaxis(-5, 20, -5, 220, 20,220,510,"");
axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
axis->SetLabelSize(15);
axis->Draw();//lower plot in pad
*/
c->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
//pad2->SetRightMargin(4.);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();  
TH1F *h3 = (TH1F*)h1->Clone("h3");
h3->SetLineColor(kBlack);
h3->GetXaxis()->SetRangeUser(0.,5);
h3->SetStats(0);      // No statistics on lower plot
//Error Calculation/Propagation for lower plot
Int_t NBins = h1->GetNbinsX();
Double_t arr1[NBins];
Double_t arr2[NBins];
Double_t arr3[NBins];
Double_t ErrCorr[NBins];
Double_t Err1[NBins];
Double_t Err2[NBins];
hErr->Fit("func","R");//Fitting errors
TF1 *func_resErr = hErr->GetFunction("func");
Double_t chi2_hErr = func_resErr->GetChisquare();
Double_t NDFErr = func_resErr ->GetNDF();
cout<<"hErr ALICE: chi2/ndf "<<chi2_hErr/NDFErr<<endl;
for(Int_t i = 0; i<NBins; i++){
	arr1[i] = h1->GetBinContent(i+1);
	Err1[i] = h1->GetBinError(i+1);
	arr2[i] = func_res2->Eval(h1->GetBinCenter(i+1));
	Err2[i] = func_resErr->Eval(hErr->GetBinCenter(i+1));
	//Abs uncorrelated
	ErrCorr[i] = sqrt(pow(Err1[i],2)+pow(Err2[i],2));
	   if(abs(arr1[i])>0 && abs(arr2[i])>0){
       	   h3->SetBinContent(i+1, arr1[i]-arr2[i]);
	   h3->SetBinError(i+1,ErrCorr[i]);
	   //cout<<"Err1"<<Err1[i]<<" Err2 "<<Err2[i]<<" Err "<<ErrCorr[i]<<endl;
	   }
}
h3->SetMarkerStyle(21);
h3->Draw("ep");       // Draw the ratio plot//h1 settings
h1->SetLineColor(kBlue+1);
h1->SetLineWidth(2);// Y axis h1 plot settings
h1->GetYaxis()->SetTitleSize(20);
h1->GetYaxis()->SetTitleFont(43);
h1->GetYaxis()->SetTitleOffset(1.55);//h2 setting
h2->SetLineColor(kRed);
h2->SetLineWidth(2);// Ratio plot (h3) settings
h3->SetTitle(""); // Remove the ratio title
h3->GetYaxis()->SetTitle("AVFD - Fit");
h3->GetYaxis()->SetRangeUser(-0.2,0.2);
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

c->SaveAs("Fit_compare.pdf");
}
// 
