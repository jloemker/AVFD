//add ratio & difference plot (and the subfile check error check here)

void RatioPlot(){
TFile *Published = new TFile("/project/alice/users/jlomker/AVFD/plots/published.root");
//the vocab. is a bit unhandy ...
//
//producing new file for 5.02TeV with corresponding centrality
//TFile* pT0eta0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/new/Result_5.02TeV_pT_0_eta_0_Cent20_30.root");//
TFile *Pb = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.3/BField0.6/Result_5.02TeV_pT_0_eta_0_Cent30_40.root");//

TKey *p = (TKey*) Published->Get("Table 5;1");
///TList *eta0pT0 = (TList*) pT0eta0->Get("FlowQCList;1");
TList *pb = (TList*) Pb->Get("FlowQCList;1");

//20 to 30 % [positive = 0][centrality = 2][harmonic = 1 = (v2)][0]
TH1F *h1 = (TH1F*) pb->FindObject("fFlowQCFinalPtDifHist[1][3][1][0]");//check 1 in full range->has 36 Entries 
//TH1F *h2 =  (TH1F*) pb->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
/*
// The Xe-Xe data 5.44TeV from ALICE 20-30%
Double_t h2_x[18] = {0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.125,1.375,1.625,1.875,2.125,2.375,2.75,3.25,3.75,4.5};
Double_t h2_y[18] = {0.0336282, 0.0442939, 0.0563404,0.067766, 0.0780708, 0.0879414, 0.0966931, 0.105392, 0.120248, 0.138635,0.154528, 0.166843, 0.179147,0.187231,0.194451 ,0.192949,0.18926,0.170559};
Double_t h2_err[18] = {0.00149402, 0.001888284, 0.002114885,0.002380251, 0.002362924, 0.001851419, 0.002107053, 0.002307215, 0.002392648, 0.00281213, 0.003222799, 0.003611608,0.00404789, 0.00447, 0.00454807,0.00527164, 0.00623852, 0.00644067};
*/
// The Pb-Pb data 5.02TeV from ALICE for 30-40%
Double_t h2_x[13] = {0.3,0.5,0.7,0.9,1.25,1.375,1.625,1.875,2.25,2.75,3.25,3.75,4.5};
Double_t h2_y[13] = {0.0477,0.0771,0.1024,0.1249,0.1695,0.192,0.2068,0.2225,0.2395,0.2494,0.2448,0.198};
Double_t h2_err[13] = {0.0021,0.003,0.0038,0.0053,0.0032,0.0039,0.0058,0.0047,0.0049,0.007,0.0057,0.0126,0.016};
//adjust fitting range and try lower order pol !
TH1F *h2 = new TH1F("ALICE Fit: Pol 6","ALICE vs AVFD",12,0.3,4.5);
TH1F *hErr = new TH1F("ErrFit","ErrFit",12, 0.3,4.5);
for(int i=1; i < 13; ++i) {
   h2->SetBinContent(i,h2_y[i-1]);
   h2->SetBinError(i,h2_err[i-1]);
   hErr->SetBinContent(i,h2_y[i-1]+h2_err[i-1]);
}

//fit hist 2 with poln and then difference from fit (poln) to my data 
TF1 *func = new TF1("func","pol 6",0,5);
//TF1 *funcErr = new TF1("funcErr","pol 3",0.,5);

h2->Fit("func","R");
TF1 *func_res2 = h2->GetFunction("func");
Double_t chi2_h2 = func_res2->GetChisquare();
Double_t NDF = func_res2 ->GetNDF();
cout<<"h2 ALICE: chi2/ndf "<<chi2_h2/NDF<<endl;

//gStyle->SetOptFit(100);
//gStyle->SetStatX(0.6);
//gStyle->SetStatW(0.2);
//gStyle->SetLabelOffset(0.1);
//gStyle->SetLabelFont(1);

TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1., 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetLeftMargin(0.2);
pad1->SetGridx();         // Vertical grid
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pad//h1->SetStats(0);
h1->GetXaxis()->SetRangeUser(0.,4.5);
h1->GetYaxis()->SetRangeUser(0.,0.4);
h1->SetLineColor(kGreen-1);
h1->SetLineWidth(2);
//h1->GetYaxis()->SetTitleOffset(2);
h1->GetYaxis()->SetTitle("differential flow v_{2}");
h1->SetStats(0);
h1->Draw();
h2->SetLineColor(kBlue-1);
h2->SetLineWidth(2);
h2->Draw("same");  
auto L_Ratio = new TLegend(0.2,0.7,0.48,0.9);
L_Ratio->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
L_Ratio->AddEntry(h1,"Neg, Centrality 30-40%", "l");
//L_Ratio->AddEntry(h2,"Pb-Pb 5.02TeV pT: 0.2 - 5 GeV, #eta: 0.8 ","l");
L_Ratio->AddEntry(h2,"ALICE,charged, Centrality 30-40%", "l");
L_Ratio->Draw();
/*TGaxis *axis = new TGaxis(-5, 20, -5, 220, 20,220,510,"");
axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
axis->SetLabelSize(15);
axis->Draw();//lower plot in pad
*/
c->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
pad2->SetLeftMargin(0.2);
//pad2->SetRightMargin(4.);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();  
TH1F *h3 = (TH1F*)h1->Clone("h3");
h3->SetLineColor(kBlack);
h3->GetXaxis()->SetRangeUser(0.,4.5);
h3->SetStats(0);      // No statistics on lower plot
//Error Calculation/Propagation for lower plot

Int_t NBins = h1->GetNbinsX();
Double_t arr1[NBins];
Double_t arr2[NBins];
Double_t arr3[NBins];
Double_t ErrCorr[NBins];
Double_t Err1[NBins];
Double_t Err2[NBins];
hErr->Fit("func","WLR");//Fitting errors
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
//h1->GetYaxis()->SetTitleOffset(0.2);//h2 setting
h2->SetLineColor(kRed);
h2->SetLineWidth(2);// Ratio plot (h3) settings
h3->SetTitle(""); // Remove the ratio title
h3->GetYaxis()->SetTitle("AVFD - Fit");
h3->GetXaxis()->SetTitle("p_{T} [GeV]");
h3->GetYaxis()->SetRangeUser(-0.2,0.25);
h3->GetYaxis()->SetNdivisions(505);
h3->GetYaxis()->SetTitleSize(18);
h3->GetYaxis()->SetTitleFont(43);
h3->GetYaxis()->SetTitleOffset(2.);
h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetYaxis()->SetLabelSize(20);// X axis ratio plot settings
h3->GetXaxis()->SetTitleSize(20);
h3->GetXaxis()->SetTitleFont(43);
h3->GetXaxis()->SetTitleOffset(3.1);//otherwise 2.9 was fine
h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
h3->GetXaxis()->SetLabelSize(18);

c->SaveAs("difference/Fit_NegcompareB0.6_tau_init0.3.pdf");
}
// 
