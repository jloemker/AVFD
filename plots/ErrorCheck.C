
void ErrorCheck(){
TFile* Split0 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_0.root");//
TFile* Split1 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_1.root");//
TFile* Split2 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_2.root");//
TFile* Split3 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_3.root");//
TFile* Split4 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_4.root");//
TFile* Split5 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_5.root");//
TFile* Split6 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_6.root");//
TFile* Split7 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_7.root");//
TFile* Split8 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_8.root");//
TFile* Split9 = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.44TeV_pTrange_0_eta_0_Cent20_30_split_9.root");//

TFile* Full = new TFile("/data/alice/jlomker/AVFD/result/dirID-0/full/Result_5.44TeV_pT_0_eta_0_Cent20_30.root");

TList *split0 = (TList*) Split0->Get("FlowQCList;1");
TList *split1 = (TList*) Split1->Get("FlowQCList;1");
TList *split2 = (TList*) Split2->Get("FlowQCList;1");
TList *split3 = (TList*) Split3->Get("FlowQCList;1");
TList *split4 = (TList*) Split4->Get("FlowQCList;1");
TList *split5 = (TList*) Split5->Get("FlowQCList;1");
TList *split6 = (TList*) Split6->Get("FlowQCList;1");
TList *split7 = (TList*) Split7->Get("FlowQCList;1");
TList *split8 = (TList*) Split8->Get("FlowQCList;1");
TList *split9 = (TList*) Split9->Get("FlowQCList;1");


TList *full = (TList*) Full->Get("FlowQCList;1");
TH1F *hfull = (TH1F*) full->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");

TH1F *h0 = (TH1F*) split0->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h1 = (TH1F*) split1->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h2 = (TH1F*) split2->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h3 = (TH1F*) split3->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h4 = (TH1F*) split4->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h5 = (TH1F*) split5->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h6 = (TH1F*) split6->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h7 = (TH1F*) split7->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h8 = (TH1F*) split8->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");
TH1F *h9 = (TH1F*) split9->FindObject("fFlowQCFinalPtDifHist[0][2][1][0]");

Int_t Nbinsx = h0->GetNbinsX();
TH1F *result = (TH1F*)h0->Clone("result");
TH1F *co;

Float_t ErrCorr;
Double_t arr1[Nbinsx];
Double_t arr2[Nbinsx];
Double_t Err1[Nbinsx];
Double_t Err2[Nbinsx];
Double_t val[Nbinsx];

for(Int_t j=1; j<9;j++){
if(j==1){co = (TH1F*)h1->Clone("co");}
if(j==2){co = (TH1F*)h2->Clone("co");}
if(j==3){co = (TH1F*)h3->Clone("co");}
if(j==4){co = (TH1F*)h4->Clone("co");}
if(j==5){co = (TH1F*)h5->Clone("co");}
if(j==6){co = (TH1F*)h6->Clone("co");}
if(j==7){co = (TH1F*)h7->Clone("co");}
if(j==8){co = (TH1F*)h8->Clone("co");}
if(j==9){co = (TH1F*)h9->Clone("co");}
	for(Int_t i=0; i<Nbinsx; i++){
	arr1[i] = result->GetBinContent(i+1);
	arr2[i] = co->GetBinContent(i+1);
	Err1[i] = result->GetBinError(i+1);
	Err2[i] = co->GetBinError(i+1);
	val[i] = (arr1[i]*(1/Err1[i])+arr2[i]*(1/Err2[i]))/2;
//cout<<"Subsample: "<< j<<" Err1 "<< (Err1[i] + Err2[i])<<endl;
	//ErrCorr = (1/2)*sqrt( Err1[i]*Err1[i] + Err2[i]*Err2[i]);
        ErrCorr = (Err1[i]+Err2[i]);//otherwise too small = 0
	//cout <<"BinNumber: "<<i<<" arr1: "<< arr1[i]<< " arr2: "<<arr2[i]<< " average: "<<val[i]<<" Err1: "<<Err1[i]<< " Err2: "<<Err2[i]<<" Error average : "<< ErrCorr<<endl;
	result->SetBinContent(i+1,val[i]);
	result->SetBinError(i+1,ErrCorr);
	/*
	Double_t f = hfull->GetBinContent(i+1);
	Double_t fErr = hfull->GetBinContent(i+1);
	hfull->SetBinContent(i+1,f/fErr);
	hfull->SetBinError(i+1,fErr);	
	*/
	}
}
//comparison from full sample vs results from subsampling method

TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetGridx();         // Vertical grid
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();               // pad1 becomes the current pad//h1->SetStats(0);
result->GetXaxis()->SetRangeUser(0.1,5.);
//hfull->SetStats(0);
result->Draw();               // Draw h1
//->SetLineWidth(4);
hfull->Draw("same");         // Draw h2 on top of h1

auto L_Ratio = new TLegend();
L_Ratio->SetHeader("Xe-Xe: 5.44TeV, 20-30%","C");
L_Ratio->AddEntry(hfull,"Full sample, pT: 0.2-5 GeV, #eta: 0.8", "l");
L_Ratio->AddEntry(result,"Subsampling, pT: 0.2-5 GeV, #eta: 0.8", "l");
L_Ratio->Draw();

hfull->GetYaxis()->SetTitle("differential flow v_{2}");
TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
axis->SetLabelSize(15);
axis->Draw();
c->cd();        
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.2);
pad2->SetGridx(); // vertical grid
pad2->Draw();
pad2->cd();     
TH1F *h = (TH1F*)hfull->Clone("h");
h->SetLineColor(kBlack);
h->GetXaxis()->SetRangeUser(0.1,5.);
Double_t Err;
for(Int_t i = 0; i<=Nbinsx; i++){
        arr1[i] = hfull->GetBinContent(i+1);
        arr2[i] = result->GetBinContent(i+1);
        Err1[i] = hfull->GetBinError(i+1);
        Err2[i] = result->GetBinError(i+1);
        val[i] = (arr1[i]-arr2[i]);
	Err = (1/2)*sqrt( Err1[i]*Err1[i] + Err2[i]*Err2[i]);
        Err = (Err1[i]+Err2[i]);//otherwise too small
	h->SetBinContent(i+1,val[i]);
	h->SetBinError(i+1,Err);
}
h->SetStats(0);
h->SetMarkerStyle(21);
h->Draw("ep");      
hfull->SetLineColor(kBlue+1);
hfull->SetLineWidth(1);
hfull->GetYaxis()->SetTitleSize(20);
hfull->GetYaxis()->SetTitleFont(43);
hfull->GetYaxis()->SetTitleOffset(1.55);
result->SetLineColor(kRed);
result->GetXaxis()->SetRangeUser(0.1,5.);
result->SetLineWidth(2);
h->SetTitle("");  
h->GetYaxis()->SetTitle("Full - Subsample");
h->GetYaxis()->SetNdivisions(505);
h->GetYaxis()->SetTitleSize(20);
h->GetYaxis()->SetTitleFont(43);
h->GetYaxis()->SetTitleOffset(1.55);
h->GetYaxis()->SetLabelFont(43); 
h->GetYaxis()->SetLabelSize(15);
h->GetXaxis()->SetTitleSize(20);
h->GetXaxis()->SetTitleFont(43);
h->GetXaxis()->SetTitleOffset(4.);
h->GetXaxis()->SetLabelFont(43);h3->GetXaxis()->SetLabelSize(15);
h->GetXaxis()->SetRangeUser(0.1,5.);
c->SaveAs("Error_Subsample_full.pdf");

}
