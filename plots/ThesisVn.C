#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>
//======================================================================================
//
//Little helps
//
//======================================================================================
void findSpread(TString input, Int_t c, Int_t harm,TString spectrum){
  TFile *f[20];
  TList *list[20];
  TH1D *fHist[20];
  //Int_t objID = 0;
  Int_t Nbins = 0;
  if(spectrum == "Pt"){Nbins = 30;}
  else if(spectrum == "Eta"){Nbins = 50;}
  TH1D *fHistPos[20];
  TH1D *fHistNeg[20];
  TH1F *fHistSpreadPos[Nbins];
  TH1F *fHistSpreadNeg[Nbins];
  TH1F *fHistSpreadDiff[Nbins];
 for(Int_t i = 0; i< Nbins; i++){
	fHistSpreadPos[i] = new TH1F(Form("fHistSpreadPos"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
	fHistSpreadNeg[i] = new TH1F(Form("fHistSpreadNeg"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
 	fHistSpreadDiff[i] = new TH1F(Form("fHistSpreadDiff"+spectrum+"Bin%d",i+1),"",1000,-0.3,0.3);
  }
  for(Int_t iFile = 0; iFile < 20; iFile++) {    //f[iFile] = TFile::Open(Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    f[iFile] = TFile::Open(input+Form("/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    if((!f[iFile])||(!f[iFile]->IsOpen())) {
    cout<<"File "<<iFile<<" not found..."<<endl; //return;
    }

    list[iFile] = dynamic_cast<TList *>(f[iFile]->Get("FlowQCList"));
    if(!list[iFile]) {
      cout<<"Input list of file "<<iFile<<" not found..."<<endl;
    }
    fHistPos[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[0][%d][%d][0]",c,harm)));
    if(!fHistPos[iFile]) {
      cout<<"Histogram of positive particles from file "<<iFile<<" not found..."<<endl;
    }
    fHistNeg[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[1][%d][%d][0]",c,harm)));
    if(!fHistNeg[iFile]) {
      cout<<"Histogram of negative particles from file "<<iFile<<" not found..."<<endl;
    }

    for(Int_t iBin = 1; iBin <= fHistPos[iFile]->GetNbinsX(); iBin++){
      fHistSpreadPos[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin));}
    for(Int_t iBin = 1; iBin <= fHistNeg[iFile]->GetNbinsX(); iBin++){
      fHistSpreadNeg[iBin-1]->Fill(fHistNeg[iFile]->GetBinContent(iBin));
      if(fabs(fHistPos[iFile]->GetBinContent(iBin))>0 && fabs(fHistNeg[iFile]->GetBinContent(iBin))>0){
      fHistSpreadDiff[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin) - fHistNeg[iFile]->GetBinContent(iBin));}
    }//loop to fill SpreadHists

  }//loop over files
  TFile *fOutput = new TFile(Form("output_findSpread/spread_"+spectrum+"_c%d0_%d0.root",c,c+1),"recreate");
  for(Int_t iBin = 1; iBin <= fHistPos[0]->GetNbinsX(); iBin++){fHistSpreadPos[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadNeg[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadDiff[iBin-1]->Write();}
  fOutput->Close();
  
 for(Int_t i = 0; i< Nbins; i++){
        delete fHistSpreadPos[i];
        delete fHistSpreadNeg[i];
        delete fHistSpreadDiff[i];
 }
//cout<<"end of find spread"<<endl;
}//end find spread

TH1F *Subsampling(TString spectrum, TString charge, TH1F *result, Int_t c, Int_t harm){
	Int_t Nbins = result->GetNbinsX();
	TH1F *s;
	Double_t sigma = 0.0;
	Double_t mean = 0.0;
    
    	TFile * of = new TFile(Form("output_findSpread/spread_"+spectrum+"_c%d0_%d0.root",c,c+1),"READ");
	for(Int_t i=1; i<=Nbins; i++){
		s =dynamic_cast<TH1F *>(of->GetKey(Form("fHistSpread"+charge+spectrum+"Bin%d",i))->ReadObj());	
		//cout<<s<<endl;
		TH1F *o = (TH1F *)s->Clone();
		sigma = o->GetRMS();
		mean = o->GetMean();//take the file check the values fro 3*sigma
		//cout<<"mean "<<mean<<" sigma "<<sigma<<endl;
		result->SetBinContent(i, mean);
		result->SetBinError(i, sigma);
	}//end Nbins

return result;
}//end subsampling


//=================================================================================================
//
//The Main from this macro
//
//
//Trouble notes:
//______________
//
// 3) Maybe improve errors on eta differential results through smaller pT integration range -> verify with 3 GeV in Analysis once the plots are prepared
//=================================================================================================
void ThesisVn(){

	gROOT->SetBatch();//to avoid opening the plots ak bad wifi struggle
	TH1F *copy[2];//Multi centrality plots
	
	//Bins from AVFD analysis
	Double_t fPtDiffNBins = 29;
        Double_t fCRCPtBins[30];
        Double_t PtBins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5.};
	for(Int_t r=0; r<30; r++){fCRCPtBins[r] = PtBins[r];}

	Double_t fEtaDiffNBins = 50;
	Double_t fCRCEtaBins[51]={0};
	Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
	for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}
	
	//Subsampling final results
	for(Int_t harm=0; harm<=3;harm++){
		//TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
		TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.2/split");
		// !!! For harm = 0 we could also try to use the fFlowQCQv1[pos=0, neg=1][cen][pt = 0, eta = 1] !!!
		cout<<"Harmonic: "<<(harm+1)<<endl;
		//apply subsampling to all histograms from pos/neg particles=======================
		findSpread(input,1,harm,"Pt");
		findSpread(input,4,harm,"Pt");
		TH1F *PvpT1 = new TH1F("pvpT1", "pvpT1", fPtDiffNBins, fCRCPtBins);
		PvpT1 = Subsampling("Pt", "Pos", PvpT1, 1, harm);
		TH1F *NvpT1 = new TH1F("nvpT1","nvpT1",fPtDiffNBins,fCRCPtBins);
		NvpT1 = Subsampling("Pt","Neg",NvpT1, 1 ,harm);
		TH1F *DvpT1 = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
		DvpT1 = Subsampling("Pt","Diff", DvpT1,1,harm);

                TH1F *PvpT4 = new TH1F("pvpT4", "pvpT4", fPtDiffNBins, fCRCPtBins);
                PvpT4 = Subsampling("Pt", "Pos", PvpT4, 4, harm);
                TH1F *NvpT4 = new TH1F("nvpT4","nvpT4",fPtDiffNBins,fCRCPtBins);
                NvpT4 = Subsampling("Pt","Neg",NvpT4, 4 ,harm);
                TH1F *DvpT4 = new TH1F("dvnpT4","dvnpT4",fPtDiffNBins,fCRCPtBins);
                DvpT4 = Subsampling("Pt","Diff", DvpT4,4,harm);

		//Multi centrality plot
		copy[0] = (TH1F*) DvpT1->Clone("copy dvpt1");
		copy[1] = (TH1F*) DvpT4->Clone("copy dvpt4");
		//The final pT differential plots================================================================================

        	short Col1,Col2;
        	double ymin, ymax,dmin,dmax;
        	double x1,x11,x2,x22,y1,y11,y2,y22;
        	TString flow;
        	if(harm == 0){
        	Col1 = kOrange;
        	Col2 = kRed;
        	ymin = -0.001;
        	ymax = 0.024;
        	dmin = -0.011;
        	dmax = 0.013;
        	x1 = 0.22;
        	x2 = 0.68;
        	y1 = 0.59;
        	y2 = 0.83;
		x11 = 0.53;
		x22 = 0.99;
		y11 = 0.59;
		y22 = 0.83;
        	flow = "directed";}
        	if(harm == 1){
        	Col1 = kMagenta;
        	Col2 = kViolet;
        	ymin = -0.01;
        	ymax = 0.31;
        	dmin = -0.031;
        	dmax = 0.023;
        	x1 = 0.22;
        	x2 = 0.58;
        	y1 = 0.59;
        	y2 = 0.84;
                x11 = 0.63;
                x22 = 0.99;
                y11 = 0.01;
                y22 = 0.26;
        	flow = "elliptic";}
        	if(harm == 2){
        	Col1 = kGreen;
        	Col2 = kBlue;
        	ymin = -0.01;
        	ymax = 0.18;
        	dmin = -0.044;
        	dmax = 0.034;
                x1 = 0.22;
                x2 = 0.58;
                y1 = 0.59;
                y2 = 0.84;
                x11 = 0.63;
                x22 = 0.99;
                y11 = 0.01;
                y22 = 0.26;
        	flow = "triangular";}

       		TCanvas *c1= new TCanvas("c1", "canvas", 800, 800);//upper plot in pad1
        	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1., 1.0);
        	pad1->SetBottomMargin(0); // Upper and lower plot are joined
        	pad1->SetTopMargin(0.16);
		pad1->SetLeftMargin(0.2);
		pad1->SetRightMargin(0.05);
        	pad1->Draw();             // Draw the upper pad: pad1
        	pad1->cd();               // pad1 becomes the current pad//and filled with pos
        	double f1 = (pad1->GetWNDC())*(pad1->GetHNDC());
		PvpT1->GetYaxis()->SetRangeUser(ymin,ymax);
		PvpT1->GetXaxis()->SetRangeUser(0.1,3);
        	PvpT1->SetLineColor(Col1+1);
        	PvpT1->SetLineWidth(1);
        	PvpT1->GetYaxis()->SetTitle("differential "+flow+Form(" flow #nu_{%d}",harm+1));
		PvpT1->GetYaxis()->CenterTitle();
       		PvpT1->GetYaxis()->SetTitleSize(0.18*(1-f1));
        	PvpT1->GetYaxis()->SetLabelSize(0.18*(1-f1));
        	PvpT1->GetYaxis()->SetTitleOffset((1.75-f1));
		PvpT1->GetYaxis()->SetNdivisions(6);
		PvpT1->SetStats(0);
        	PvpT1->SetMarkerStyle(22);
        	PvpT1->SetMarkerColor(Col1+1);
        	PvpT1->SetTitle(" ");
        	PvpT1->DrawCopy();
        	//NvpT1->GetXaxis()->SetRangeUser(0.,3);//range in bins
        	NvpT1->SetLineColor(Col1+3);
        	NvpT1->SetLineWidth(1);
        	NvpT1->SetMarkerStyle(23);
        	NvpT1->SetMarkerColor(Col1+3);
        	NvpT1->SetTitle(" ");
        	NvpT1->DrawCopy("same");

        	PvpT4->SetLineColor(Col2+1);
        	PvpT4->SetLineWidth(1);
        	PvpT4->SetMarkerStyle(26);
        	PvpT4->SetMarkerColor(Col2+1);
        	PvpT4->SetTitle(" ");
        	PvpT4->DrawCopy("same");
        	//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins
        	NvpT4->SetLineColor(Col2+3);
        	NvpT4->SetLineWidth(1);
        	NvpT4->SetMarkerStyle(32);
        	NvpT4->SetMarkerColor(Col2+3);
        	NvpT4->SetTitle(" ");
        	NvpT4->DrawCopy("same");

		auto L0 = new TLegend(0.05,0.9,0.95,1);
		L0->SetTextSize(0.035);
		L0->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, #tau_{0} = 0.4 fm/c, #tau_{B} = 0.2 fm/c, |#eta| #leq 0.8 ","C");
		L0->SetFillStyle(0);
		L0->SetBorderSize(0);

		auto L01 = new TLegend(x1+0.006,0.3,x2+0.006,0.42);
		L01->SetNColumns(2);
		L01->SetTextSize(0.033);//L01->AddEntry(DvpT1,Form("#Delta q, Cent. 10%%-20%%"), "lep");
                L01->AddEntry(DvpT4,Form("#bf{#Delta q, 40%%-50%%}"), "lep");
                L01->SetBorderSize(0);
                L01->SetFillStyle(0);

                auto L011 = new TLegend(x1+0.006,0.06,x2+0.006,0.18);
                L011->SetNColumns(2);
                L011->SetTextSize(0.033);
                L011->AddEntry(DvpT1,Form("#bf{#Delta q, 10%%-20%%}"), "lep");
                //L011->AddEntry(DvpT4,Form("#Delta q, Cent. 40%%-50%%"), "lep");
                L011->SetBorderSize(0);
                L011->SetFillStyle(0);

        	auto L1 = new TLegend(x1,y1,x2,y2);//L1->SetTextSize(0.04);//L1->SetHeader("#bf{AVFD Simulation Pb-Pb} @ #sqrt{s} = 5.02TeV","C");
        	L1->SetTextSize(0.055);//L2->SetTextAlign(11);//or13 //L1->SetEntrySeparation(0.06);//L1->SetNColumns(2); 
	        //L1->AddEntry(PvpT1,Form("+h, Cent. 10%%-20%%"), "lep");
        	L1->AddEntry(PvpT4,Form("#bf{+h, 40%%-50%%}"), "lep");
		//L1->AddEntry(NvpT1,Form("-h, Cent. 10%%-20%%"), "lep");
        	L1->AddEntry(NvpT4,Form("#bf{-h, 40%%-50%%}"), "lep");//L1->AddEntry(DvpT1,Form("#Delta q, Centrality 10%%-20%%"), "lep");//L1->AddEntry(DvpT4,Form("#Delta q, Centrality 40%%-50%%"), "lep");
		L1->SetBorderSize(0);
	        L1->SetFillStyle(0);
		L1->Draw();

                auto L11 = new TLegend(x11,y11,x22,y22);//L1->SetTextSize(0.04);//L1->SetHeader("#bf{AVFD Simulation Pb-Pb} @ #sqrt{s} = 5.02TeV","C");
                L11->SetTextSize(0.055);//L2->SetTextAlign(11);//or13 //L1->SetEntrySeparation(0.06);//L1->SetNColumns(2);
                L11->AddEntry(PvpT1,Form("#bf{+h, 10%%-20%%}"), "lep");
                //L11->AddEntry(PvpT4,Form("+h, Cent. 40%%-50%%"), "lep");
                L11->AddEntry(NvpT1,Form("#bf{-h, 10%%-20%%}"), "lep");
                //L11->AddEntry(NvpT4,Form("-h, Cent. 40%%-50%%"), "lep");//L1->AddEntry(DvpT1,Form("#Delta q, Centrality 10%%-20%%"), "lep");//L1->AddEntry(DvpT4,Form("#Delta q, Centrality 40%%-50%%"), "lep");
                L11->SetBorderSize(0);
                L11->SetFillStyle(0);
                L11->Draw();

        	c1->cd();          // Go back to the main canvas before defining pad2
        	L0->Draw();
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.4);
        	pad2->SetTopMargin(0);
        	pad2->SetLeftMargin(0.2);      //pad2->SetRightMargin(4.);
        	pad2->SetBottomMargin(0.23);
		pad2->SetRightMargin(0.05);
        	//pad2->SetGridx(); // vertical grid
        	pad2->Draw();
        	pad2->cd();
		double f2 = (pad2->GetWNDC())*(pad2->GetHNDC());
        	DvpT1->SetTitle(" ");
        	DvpT1->SetStats(0);
		DvpT1->GetYaxis()->SetRangeUser(dmin,dmax);
        	DvpT1->SetMarkerStyle(20);
        	DvpT1->SetMarkerColor(Col1+2);
        	DvpT1->GetXaxis()->SetRangeUser(0.1,3);
        	DvpT1->SetLineColor(Col1+2);
        	DvpT1->SetLineWidth(1);
        	DvpT1->GetXaxis()->SetTitle("p_{T} [GeV]");
        	DvpT1->GetYaxis()->SetTitle(Form("#nu_{%d}(+h) - #nu_{%d}(-h)",harm+1,harm+1));
        	DvpT1->GetYaxis()->CenterTitle();
		DvpT1->GetYaxis()->SetTitleSize(0.17*(1-f2));
        	DvpT1->GetYaxis()->SetLabelSize(0.18*(1-f2));
	        DvpT1->GetYaxis()->SetNdivisions(4);
        	DvpT1->GetXaxis()->SetTitleSize(0.16*(1-f2));
	        DvpT1->GetXaxis()->SetLabelSize(0.18*(1-f2));
	        DvpT1->GetYaxis()->SetTitleOffset((1.17-f2));
		DvpT1->DrawCopy();
        
		DvpT4->SetMarkerStyle(24);
        	DvpT4->SetMarkerColor(Col2+2);
        	DvpT4->GetXaxis()->SetRangeUser(0.1,3);
        	DvpT4->SetLineColor(Col2+2);
        	DvpT4->SetLineWidth(1);
        	DvpT4->DrawCopy("same");
		c1->cd();
		L01->Draw();
		L011->Draw();
        	c1->SaveAs(Form("vn/v%d_pT.pdf",harm+1));
       		delete c1;

                PvpT1->Delete();
                NvpT1->Delete();
                DvpT1->Delete();
        	PvpT4->Delete();
        	NvpT4->Delete();
        	DvpT4->Delete();
	
		//Now for eta differential flow====================================================================================
                findSpread(input,1,harm,"Eta");
                findSpread(input,4,harm,"Eta");
                
		TH1F *PvEta1 = new TH1F("pvEta", "pvEta", fEtaDiffNBins, fCRCEtaBins);
                PvEta1 = Subsampling("Eta", "Pos", PvEta1, 1, harm);
                TH1F *NvEta1 = new TH1F("nvEta1", "nvEta1", fEtaDiffNBins, fCRCEtaBins);
                NvEta1 = Subsampling("Eta", "Neg", NvEta1, 1, harm);
                TH1F *DvEta1 = new TH1F("deltavneta1","deltavneta1",fEtaDiffNBins,fCRCEtaBins);
                DvEta1 = Subsampling("Eta", "Diff", DvEta1, 1, harm);
	
		TH1F *PvEta4 = new TH1F("pvEta4", "pvEta4", fEtaDiffNBins, fCRCEtaBins);
		PvEta4 = Subsampling("Eta", "Pos", PvEta4, 4, harm);		
		TH1F *NvEta4 = new TH1F("nvEta4", "nvEta4", fEtaDiffNBins, fCRCEtaBins);
		NvEta4 = Subsampling("Eta", "Neg", NvEta4, 4, harm);
		TH1F *DvEta4 = new TH1F("deltavneta4","deltavneta4",fEtaDiffNBins,fCRCEtaBins);
		DvEta4 = Subsampling("Eta", "Diff", DvEta4, 4, harm);
		//Multi centrality plot
		copy[0] = (TH1F*) DvEta1->Clone("copy dvpt1");
                copy[1] = (TH1F*) DvEta4->Clone("copy dvpt4");
		//Plot the pos/neg with difference=================================================
		//PlotEta(PvEta1, NvEta1, DvEta1, 1, harm);
                //PlotEta(PvEta4, NvEta4, DvEta4, 4, harm);

		//The final Eta differential plots================================================================================== 	
		if(harm == 0){
        	Col1 = kOrange;
        	Col2 = kRed;
        	ymin = -0.0005;
        	ymax = 0.0055;
        	dmin = -0.0081;
	        dmax = 0.0081;
        	x1 = 0.22;
	        x2 = 0.72;
        	y1 = 0.59;
	        y2 = 0.83;
		flow = "directed";}
	        if(harm == 1){
        	Col1 = kMagenta;
        	Col2 = kViolet;
        	ymin = 0.042;
        	ymax = 0.092;
        	dmin = -0.011;
        	dmax = 0.011;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.26;
        	y2 = 0.50;
		flow = "elliptic";}
        	if(harm == 2){
        	Col1 = kGreen;
        	Col2 = kBlue;
        	ymin = 0.015;
        	ymax = 0.055;
        	dmin = -0.022;
        	dmax = 0.023;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.60;
        	y2 = 0.84;
		flow = "triangular";}

        	TCanvas *cEta = new TCanvas("ceta", "etacanvas", 800, 800);
        	TPad *padeta1 = new TPad("padeta1", "padeta1", 0, 0.4, 1., 1.);
        	padeta1->SetBottomMargin(0);
		padeta1->SetTopMargin(0.16);
        	padeta1->SetLeftMargin(0.2);
		padeta1->SetRightMargin(0.05);
        	padeta1->Draw();
        	padeta1->cd();
		double f3 = (padeta1->GetWNDC())*(padeta1->GetHNDC());
		//cout<<f3<<endl;
		PvEta1->GetYaxis()->SetRangeUser(ymin,ymax);
        	PvEta1->GetXaxis()->SetRangeUser(-0.8,0.8);//PvEta->GetYaxis()->SetNdivisions(4) //PvEta->GetYaxis()->SetRangeUser(y-y*0.35,y+y*0.35);
        	PvEta1->SetLineColor(Col1+1);
        	PvEta1->SetLineWidth(1);
        	PvEta1->GetYaxis()->SetTitle("differential "+flow+" flow "+Form("#nu_{%d}",harm+1));
		PvEta1->GetYaxis()->CenterTitle();
        	PvEta1->GetYaxis()->SetTitleSize(0.18*(1-f3));
		PvEta1->GetYaxis()->SetLabelSize(0.18*(1-f3));
		PvEta1->GetYaxis()->SetTitleOffset((1.85-f3));
		PvEta1->GetYaxis()->SetNdivisions(6);
		//PvEta1->GetYaxis()->CenterTitle(true);
		PvEta1->SetStats(0);
		//PvEta1->GetYaxis()->SetNdivisions(6);
        	PvEta1->SetMarkerStyle(22);
        	PvEta1->SetMarkerColor(Col1+1);
        	PvEta1->SetTitle(" ");
        	PvEta1->DrawCopy();      //NvEta->GetYaxis()->SetNdivisions(4);
        	NvEta1->SetLineColor(Col1+3);
        	NvEta1->SetLineWidth(1);
        	NvEta1->SetTitle(" ");
        	NvEta1->SetMarkerStyle(23);
        	NvEta1->SetMarkerColor(Col1+3);
        	NvEta1->DrawCopy("same");

        	PvEta4->SetLineColor(Col2+1);
        	PvEta4->SetLineWidth(1);
        	PvEta4->SetMarkerStyle(26);
        	PvEta4->SetMarkerColor(Col2+1);
        	PvEta4->SetTitle(" ");
        	PvEta4->DrawCopy("same");      //NvEta->GetYaxis()->SetNdivisions(4);
        	NvEta4->SetLineColor(Col2+3);
        	NvEta4->SetLineWidth(1);
        	NvEta4->SetTitle(" ");
        	NvEta4->SetMarkerStyle(32);
        	NvEta4->SetMarkerColor(Col2+3);
        	NvEta4->DrawCopy("same");

                auto L0E = new TLegend(0.05,0.9,0.95,1);
                L0E->SetTextSize(0.032);
                L0E->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, #tau_{0} = 0.4 fm/c, #tau_{B} = 0.2 fm/c, p_{T} #in [0.2,5]GeV ","C");
                L0E->SetFillStyle(0);
                L0E->SetBorderSize(0);

                auto L02 = new TLegend(x1,0.3,x2,0.4);
                L02->SetNColumns(2);
		L02->SetColumnSeparation(0.75);
		L02->SetTextSize(0.035);
                L02->AddEntry(DvEta1,Form("#bf{#Delta q, 10%%-20%%}"), "lep");
                L02->AddEntry(DvEta4,Form("#bf{#Delta q, 40%%-50%%}"), "lep");
                L02->SetBorderSize(0);
                L02->SetFillStyle(0);

        	auto L2 = new TLegend(x1,y1,x2,y2);
		//L2->SetTextSize(0.04);
        	//L2->SetHeader("#bf{AVFD Simulation Pb-Pb} @ #sqrt{s} = 5.02TeV","C");
		L2->SetTextSize(0.056);//L2->SetTextAlign(11);//or13
		//L2->AddEntry((TObject*)0, "p_{T} #in {0.2,5} GeV, |#eta| #leq 0.8"," ");
        	//L2->AddEntry((TObject*)0, "#tau_{0} = 0.4 fm/c, #tau_{B} = 0.2 fm/c"," ");
		//L2->SetColumnSeparation(0.009);
		L2->SetNColumns(2);
                L2->SetColumnSeparation(0.75);
		if(harm == 0){L2->SetColumnSeparation(0.1);}
		L2->AddEntry(PvEta1,Form("#bf{+h, 10%%-20%%}"), "lep");
        	L2->AddEntry(PvEta4,Form("#bf{+h, 40%%-50%%}"), "lep");
        	L2->AddEntry(NvEta1,Form("#bf{-h, 10%%-20%%}"), "lep");
		L2->AddEntry(NvEta4,Form("#bf{-h, 40%%-50%%}"), "lep");
		//L2->AddEntry(DvEta1,Form("#Delta q, Centrality 10%%-20%%"),"lep");
		//L2->AddEntry(DvEta4,Form("#Delta q, Centrality 40%%-50%%"), "lep");
		L2->SetBorderSize(0);
		L2->SetFillStyle(0);
		L2->Draw();
        	cEta->cd();
		L0E->Draw();

        	TPad *padeta2 = new TPad("padeta2", "padeta2", 0, 0., 1, 0.4);
        	padeta2->SetTopMargin(0);
        	padeta2->SetLeftMargin(0.2);       //padeta2->SetRightMargin(4.);
        	padeta2->SetBottomMargin(0.23);
        	padeta2->SetRightMargin(0.05);
        	padeta2->Draw();
        	padeta2->cd();
		double f4 = (padeta2->GetWNDC())*(padeta2->GetHNDC());
        	//cout<<f4<<endl;;
		DvEta1->SetStats(0);
		DvEta1->GetYaxis()->SetRangeUser(dmin,dmax);
        	DvEta1->GetXaxis()->SetRangeUser(-0.8,0.8);
        	DvEta1->GetXaxis()->SetTitle("#eta");
        	DvEta1->SetTitle(" ");
        	DvEta1->SetMarkerStyle(20);
        	DvEta1->SetMarkerColor(Col1+2);
        	DvEta1->SetLineColor(Col1+2);
        	DvEta1->SetLineWidth(1);
        	DvEta1->GetYaxis()->SetTitle(Form("#nu_{%d}(+h) - #nu_{%d}(-h)",harm+1,harm+1));
		DvEta1->GetYaxis()->SetTitleSize(0.17*(1-f4));
        	DvEta1->GetYaxis()->CenterTitle();
		DvEta1->GetYaxis()->SetLabelSize(0.18*(1-f4));
		DvEta1->GetYaxis()->SetNdivisions(4);
		DvEta1->GetXaxis()->SetTitleSize(0.16*(1-f4));
		DvEta1->GetXaxis()->SetLabelSize(0.17*(1-f4));
		DvEta1->GetYaxis()->SetTitleOffset((1.2-f4));
		DvEta1->Draw();
		DvEta4->SetMarkerStyle(24);
		DvEta4->SetMarkerColor(Col2+2);
		DvEta4->SetLineColor(Col2+2);
		DvEta4->SetLineWidth(1);
		DvEta4->DrawCopy("same");
        	cEta->cd();
		L02->Draw();
                cEta->SaveAs(Form("vn/v%d_eta.pdf",harm+1));
		delete cEta;

	
        	PvEta1->Delete();
	       	NvEta1->Delete();
	        DvEta1->Delete();
                PvEta4->Delete();
                NvEta4->Delete();
                DvEta4->Delete();

	}//harmonic loop
/*
	//Multi Centrality plot=========================================================
	TCanvas *Cen = new TCanvas("Cen","Cen",400,400);
	Cen->SetLeftMargin(0.2);
	auto C = new TLegend();
	C->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
	copy[cent]->GetYaxis()->SetTitleOffset(2.0);
	copy[cent]->GetYaxis()->SetTitle(Form("differential flow #Delta v_{%d}",harm+1));
	if(pT==true){
		copy[cent]->GetXaxis()->SetTitle("p_{T} [GeV]");
		copy[cent]->GetXaxis()->SetRangeUser(0,5);
	}else if(pT==false){
        	copy[cent]->GetXaxis()->SetRangeUser(-0.9,0.9);
        	copy[cent]->GetXaxis()->SetTitle("#eta");
	}
	copy[cent]->SetStats(0);
	copy[cent]->SetTitle(" ");
	//copy[cent]->SetMarkerStyle(20);
	copy[cent]->SetLineColor(cent);
	C->AddEntry(copy[cent],Form("AVFD, Centrality %d0-%d0",cent,cent+1), "l");
	copy[cent]->DrawCopy();
	for(Int_t c = cent+1; c<=cmax; c++){
		copy[c]->SetLineColor(c);
                C->AddEntry(copy[c],Form("AVFD, Centrality %d0-%d0",c,c+1), "l");
		copy[c]->DrawCopy("same");
	}
	C->Draw();
	if(pT==true){Cen->SaveAs(Form("vn/tau_init0.4/BField0.4/v%d_Pt_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
	else if(pT == false){Cen->SaveAs(Form("vn/tau_init0.4/BField0.4/v%d_Eta_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
	delete Cen;

	//}//end of hr<fkFlowNHarm
*/
}//end void


