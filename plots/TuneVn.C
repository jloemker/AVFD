#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>
//======================================================================================
//
//Uncertainty handling
//
//======================================================================================
void findSpread(TString input, Int_t c, Int_t harm,TString spectrum, Int_t tau, Int_t tauB){
  TFile *f[20];
  TList *list[20];
  TH1D *fHist[20];
  Int_t Nbins = 0;
  if(spectrum == "Pt"){Nbins = 30;}
  else if(spectrum == "Eta"){Nbins = 50;}
  TH1D *fHistPos[20];
  TH1D *fHistNeg[20];
  TH1F *fHistSpreadPos[Nbins];
  TH1F *fHistSpreadNeg[Nbins];
  TH1F *fHistSpreadDiff[Nbins];
  TH1F *fHistSpreadWavg[Nbins];
 for(Int_t i = 0; i< Nbins; i++){
	fHistSpreadPos[i] = new TH1F(Form("fHistSpreadPos"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
	fHistSpreadNeg[i] = new TH1F(Form("fHistSpreadNeg"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
 	fHistSpreadDiff[i] = new TH1F(Form("fHistSpreadDiff"+spectrum+"Bin%d",i+1),"",1000,-0.3,0.3);
        fHistSpreadWavg[i] = new TH1F(Form("fHistSpreadWavg"+spectrum+"Bin%d",i+1),"",1000,-0.3,0.3);
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
    double p, pErr, wp, n, nErr, wn; 
    for(Int_t iBin = 1; iBin <= fHistPos[iFile]->GetNbinsX(); iBin++){
      fHistSpreadPos[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin));
	/*Weighted average of Pos and neg
      p = fHistPos[iFile-1]->GetBinContent(iBin);
      pErr = fHistPos[iBin-1]->GetBinError(iBin);
      wp = 1/(pErr*pErr);
      n = fHistNeg[iFile]->GetBinContent(iBin);
      nErr = fHistNeg[iFile]->GetBinError(iBin);
      wn = 1/(nErr*nErr);          
      fHistSpreadWavg[iBin-1]->Fill((p*wp+ n*wn)/(wp+wn));*/
    }
    for(Int_t iBin = 1; iBin <= fHistNeg[iFile]->GetNbinsX(); iBin++){
      fHistSpreadNeg[iBin-1]->Fill(fHistNeg[iFile]->GetBinContent(iBin));
      if(fabs(fHistPos[iFile]->GetBinContent(iBin))>0 && fabs(fHistNeg[iFile]->GetBinContent(iBin))>0){
      fHistSpreadDiff[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin) - fHistNeg[iFile]->GetBinContent(iBin));}
    }//loop to fill SpreadHists

  }//loop over files
  TFile *fOutput = new TFile(Form("output_findSpread/spread_"+spectrum+"_0.%d_B0.%d.root",tau,tauB),"recreate");
  for(Int_t iBin = 1; iBin <= fHistPos[0]->GetNbinsX(); iBin++){fHistSpreadPos[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadNeg[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadDiff[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistPos[0]->GetNbinsX(); iBin++){fHistSpreadWavg[iBin-1]->Write();}
  fOutput->Close();
  
 for(Int_t i = 0; i< Nbins; i++){
        delete fHistSpreadPos[i];
        delete fHistSpreadNeg[i];
        delete fHistSpreadDiff[i];
 	delete fHistSpreadWavg[i];
 }
cout<<"end of find spread"<<endl;
}//end find spread

TH1F *Subsampling(TString spectrum, TString charge, TH1F *result, Int_t c, Int_t harm, Int_t tau, Int_t tauB){
	Int_t Nbins = result->GetNbinsX();
	TH1F *s;
	Double_t sigma = 0.0;
	Double_t mean = 0.0;
    
    	TFile * of = new TFile(Form("output_findSpread/spread_"+spectrum+"_0.%d_B0.%d.root",tau,tauB),"READ");
	for(Int_t i=1; i<=Nbins; i++){
		s =dynamic_cast<TH1F *>(of->GetKey(Form("fHistSpread"+charge+spectrum+"Bin%d",i))->ReadObj());
		TH1F *o = (TH1F *)s->Clone();
		sigma = o->GetRMS();
		mean = o->GetMean();
		result->SetBinContent(i, mean);
		result->SetBinError(i, sigma);
	}//end Nbins
	of->Close();
return result;
}//end subsampling

TH1F *Result(TString spectrum, TH1F *res, Int_t c, Int_t harm, Int_t tau, Int_t tauB){
        Int_t Nbins = res->GetNbinsX();
        TH1F *s1, *s2;
        Double_t sigma = 0.0;
        Double_t mean = 0.0;
        double p, pErr, wp, n, nErr, wn;
	TH1F *pos = (TH1F *)res->Clone();
        TH1F *neg = (TH1F *)res->Clone();
        TFile * of = new TFile(Form("output_findSpread/spread_"+spectrum+"_0.%d_B0.%d.root",tau,tauB),"READ");
        for(Int_t i=1; i<=Nbins; i++){
                s1 =dynamic_cast<TH1F *>(of->Get(Form("fHistSpreadPos"+spectrum+"Bin%d",i))); 
		TH1F *o1 = (TH1F *)s1->Clone();
                p = o1->GetMean();
                pErr = o1->GetRMS();
		pos->SetBinContent(i, p);
		pos->SetBinError(i,pErr);
                s2 =dynamic_cast<TH1F *>(of->Get(Form("fHistSpreadNeg"+spectrum+"Bin%d",i))); 
                TH1F *o2 = (TH1F *)s2->Clone();
                n = o2->GetMean();
                nErr = o2->GetRMS();
                neg->SetBinContent(i,n);
		neg->SetBinError(i,nErr);
		}		
	for(Int_t k=1; k<=Nbins; k++){
		p = pos->GetBinContent(k);
		pErr = pos->GetBinError(k);
		if(pErr != 0){wp = 1/(pErr*pErr);
		}else{wp = 1;} 
		n = neg->GetBinContent(k);
		nErr = neg->GetBinError(k);
		if(nErr!= 0){wn = 1/(nErr*nErr);
		}else{wn = 1;}
                res->SetBinContent(k, (p*wp+ n*wn)/(wp+wn));//
		res->SetBinError(k, sqrt(1/(wp+wn)));
        cout<<"weighted avg result"<<(p*wp+ n*wn)/(wp+wn)<<" with error: 1/(wp+wn) = "<<1/(wp+wn)<<" vs. the sqrt(error): "<<sqrt(1/(wp+wn))<<endl;
	}//end Nbins
	of->Close();
	return res;
}

void FitResult(TH1F *result, TString save){
        TH1F *v2 = (TH1F*)result->Clone("v2");
        TH1F *h3 = (TH1F*)v2->Clone("h3");
        Int_t NBins = v2->GetNbinsX();
        Double_t arr1[NBins];
        Double_t arr2[NBins];
        Double_t ErrCorr[NBins];
        Double_t Err1[NBins];
        Double_t Err2[NBins];
	//The binnig is perfect:version2 tab 8 from https://www.hepdata.net/record/ins1419244
	//Double_t fPtBins[14] = {0.3,0.5,0.7,0.9,1.25,1.375,1.625,1.875,2.25,2.75,3.25,3.75,4.5,5.};//original
	//edges
	Double_t fPtBins[13] = {0.2,0.4,0.6,0.8,1.,1.25,1.5,1.75,2,2.5,3.,3.5,4};
	Double_t h2_x[13];
	for(int r = 0; r<13; r++){h2_x[r] = fPtBins[r];}
	Double_t h2_y[12] = {0.0477,0.0771,0.1024,0.1249,0.1695,0.192,0.2068,0.2225,0.2395,0.2494,0.2448,0.198};
	Double_t h2_err[12] = {0.0021,0.003,0.0038,0.0053,0.0032,0.0039,0.0058,0.0047,0.0049,0.007,0.0057,0.0126};//,0.016};
	TH1F *h2 = new TH1F("ALICE Fit: Pol 6","ALICE vs AVFD",12,h2_x);//0.3,4.5);
	TH1F *hErr = new TH1F("ErrFit","ErrFit",12, h2_x);

	for(int i=1; i <= 12; ++i) {
   	h2->SetBinContent(i,h2_y[i-1]);
   	h2->SetBinError(i,h2_err[i-1]);
   	hErr->SetBinContent(i,(h2_err[i-1]));//+h2_y[i-1]));//+the central value
	cout<<"Set Bin :"<<i<<" at "<<h2->GetBinCenter(i)<<" with "<<h2->GetBinContent(i)<<endl;
	}

	//error function for v2 fit
	TH1F *v2Err = (TH1F*)v2->Clone("v2Err");
	for(int t = 1; t<=v2->GetNbinsX();t++){
	v2Err->SetBinContent(t,v2->GetBinError(t));
	}

	TF1 *func = new TF1("func","pol 2",0,3.8);
	
	h2->Fit("func","R");//Fit ALICE
	TF1 *func_res = h2->GetFunction("func");
	Double_t chi2_h2 = func_res->GetChisquare();
	Double_t NDF = func_res ->GetNDF();
	cout<<"h2 ALICE: chi2/ndf "<<chi2_h2/NDF<<endl;

        hErr->Fit("func","R");//Fitting errors from ALICE
        TF1 *func_resErr = hErr->GetFunction("func");
        Double_t chi2_hErr = func_resErr->GetChisquare();
        Double_t NDFErr = func_resErr ->GetNDF();
        cout<<"hErr ALICE: chi2/ndf "<<chi2_hErr/NDFErr<<endl;
	double p00 = func_res->GetParameter(0);
	double p11 = func_res->GetParameter(1);
	double p22 = func_res->GetParameter(2);
        double pe00 = func_res->GetParError(0);
        double pe11 = func_res->GetParError(1);
        double pe22 = func_res->GetParError(2);

	//double p33 = func_res->GetParameter(3);

	cout<<"Alice: p0 "<<p00<<"+-"<<pe00<<" p1 "<<p11<<"+-"<<pe11<<" p2 "<<p22<<"+-"<<pe22<<endl;

        v2->Fit("func","R");//FIt AVFD
        TF1 *func_v2 = v2->GetFunction("func");
        Double_t chi2_v2 = func_v2->GetChisquare();
        Double_t NDF_v2 = func_v2->GetNDF();
        cout<<"v2 AVFD: ch2/ndf "<<chi2_v2/NDF_v2<<endl;

        v2Err->Fit("func","R");//Fitting errors from AVFD
        TF1 *func_v2Err = v2Err->GetFunction("func");
        Double_t chi2_v2Err = func_v2Err->GetChisquare();
        Double_t NDF_v2Err = func_v2Err->GetNDF();
        cout<<"v2 AVFD: ch2/ndf "<<chi2_v2Err/NDF_v2Err<<endl;
        double p0 = func_v2->GetParameter(0);
        double p1 = func_v2->GetParameter(1);
        double p2 = func_v2->GetParameter(2);
        double pe0 = func_v2->GetParError(0);
        double pe1 = func_v2->GetParError(1);
        double pe2 = func_v2->GetParError(2);

        //double p3 = func_v2->GetParameter(3);
	 
	cout<<"AVFD: p0 "<<p0<<"+-"<<pe0<<" p1 "<<p1<<"+-"<<pe1<<" p2 "<<p2<<"+-"<<pe2<<" AL/AV "<<(p00/p0) + (p11/p1) + (p22/p2) <<endl;
	for(Int_t i = 0; i<NBins; i++){
        	arr1[i] = v2->GetBinContent(i+1);
        	Err1[i] = v2->GetBinError(i+1);
        	arr2[i] = func_res->Eval(v2->GetBinCenter(i+1));
        	Err2[i] = func_resErr->Eval(v2->GetBinCenter(i+1));
        	ErrCorr[i] = sqrt(pow(Err1[i],2)+pow(Err2[i],2));
           	if(abs(arr1[i])>0 && abs(arr2[i])>0){
           		h3->SetBinContent(i+1, arr1[i]-arr2[i]);
           		h3->SetBinError(i+1,ErrCorr[i]);
           	}
	}

	TCanvas *c = new TCanvas("c", "canvas", 800, 800);//upper plot in pad1
	TPad *pad1 = new TPad("padtest", "padtest", 0, 0.3, 1,0.9);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetLeftMargin(0.2);
	pad1->SetRightMargin(0.1);
	//pad1->SetGridx();         // Vertical grid
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad//h1->SetStats(0);
	v2->GetXaxis()->SetRangeUser(0.2,3.5);
	v2->GetYaxis()->SetRangeUser(0.001,0.26);
	v2->SetFillStyle(2);
	v2->SetFillColor(kMagenta+2);
	v2->SetMarkerStyle(8);
	v2->SetMarkerColor(kMagenta+2);
	v2->SetLineColor(kMagenta+2);
	v2->GetYaxis()->SetTitle("differential ellipcit flow #nu_{2}");
	v2->SetStats(0);
        v2->GetYaxis()->SetTitleSize(0.07);
        v2->GetYaxis()->SetTitleOffset(0.9);
        v2->GetYaxis()->SetLabelSize(0.05);
	v2->SetTitle(" ");
	v2->Draw("h");
        h2->SetFillStyle(1);
        h2->SetMarkerStyle(8);
	h2->SetLineColor(kBlack);
        h2->SetMarkerColor(kBlack);
	h2->SetFillColor(kBlack);
	func_v2->SetLineStyle(2);
	h2->DrawCopy("hsame");//E3,E4 also nice
	c->cd();
	auto L_Ratio = new TLegend(0.2,0.9,0.8,1);
	L_Ratio->SetTextSize(0.04);
	L_Ratio->SetHeader("#bf{Pb-Pb @ #sqrt{s} = 5.02TeV}","C");	
        L_Ratio->SetFillStyle(0);
        L_Ratio->SetBorderSize(0);
        L_Ratio->Draw();
/*
        auto L_Ro = new TLegend(0.4,0.3,0.73,0.53);
        L_Ro->SetTextSize(0.028);
	L_Ro->SetHeader("#bf{f(x) = p_{0} + p_{1}*x + p_{2}*x^{2} } ","C");
        L_Ro->AddEntry(func_v2,Form("%f + %f x %f x^{2}",p0,p1,p2), "l");
        L_Ro->AddEntry(func_res,Form("%f + %f x %f x^{2}",p00,p11,p22), "l");
        L_Ro->SetFillStyle(0);
        L_Ro->SetBorderSize(0);
        L_Ro->Draw();
*/
        auto L_Ra = new TLegend(0.2,0.83,0.8,0.92);
        L_Ra->SetTextSize(0.03);
        L_Ra->SetHeader("Cent. 30%-40%, p_{T} #in [0.2,5] GeV, |#eta| #leq  0.8","C");
        L_Ra->SetFillStyle(0);
        L_Ra->SetBorderSize(0);
        L_Ra->Draw();

	auto L_R = new TLegend(0.38,0.32,0.83,0.55);
        L_R->SetHeader("#bf{f(x) = p_{0} + p_{1}*x + p_{2}*x^{2} } ","C");
	L_R->SetTextSize(0.025);
	//L_R->SetNColumns(2);
	L_R->AddEntry(v2,"AVFD: #tau_{0} = 0.3 fm/c, #tau_{B} = 0.2 fm/c", "lep");
	L_R->AddEntry(func_v2, Form("%f + %f x %f x^{2}",p0,p1,p2), "l");
	L_R->AddEntry(h2,"ALICE: |#Delta #eta| > 1", "lep");
	L_R->AddEntry(func_res, Form("%f + %f x %f x^{2}",p00,p11,p22),"l");
	L_R->SetFillStyle(0);
	L_R->SetBorderSize(0);
	L_R->Draw();

	c->cd();          // Go back to the main canvas before defining pad2
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetLeftMargin(0.2);
	pad2->SetRightMargin(0.1);
	pad2->SetBottomMargin(0.3);
	//pad2->SetGridx(); // vertical grid
	pad2->Draw();
	pad2->cd();
        h3->GetXaxis()->SetRangeUser(0.2,3.5);
	double max = h3->GetMaximum();
	double min = h3->GetMinimum();
	h3->SetLineColor(kBlack);
	//h3->GetXaxis()->SetRangeUser(0.2,3.5);
	h3->SetStats(0);      // No statistics on lower plot
	h3->SetMarkerStyle(21);
	h3->SetTitle(" "); // Remove the ratio title
	h3->GetYaxis()->SetTitle("AVFD - Fit");
	h3->GetXaxis()->SetTitle("p_{T} [GeV]");
	h3->GetYaxis()->SetNdivisions(5);
	h3->GetYaxis()->SetTitleSize(0.12);
	h3->GetYaxis()->SetTitleOffset(0.52);
	h3->GetYaxis()->SetLabelSize(0.11);
	h3->GetXaxis()->SetTitleSize(0.11);
	h3->GetXaxis()->SetTitleOffset(1.1);
	h3->GetXaxis()->SetLabelSize(0.1);
	h3->Draw("");
	c->cd();
        auto L_Ro = new TLegend(0.2,0.18,0.8,0.32);
        L_Ro->SetTextSize(0.028);
        L_Ro->SetHeader(Form("#bf{ #Delta_{min} = %f, #Delta_{max} = %f } ",min,max),"C");
        //L_Ro->AddEntry(func_v2,Form("%f + %f x %f x^{2}",p0,p1,p2), "l");
        //L_Ro->AddEntry(func_res,Form("%f + %f x %f x^{2}",p00,p11,p22), "l");
        L_Ro->SetFillStyle(0);
        L_Ro->SetBorderSize(0);
        //L_Ro->Draw();

	c->SaveAs(save);
	h3->Delete();
	h2->Delete();
	v2->Delete();
	delete c;
	func->Delete();
	func_resErr->Delete();	
	//func_v2->Delete();
	//func_v2Err->Delete();
}

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
void TuneVn(){

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
	//multi tauB plots for tau init 0.3
	TString input32 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.3/BField0.2/split");
        TString input34 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.3/BField0.4/split");
        TString input36 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.3/BField0.6/split");
cout<<"start?"<<endl;
	TString input42 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.2/split");
	TString input44 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.4/split");
        TString input40 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.0/split");
        TString input12 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.1/BField0.2/split");
        //TString input10 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.1/BField0.0/split");
        TString input92 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.9/BField0.2/split");



        //apply subsampling to all histograms from pos/neg particles======================= Check that findSpread and subsampling is associated with the t0 and not c anymore !
        findSpread(input32,3,1,"Pt",3,2);
        findSpread(input34,3,1,"Pt",3,4);
        findSpread(input36,3,1,"Pt",3,6);

	//test        
        findSpread(input42,3,1,"Pt",4,2);       
        TH1F *v42 = new TH1F("avgV42","avgV42", fPtDiffNBins, fCRCPtBins);
        v42 = Result("Pt", v42, 3, 1, 4, 2);
        FitResult(v42, "vn/v2_fit42.pdf");
/*
        findSpread(input44,3,1,"Pt",4,4);
        TH1F *v44 = new TH1F("avgV44","avgV44", fPtDiffNBins, fCRCPtBins);
        v44 = Result("Pt", v44, 3, 1, 4, 4);
        FitResult(v44, "vn/v2_fit44.pdf");
*/	
        findSpread(input40,3,1,"Pt",4,0);
        TH1F *v40 = new TH1F("avgV40","avgV40", fPtDiffNBins, fCRCPtBins);
        v40 = Result("Pt", v40, 3, 1, 4, 0);
        FitResult(v40, "vn/v2_fit40.pdf");

        findSpread(input12,3,1,"Pt",1,2);
        TH1F *v12 = new TH1F("avgV12","avgV12", fPtDiffNBins, fCRCPtBins);
        v12 = Result("Pt", v12, 3, 1, 1, 2);
        FitResult(v12, "vn/v2_fit12.pdf");

        /*findSpread(input10,3,1,"Pt",1,0);
        TH1F *v10 = new TH1F("avgV10","avgV10", fPtDiffNBins, fCRCPtBins);
        v10 = Result("Pt", v10, 3, 1, 1, 0);
        FitResult(v10, "vn/v2_fit10.pdf");
 	*/
        findSpread(input92,3,1,"Pt",9,2);
        TH1F *v92 = new TH1F("avgV92","avgV92", fPtDiffNBins, fCRCPtBins);
        v92 = Result("Pt", v92, 3, 1, 9, 2);
        FitResult(v92, "vn/v2_fit92.pdf");

        TH1F *v22 = new TH1F("avgV22","avgV22", fPtDiffNBins, fCRCPtBins);
        v22 = Result("Pt", v22, 3, 1, 3, 2);
 	FitResult(v22, "vn/v2_fit32.pdf");                                                  
        TH1F *PvpT32 = new TH1F("pvpT32", "pvpT32", fPtDiffNBins, fCRCPtBins);
        PvpT32 = Subsampling("Pt", "Pos", PvpT32, 3, 1,3,2);//cent,harm-1,tauInit,tauB
        TH1F *NvpT32 = new TH1F("nvpT32","nvpT32",fPtDiffNBins,fCRCPtBins);
        NvpT32 = Subsampling("Pt","Neg",NvpT32, 3, 1,3,2);
        TH1F *DvpT32 = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
        DvpT32 = Subsampling("Pt","Diff", DvpT32,3, 1,3,2);

        TH1F *v24 = new TH1F("v24","v24", fPtDiffNBins, fCRCPtBins);
        v24 = Result("Pt", v24, 3,1,3,4);//actually 34 !..I believe
        FitResult(v24, "vn/v2_fit34.pdf");
        TH1F *PvpT34 = new TH1F("pvpT34", "pvpT34", fPtDiffNBins, fCRCPtBins);
        PvpT34 = Subsampling("Pt", "Pos", PvpT34, 3, 1,3,4);
        TH1F *NvpT34 = new TH1F("nvpT34","nvpT34",fPtDiffNBins,fCRCPtBins);
        NvpT34 = Subsampling("Pt","Neg",NvpT34, 3 ,1,3,4);
        TH1F *DvpT34 = new TH1F("dvnpT34","dvnpT34",fPtDiffNBins,fCRCPtBins);
        DvpT34 = Subsampling("Pt","Diff", DvpT34,3,1,3,4);
        
	TH1F *v26 = new TH1F("v26","v26", fPtDiffNBins, fCRCPtBins);
        v26 = Result("Pt", v26,3,1,3,6);//actually 34 !..I believe
        FitResult(v26, "vn/v2_fit36.pdf");
        TH1F *PvpT36 = new TH1F("pvpT36", "pvpT36", fPtDiffNBins, fCRCPtBins);
        PvpT36 = Subsampling("Pt", "Pos", PvpT36, 3, 1,3,6);
        TH1F *NvpT36 = new TH1F("nvpT36","nvpT36",fPtDiffNBins,fCRCPtBins);
        NvpT36 = Subsampling("Pt","Neg",NvpT36, 3, 1,3,6);
        TH1F *DvpT36 = new TH1F("dvnpT36","dvnpT36",fPtDiffNBins,fCRCPtBins);
        DvpT36 = Subsampling("Pt","Diff", DvpT36,3, 1,3,6);
        
        auto L0 = new TLegend(0.1,0.92,0.81,1);
        L0->SetTextSize(0.025);
        L0->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 30%-40%, #tau_{0} = 0.3 fm/c, |#eta| #leq 0.8 ","C");
        L0->SetFillStyle(0);
        L0->SetBorderSize(0);

        //        TCanvas *ctest = new TCanvas("ctest", "canvastest", 800, 800);//upper plot in pad1
        //        ctest->cd();
	//	v22->Draw();
        //        ctest->SaveAs("vn/test.pdf");
	
        auto L00 = new TLegend(0.28,0.52,0.46,0.72);
        L00->SetTextSize(0.026);
        L00->AddEntry(PvpT32,Form("+h, #tau_{B} = 0.2 fm/c"), "lep");
        L00->AddEntry(PvpT34,Form("+h, #tau_{B} = 0.4 fm/c"), "lep");
        L00->AddEntry(PvpT36,Form("+h, #tau_{B} = 0.6 fm/c"), "lep");
        L00->SetBorderSize(0);
        L00->SetFillStyle(0);

        auto L11 = new TLegend(0.68,0.52,0.86,0.72);
        L11->SetTextSize(0.026);
        L11->AddEntry(NvpT32,Form("-h, #tau_{B} = 0.2 fm/c"), "lep");
        L11->AddEntry(NvpT34,Form("-h, #tau_{B} = 0.4 fm/c"), "lep");
        L11->AddEntry(NvpT36,Form("-h, #tau_{B} = 0.6 fm/c"), "lep");
        L11->SetBorderSize(0);
        L11->SetFillStyle(0);


        auto L01 = new TLegend(0.12,0.1,0.72,0.2);
        L01->SetNColumns(3);
        L01->SetTextSize(0.026);
	L01->AddEntry(DvpT32,Form("#tau_{B} = 0.2 fm/c"), "lep");
	L01->AddEntry(DvpT34,Form("#tau_{B} = 0.4 fm/c"), "lep");
	L01->AddEntry(DvpT36,Form("#tau_{B} = 0.6 fm/c"), "lep");
        L01->SetBorderSize(0);
        L01->SetFillStyle(0);


        TCanvas *m1= new TCanvas("c1", "canvas", 800, 800);//upper plot in pad1
        TPad *pad = new TPad("pad1", "pad1", 0, 0.4, 0.5, 0.99);
        pad->SetBottomMargin(0.15); // Upper and lower plot are joined
        pad->SetLeftMargin(0.21);
        pad->SetRightMargin(0);
	pad->SetGridy();
	pad->Draw();             // Draw the upper pad: pad1
        pad->cd();               // pad1 becomes the current pad//and filled with pos
        double b1 = (pad->GetWNDC())*(pad->GetHNDC());
        PvpT32->GetYaxis()->SetRangeUser(0,0.25);
        PvpT32->GetXaxis()->SetRangeUser(0.1,3.1);
	PvpT32->SetLineColor(kMagenta+1);
        PvpT32->SetLineWidth(1);
        PvpT32->GetYaxis()->SetTitle("elliptic differential flow #nu_{2}");
        PvpT32->GetYaxis()->SetTitleSize(0.1*(0.92-b1));
        PvpT32->GetYaxis()->SetLabelSize(0.1*(0.8-b1));
        PvpT32->GetYaxis()->SetTitleOffset((6*b1));
        PvpT32->GetXaxis()->SetTitle("p_{T} [GeV]");
        PvpT32->GetXaxis()->SetTitleSize(0.1*(0.8-b1));
        PvpT32->GetXaxis()->SetLabelSize(0.1*(0.8-b1));
        PvpT32->GetXaxis()->SetTitleOffset((2.9*b1));
        PvpT32->SetStats(0);
        PvpT32->SetMarkerStyle(22);
	PvpT32->SetMarkerColor(kMagenta+1);
        PvpT32->SetTitle(" ");
        PvpT32->DrawCopy();
        PvpT34->SetLineWidth(1);
	PvpT34->SetLineColor(kViolet+1);
        PvpT34->SetMarkerStyle(26);
	PvpT34->SetMarkerColor(kViolet+1);
        PvpT34->SetTitle(" ");
        PvpT34->DrawCopy("same");//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins
        PvpT36->SetLineWidth(1);
	PvpT36->SetLineColor(kMagenta+4);
        PvpT36->SetMarkerStyle(26);
	PvpT36->SetMarkerColor(kMagenta+4);
        PvpT36->SetTitle(" ");
        PvpT36->DrawCopy("same");//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins

	m1->cd();
        TPad *pa = new TPad("pa", "pa", 0.5, 0.4, 1.0, 0.99);
        pa->SetBottomMargin(0.15); // Upper and lower plot are joined
        pa->SetLeftMargin(0.);
	pa->SetRightMargin(0.21);
	pa->SetGridy();
        pa->Draw();             // Draw the upper pad: pad1
        pa->cd();               // pad1 becomes the current pad//and filled with pos
        NvpT32->GetXaxis()->SetRangeUser(0.1,3.1);
	NvpT32->GetYaxis()->SetRangeUser(0,0.25);
	NvpT32->SetLineColor(kMagenta+1);
        NvpT32->SetLineWidth(1);
        NvpT32->GetYaxis()->SetTitle(" ");
        NvpT32->GetXaxis()->SetTitle("p_{T} [GeV]");
        NvpT32->GetXaxis()->SetTitleSize(0.1*(0.8-b1));
        NvpT32->GetXaxis()->SetLabelSize(0.1*(0.8-b1));
        NvpT32->GetXaxis()->SetTitleOffset((2.8*b1));
        NvpT32->SetStats(0);
        NvpT32->SetMarkerStyle(23);
	NvpT32->SetMarkerColor(kMagenta+1);
        NvpT32->SetTitle(" ");
        NvpT32->DrawCopy();

        NvpT34->SetLineColor(kViolet+1);
        NvpT34->SetMarkerColor(kViolet+1);
        NvpT34->SetLineWidth(1);
        NvpT34->SetMarkerStyle(32);
        NvpT34->SetTitle(" ");
        NvpT34->DrawCopy("same");
        NvpT36->SetLineColor(kMagenta+4);
        NvpT36->SetMarkerColor(kMagenta+4);
        NvpT36->SetLineWidth(1);
        NvpT36->SetMarkerStyle(32);
        NvpT36->SetTitle(" ");
        NvpT36->DrawCopy("same");
	m1->cd();          // Go back to the main canvas before defining pad2
	L0->Draw();
	L00->Draw();
	L11->Draw();
        TPad *pado = new TPad("pad2", "pad2", 0, 0., 1, 0.4);
        pado->SetTopMargin(0.);
        pado->SetLeftMargin(0.1);      //pad2->SetRightMargin(4.);
        pado->SetBottomMargin(0.2);
        //pado->SetGridx(); // vertical grid
        pado->Draw();
        pado->cd();
        double b2 = (pado->GetWNDC())*(pado->GetHNDC());
        DvpT32->SetTitle(" ");
        DvpT32->SetStats(0);
	    //DvpT32->GetYaxis()->SetRangeUser(dmin,dmax);
        DvpT32->SetMarkerStyle(20);
        DvpT32->SetMarkerColor(kMagenta+1);
        DvpT32->GetXaxis()->SetRangeUser(0.1,3.1);
        DvpT32->SetLineColor(kMagenta);
        DvpT32->SetLineWidth(1);
        DvpT32->GetXaxis()->SetTitle("p_{T} [GeV]");
        DvpT32->GetYaxis()->SetTitle(" #nu_{2}(+h) - #nu_{2}(-h)");
        DvpT32->GetYaxis()->SetTitleSize(0.12*(1-b2));
        DvpT32->GetYaxis()->SetLabelSize(0.1*(1.1-b2));
        DvpT32->GetYaxis()->SetNdivisions(5);
        DvpT32->GetXaxis()->SetTitleSize(0.11*(1-b2));
        DvpT32->GetXaxis()->SetLabelSize(0.1*(1.1-b2));
        DvpT32->GetYaxis()->SetTitleOffset((1.15-b2));
        DvpT32->DrawCopy();
        DvpT34->SetMarkerStyle(24);
        DvpT34->SetMarkerColor(kViolet+1);
        DvpT34->GetXaxis()->SetRangeUser(0.1,3);
        DvpT34->SetLineColor(kViolet+1);
        DvpT34->SetLineWidth(1);
        DvpT34->DrawCopy("same");
        DvpT36->SetMarkerStyle(24);
        DvpT36->SetMarkerColor(kMagenta+4);
        DvpT36->GetXaxis()->SetRangeUser(0.1,3);
        DvpT36->SetLineColor(kMagenta+4);
        DvpT36->SetLineWidth(1);
        DvpT36->DrawCopy("same");
	    m1->cd();
	    L01->Draw();

        m1->SaveAs("vn/v2_multiB_tau0.3_pT.pdf");
        delete m1;
        PvpT32->Delete();
        NvpT32->Delete();
        DvpT32->Delete();
        PvpT34->Delete();
        NvpT34->Delete();
        DvpT34->Delete();
        PvpT36->Delete();
        NvpT36->Delete();
        DvpT36->Delete();

        findSpread(input32,3,1,"Eta",3,2);
        findSpread(input34,3,1,"Eta",3,4);
        findSpread(input36,3,1,"Eta",3,6);//c, harm-1,spec,tauInit,tauB
    
        TH1F *PvEta32 = new TH1F("pvEta32", "pvEta32", fEtaDiffNBins, fCRCEtaBins);
        PvEta32 = Subsampling("Eta", "Pos", PvEta32, 3, 1,3,2);
        TH1F *NvEta32 = new TH1F("nvEta32", "nvEta32", fEtaDiffNBins, fCRCEtaBins);
        NvEta32 = Subsampling("Eta", "Neg", NvEta32, 3, 1,3,2);
        TH1F *DvEta32 = new TH1F("deltavneta32","deltavneta32",fEtaDiffNBins,fCRCEtaBins);
        DvEta32 = Subsampling("Eta", "Diff", DvEta32, 3, 1,3,2);

        TH1F *PvEta34 = new TH1F("pvEta34", "pvEta34", fEtaDiffNBins, fCRCEtaBins);
        PvEta34 = Subsampling("Eta", "Pos", PvEta34, 3, 1,3,4);
        TH1F *NvEta34 = new TH1F("nvEta34", "nvEta34", fEtaDiffNBins, fCRCEtaBins);
        NvEta34 = Subsampling("Eta", "Neg", NvEta34, 3, 1,3,4);
        TH1F *DvEta34 = new TH1F("deltavneta34","deltavneta34",fEtaDiffNBins,fCRCEtaBins);
        DvEta34 = Subsampling("Eta", "Diff", DvEta34, 3, 1,3,4);

        TH1F *PvEta36 = new TH1F("pvEta36", "pvEta36", fEtaDiffNBins, fCRCEtaBins);
        PvEta36 = Subsampling("Eta", "Pos", PvEta36, 3, 1,3,6);
        TH1F *NvEta36 = new TH1F("nvEta36", "nvEta36", fEtaDiffNBins, fCRCEtaBins);
        NvEta36 = Subsampling("Eta", "Neg", NvEta36, 3, 1,3,6);
        TH1F *DvEta36 = new TH1F("deltavneta36","deltavneta36",fEtaDiffNBins,fCRCEtaBins);
        DvEta36 = Subsampling("Eta", "Diff", DvEta36, 3, 1,3,6);

        auto L = new TLegend(0.1, 0.92, 0.81, 1);
        L->SetTextSize(0.025);
        L->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 30%-40%, #tau_{0} = 0.3 fm/c, pT #in [0.2,5]GeV ", "C");
        L->SetFillStyle(0);
        L->SetBorderSize(0);

        auto L2 = new TLegend(0.28, 0.52, 0.46, 0.72);
        L2->SetTextSize(0.026);
        L2->AddEntry(PvEta32, Form("+h, #tau_{B} = 0.2 fm/c"), "lep");
        L2->AddEntry(PvEta34, Form("+h, #tau_{B} = 0.4 fm/c"), "lep");
        L2->AddEntry(PvEta36, Form("+h, #tau_{B} = 0.6 fm/c"), "lep");
        L2->SetBorderSize(0);
        L2->SetFillStyle(0);

        auto L1 = new TLegend(0.68, 0.52, 0.86, 0.72);
        L1->SetTextSize(0.026);
        L1->AddEntry(NvEta32, Form("-h, #tau_{B} = 0.2 fm/c"), "lep");
        L1->AddEntry(NvEta34, Form("-h, #tau_{B} = 0.4 fm/c"), "lep");
        L1->AddEntry(NvEta36, Form("-h, #tau_{B} = 0.6 fm/c"), "lep");
        L1->SetBorderSize(0);
        L1->SetFillStyle(0);

        auto L21 = new TLegend(0.12, 0.1, 0.72, 0.2);
        L21->SetNColumns(3);
        L21->SetTextSize(0.026);
        L21->AddEntry(DvEta32, Form("#tau_{B} = 0.2 fm/c"), "lep");
        L21->AddEntry(DvEta34, Form("#tau_{B} = 0.4 fm/c"), "lep");
        L21->AddEntry(DvEta36, Form("#tau_{B} = 0.6 fm/c"), "lep");
        L21->SetBorderSize(0);
        L21->SetFillStyle(0);

 	    TCanvas *m2 = new TCanvas("ceta", "etacanvas", 800, 800);
        TPad *padeta = new TPad("padeta1", "padeta1", 0, 0.4, 0.5, 0.99);
        padeta->SetBottomMargin(0.15);
        padeta->SetLeftMargin(0.21);
	padeta->SetRightMargin(0);
	padeta->SetGridy();
        padeta->Draw();
        padeta->cd();
        PvEta32->SetStats(0);
        PvEta32->GetYaxis()->SetRangeUser(0.03,0.105);
        PvEta32->GetXaxis()->SetRangeUser(-0.81,0.81);//PvEta->GetYaxis()->SetNdivisions(4) //PvEta->GetYaxis()->SetRangeUser(y-y*0.35,y+y*0.35);
        PvEta32->SetLineColor(kMagenta+1);
        PvEta32->SetLineWidth(1);
        PvEta32->GetYaxis()->SetTitle("elliptic differential flow #nu_{2}");
        PvEta32->GetYaxis()->SetTitleSize(0.1 * (0.92 - b1));
        PvEta32->GetYaxis()->SetLabelSize(0.1 * (0.8 - b1));
        PvEta32->GetYaxis()->SetTitleOffset((6 * b1));
	    PvEta32->GetXaxis()->SetTitle("#eta");
	    PvEta32->GetXaxis()->SetTitleSize(0.1 * (0.8 - b1));
	    PvEta32->GetXaxis()->SetLabelSize(0.1 * (0.8 - b1));
	    PvEta32->GetXaxis()->SetTitleOffset((2.99 * b1));
        PvEta32->SetStats(0);
	PvEta32->GetXaxis()->SetNdivisions(4);
        PvEta32->SetMarkerStyle(22);
        PvEta32->SetMarkerColor(kMagenta+1);
        PvEta32->SetTitle(" ");
        PvEta32->DrawCopy();      //NvEta->GetYaxis()->SetNdivisions(4);
        PvEta34->SetLineColor(kViolet+1);
        PvEta34->SetLineWidth(1);
        PvEta34->SetMarkerStyle(26);
        PvEta34->SetMarkerColor(kViolet+1);
        PvEta34->SetTitle(" ");
        PvEta34->DrawCopy("same"); 
        PvEta36->SetLineColor(kMagenta+4);
        PvEta36->SetLineWidth(1);
        PvEta36->SetMarkerStyle(26);
        PvEta36->SetMarkerColor(kMagenta+4);
        PvEta36->SetTitle(" ");
        PvEta36->DrawCopy("same");      //NvEta->GetYaxis()->SetNdivisions(4);
        m2->cd();

        TPad *padeta2 = new TPad("padeta2", "padeta2", 0.5, 0.4, 1, 0.99);
        padeta2->SetBottomMargin(0.15);
        padeta2->SetLeftMargin(0);
        padeta2->SetRightMargin(0.21);
        padeta2->SetGridy();
        padeta2->Draw();
        padeta2->cd();
	NvEta32->SetStats(0);
        NvEta32->GetXaxis()->SetRangeUser(-0.81, 0.81);
        NvEta32->GetYaxis()->SetRangeUser(0.03, 0.105);
        NvEta32->GetYaxis()->SetTitle(" ");
        NvEta32->GetXaxis()->SetTitle("#eta");
        NvEta32->GetXaxis()->SetTitleSize(0.1 * (0.8 - b1));
        NvEta32->GetXaxis()->SetLabelSize(0.1 * (0.8 - b1));
        NvEta32->GetXaxis()->SetTitleOffset((2.99 * b1));
	NvEta32->GetXaxis()->SetNdivisions(4);
        NvEta32->SetLineColor(kMagenta+1);
        NvEta32->SetLineWidth(1);
        NvEta32->SetTitle(" ");
        NvEta32->SetMarkerStyle(23);
        NvEta32->SetMarkerColor(kMagenta+1);
        NvEta32->DrawCopy();
        NvEta34->SetLineWidth(1);
	    NvEta34->SetLineColor(kViolet+1);
        NvEta34->SetTitle(" ");
        NvEta34->SetMarkerStyle(32); 
	    NvEta34->SetMarkerColor(kViolet+1);
        NvEta34->DrawCopy("same"); 
        NvEta36->SetLineColor(kMagenta+4);
        NvEta36->SetLineWidth(1);
        NvEta36->SetTitle(" ");
        NvEta36->SetMarkerStyle(32);
        NvEta36->SetMarkerColor(kMagenta+4);
        NvEta36->DrawCopy("same");
        m2->cd();
        L->Draw();
        L1->Draw();
        L2->Draw();
        TPad *padetao = new TPad("padeta2", "padeta2", 0, 0., 1, 0.4);
        padetao->SetTopMargin(0);
        padetao->SetLeftMargin(0.1);       //padeta2->SetRightMargin(4.);
        padetao->SetBottomMargin(0.2);
        padetao->Draw();
        padetao->cd();
        DvEta32->SetStats(0);
        DvEta32->GetYaxis()->SetRangeUser(-0.025,0.0081);
        DvEta32->GetXaxis()->SetRangeUser(-0.8,0.8);
        DvEta32->GetXaxis()->SetTitle("#eta");
        DvEta32->SetTitle(" ");
        DvEta32->SetMarkerStyle(20);
        DvEta32->SetMarkerColor(kMagenta+1);
        DvEta32->SetLineColor(kMagenta+1);
        DvEta32->SetLineWidth(1);
        DvEta32->GetYaxis()->SetTitle("#nu_{2}(+h) - #nu_{2}(-h)");
        DvEta32->GetYaxis()->SetTitleSize(0.12 * (1 - b2));
        DvEta32->GetYaxis()->SetLabelSize(0.1 * (1.1 - b2));
        DvEta32->GetYaxis()->SetNdivisions(5);
        DvEta32->GetXaxis()->SetTitleSize(0.11 * (1 - b2));
        DvEta32->GetXaxis()->SetLabelSize(0.1 * (1.1 - b2));
        DvEta32->GetYaxis()->SetTitleOffset((1.15 - b2));
        DvEta32->GetYaxis()->SetNdivisions(5);
        DvEta32->Draw();
        DvEta34->SetMarkerStyle(24);
        DvEta34->SetMarkerColor(kViolet+1);
        DvEta34->SetLineColor(kViolet+1);
        DvEta34->SetLineWidth(1);
        DvEta34->DrawCopy("same");
        DvEta36->SetMarkerStyle(24);
        DvEta36->SetMarkerColor(kMagenta+4);
        DvEta36->SetLineColor(kMagenta+4);
        DvEta36->SetLineWidth(1);
        DvEta36->DrawCopy("same");
        m2->cd();
        L21->Draw();
        m2->SaveAs("vn/v2_multiB_tau0.3_eta.pdf");
        delete m2;
        PvEta32->Delete();
        NvEta32->Delete();
        DvEta32->Delete();
        PvEta34->Delete();
        NvEta34->Delete();
        DvEta34->Delete();
        PvEta36->Delete();
        NvEta36->Delete();
        DvEta36->Delete();


//Subsampling final results for the multi tauInit plot
	for(Int_t harm=0; harm<=2;harm++){
		//TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
		
		TString input1 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.1/BField0.2/split");
		TString input3 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.3/BField0.2/split");
		TString input6 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.6/BField0.2/split");
		TString input9 = Form("/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.9/BField0.2/split");
		
		// !!! For harm = 0 we could also try to use the fFlowQCQv1[pos=0, neg=1][cen][pt = 0, eta = 1] !!!
		cout<<"Harmonic: "<<(harm+1)<<endl;
		//apply subsampling to all histograms from pos/neg particles======================= Check that findSpread and subsampling is associated with the t0 and not c anymore !
		findSpread(input1,3,harm,"Pt",1,2);
		findSpread(input3,3,harm,"Pt",3,2);
                findSpread(input6,3,harm,"Pt",6,2);
                findSpread(input9,3,harm,"Pt",9,2);
		cout<<"spread found"<<endl;
		TH1F *PvpT1 = new TH1F("pvpT1", "pvpT1", fPtDiffNBins, fCRCPtBins);
		PvpT1 = Subsampling("Pt", "Pos", PvpT1, 3, harm,1,2);
		TH1F *NvpT1 = new TH1F("nvpT1","nvpT1",fPtDiffNBins,fCRCPtBins);
		NvpT1 = Subsampling("Pt","Neg",NvpT1, 3, harm,1,2);
		TH1F *DvpT1 = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
		DvpT1 = Subsampling("Pt","Diff", DvpT1,3, harm,1,2);

        TH1F *PvpT3 = new TH1F("pvpT3", "pvpT3", fPtDiffNBins, fCRCPtBins);
        PvpT3 = Subsampling("Pt", "Pos", PvpT3, 3, harm,3,2);
        TH1F *NvpT3 = new TH1F("nvpT3","nvpT3",fPtDiffNBins,fCRCPtBins);
        NvpT3 = Subsampling("Pt","Neg",NvpT3, 3 ,harm,3,2);
        TH1F *DvpT3 = new TH1F("dvnpT3","dvnpT3",fPtDiffNBins,fCRCPtBins);
        DvpT3 = Subsampling("Pt","Diff", DvpT3,3,harm,3,2);
	cout<<"tau init 0.3 subsampled"<<endl;
        TH1F *PvpT6 = new TH1F("pvpT6", "pvpT6", fPtDiffNBins, fCRCPtBins);
        PvpT6 = Subsampling("Pt", "Pos", PvpT6, 3, harm,6,2);
        TH1F *NvpT6 = new TH1F("nvpT6","nvpT6",fPtDiffNBins,fCRCPtBins);
        NvpT6 = Subsampling("Pt","Neg",NvpT6, 3, harm,6,2);
        TH1F *DvpT6 = new TH1F("dvnpT6","dvnpT6",fPtDiffNBins,fCRCPtBins);
        DvpT6 = Subsampling("Pt","Diff", DvpT6,3, harm,6,2);
	cout<<"tau 0.6"<<endl;

        TH1F *PvpT9 = new TH1F("pvpT9", "pvpT9", fPtDiffNBins, fCRCPtBins);
        PvpT9 = Subsampling("Pt", "Pos", PvpT9, 3, harm,9,2);
        TH1F *NvpT9 = new TH1F("nvpT9","nvpT9",fPtDiffNBins,fCRCPtBins);
        NvpT9 = Subsampling("Pt","Neg",NvpT9, 3 ,harm,9,2);
        TH1F *DvpT9 = new TH1F("dvnpT9","dvnpT9",fPtDiffNBins,fCRCPtBins);
        DvpT9 = Subsampling("Pt","Diff", DvpT9,3,harm,9,2);
        
		//Multi centrality plot
		//copy[0] = (TH1F*) DvpT1->Clone("copy dvpt1");
		//copy[1] = (TH1F*) DvpT4->Clone("copy dvpt4");
		//The final pT differential plots================================================================================

        	short Col1,Col2;
        	double ymin, ymax,dmin,dmax;
        	double x1,x2,y1,y2;
        	TString flow;
        	if(harm == 0){
        	Col1 = kOrange;
        	Col2 = kRed;
        	ymin = -0.005;
        	ymax = 0.03;
        	dmin = -0.012;
        	dmax = 0.046;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.65;
        	y2 = 0.89;
        	flow = "directed";}
        	if(harm == 1){
        	Col1 = kMagenta;
        	Col2 = kViolet;
        	ymin = -0.001;
        	ymax = 0.27;
        	dmin = -0.061;
        	dmax = 0.041;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.65;
        	y2 = 0.89;
        	flow = "elliptic";}
        	if(harm == 2){
        	Col1 = kGreen;
        	Col2 = kBlue;
        	ymin = -0.01;
        	ymax = 0.25;
        	dmin = -0.074;
        	dmax = 0.054;
        	x1 = 0.39;
        	x2 = 0.89;
        	y1 = 0.01;
        	y2 = 0.25;
        	flow = "triangular";}

        	auto l = new TLegend(0.1, 0.92, 0.81, 1);
        	l->SetTextSize(0.025);
        	l->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 30%-40%, #tau_{B} = 0.2 fm/c, |#eta| #leq 0.8 ", "C");
        	l->SetFillStyle(0);
        	l->SetBorderSize(0);
                auto l1 = new TLegend(0.1, 0.8, 0.3, 0.95);
                l1->SetTextSize(0.026);
                l1->SetHeader("#tau_{0} = 0.1 fm/c", "C");
                l1->AddEntry(PvpT1, Form("+h"), "lep");
                l1->AddEntry(NvpT1, Form("-h"), "lep");
                l1->SetBorderSize(0);
                l1->SetFillStyle(0);
                auto l6 = new TLegend(0.38, 0.8, 0.58, 0.95);
                l6->SetTextSize(0.026);
                l6->SetHeader("#tau_{0} = 0.6 fm/c", "C");
                l6->AddEntry(PvpT6, Form("+h"), "lep");
                l6->AddEntry(NvpT6, Form("-h"), "lep");
                l6->SetBorderSize(0);
                l6->SetFillStyle(0);
                auto l9 = new TLegend(0.68, 0.8, 0.88, 0.95);
                l9->SetTextSize(0.026);
                l9->SetHeader("#tau_{0} = 0.9 fm/c", "C");
                l9->AddEntry(PvpT9, Form("+h"), "lep");
                l9->AddEntry(NvpT9, Form("-h"), "lep");
                l9->SetBorderSize(0);
                l9->SetFillStyle(0);
       		
		auto d = new TLegend(0.1,0.3,0.7,0.4);
                d->SetTextSize(0.026);
		d->SetNColumns(3);
                d->AddEntry(DvpT1, "#tau_{0} = 0.1 fm/c", "lep");
                d->AddEntry(DvpT6, "#tau_{0} = 0.6 fm/c", "lep");
                d->AddEntry(DvpT9, "#tau_{0} = 0.9 fm/c", "lep");
                d->SetBorderSize(0);
                d->SetFillStyle(0);

		TCanvas *c1= new TCanvas("c1", "canvas", 800, 800);//upper plot in pad1
        	TPad *pad1 = new TPad("pad1", "pad1", 0., 0.4, 0.38, 0.95);
        	pad1->SetLeftMargin(0.25);
		pad1->SetRightMargin(0);
            	pad1->SetBottomMargin(0.15);
		pad1->SetGridy();
        	pad1->Draw();             // Draw the upper pad: pad1
        	pad1->cd();               // pad1 becomes the current pad//and filled with pos
        	double f1 = (pad1->GetWNDC())*(pad1->GetHNDC());
		PvpT1->SetTitle(" ");
		PvpT1->GetYaxis()->SetRangeUser(ymin,ymax);
		PvpT1->GetXaxis()->SetRangeUser(0.1,3.1);
        	PvpT1->SetLineColor(Col1+1);
        	PvpT1->SetLineWidth(1);
        	PvpT1->GetYaxis()->SetTitle("differential "+flow+Form(" flow #nu_{%d}",harm+1));
       		PvpT1->GetYaxis()->SetTitleSize(0.12*(0.8-f1));
        	PvpT1->GetYaxis()->SetLabelSize(0.1*(0.8-f1));
        	PvpT1->GetYaxis()->SetTitleOffset((8*f1));
                PvpT1->GetXaxis()->SetLabelSize(0.1 * (0.9 - f1));
                PvpT1->GetXaxis()->SetTitleSize(0.1 * (0.8 - f1));
                PvpT1->GetXaxis()->SetTitleOffset((6 * f1));
		PvpT1->GetXaxis()->SetTitle("p_{T} [GeV]");
                PvpT1->SetStats(0);
        	PvpT1->SetMarkerStyle(26);
        	PvpT1->SetMarkerColor(Col1+1);
        	PvpT1->DrawCopy();//NvpT1->GetXaxis()->SetRangeUser(0.,3);//range in bins
        	NvpT1->SetLineColor(Col1+1);
        	NvpT1->SetLineWidth(1);
        	NvpT1->SetMarkerStyle(23);
        	NvpT1->SetMarkerColor(Col1+1);
          	NvpT1->DrawCopy("same");
                c1->cd();
                TPad* pad11 = new TPad("pad11", "pad11", 0.38, 0.4, 0.68, 0.95);
                pad11->SetLeftMargin(0);
	 	pad11->SetRightMargin(0);
                pad11->SetBottomMargin(0.15);
		pad11->SetGridy();
                pad11->Draw();             // Draw the upper pad: pad1
                pad11->cd();               // pad1 becomes the current pad//and filled with po
		f1 =(pad11->GetWNDC())*(pad11->GetHNDC());
                PvpT6->SetStats(0);
                PvpT6->GetYaxis()->SetRangeUser(ymin, ymax);
                PvpT6->GetXaxis()->SetRangeUser(0.1, 3.1);
                PvpT6->GetXaxis()->SetLabelSize(0.1 * (1 - f1));
                PvpT6->GetXaxis()->SetTitleSize(0.1 * (0.9 - f1));
                PvpT6->GetXaxis()->SetTitleOffset((6 * f1));
                PvpT6->GetXaxis()->SetTitle("p_{T} [GeV]");
                PvpT6->SetLineColor(Col2+1);
                PvpT6->SetLineWidth(1);
                PvpT6->SetMarkerStyle(26);
                PvpT6->SetMarkerColor(Col2+1);
                PvpT6->SetTitle(" ");
                PvpT6->DrawCopy();//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins
                NvpT6->SetLineColor(Col2+1);
                NvpT6->SetLineWidth(1);
                NvpT6->SetMarkerStyle(23);
                NvpT6->SetMarkerColor(Col2+1);
                NvpT6->SetTitle(" ");
                NvpT6->DrawCopy("same");
		c1->cd();
                TPad* pad12 = new TPad("pad12", "pad12", 0.68, 0.4, 1, 0.95);
                pad12->SetLeftMargin(0);
                pad12->SetRightMargin(0.05);
                pad12->SetBottomMargin(0.15);
		pad12->SetGridy();
                pad12->Draw();             // Draw the upper pad: pad1
                pad12->cd();
		PvpT9->SetStats(0);
                PvpT9->GetYaxis()->SetRangeUser(ymin, ymax);
                PvpT9->GetXaxis()->SetRangeUser(0.1, 3.1);
                PvpT9->GetXaxis()->SetLabelSize(0.1 * (1 - f1));
                PvpT9->GetXaxis()->SetTitleSize(0.1 * (0.9 - f1));
                PvpT9->GetXaxis()->SetTitleOffset((6 * f1));
                PvpT9->GetXaxis()->SetTitle("p_{T} [GeV]");
            	PvpT9->SetLineColor(Col1+4);
            	PvpT9->SetLineWidth(1);
            	PvpT9->SetMarkerStyle(26);
            	PvpT9->SetMarkerColor(Col1+4);
            	PvpT9->SetTitle(" ");
            	PvpT9->DrawCopy();//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins
            	NvpT9->SetLineColor(Col1+4);
            	NvpT9->SetLineWidth(1);
            	NvpT9->SetMarkerStyle(23);
            	NvpT9->SetMarkerColor(Col1+4);
            	NvpT9->SetTitle(" ");
            	NvpT9->DrawCopy("same");
                /*
		PvpT9->SetLineColor(Col2+4);
                PvpT9->SetLineWidth(1);
                PvpT9->SetMarkerStyle(26);
                PvpT9->SetMarkerColor(Col2+4);
                PvpT9->SetTitle(" ");
                PvpT9->DrawCopy("same");//NvpT4->GetXaxis()->SetRangeUser(0.,3);//range in bins
                NvpT9->SetLineColor(Col2+4);
                NvpT9->SetLineWidth(1);
                NvpT9->SetMarkerStyle(32);
                NvpT9->SetMarkerColor(Col2+4);
                NvpT9->SetTitle(" ");
                NvpT9->DrawCopy("same");*/

        	c1->cd();          // Go back to the main canvas before defining pad2
		l->Draw();
		l1->Draw();
		l6->Draw();
		l9->Draw();
        	TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0., 1, 0.4);
        	pad2->SetTopMargin(0);
        	pad2->SetLeftMargin(0.98);      //pad2->SetRightMargin(4.);
        	pad2->SetBottomMargin(0.2);
		pad2->SetRightMargin(0.01);
        	pad2->Draw();
        	pad2->cd();
		double f2 = (pad2->GetWNDC())*(pad2->GetHNDC());
        	DvpT1->SetTitle(" ");
        	DvpT1->SetStats(0);
		DvpT1->GetYaxis()->SetRangeUser(dmin,dmax);
        	DvpT1->SetMarkerStyle(20);
        	DvpT1->SetMarkerColor(Col1+1);
        	DvpT1->GetXaxis()->SetRangeUser(0.1,3.1);
        	DvpT1->SetLineColor(Col1+1);
        	DvpT1->SetLineWidth(1);
        	DvpT1->GetXaxis()->SetTitle("p_{T} [GeV]");
        	DvpT1->GetYaxis()->SetTitle(Form("difference #nu_{%d}(+h) - #nu_{%d}(-h)",harm+1,harm+1));
        	DvpT1->GetYaxis()->SetTitleSize(0.12*(1-f2));
        	DvpT1->GetYaxis()->SetLabelSize(0.1*(1-f2));
	        DvpT1->GetYaxis()->SetNdivisions(5);
        	DvpT1->GetXaxis()->SetTitleSize(0.12*(0.9-f2));
	        DvpT1->GetXaxis()->SetLabelSize(0.1*(1-f2));
	        DvpT1->GetYaxis()->SetTitleOffset((1.1-f2));
		DvpT1->DrawCopy();  
                DvpT6->SetMarkerStyle(24);
                DvpT6->SetMarkerColor(Col2+1);
                DvpT6->GetXaxis()->SetRangeUser(0.1,31);
                DvpT6->SetLineColor(Col2+1);
                DvpT6->SetLineWidth(1);
                DvpT6->DrawCopy("same");
                DvpT9->SetMarkerStyle(20);
                DvpT9->SetMarkerColor(Col1+4);
                DvpT9->GetXaxis()->SetRangeUser(0.1,31);
                DvpT9->SetLineColor(Col1+4);
                DvpT9->SetLineWidth(1);
                DvpT9->DrawCopy("same");
		c1->cd();
		d->Draw();
        	c1->SaveAs(Form("vn/v%d_multiTauIn_pT.pdf",harm+1));
       		delete c1;

                PvpT1->Delete();
                NvpT1->Delete();
                DvpT1->Delete();
        	PvpT3->Delete();
        	NvpT3->Delete();
        	DvpT3->Delete();
	        PvpT6->Delete();
                NvpT6->Delete();
                DvpT6->Delete();
                PvpT9->Delete();
                NvpT9->Delete();
                DvpT9->Delete();

		//Now for eta differential flow====================================================================================
                findSpread(input1,3,harm,"Eta",1,2);
                findSpread(input3,3,harm,"Eta",3,2);
                findSpread(input6,3,harm,"Eta",6,2);
                findSpread(input9,3,harm,"Eta",9,2);

		TH1F *PvEta1 = new TH1F("pvEta", "pvEta", fEtaDiffNBins, fCRCEtaBins);
                PvEta1 = Subsampling("Eta", "Pos", PvEta1, 3, harm,1,2);
                TH1F *NvEta1 = new TH1F("nvEta1", "nvEta1", fEtaDiffNBins, fCRCEtaBins);
                NvEta1 = Subsampling("Eta", "Neg", NvEta1, 3, harm,1,2);
                TH1F *DvEta1 = new TH1F("deltavneta1","deltavneta1",fEtaDiffNBins,fCRCEtaBins);
                DvEta1 = Subsampling("Eta", "Diff", DvEta1, 3, harm,1,2);
	
		TH1F *PvEta3 = new TH1F("pvEta3", "pvEta3", fEtaDiffNBins, fCRCEtaBins);
		PvEta3 = Subsampling("Eta", "Pos", PvEta3, 3, harm,3,2);		
		TH1F *NvEta3 = new TH1F("nvEta3", "nvEta3", fEtaDiffNBins, fCRCEtaBins);
		NvEta3 = Subsampling("Eta", "Neg", NvEta3, 3, harm,3,2);
		TH1F *DvEta3 = new TH1F("deltavneta3","deltavneta3",fEtaDiffNBins,fCRCEtaBins);
		DvEta3 = Subsampling("Eta", "Diff", DvEta3, 3, harm,3,2);

                TH1F *PvEta6 = new TH1F("pvEta6", "pvEta6", fEtaDiffNBins, fCRCEtaBins);
                PvEta6 = Subsampling("Eta", "Pos", PvEta6, 3, harm,6,2);
                TH1F *NvEta6 = new TH1F("nvEta6", "nvEta1", fEtaDiffNBins, fCRCEtaBins);
                NvEta6 = Subsampling("Eta", "Neg", NvEta6, 3, harm,6,2);
                TH1F *DvEta6 = new TH1F("deltavneta6","deltavneta6",fEtaDiffNBins,fCRCEtaBins);
                DvEta6 = Subsampling("Eta", "Diff", DvEta6, 3, harm,6,2);

                TH1F *PvEta9 = new TH1F("pvEta9", "pvEta9", fEtaDiffNBins, fCRCEtaBins);
                PvEta9 = Subsampling("Eta", "Pos", PvEta9, 3, harm,9,2);
                TH1F *NvEta9 = new TH1F("nvEta9", "nvEta9", fEtaDiffNBins, fCRCEtaBins);
                NvEta9 = Subsampling("Eta", "Neg", NvEta9, 3, harm,9,2);
                TH1F *DvEta9 = new TH1F("deltavneta9","deltavneta9",fEtaDiffNBins,fCRCEtaBins);
                DvEta9 = Subsampling("Eta", "Diff", DvEta9, 3, harm,9,2);
                
		//Plot the pos/neg with difference=================================================
		//PlotEta(PvEta1, NvEta1, DvEta1, 1, harm);
                //PlotEta(PvEta4, NvEta4, DvEta4, 4, harm);

		//The final Eta differential plots================================================================================== 	
		if(harm == 0){
        	Col1 = kOrange;
        	Col2 = kRed;
        	ymin = -0.001;
        	ymax = 0.0089;
        	dmin = -0.011;
	        dmax = 0.011;
        	x1 = 0.22;
	        x2 = 0.72;
        	y1 = 0.65;
	        y2 = 0.89;
		flow = "directed";}
	        if(harm == 1){
        	Col1 = kMagenta;
        	Col2 = kViolet;
        	ymin = 0.05;
        	ymax = 0.095;
        	dmin = -0.032;
        	dmax = 0.022;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.36;
        	y2 = 0.60;
		flow = "elliptic";}
        	if(harm == 2){
        	Col1 = kGreen;
        	Col2 = kBlue;
        	ymin = 0.02;
        	ymax = 0.15;
        	dmin = -0.082;
        	dmax = 0.022;
        	x1 = 0.22;
        	x2 = 0.72;
        	y1 = 0.64;
        	y2 = 0.88;
		flow = "triangular";}

        	auto le = new TLegend(0.1, 0.92, 0.81, 1);
        	le->SetTextSize(0.025);
        	le->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 30%-40%, #tau_{B} = 0.2 fm/c, p_{T} #in [0.2,5]GeV ", "C");
        	le->SetFillStyle(0);
        	le->SetBorderSize(0);
                auto le1 = new TLegend(0.1, 0.85, 0.3, 0.95);
                le1->SetTextSize(0.026);
		le1->SetNColumns(2);
                le1->SetHeader("#tau_{0} = 0.1 fm/c", "C");
                le1->AddEntry(PvEta1, Form("+h"), "lep");
                le1->AddEntry(NvEta1, Form("-h"), "lep");
                le1->SetBorderSize(0);
                le1->SetFillStyle(0);
                auto le6 = new TLegend(0.38, 0.85, 0.58, 0.95);
                le6->SetTextSize(0.026);
		le6->SetNColumns(2);
                le6->SetHeader("#tau_{0} = 0.6 fm/c", "C");
                le6->AddEntry(PvEta6, Form("+h"), "lep");
                le6->AddEntry(NvEta6, Form("-h"), "lep");
                le6->SetBorderSize(0);
                le6->SetFillStyle(0);
                auto le9 = new TLegend(0.68, 0.85, 0.88, 0.95);
                le9->SetTextSize(0.026);
		le9->SetNColumns(2);
                le9->SetHeader("#tau_{0} = 0.9 fm/c", "C");
                le9->AddEntry(PvEta9, Form("+h"), "lep");
                le9->AddEntry(NvEta9, Form("-h"), "lep");
                le9->SetBorderSize(0);
                le9->SetFillStyle(0);
                auto de = new TLegend(0.1,0.08,0.7,0.18);
                de->SetTextSize(0.026);
                de->SetNColumns(3);
                de->AddEntry(DvEta1, "#tau_{0} = 0.1 fm/c", "lep");
                de->AddEntry(DvEta6, "#tau_{0} = 0.6 fm/c", "lep");
                de->AddEntry(DvEta9, "#tau_{0} = 0.9 fm/c", "lep");
                de->SetBorderSize(0);
                de->SetFillStyle(0);

        	TCanvas *cEta = new TCanvas("ceta", "etacanvas", 800, 800);
        	TPad *padeta1 = new TPad("padeta1", "padeta1",  0., 0.4, 0.38, 0.95);
        	padeta1->SetLeftMargin(0.25);
		padeta1->SetRightMargin(0);
            	padeta1->SetBottomMargin(0.15);
		padeta1->SetGridy();
        	padeta1->Draw();
        	padeta1->cd();
		double f3 = (padeta1->GetWNDC())*(padeta1->GetHNDC());
		PvEta1->GetYaxis()->SetRangeUser(ymin,ymax);
        	PvEta1->GetXaxis()->SetRangeUser(-0.8,0.8);//PvEta->GetYaxis()->SetNdivisions(4) //PvEta->GetYaxis()->SetRangeUser(y-y*0.35,y+y*0.35);
        	PvEta1->SetLineColor(Col1+1);
        	PvEta1->SetLineWidth(1);
        	PvEta1->GetYaxis()->SetTitle("differential "+flow+" flow "+Form("#nu_{%d}",harm+1));
        	PvEta1->GetYaxis()->SetTitleSize(0.12*(0.8-f3));
		PvEta1->GetYaxis()->SetLabelSize(0.1*(0.8-f3));
		PvEta1->GetYaxis()->SetTitleOffset((8*f3));//PvEta1->GetYaxis()->CenterTitle(true);
		PvEta1->SetStats(0);//PvEta1->GetYaxis()->SetNdivisions(6);
                PvEta1->GetXaxis()->SetLabelSize(0.1 * (0.9 - f3));
                PvEta1->GetXaxis()->SetTitleSize(0.1 * (0.8 - f3));
                PvEta1->GetXaxis()->SetTitleOffset((3 * f3));
		PvEta1->GetXaxis()->SetTitle("#eta");
		PvEta1->GetXaxis()->SetNdivisions(4);
        	PvEta1->SetMarkerStyle(26);
        	PvEta1->SetMarkerColor(Col1+1);
        	PvEta1->SetTitle(" ");
        	PvEta1->DrawCopy();      //NvEta->GetYaxis()->SetNdivisions(4);
        	NvEta1->SetLineColor(Col1+1);
        	NvEta1->SetLineWidth(1);
        	NvEta1->SetTitle(" ");
        	NvEta1->SetMarkerStyle(23);
        	NvEta1->SetMarkerColor(Col1+1);
        	NvEta1->DrawCopy("same");
                cEta->cd();
                TPad* padE11 = new TPad("padE11", "padE11", 0.38, 0.4, 0.68, 0.95);
                padE11->SetLeftMargin(0);
	 	padE11->SetRightMargin(0);
                padE11->SetBottomMargin(0.15);
		padE11->SetGridy();
                padE11->Draw();             // Draw the upper pad: pad1
                padE11->cd();               // pad1 becomes the current pad//and filled with po
		f3 =(padE11->GetWNDC())*(padE11->GetHNDC());
                PvEta6->SetStats(0);
                PvEta6->GetYaxis()->SetRangeUser(ymin, ymax);
                PvEta6->GetXaxis()->SetRangeUser(-0.8, 0.8);
                PvEta6->GetXaxis()->SetLabelSize(0.1 * (1 - f3));
                PvEta6->GetXaxis()->SetTitleSize(0.1 * (0.9 - f3));
                PvEta6->GetXaxis()->SetTitleOffset((3 * f3));
                PvEta6->GetXaxis()->SetTitle("#eta");
		PvEta6->GetXaxis()->SetNdivisions(4);
                PvEta6->SetLineColor(Col2+1);
                PvEta6->SetLineWidth(1);
                PvEta6->SetMarkerStyle(26);
                PvEta6->SetMarkerColor(Col2+1);
                PvEta6->SetTitle(" ");
                PvEta6->DrawCopy();      //NvEta->GetYaxis()->SetNdivisions(4);
                NvEta6->SetLineColor(Col2+1);
                NvEta6->SetLineWidth(1);
                NvEta6->SetTitle(" ");
                NvEta6->SetMarkerStyle(23);
                NvEta6->SetMarkerColor(Col2+1);
                NvEta6->DrawCopy("same");
		cEta->cd();
                TPad* padE12 = new TPad("padE12", "padE12", 0.68, 0.4, 1, 0.95);
                padE12->SetLeftMargin(0);
                padE12->SetRightMargin(0.05);
                padE12->SetBottomMargin(0.15);
		padE12->SetGridy();
                padE12->Draw();             // Draw the upper pad: pad1
                padE12->cd();
		PvEta9->SetStats(0);
                PvEta9->GetYaxis()->SetRangeUser(ymin, ymax);
                PvEta9->GetXaxis()->SetRangeUser(-0.8, 0.8);
                PvEta9->GetXaxis()->SetLabelSize(0.1 * (1 - f3));
                PvEta9->GetXaxis()->SetTitleSize(0.1 * (0.9 - f3));
                PvEta9->GetXaxis()->SetTitleOffset((3 * f3));
                PvEta9->GetXaxis()->SetTitle("#eta");
		PvEta9->GetXaxis()->SetNdivisions(4);
                PvEta9->SetLineColor(Col1+4);
                PvEta9->SetLineWidth(1);
                PvEta9->SetMarkerStyle(26);
                PvEta9->SetMarkerColor(Col1+4);
                PvEta9->SetTitle(" ");
                PvEta9->DrawCopy();      //NvEta->GetYaxis()->SetNdivisions(4);
                NvEta9->SetLineColor(Col1+4);
                NvEta9->SetLineWidth(1);
                NvEta9->SetTitle(" ");
                NvEta9->SetMarkerStyle(23);
                NvEta9->SetMarkerColor(Col1+4);
                NvEta9->DrawCopy("same");

        	cEta->cd();
                le->Draw();
		le1->Draw();
		le6->Draw();
		le9->Draw();
        	TPad *padeta2 = new TPad("padeta2", "padeta2", 0.01, 0., 1, 0.4);
        	padeta2->SetTopMargin(0);
        	padeta2->SetLeftMargin(0.97);       //padeta2->SetRightMargin(4.);
        	padeta2->SetBottomMargin(0.2);
		padeta2->SetRightMargin(0.015);
        	padeta2->Draw();
        	padeta2->cd();
		double f4 = (padeta2->GetWNDC())*(padeta2->GetHNDC());
        	DvEta1->SetTitle(" ");
        	DvEta1->SetStats(0);
		DvEta1->GetYaxis()->SetRangeUser(dmin,dmax);
        	DvEta1->SetMarkerStyle(20);
        	DvEta1->SetMarkerColor(Col1+1);
        	DvEta1->GetXaxis()->SetRangeUser(-0.8,0.8);
        	DvEta1->SetLineColor(Col1+1);
        	DvEta1->SetLineWidth(1);
        	DvEta1->GetXaxis()->SetTitle("#eta");
        	DvEta1->GetYaxis()->SetTitle(Form("difference #nu_{%d}(+h) - #nu_{%d}(-h)",harm+1,harm+1));
        	DvEta1->GetYaxis()->SetTitleSize(0.12*(1-f4));
        	DvEta1->GetYaxis()->SetLabelSize(0.1*(1-f4));
	        DvEta1->GetYaxis()->SetNdivisions(5);
        	DvEta1->GetXaxis()->SetTitleSize(0.12*(0.9-f4));
	        DvEta1->GetXaxis()->SetLabelSize(0.1*(1-f4));
	        DvEta1->GetYaxis()->SetTitleOffset((1.1-f4));
		DvEta1->DrawCopy();
                DvEta6->SetMarkerStyle(24);
                DvEta6->SetMarkerColor(Col2+1);
                DvEta6->GetXaxis()->SetRangeUser(-0.8,0.8);
                DvEta6->SetLineColor(Col2+1);
                DvEta6->SetLineWidth(1);
                DvEta6->DrawCopy("same");
                DvEta9->SetMarkerStyle(20);
                DvEta9->SetMarkerColor(Col1+4);
                DvEta9->GetXaxis()->SetRangeUser(-0.8,0.8);
                DvEta9->SetLineColor(Col1+4);
                DvEta9->SetLineWidth(1);
                DvEta9->DrawCopy("same");
		cEta->cd();
		de->Draw();
        	cEta->SaveAs(Form("vn/v%d_mutliTauIn_eta.pdf",harm+1));
        	delete cEta;

	
        	PvEta1->Delete();
	       	NvEta1->Delete();
	        DvEta1->Delete();
                PvEta3->Delete();
                NvEta3->Delete();
                DvEta3->Delete();
                PvEta6->Delete();
                NvEta6->Delete();
                DvEta6->Delete();
                PvEta9->Delete();
                NvEta9->Delete();
                DvEta9->Delete();


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


