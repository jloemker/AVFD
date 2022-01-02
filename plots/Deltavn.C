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
//
//======================================================================================
TH1F *Subsampling(TString file, TString obj, TH1F *result, Int_t c, Int_t harm){
TFile *f;
TH1F *co, *h;
TH1D *Mean = new TH1D("Unweighted Bin Content","Unweighted Bin Content",11,0,10);
//TH1F *Mean = (TH1F*)result;
Int_t Nbins = result->GetNbinsX();
cout<<"subsampling Nbins"<<Nbins<<endl;
for(Int_t i=1; i<=Nbins; i++){
        Double_t w[100];
        Double_t wTotal = 0.0;
        Double_t gCombinedValue = 0.0;
        Double_t gCombinedError = 0.0;
        for(Int_t split = 0; split<10;split ++){
	        f = new TFile(Form(file+"_split_%d.root",split));
        	TList *l = (TList*) f->Get("FlowQCList;1");
        	h = (TH1F*) l->FindObject(obj);

                co = (TH1F*)h->Clone("co");
		Mean->Fill(split+1,co->GetBinContent(i));

                w[i] = 1./TMath::Power(co->GetBinError(i),2);
                wTotal +=w [i];

                gCombinedValue += co->GetBinContent(i)*w[i];
                gCombinedError += TMath::Power(co->GetBinError(i)*w[i],2);
		
		delete f;
	}
        gCombinedValue /= wTotal;
        gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
        if(abs(gCombinedValue)>0){
                result->SetBinContent(i, gCombinedValue);
                result->SetBinError(i, gCombinedError);
        }
}//end Nbins

TCanvas *M = new TCanvas("M","M",400,400);
M->SetLeftMargin(0.2);
Mean->GetYaxis()->SetTitleOffset(2.0);
Mean->GetYaxis()->SetTitle(Form("All values from differential flow v_{%d}",harm+1));
Mean->GetXaxis()->SetRangeUser(0,10);
Mean->GetXaxis()->SetTitle("Split sample");
Mean->SetLineColor(kRed);
Mean->DrawCopy();
M ->SaveAs(Form("v%d/Split_c%d_"+obj+".pdf",harm+1,c));
delete M;
delete Mean;
return result;
}//end Subsampling

TH1F *Subsample_Dv(TString file, TString obj1,TString obj2, TH1F *result, Int_t c, Int_t harm){
TFile *f;
TH1F *Dv = (TH1F*)result;
TH1D *MeanD = new TH1D("Unweighted Bin Difference","Unweighted Bin Difference",11,0,10);//(TH1F*)result;
TH1F *P, *N;
Int_t Nbins = result->GetNbinsX();

cout<<"subsampling Nbins"<<Nbins<<endl;
for(Int_t i=1; i<Nbins; i++){// <= to check if pT plot improves
	Double_t w[100];
	Double_t wTotal = 0.0;
	Double_t gCombinedValue = 0.0;
	Double_t gCombinedError = 0.0;

        Double_t vpos = 0.;
        Double_t vneg = 0.;
        Double_t deltav = 0.;
        Double_t vposError = 0.;
        Double_t vnegError = 0.;
        Double_t deltavError = 0.;
	for(Int_t split = 0; split<10;split ++){
        	f = new TFile(Form(file+"_split_%d.root",split));
        	TList *l = (TList*) f->Get("FlowQCList;1");
        	P = (TH1F*) l->FindObject(obj1);
        	N = (TH1F*) l->FindObject(obj2);
       
        	vpos = P->GetBinContent(i);
        	vneg = N->GetBinContent(i);
        	vposError = P->GetBinError(i);
        	vnegError = N->GetBinError(i);
        	deltav = vpos - vneg;
        	if(fabs(vpos)>0.&&fabs(vneg)>0.){
                         deltavError = sqrt(abs(pow(vposError,2)+pow(vnegError,2)-2*vnegError*vposError));
                         Dv->SetBinContent(i,deltav);
                         Dv->SetBinError(i, deltavError);
			
			MeanD->Fill(split+1,Dv->GetBinContent(i));
        	}
	w[i] = 1./TMath::Power(Dv->GetBinError(i),2);//maybe i+1
	wTotal +=w [i];

	gCombinedValue += Dv->GetBinContent(i)*w[i];
	gCombinedError += TMath::Power(Dv->GetBinError(i)*w[i],2);
	//MeanD->Fill(split+1,Dv->GetBinContent(i)*w[i]);

	//evtl fill for gaussian here ?
	delete f;
	}//end of split		
	gCombinedValue /= wTotal;
	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
	if(abs(gCombinedValue)>0){
		result->SetBinContent(i, gCombinedValue);
		result->SetBinError(i, gCombinedError);
	}	
}//end Nbins
TCanvas *M = new TCanvas("M","M",400,400);
M->SetLeftMargin(0.2);
MeanD->GetYaxis()->SetTitleOffset(2.0);
MeanD->GetYaxis()->SetTitle(Form("All values from differential flow v_{%d}",harm+1));
MeanD->GetXaxis()->SetRangeUser(0,10);
MeanD->GetXaxis()->SetTitle("split");
MeanD->SetLineColor(kRed);
MeanD->DrawCopy();
M ->SaveAs(Form("v%d/SplitDelta_c%d_"+obj1+".pdf",harm+1,c));
delete M;
delete MeanD;
return result;
}//end Subsample

void PlotpT(TH1F *PvpT, TH1F *NvpT, TH1F *DvpT, Int_t c, Int_t harm){
        //Plotting the subsampled pos/neg and their difference as in Multiplot.C for every harmonic==========
        TCanvas *c1= new TCanvas("c1", "canvas", 800, 800);//upper plot in pad1
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
	PvpT->SetTitle(" ");
        PvpT->DrawCopy();
        NvpT->SetLineColor(kBlue-1);
        NvpT->SetLineWidth(2);
	NvpT->SetTitle(" ");
        NvpT->DrawCopy("same");
        auto L1 = new TLegend();
        L1->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
        L1->AddEntry(PvpT,Form("AVFD, pos, Cent %d0-%d0",c,c+1), "l");
        L1->AddEntry(NvpT,Form("AVFD, neg, Cent %d0-%d0",c,c+1), "l");
        L1->Draw();
        c1->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.5);
        pad2->SetTopMargin(0);
        pad2->SetLeftMargin(0.2);
        pad2->SetRightMargin(4.);
        pad2->SetBottomMargin(0.2);
        pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();
	DvpT->SetTitle(" ");
	DvpT->SetStats(0);
        DvpT->GetXaxis()->SetRangeUser(0.,5.);
        DvpT->SetLineColor(kGreen-1);
        DvpT->SetLineWidth(2);
        DvpT->GetXaxis()->SetTitle("p_{T} [GeV]");
        DvpT->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
        DvpT->DrawCopy();
        c1->SaveAs(Form("v%d/pT_c%d.pdf",harm+1,c));
/*
        PvpT->Delete();
        NvpT->Delete();
        DvpT->Delete();
*/ 
       delete c1;
}

void PlotEta(TH1F *PvEta, TH1F *NvEta, TH1F *DvEta,Int_t c, Int_t harm){
        //Plotting the subsampled pos/neg and their difference as in Multiplot.C for every harmonic==========
        TCanvas *cEta = new TCanvas("ceta", "etacanvas", 800, 800);
        TPad *padeta1 = new TPad("padeta1", "padeta1", 0, 0.5, 1., 1.);
        padeta1->SetBottomMargin(0);
        padeta1->SetLeftMargin(0.2);
        padeta1->Draw();
        padeta1->cd();
        PvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
        PvEta->SetLineColor(kRed+1);
        PvEta->SetLineWidth(2);
        PvEta->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
        PvEta->SetStats(0);
	PvEta->SetTitle(" ");
        PvEta->DrawCopy();
        NvEta->SetLineColor(kBlue+1);
        NvEta->SetLineWidth(2);
	NvEta->SetTitle(" ");
        NvEta->DrawCopy("same");
        auto L2 = new TLegend();
        L2->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
        L2->AddEntry(PvEta,Form("AVFD, pos, Cent %d0-%d0",c,c+1), "l");
        L2->AddEntry(NvEta,Form("AVFD, neg, Cent %d0-%d0",c,c+1), "l");
        L2->Draw();
        cEta->cd();
        TPad *padeta2 = new TPad("padeta2", "padeta2", 0, 0., 1, 0.5);
        padeta2->SetTopMargin(0);
        padeta2->SetLeftMargin(0.2);
        padeta2->SetRightMargin(4.);
        padeta2->SetBottomMargin(0.2);
        padeta2->SetGridx();
        padeta2->Draw();
        padeta2->cd();
	DvEta->SetStats(0);
        DvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
        DvEta->GetXaxis()->SetTitle("#eta");
	DvEta->SetTitle(" ");
        DvEta->SetLineColor(kGreen+1);
        DvEta->SetLineWidth(2);
        DvEta->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
        DvEta->DrawCopy();
        cEta->SaveAs(Form("v%d/eta_c%d.pdf",harm+1,c));
/*
        PvEta->Delete();
        NvEta->Delete();
        DvEta->Delete();
*/
        delete cEta;

}

//from subsampled histos -> calculate diff pos - neg and propagate final error magically without knowledge of the degree of correlation ...
//took the saem equation for the \Delta vn Error as  in the CalculateFlowCME 

//=================================================================================================
//
//The Main from this macro
//
//
//Trouble notes:
//______________
//
// 1) Error Propagation on \Deltav1 - here and in CalculateFlowCME.cxx
// 2) Filling of hists at right place in the ClaculateFlowCME.cxx - maybe a loop earlier ? 
// 3) Automate things and make centrality evolution plots (via mean or multiple in 1 canvas) ?
// 4) Change limits/ranges and make things pretty
//=================================================================================================
void Deltavn(bool pT, Int_t cent, Int_t cmax, Int_t harm){
	//Bins from CalculateFlowCME.cxx, for the not yet initialized histograms (\Delta vn)=======
	Double_t fPtDiffNBins = 36;
	Double_t fCRCPtBins[36];//actually 37 before
	Double_t PtBins[37] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,17.,20.,25.,30.,40.,50.};
	for(Int_t r=0; r<37; r++){fCRCPtBins[r] = PtBins[r];}

	Double_t fEtaDiffNBins = 50;
	Double_t fCRCEtaBins[51];
	Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
	for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}

//=========================================================================================
	TH1F *cen;
	if(pT == true){cen = new TH1F("Multiple Centralities","Multiple Centralities",fPtDiffNBins, fCRCPtBins);}
	else if(pT == false){cen = new TH1F("Multiple Centralities", "Multiple Centralities", fEtaDiffNBins, fCRCEtaBins);}
	//generate empty hist with bins from above !
        TCanvas *Cen = new TCanvas("Cen","Cen",400,400);//one for pT's and one for eta's
        Cen->SetLeftMargin(0.2);
        cen->GetYaxis()->SetTitleOffset(2.0);
        cen->GetYaxis()->SetTitle(Form("differential flow #Delta v_{%d}",harm+1));//cen->GetXaxis()->SetRangeUser(0,5);
	if(pT == true){
	cen->GetXaxis()->SetRangeUser(0.1,5);
	cen->GetXaxis()->SetTitle("p_{T} [GeV]");
	cen->GetYaxis()->SetRangeUser(-0.1,0.1);
	}
	else{
	cen->GetXaxis()->SetRangeUser(-1,1);
	cen->GetXaxis()->SetTitle("#eta");
	}
	cen->SetStats(0);        
	cen->DrawCopy();
        auto C = new TLegend();
        C->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");

	for(Int_t c=cent; c<=cmax;c++){
	cout<<"cen "<<c<<endl;
	TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
	// !!! For harm = 0 we could also try to use the fFlowQCQv1[pos=0, neg=1][cen][pt = 0, eta = 1] !!!
	//for(Int_t harm=0; harm<3; harm++){
		cout<<"harm "<<harm<<endl;;
		//apply subsampling to all histograms from pos/neg particles=======================
		//mean from subsamples:
		if(pT == true){
		TH1F *PvpT = new TH1F("pvpT", "pvpT", fPtDiffNBins, fCRCPtBins);
		PvpT = Subsampling(input,Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm), PvpT, c, harm);
		TH1F *NvpT = new TH1F("nvpT","nvpT",fPtDiffNBins,fCRCPtBins);
		NvpT = Subsampling(input,Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm), NvpT, c ,harm);
		TH1F *DvpT = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
		DvpT = Subsample_Dv(input, Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm),Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm), DvpT, c, harm);
	 	//Plot the difference for multiple centralities====================================	
		DvpT->SetLineColor(c);
		C->AddEntry(DvpT,Form("AVFD, Cent %d0-%d0",c,c+1), "l");
		DvpT->DrawCopy("same");//draws hoefully in Cen canvas
		//Plot the pos/neg with difference================================================= 
		PlotpT(PvpT, NvpT, DvpT, c, harm);// plot will delete all hists!
		
        	PvpT->Delete();
        	NvpT->Delete();
        	//DvpT->Delete();
		}//end of pT if
		else{
		TH1F *PvEta = new TH1F("pvEta","pvEta",fEtaDiffNBins,fCRCEtaBins);
		PvEta = Subsampling(input,Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm), PvEta, c, harm);
		TH1F *NvEta = new TH1F("nvEta", "nvEta", fEtaDiffNBins, fCRCEtaBins);
		NvEta = Subsampling(input,Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm), NvEta, c, harm);
		TH1F *DvEta = new TH1F("deltavneta","deltavneta",fEtaDiffNBins,fCRCEtaBins);
		DvEta = Subsample_Dv(input, Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm), Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm), DvEta, c, harm);
		//Plot the difference for multiple centralities====================================
		DvEta->SetLineColor(c);
        	C->AddEntry(DvEta,Form("AVFD, Cent %d0-%d0",c,c+1), "l");
		DvEta->DrawCopy("same");
		//Plot the pos/neg with difference=================================================
		PlotEta(PvEta, NvEta, DvEta, c, harm);
		
        	PvEta->Delete();
        	NvEta->Delete();
	        DvEta->Delete();
		}//end of eta if
		/*
      		else if(subsample_deltavn == false){
		Double_t vpos = 0.;
        	Double_t vneg = 0.;
      		Double_t deltav = 0.;
        	Double_t vposError = 0.;
        	Double_t vnegError = 0.;
       	 	Double_t deltavError = 0.;
		for(Int_t pt=1; pt<=fPtDiffNBins; pt++){
                    vpos = PvpT->GetBinContent(pt);
                    vneg = NvpT->GetBinContent(pt);
                    vposError = PvpT->GetBinError(pt);
                    vnegError = NvpT->GetBinError(pt);
                    deltav = vpos - vneg;
                   	if(fabs(vpos)>0.&&fabs(vneg)>0.){
				deltavError = sqrt(abs(pow(vposError,2)+pow(vnegError,2)-2*vnegError*vposError));
                		DvpT->SetBinContent(pt,deltav);
                		DvpT->SetBinError(pt, deltavError);   
			}
		}//end of pt bin
		for(Int_t eta=1; eta<=fEtaDiffNBins; eta++){
                    vpos = PvEta->GetBinContent(eta);
                    vneg = NvEta->GetBinContent(eta);
                    vposError = PvEta->GetBinError(eta);
                    vnegError = NvEta->GetBinError(eta);
                    deltav = vpos - vneg;
                        if(fabs(vpos)>0.&&fabs(vneg)>0.){
				deltavError = sqrt(abs(pow(vposError,2)+pow(vnegError,2)-2*vnegError*vposError));
                                DvEta->SetBinContent(eta,deltav);
                                DvEta->SetBinError(eta, deltavError);
			}
        	}//end of eta bin
		}//end if false
		*/
        //===================================================================================================

}//centrality loop
C->Draw();
Cen->SaveAs(Form("v%d/Multi_cen%d_%d.pdf",harm+1,cent,cmax));
//delete C;
delete Cen;
delete cen;

//}//end of hr<fkFlowNHarm
//}//end of centrality loop

}//end void


