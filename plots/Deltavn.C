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
//let’s try to confirm that the spread we see on the left plot is correct then as a cross-check; 
//let’s pick three pT bins and make an 1D histogram for each of these bins that will be filled with the v1 values of all samples; 
//then compare the std deviation of the 1D histo with the spread of the TProfile for this pT bin
//
//======================================================================================
void CrossCheck(TProfile *Prof, TH1D *hist, Int_t bin, Double_t binValue){
//calculate spread for Profile and std from TH1D
//for(Int_t bin=1; bin<36; bin++){
	Double_t stats[6] = {0.};
	Prof->GetXaxis()->SetRange(bin,bin);
	Prof->GetStats(stats);
	Double_t SumWeig   = stats[0];
        Double_t SumWeigSq  = stats[1];
        Double_t SumTwo  = stats[4];
        Double_t SumTwoSq = stats[5];
//wtf one is 0 the other has a value in hist->(bin ?!)
	cout<<"bin Nr:"<<bin<<" binValue "<<binValue<<" contentProfile "<<Prof->GetBinContent(bin)<<" contentHist "<< hist->GetBinContent(bin)<<endl;
	Double_t spread = pow( fabs(pow(SumTwo/SumWeig,2)- SumTwoSq/SumWeigSq) ,0.5);
	Double_t histErr = hist->GetRMS();
	Double_t histErr1 = hist->GetStdDev();
	Double_t histErr2 = hist->GetBinError(bin);
	cout<<" TProfile spread "<<spread<<"; hist->GetBinError "<<histErr2<<"; (hist->GetRMS)"<<histErr<<"; (hist->GetStdDev) "<<histErr1<<endl;
//}
}//end void

TH1F *Subsampling(TString file, TString obj, Double_t binArr[], TH1F *result, Int_t c, Int_t harm){
Int_t binMin=15;
Double_t binMinArr = binArr[binMin];
TFile *f;
TH1F *co, *h;
Int_t Nbins = result->GetNbinsX();
TH1D *hist = new TH1D("hist","hist",Nbins,binArr);
TProfile *Mint = new TProfile("Weighted Bin Content","Weighted Bin Content",11,0,10,"s");
TProfile *Mdiff = new TProfile("Unweighted Bin Content","Unweighted Bin Content",Nbins,binArr);
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

		Mdiff->Fill(binArr[i-1],co->GetBinContent(i));//differential value for all subsamples without weight
		Mint->Fill(split+1,co->GetBinContent(i),Nbins);//integrated value over all bins, weighted by #bins
		if(i==binMin){hist->Fill(binArr[i-1],co->GetBinContent(i));}

                w[i] = 1./TMath::Power(co->GetBinError(i),2);
                wTotal +=w [i];
                gCombinedValue += co->GetBinContent(i)*w[i];
                gCombinedError += TMath::Power(co->GetBinError(i)*w[i],2);
		
		delete f;
	}
        gCombinedValue /= wTotal;
        gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
        if(fabs(gCombinedValue)>0){
                result->SetBinContent(i, gCombinedValue);
                result->SetBinError(i, gCombinedError);
        }
}//end Nbins
CrossCheck(Mdiff,hist,binMin,binMinArr);//Profile, hist, bin range is integer -> relevant for spread calculation
TCanvas *M = new TCanvas("M","M",400,400);
M->SetLeftMargin(0.2);
Mdiff->GetYaxis()->SetTitleOffset(2.0);
Mdiff->GetYaxis()->SetTitle(Form("differential flow v_{%d} from all splitfiles",harm+1));
if(Nbins<40){//pT bins
Mdiff->GetXaxis()->SetTitle("p_{T} [GeV]");
Mdiff->GetXaxis()->SetRange(1,22);//range in bins !
}else if(Nbins>40){
Mdiff->GetXaxis()->SetTitle("#eta");
Mdiff->GetXaxis()->SetRangeUser(-0.9,0.9);
}
Mdiff->SetLineColor(c+harm);
Mdiff->DrawCopy();
M ->SaveAs(Form("v%d/Mdiff_c%d_"+obj+".pdf",harm+1,c));
delete M;
delete Mdiff;

TCanvas *Mi = new TCanvas("Mi","Mi",400,400);
Mi->SetLeftMargin(0.2);
Mint->GetYaxis()->SetTitleOffset(2.0);
Mint->GetYaxis()->SetTitle(Form("differential flow v_{%d} from all bins",harm+1));
Mint->GetXaxis()->SetTitle("Split file");
Mint->GetXaxis()->SetRangeUser(0,10);
Mint->SetLineColor(c+harm);
Mint->DrawCopy();
Mi ->SaveAs(Form("v%d/Mint_c%d_"+obj+".pdf",harm+1,c));
delete Mi;
delete Mint;
return result;
}//end Subsampling

TH1F *Subsample_Dv(TString file, TString obj1,TString obj2, Double_t binArr[], TH1F *result, Int_t c, Int_t harm){
TFile *f;
TH1F *Dv = (TH1F*)result;
Int_t Nbins = result->GetNbinsX();
TH1D *MeanD = new TH1D("Unweighted Bin Difference","Unweighted Bin Difference",11,0,10);
TProfile *Dint = new TProfile("Weighted Bin Content","Weighted Bin Content",11,0,10,"s");
TProfile *Diff = new TProfile("Unweighted Bin Content","Unweighted Bin Content",Nbins,binArr);
TH1F *P, *N;

cout<<"subsampling Nbins"<<Nbins<<endl;
for(Int_t i=1; i<=Nbins; i++){
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
	
	                 Diff->Fill(binArr[i],Dv->GetBinContent(i));//differential value for all subsamples without weight
        	         Dint->Fill(split+1,Dv->GetBinContent(i),Nbins);//integrated value over all bins, weighted by #bins
        	}
	w[i] = 1./TMath::Power(Dv->GetBinError(i),2);//maybe i+1
	wTotal +=w [i];
	gCombinedValue += Dv->GetBinContent(i)*w[i];
	gCombinedError += TMath::Power(Dv->GetBinError(i)*w[i],2);

	delete f;
	}//end of split		
	gCombinedValue /= wTotal;
	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
	if(fabs(gCombinedValue)>0){
		result->SetBinContent(i, gCombinedValue);
		result->SetBinError(i, gCombinedError);
	}	
}//end Nbins
TCanvas *M = new TCanvas("M","M",400,400);
M->SetLeftMargin(0.2);
Diff->GetYaxis()->SetTitleOffset(2.0);
Diff->GetYaxis()->SetTitle(Form("#Delta v_{%d} from all splitfiles",harm+1));
if(Nbins<40){//pT bins
Diff->GetXaxis()->SetTitle("p_{T} [GeV]");
Diff->GetXaxis()->SetRange(1,22);//range in bins !
}else if(Nbins>40){
Diff->GetXaxis()->SetTitle("#eta");
Diff->GetXaxis()->SetRangeUser(-0.9,0.9);
}
Diff->SetLineColor(c+harm);
Diff->DrawCopy();
M ->SaveAs(Form("v%d/Diff_c%d_"+obj1+".pdf",harm+1,c));
delete M;
delete Diff;

TCanvas *Mi = new TCanvas("Mi","Mi",400,400);
Mi->SetLeftMargin(0.2);
Dint->GetYaxis()->SetTitleOffset(2.0);
Dint->GetYaxis()->SetTitle(Form("#Delta v_{%d} from all bins",harm+1));
Dint->GetXaxis()->SetTitle("Split file");
Dint->GetXaxis()->SetRangeUser(0,10);
Dint->SetLineColor(c+harm);
Dint->DrawCopy();
Mi ->SaveAs(Form("v%d/Dint_c%d_"+obj1+".pdf",harm+1,c));
delete Mi;
delete Dint;
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
        PvpT->GetXaxis()->SetRange(1,22);
        PvpT->SetLineColor(c+1);
        PvpT->SetLineWidth(2);
        PvpT->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
        PvpT->SetStats(0);
	PvpT->SetTitle(" ");
        PvpT->DrawCopy();
	NvpT->GetXaxis()->SetRange(1,22);//range in bins
        NvpT->SetLineColor(c+2);
        NvpT->SetLineWidth(2);
	NvpT->SetTitle(" ");
        NvpT->DrawCopy("same");
        auto L1 = new TLegend();
        L1->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
        L1->AddEntry(PvpT,Form("+h, Centrality %d0-%d0",c,c+1), "l");
        L1->AddEntry(NvpT,Form("-h, Centrality %d0-%d0",c,c+1), "l");
        L1->Draw();
        c1->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.5);
        pad2->SetTopMargin(0);
        pad2->SetLeftMargin(0.2);
        //pad2->SetRightMargin(4.);
        pad2->SetBottomMargin(0.2);
        pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();
	DvpT->SetTitle(" ");
	DvpT->SetStats(0);
        DvpT->GetXaxis()->SetRange(1,22);
        DvpT->SetLineColor(c);
        DvpT->SetLineWidth(2);
        DvpT->GetXaxis()->SetTitle("p_{T} [GeV]");
        DvpT->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
        DvpT->DrawCopy();
        c1->SaveAs(Form("v%d/pT_c%d.pdf",harm+1,c));
 
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
	PvEta->GetXaxis()->SetNdivisions(4);
        PvEta->SetLineColor(c+1);
        PvEta->SetLineWidth(2);
        PvEta->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
        PvEta->SetStats(0);
	PvEta->SetTitle(" ");
        PvEta->DrawCopy();
	NvEta->GetXaxis()->SetNdivisions(4);
        NvEta->SetLineColor(c+2);
        NvEta->SetLineWidth(2);
	NvEta->SetTitle(" ");
        NvEta->DrawCopy("same");
        auto L2 = new TLegend();
        L2->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");
        L2->AddEntry(PvEta,Form("+h, Centrality %d0-%d0",c,c+1), "l");
        L2->AddEntry(NvEta,Form("-h, Centrality %d0-%d0",c,c+1), "l");
        L2->Draw();
        cEta->cd();
        TPad *padeta2 = new TPad("padeta2", "padeta2", 0, 0., 1, 0.5);
        padeta2->SetTopMargin(0);
        padeta2->SetLeftMargin(0.2);
        //padeta2->SetRightMargin(4.);
        padeta2->SetBottomMargin(0.2);
        padeta2->SetGridx();
        padeta2->Draw();
        padeta2->cd();
	DvEta->SetStats(0);
        DvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
        DvEta->GetXaxis()->SetTitle("#eta");
	DvEta->SetTitle(" ");
        DvEta->SetLineColor(c);
        DvEta->SetLineWidth(2);
        DvEta->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
        DvEta->DrawCopy();
        cEta->SaveAs(Form("v%d/eta_c%d.pdf",harm+1,c));

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
// 3) Change limits/ranges and make things pretty
//=================================================================================================
void Deltavn(bool pT, Int_t cent, Int_t cmax, Int_t harm){
	gROOT->SetBatch();//to avoid opening the plots ak bad wifi struggle
	//Bins from CalculateFlowCME.cxx, for the not yet initialized histograms (\Delta vn)=======
	Double_t fPtDiffNBins = 36;
	Double_t fCRCPtBins[37]={0};
	Double_t PtBins[37] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,17.,20.,25.,30.,40.,50.};
	for(Int_t r=0; r<37; r++){fCRCPtBins[r] = PtBins[r];}

	Double_t fEtaDiffNBins = 50;
	Double_t fCRCEtaBins[51]={0};
	Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
	for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}

//================================================================================================
	//generate empty canvas for multiple centrality plots
        TCanvas *Cen = new TCanvas("Cen","Cen",400,400);//one for pT's and one for eta's
        Cen->SetLeftMargin(0.2);
	auto C = new TLegend();
        C->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #eta: 0.8","C");

	for(Int_t c=cent; c<=cmax;c++){
	cout<<"centrality "<<c<<endl;
	TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
	// !!! For harm = 0 we could also try to use the fFlowQCQv1[pos=0, neg=1][cen][pt = 0, eta = 1] !!!
	//for(Int_t harm=0; harm<3; harm++){
		cout<<"harmonic - 1: "<<harm<<endl;;
		//apply subsampling to all histograms from pos/neg particles=======================
		//mean from subsamples:
		if(pT == true){
		TH1F *PvpT = new TH1F("pvpT", "pvpT", fPtDiffNBins, fCRCPtBins);
		PvpT = Subsampling(input,Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm),fCRCPtBins, PvpT, c, harm);
		TH1F *NvpT = new TH1F("nvpT","nvpT",fPtDiffNBins,fCRCPtBins);
		NvpT = Subsampling(input,Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm),fCRCPtBins, NvpT, c ,harm);
		TH1F *DvpT = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
		DvpT = Subsample_Dv(input, Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm),Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm), fCRCPtBins, DvpT, c, harm);
	 	//Plot the difference for multiple centralities====================================	
		if(c == cent){//draw vn difference for multiple centralities       
			DvpT->GetYaxis()->SetTitleOffset(2.0);
        		DvpT->GetYaxis()->SetTitle(Form("differential flow #Delta v_{%d}",harm+1));//cen->GetXaxis()->SetRangeUser(0,5);
        		DvpT->SetStats(0);
			DvpT->SetTitle(" ");
			DvpT->GetXaxis()->SetRange(1,22);
        		DvpT->GetXaxis()->SetTitle("p_{T} [GeV]");
        		//DvpT->SetAxisRange(0,5,"X");
			DvpT->SetLineColor(c);
			C->AddEntry(DvpT,Form("AVFD, Centrality %d0-%d0",c,c+1),"l");
			DvpT->DrawCopy();
		   }else if(c>cent){
			DvpT->GetXaxis()->SetRange(1,22);
			//DvpT->SetAxisRange(0,5,"X");
			DvpT->SetLineColor(c);
			C->AddEntry(DvpT,Form("AVFD, Centrylity %d0-%d0",c,c+1), "l");
			DvpT->DrawCopy("same");
		}
		//Plot the pos/neg with difference================================================= 
		PlotpT(PvpT, NvpT, DvpT, c, harm);// plot will delete all hists!
		
        	PvpT->Delete();
        	NvpT->Delete();
        	//DvpT->Delete();
		}//end of pT if
		else if(pT == false){
		TH1F *PvEta = new TH1F("pvEta","pvEta",fEtaDiffNBins,fCRCEtaBins);
		PvEta = Subsampling(input,Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm),fCRCEtaBins, PvEta, c, harm);
		TH1F *NvEta = new TH1F("nvEta", "nvEta", fEtaDiffNBins, fCRCEtaBins);
		NvEta = Subsampling(input,Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm),fCRCEtaBins, NvEta, c, harm);
		TH1F *DvEta = new TH1F("deltavneta","deltavneta",fEtaDiffNBins,fCRCEtaBins);
		DvEta = Subsample_Dv(input, Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm), Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm), fCRCEtaBins, DvEta, c, harm);
		//Plot the difference for multiple centralities====================================
		if(c==cent){ 
			DvEta->GetYaxis()->SetTitleOffset(2.0);
			DvEta->GetYaxis()->SetTitle(Form("differential flow #Delta v_{%d}",harm+1));
        		DvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
        		DvEta->GetXaxis()->SetTitle("#eta");
        		DvEta->SetStats(0);
			DvEta->SetTitle(" ");
			DvEta->SetLineColor(c);
			C->AddEntry(DvEta,Form("AVFD, Centrality %d0-%d0",c,c+1),"l");
	        	DvEta->DrawCopy();
		   }else if(c>cent){
                        DvEta->SetLineColor(c);
			DvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
                        C->AddEntry(DvEta,Form("AVFD, Centrylity %d0-%d0",c,c+1), "l");
                        DvEta->DrawCopy("same");
		}

		//Plot the pos/neg with difference=================================================
		PlotEta(PvEta, NvEta, DvEta, c, harm);
		
        	PvEta->Delete();
        	NvEta->Delete();
	        //DvEta->Delete();
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
if(pT==true){Cen->SaveAs(Form("v%d/Pt_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
else if(pT == false){Cen->SaveAs(Form("v%d/Eta_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
//delete C;
delete Cen;

//}//end of hr<fkFlowNHarm

}//end void


