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
void findSpread(Int_t c, Int_t harm,TString spectrum){
  TFile *f[10];
  TList *list[10];
  TH1D *fHist[10];
  //Int_t objID = 0;
  Int_t Nbins = 0;
  if(spectrum == "Pt"){Nbins = 30;}
  else if(spectrum == "Eta"){Nbins = 50;}
  TH1D *fHistPos[10];
  TH1D *fHistNeg[10];
  TH1F *fHistSpreadPos[Nbins];
  TH1F *fHistSpreadNeg[Nbins];
  TH1F *fHistSpreadDiff[Nbins];
/*
  if(obj == Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm)){
	objID = 1;//pos pT
  }else if(obj == Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm)){
  	objID = 2;//neg pT
  }else if(obj == Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",c,harm)){
	objID = 3;//del pT
  }else if(obj == Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm)){
	objID = 4;//pos eta
  }else if(obj == Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm)){
        objID = 5;//neg eta
  }else if(obj == Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",c,harm)){
        objID = 6;//del eta
  }

  TH1F *fHistSpread[Nbins];

  for(Int_t i = 0; i < Nbins; i++){
	TString histo_string;
                if(objID== 0){
                        cout<<"No obj to findSpread"<<endl;
                }else if(objID==1){
                        histo_string = Form("fHistSpreadPosPtBin%d",i+1);
		}else if(objID == 2){
                        histo_string = Form("fHistSpreadNegPtBin%d",i+1);
                }else if(objID==3){
                        histo_string = Form("fHistSpreadDPtBin%d",i+1);
                }else if(objID==4){
                        histo_string = Form("fHistSpreadPosEtaBin%d",i+1);
                }else if(objID==5){
                        histo_string = Form("fHistSpreadNegEtaBin%d",i+1);
                }else if(objID == 6){
                        histo_string = Form("fHistSpreadDEtaBin%d",i+1);
        }
	fHistSpread[i] = new TH1F(histo_string,"",1000,-0.3,0.3);
   }
*/
 for(Int_t i = 0; i< Nbins; i++){
	fHistSpreadPos[i] = new TH1F(Form("fHistSpreadPos"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
	fHistSpreadNeg[i] = new TH1F(Form("fHistSpreadNeg"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
 	fHistSpreadDiff[i] = new TH1F(Form("fHistSpreadDiff"+spectrum+"Bin%d",i+1),"",1000,-0.3,0.3);
 }
  for(Int_t iFile = 0; iFile < 10; iFile++) {
    //f[iFile] = TFile::Open(Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    f[iFile] = TFile::Open(Form("/data/alice/jlomker/AVFD/result/dirID-0/NoBField/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    if((!f[iFile])||(!f[iFile]->IsOpen())) {
    cout<<"File "<<iFile<<" not found..."<<endl;
    return;
    }

    list[iFile] = dynamic_cast<TList *>(f[iFile]->Get("FlowQCList"));
    if(!list[iFile]) {
      cout<<"Input list of file "<<iFile<<" not found..."<<endl;
      return;
    }
    fHistPos[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[0][%d][%d][0]",c,harm)));
    if(!fHistPos[iFile]) {
      cout<<"Histogram of positive particles from file "<<iFile<<" not found..."<<endl;
      return;
    }
    fHistNeg[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[1][%d][%d][0]",c,harm)));
    if(!fHistNeg[iFile]) {
      cout<<"Histogram of negative particles from file "<<iFile<<" not found..."<<endl;
      return;
    }

/*
    if(Nbins<40){
    fHist[iFile]->GetXaxis()->SetRangeUser(0.2,5);}//if
    if(Nbins>40){
    fHist[iFile]->GetXaxis()->SetRangeUser(-0.9,0.9);}//if
    fHist[iFile]->SetMarkerColor(iFile+1);
    fHist[iFile]->SetLineColor(iFile+1);
    fHist[iFile]->Draw("ESAME");
*/
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
}//end find spread

double CrossCheck(TProfile *Prof, TH1D *hist, Int_t bin){//calculate spread for Profile and std from TH1D
	Double_t N = Prof->GetNbinsX();
	Double_t RMS = hist->GetRMS();
	//Double_t RMS = Prof->GetBinError(bin);
	//for(Int_t i = 1; i<=hist->GetNbinsX(); i++){cout<<"Content 1D hist "<<hist->GetBinContent(i)<<endl;}
	//cout<<" TProfle->GetBinError(bin_Nr) "<<Prof->GetBinError(bin)<<" hist->GetRMS() "<<hist->GetRMS()<<endl;
	return RMS;
}//end crosscheck
//this is where I have to continue after food
TH1F *Subsampling(TString file, TString spectrum, TString charge, TH1F *result, Int_t c, Int_t harm){
	TFile *f, *t;
	TFile *spread, *output;
	TList *l;
	TH1F *h, *o;

	Int_t objID = 0;
	Int_t Nbins = result->GetNbinsX();
	Int_t MinBin, MaxBin;
	//Double_t RMS[51] = {0.};
	//TProfile *Mint = (TProfile*)mint->Clone("Mint");
	//TProfile *Mdiff = (TProfile*)mdiff->Clone("Mdiff");

	//TH1F *rms = (TH1F*)RMs->Clone("rms");
	//TH1D *cc = (TH1D*)CC->Clone("cc");
	//TH1D *hist = (TH1D*)Hist->Clone("hist");
	//setting for pT/eta
	/*
	if(obj == Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm) ){
		cout<<"Subsampling differential pT from positive particles"<<endl;
		objID = 1;
		MinBin = 1;
		MaxBin = Nbins;
	}else if(obj == Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm)){
                cout<<"Subsampling differential pT from negative particles"<<endl;
                objID = 2;
                MinBin = 1;
                MaxBin = Nbins;
	}else if(obj == Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",c,harm)){
                cout<<"Subsampling differential delta pT"<<endl;
                objID = 3;
                MinBin = 1;
                MaxBin = Nbins;
	}else if(obj == Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm)){
		cout<<"Subsampling differential eta for positive particles"<<endl;
		objID = 4;
		MinBin = 18;
		MaxBin = 34;
	}else if(obj == Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm)){
                cout<<"Subsampling differential eta for negative particles"<<endl;
                objID = 5;
                MinBin = 18;
                MaxBin = 34;
	}else if(obj == Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",c,harm)){
                cout<<"Subsampling differential delta eta"<<endl;
                objID = 6;
                MinBin = 18;
                MaxBin = 34;
	}
	*/
	for(Int_t i=1; i<=Nbins; i++){
		Double_t sigma = 0.;
		Double_t mean = 0.;
		/*Double_t w[100];
	        Double_t wTotal = 0.0;
        	Double_t gCombinedValue = 0.0;
	        Double_t gCombinedError = 0.0;
		*/
		TString spreadstring;
		output = spread->Open(Form("output_findSpread/spread_"+spectrum+"_c%d0_%d0.root",c,c+1));
			/*if(objID==0){
				cout<<"No obj to subsample"<<endl;
			}else if(objID==1){
				spreadstring = Form("fHistSpreadPosPtBin%d",i);
			}else if(objID ==2){
				spreadstring = Form("fHistSpreadNegPtBin%d",i);
			}else if(objID ==3){
				spreadstring = Form("fHistSpreadDPtBin%d",i);
			}else if(objID ==4){
				spreadstring = Form("fHistSpreadPosEtaBin%d",i);
			}else if(objID ==5){
				spreadstring = Form("fHistSpreadNegEtaBin%d",i);
			}else if(objID==6){
				spreadstring = Form("fHistSpreadDEtaBin%d",i);
			}*/
		spreadstring = Form("fHistSpread"+charge+spectrum+"Bin%d",i);	
		//cout<<spreadstring<<"string"<<endl;
		o = (TH1F*) output->Get(spreadstring);
		sigma = o->GetRMS();
		mean = o->GetMean();//take the file check the values fro 3*sigma
		result->SetBinContent(i, mean);
		result->SetBinError(i, sigma);
		o->Reset();
		/*
	        for(Int_t split = 0; split<10;split ++){
		        t = f->Open(Form(file+"_split_%d.root",split));//close files again !
        		l = (TList*) t->Get("FlowQCList;1");
	        	h = (TH1F*) l->FindObject(obj);
		
			Double_t val = h->GetBinContent(i);
			Double_t center = h->GetBinCenter(i);
			//cout<<" rms("<<i<<") from split histogram "<<h->GetBinError(i)<<" from find spread: 3* sigma "<<3*sigma<<" mean "<<mean<<endl;
			if(fabs(h->GetBinError(i))>fabs(3*sigma)){cout << "left out flow value: "<< val<<" and error "<<h->GetBinError(i)<<" with bin Nr. "<<i<<" from splitfile: "<<split<<endl;}
			if(fabs(h->GetBinError(i))<=fabs(3*sigma)){
				//Mdiff->Fill(center,val);//differential value for all subsamples without weight
				//Mint->Fill(split+1,h->GetBinContent(i),Nbins);//integrated value over all bins, weighted by #bins		
				//hist->Fill(split,h->GetBinContent(i));//For RMS vs split evolution//ggf without +1
				//cc->Fill(val);//For CrossCheck with Mdiff
				//Int_t len = Mdiff->GetNbinsX();
				//Mdiff->GetXaxis()->SetRange(i,i);
				//cout<<"split "<<split<<",  error 1D cc: "<<cc->GetRMS()<<", content 1D hist: "<<hist->GetRMS()<<", error TProfile(bin): "<<Mdiff->GetBinError(i)<<", center: "<<Mdiff->GetBinCenter(i)<<endl;
				//Rms[split] = CrossCheck(Mdiff,cc,i);
				//Mdiff->GetXaxis()->SetRange(0,len);
               			w[i] = 1./TMath::Power(h->GetBinError(i),2);
                		wTotal +=w [i];
	                	gCombinedValue += h->GetBinContent(i)*w[i];
        	        	gCombinedError += TMath::Power(h->GetBinError(i)*w[i],2);
			}//end if (...)
			//cc->Reset();
			h->Reset();
			f->Close();
			t->Close();
		}//end split
		o->Reset();//may fix memory problem
        	gCombinedValue /= wTotal;
        	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
		t = f->Open(Form("/data/alice/jlomker/AVFD/result/dirID-0/new/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",c,c+1));
	        
                l = (TList*) t->Get("FlowQCList;1");
                h = (TH1F*) l->FindObject(obj);
		if(fabs(gCombinedError)>0&&fabs(h->GetBinContent(i)>0)){
        	        result->SetBinContent(i, h->GetBinContent(i));
                	result->SetBinError(i, gCombinedError);
        	}
		h->Reset();
		f->Close();
		t->Close();
		*/
	}//end Nbins
return result;
}//end subsampling
/*
	cout<<"after loop"<<endl;
	RMS[i] = CrossCheck(Mdiff,cc,i);
	
	if(RMS[i]<RMS[i-1]){//Checking the drops in RMS pT evolution
	TCanvas *D = new TCanvas("Drops","Drops",400,400);
	D->SetLeftMargin(0.2);
	hist->GetYaxis()->SetTitleOffset(2.0);
	if(Nbins<40){hist->GetYaxis()->SetTitle(Form("Differential flow v_{%d}(p_{T} = %f)",harm+1,binArr[i]));}
	else if(Nbins>40){hist->GetYaxis()->SetTitle(Form("Differential flow v_{%d}(#eta = %f)",harm+1,binArr[i]));}
	hist->GetXaxis()->SetTitle("Splitfile");
	hist->SetLineColor(c+harm);
	hist->DrawCopy();
	D->SaveAs(Form("v%d/drop_bin%d_c%d_"+obj+".pdf",harm+1,i,c));
	delete D;
	}//end of single bin subsample plots for strange pT drops
	else if(RMS[i]>=RMS[i-1]){
        TCanvas *D = new TCanvas("Drops","Drops",400,400);
        D->SetLeftMargin(0.2);
        hist->GetYaxis()->SetTitleOffset(2.0);
        if(Nbins<40){hist->GetYaxis()->SetTitle(Form("Differential flow v_{%d}(p_{T} = %f)",harm+1,binArr[i]));}
        else if(Nbins>40){hist->GetYaxis()->SetTitle(Form("Differential flow v_{%d}(#eta = %f)",harm+1,binArr[i]));}
        hist->GetXaxis()->SetTitle("Splitfile");
        hist->SetLineColor(c+harm);
        hist->DrawCopy();
        D->SaveAs(Form("v%d/rise_bin%d_c%d_"+obj+".pdf",harm+1,i,c));
        delete D;
	}//end of else if -> rising RMS pT evolution
	//delete hist;
	rms->SetBinContent(i,RMS[i]);//RMS vs bin evolution
//}//end Nbins
delete hist;

TCanvas *R = new TCanvas("R","R",400,400);
R->SetLeftMargin(0.2);
rms->GetYaxis()->SetTitleOffset(2.0);
rms->GetYaxis()->SetTitle(Form("RMS of 1D Hist from all splitfiles v_{%d}",harm+1));
if(Nbins<40){//pT bins
rms->GetXaxis()->SetTitle("p_{T} [GeV]");
rms->GetXaxis()->SetRange(1,22);//range in bins !
}else if(Nbins>40){
rms->GetXaxis()->SetTitle("#eta");
rms->GetXaxis()->SetRangeUser(-0.9,0.9);
}
rms->SetLineColor(c+harm);
rms->DrawCopy();
R->SaveAs(Form("v%d/RMS_c%d_"+obj+".pdf",harm+1,c));
delete R;
delete rms;

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
return result;*/
//}//end Subsampling

TH1F *Dv(Int_t c, Int_t harm, TH1F *P, TH1F *N, TString obj, TH1F *result){//to get delta vn from subsampled vn(+/-q)
	Int_t Nbins = result->GetNbinsX();
	Int_t objID = 0;
	Double_t vpos, vneg, vposError, vnegError, deltav, deltavError;
	Int_t MinBin, MaxBin;
	TH1F *o;
	TFile *output, *spread;
 	if(obj == Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",c,harm)){
                cout<<"Calculating differential delta pT"<<endl;
                objID = 1;
                MinBin = 1;
                MaxBin = Nbins;
        }else if(obj == Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",c,harm)){
                cout<<"Calculating differential delta eta"<<endl;
                objID = 2;
                MinBin = 18;
                MaxBin = 34;
        }
        for(Int_t i=MinBin; i<=MaxBin; i++){
                TString spreadstring;
                output = spread->Open(Form("output_findSpread/spread_"+obj+"_c%d0_%d0.root",c,c+1));
                        if(objID==0){
                                cout<<"No obj to subsample"<<endl;
                        }else if(objID ==1){
                                spreadstring = Form("fHistSpreadDPtBin%d",i);
                        }else if(objID==2){
                                spreadstring = Form("fHistSpreadDEtaBin%d",i);
                        }
                o = (TH1F*) output->Get(spreadstring);
                deltavError = o->GetRMS();// average error on spread from the delta vn analysis results
		vpos = P->GetBinContent(i);// contents from subsampled results, weighted by combined (over all splits) inverse errors suqared for each bin
        	vneg = N->GetBinContent(i);
	        //vposError = P->GetBinError(i);
        	//vnegError = N->GetBinError(i);
	        deltav = fabs(vpos) - fabs(vneg);
        	if(fabs(vpos)>0.&&fabs(vneg)>0.){//switch to completely correlated error
        		//deltavError = sqrt(abs(pow(vposError,2)-pow(vnegError,2)));
               		result->SetBinContent(i,deltav);
                	result->SetBinError(i, deltavError);
        	}//end if
		o->Reset();
  	}//end for	
  return result;
}//end Dv

TH1F *Subsample_Dv(TString file, TString obj1,TString obj2, Double_t binArr[], TH1F *result, Int_t c, Int_t harm,TProfile *mint, TProfile *mdiff, TH1F *RMs, TH1D *CC, TH1D *Hist){
TFile *f, *t;
TH1F *Dv = (TH1F*)result;
Int_t Nbins = result->GetNbinsX();

//Double_t RMS[51] = {0.};
TProfile *Dint = (TProfile*)mint->Clone("Mint");
TProfile *Diff = (TProfile*)mdiff->Clone("Mdiff");

//TH1D *MeanD = new TH1D("Unweighted Bin Difference","Unweighted Bin Difference",11,0,10);
//TProfile *Dint = new TProfile("Weighted Bin Content","Weighted Bin Content",11,0,10,"s");
//TProfile *Diff = new TProfile("Unweighted Bin Content","Unweighted Bin Content",Nbins,binArr);
TH1F *P, *N;
if(Nbins<50){cout<<"Subsampling delta pT"<<endl;}
else if(Nbins>50){cout<<"Subsampling #eta"<<endl;}

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
        	t = f->Open(Form(file+"_split_%d.root",split));
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

	f->Close();
	t->Close();
	}//end of split		
	gCombinedValue /= wTotal;
	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
	if(fabs(gCombinedValue)>0){
		result->SetBinContent(i, gCombinedValue);
		result->SetBinError(i, gCombinedError);
	}	
}//end Nbins
return result;

}//end subsampling
/*
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
*/
//}//end Subsample

void PlotpT(TH1F *PvpT, TH1F *NvpT, TH1F *DvpT, Int_t c, Int_t harm){
        //Plotting the subsampled pos/neg and their difference as in Multiplot.C for every harmonic==========
        TCanvas *c1= new TCanvas("c1", "canvas", 800, 800);//upper plot in pad1
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1., 1.0);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetLeftMargin(0.2);
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad//and filled with pos
        PvpT->GetXaxis()->SetRangeUser(0,5);
        PvpT->SetLineColor(c+1);
        PvpT->SetLineWidth(1);
        PvpT->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
        PvpT->SetStats(0);
	PvpT->SetMarkerStyle(22);
	PvpT->SetMarkerColor(c+1);
	PvpT->SetTitle(" ");
        PvpT->DrawCopy();
	NvpT->GetXaxis()->SetRangeUser(0.,5);//range in bins
        NvpT->SetLineColor(c+2);
        NvpT->SetLineWidth(1);
	NvpT->SetMarkerStyle(23);
	NvpT->SetMarkerColor(c+2);
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
	DvpT->SetMarkerStyle(20);
	DvpT->SetMarkerColor(c);
        DvpT->GetXaxis()->SetRangeUser(0,5);
        DvpT->SetLineColor(c);
        DvpT->SetLineWidth(1);
        DvpT->GetXaxis()->SetTitle("p_{T} [GeV]");
        DvpT->GetYaxis()->SetTitle(Form("#Delta v_{%d} = v_{%d}(+h) - v_{%d}(-h)",harm+1,harm+1,harm+1));
        DvpT->DrawCopy();
        c1->SaveAs(Form("v%d/pT_c%d.pdf",harm+1,c));
       delete c1;
}

void PlotEta(TH1F *PvEta, TH1F *NvEta, TH1F *DvEta,Int_t c, Int_t harm){
	Double_t y = PvEta->GetBinContent(21);//to set the y range automatically 
        //Plotting the subsampled pos/neg and their difference as in Multiplot.C for every harmonic==========
        TCanvas *cEta = new TCanvas("ceta", "etacanvas", 800, 800);
        TPad *padeta1 = new TPad("padeta1", "padeta1", 0, 0.5, 1., 1.);
	padeta1->SetBottomMargin(0);
        padeta1->SetLeftMargin(0.2);
        padeta1->Draw();
        padeta1->cd();
        PvEta->GetXaxis()->SetRangeUser(-0.9,0.9);
	//PvEta->GetYaxis()->SetNdivisions(4);
        PvEta->GetYaxis()->SetRangeUser(y-y*0.35,y+y*0.35);
	PvEta->SetLineColor(c+1);
        PvEta->SetLineWidth(1);
        PvEta->GetYaxis()->SetTitle(Form("differential flow v_{%d}",harm+1));
        PvEta->SetStats(0);
	PvEta->SetMarkerStyle(22);
	PvEta->SetMarkerColor(c+1);
	PvEta->SetTitle(" ");
        PvEta->DrawCopy();
	//NvEta->GetYaxis()->SetNdivisions(4);
        NvEta->SetLineColor(c+2);
        NvEta->SetLineWidth(1);
	NvEta->SetTitle(" ");
	NvEta->SetMarkerStyle(23);
	NvEta->SetMarkerColor(c+2);
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
	DvEta->SetMarkerStyle(20);
	DvEta->SetMarkerColor(c);
        DvEta->SetLineColor(c);
        DvEta->SetLineWidth(1);
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
        Double_t fPtDiffNBins = 29;
        Double_t fCRCPtBins[30];
        Double_t PtBins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5.};
	for(Int_t r=0; r<30; r++){fCRCPtBins[r] = PtBins[r];}

	Double_t fEtaDiffNBins = 50;
	Double_t fCRCEtaBins[51]={0};
	Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
	for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}
		//hists for crosscheck
		TProfile *mdiff;
		TH1F* RMs;
		if(pT == true){
			 mdiff = new TProfile("Unweighted Bin Content","Unweighted Bin Content",fPtDiffNBins,fCRCPtBins,"s");
			// mdiff->SetDefaultSumw2(kTRUE);
        		//Mdiff->SetDefaultSumw2(kTRUE);
        		 RMs = new TH1F("rms evolution","rms evolution",fPtDiffNBins,fCRCPtBins);	 	
		}else if(pT == false){
                	 mdiff = new TProfile("Unweighted Bin Content","Unweighted Bin Content",fEtaDiffNBins,fCRCEtaBins,"s");
                	 RMs = new TH1F("rms evolution","rms evolution",fEtaDiffNBins,fCRCEtaBins);
		}
	TProfile *mint = new TProfile("Weighted Bin Content","Weighted Bin Content",10,0,9,"s");
	TH1D *CC = new TH1D();
	TH1D *Hist = new TH1D("CrossCheck","CrossCheck",10,0,9);
//================================================================================================
	TH1F *copy[cmax];//Multi centrality plots
	for(Int_t c=cent; c<=cmax;c++){
	cout<<"Centrality: "<<c<<"0-"<<(c+1)<<"0"<<endl;//actually 5.02
	//TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
	TString input = Form("/data/alice/jlomker/AVFD/result/dirID-0/NoBField/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0",c,c+1);
	// !!! For harm = 0 we could also try to use the fFlowQCQv1[pos=0, neg=1][cen][pt = 0, eta = 1] !!!
	//for(Int_t harm=0; harm<3; harm++){
		cout<<"Harmonic: "<<(harm+1)<<endl;;
		//apply subsampling to all histograms from pos/neg particles=======================
		if(pT == true){
			findSpread(c,harm,"Pt");//finds spread from the 10 splitfiles for pos, neg and delta histogram as input for the subsampling
			TH1F *PvpT = new TH1F("pvpT", "pvpT", fPtDiffNBins, fCRCPtBins);
			PvpT = Subsampling(input, "Pt", "Pos", PvpT, c, harm);
			//findSpread(c,Form("fFlowQCQv1[0][%d][0]",c),fPtDiffNBins);
			//findSpread(c,harm,Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm),fPtDiffNBins);
			//PvpT = Subsampling(input,Form("fFlowQCFinalPtDifHist[0][%d][%d][0]",c,harm),fCRCPtBins, PvpT, c, harm, mint, mdiff,RMs, CC,Hist);
			//PvpT = Subsampling(input,Form("fFlowQCQv1[0][%d][0]",c),fCRCPtBins, PvpT, c, harm, mint, mdiff,RMs, CC,Hist);
			TH1F *NvpT = new TH1F("nvpT","nvpT",fPtDiffNBins,fCRCPtBins);
			//findSpread(c,Form("fFlowQCQv1[1][%d][0]",c),fPtDiffNBins);
			//findSpread(c,harm,Form("fFlowQCFinalPtDifHist[1][%d][%d][0]",c,harm),fPtDiffNBins);
			//NvpT = Subsampling(input,Form("fFlowQCQv1[1][%d][0]",c),fCRCPtBins, NvpT, c ,harm,mint, mdiff,RMs, CC,Hist);
			NvpT = Subsampling(input,"Pt","Neg",NvpT, c ,harm);
			TH1F *DvpT = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
			//findSpread(c,harm,Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",c,harm),fPtDiffNBins);
			DvpT = Subsampling(input,"Pt","Diff", DvpT,c,harm);
			//DvpT = Subsampling(input,Form("fFlowQCFinalPtDifDeltaHist[%d][%d]",c,harm), fCRCPtBins, DvpT, c, harm,mint, mdiff,RMs, CC,Hist);
			//Multi centrality plot
			copy[c] = (TH1F*) DvpT->Clone("copy dvpt");
			//Plot the pos/neg with difference================================================= 
			PlotpT(PvpT, NvpT, DvpT, c, harm);// plot will delete all hists!
		
        		PvpT->Delete();
        		NvpT->Delete();
        		DvpT->Delete();
		}//end of pT if
		else if(pT == false){
			findSpread(c,harm,"Eta");
			TH1F *PvEta = new TH1F("pvEta", "pvEta", fEtaDiffNBins, fCRCEtaBins);
			PvEta = Subsampling(input, "Eta", "Pos", PvEta, c, harm);		
			//PvEta = Subsampling(input,Form("fFlowQCFinalEtaDifHist[0][%d][%d][0]",c,harm),fCRCEtaBins, PvEta, c, harm,mint, mdiff,RMs, CC,Hist);
			TH1F *NvEta = new TH1F("nvEta", "nvEta", fEtaDiffNBins, fCRCEtaBins);
			//findSpread(c,harm,Form("fFlowQCFinalEtaDifHist[1][%d][%d][0]",c,harm),fEtaDiffNBins);
			NvEta = Subsampling(input, "Eta", "Neg", NvEta, c, harm);
			TH1F *DvEta = new TH1F("deltavneta","deltavneta",fEtaDiffNBins,fCRCEtaBins);
			//DvEta = Dv(PvEta, NvEta, fCRCEtaBins, DvEta);
			//findSpread(c, harm,Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",c,harm),fEtaDiffNBins);
			DvEta = Subsampling(input, "Eta", "Diff", DvEta, c, harm);
			//DvEta = Subsampling(input, Form("fFlowQCFinalEtaDifDeltaHist[%d][%d]",c,harm), fCRCEtaBins, DvEta, c, harm,mint, mdiff,RMs, CC,Hist);
			//Multi centrality plot
			copy[c] = (TH1F*) DvEta->Clone("copy dvpt");
			//Plot the pos/neg with difference=================================================
			PlotEta(PvEta, NvEta, DvEta, c, harm);
			
        		PvEta->Delete();
	        	NvEta->Delete();
		        DvEta->Delete();
		}//end of eta if
	}//centrality loop

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
	if(pT==true){Cen->SaveAs(Form("v%d/Pt_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
	else if(pT == false){Cen->SaveAs(Form("v%d/Eta_Multi_cen%d_%d.pdf",harm+1,cent,cmax));}
	delete Cen;

	//}//end of hr<fkFlowNHarm

}//end void


