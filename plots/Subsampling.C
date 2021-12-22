#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>
TH1F* Subsampling(TString file, TString obj){
TH1F *result;
TH1F *co;
for(Int_t Nbinsx = 0; Nbinsx<52; Nbinsx++){
        Double_t w[100];
        Double_t wTotal = 0.0;
	Double_t gCombinedValue = 0.0;
	Double_t gCombinedError = 0.0;

	for(Int_t split = 0; split<10;split++){
		TFile *f = new TFile(Form(file+"_split_%d.root",split);
		TList *l = (TList*) f->Get("FlowQCList;1");
		TH1F *h = (TH1F*) l->FindObject(obj);
		if(split = 0){
			Nbinsx = h0->GetNbinsX();
			result = (TH1F*)h0->Clone("result");
		}
		co = (TH1F*)h9->Clone("co");
		w[i] = 1./TMath::Power(co->GetBinError(i),2);//maybe i+1
        	wTotal +=w [i];

        	gCombinedValue += co->GetBinContent(i)*w[i];
        	gCombinedError += TMath::Power(co->GetBinError(i)*w[i],2);
	}	
	gCombinedValue /= wTotal;
	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
	result->SetBinContent(i, gCombinedValue);
	result->SetBinError(i, gCombinedError);
	}//Nbins
return result;
}//end function

