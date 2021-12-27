#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>


TH1F *Subsampling(TString file, TString obj, TH1F *result){
TFile *f;
TH1F *co;
TH1F *h0,*h1,*h2,*h3,*h4,*h5,*h6,*h7,*h8,*h9;
Int_t Nbins = result->GetNbinsX();
cout<<"subsampling Nbins"<<Nbins<<endl;
for(Int_t split = 0; split<10;split ++){
	f = new TFile(Form(file+"_split_%d.root",split));
	TList *l = (TList*) f->Get("FlowQCList;1");
	if(split == 0){h0 = (TH1F*) l->FindObject(obj);}
	if(split == 1){h1 = (TH1F*) l->FindObject(obj);}
	if(split == 2){h2 = (TH1F*) l->FindObject(obj);}
	if(split == 3){h3 = (TH1F*) l->FindObject(obj);}
	if(split == 4){h4 = (TH1F*) l->FindObject(obj);}
	if(split == 5){h5 = (TH1F*) l->FindObject(obj);}
	if(split == 6){h6 = (TH1F*) l->FindObject(obj);}
	if(split == 7){h7 = (TH1F*) l->FindObject(obj);}
	if(split == 8){h8 = (TH1F*) l->FindObject(obj);}
	if(split == 9){h9 = (TH1F*) l->FindObject(obj);}
}
for(Int_t i=1; i<=Nbins; i++){
	Double_t w[100];
	Double_t wTotal = 0.0;
	Double_t gCombinedValue = 0.0;
	Double_t gCombinedError = 0.0;
	for(Int_t s=0; s<9;s++){
		if(s==0){co = (TH1F*)h0->Clone("co");}
		if(s==1){co = (TH1F*)h1->Clone("co");}
		if(s==2){co = (TH1F*)h2->Clone("co");}
		if(s==3){co = (TH1F*)h3->Clone("co");}
		if(s==4){co = (TH1F*)h4->Clone("co");}
		if(s==5){co = (TH1F*)h5->Clone("co");}
		if(s==6){co = (TH1F*)h6->Clone("co");}
		if(s==7){co = (TH1F*)h7->Clone("co");}
		if(s==8){co = (TH1F*)h8->Clone("co");}
		if(s==9){co = (TH1F*)h9->Clone("co");}	
	
		w[i] = 1./TMath::Power(co->GetBinError(i),2);//maybe i+1
		wTotal +=w [i];

		gCombinedValue += co->GetBinContent(i)*w[i];
		gCombinedError += TMath::Power(co->GetBinError(i)*w[i],2);
	}
	gCombinedValue /= wTotal;
	gCombinedError = (1./wTotal)*TMath::Sqrt(gCombinedError);
	if(abs(gCombinedValue)>0){
		result->SetBinContent(i, gCombinedValue);
		result->SetBinError(i, gCombinedError);
	}	
}//end Nbins
delete f;
return result;
}//end Subsampling
