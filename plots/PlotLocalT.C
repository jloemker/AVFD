#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotLocalT(Int_t snap1,Int_t snap2, Int_t snap3, Int_t snap4){//argument for time stamps ?
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	TFile* f = new TFile("/project/alice/users/jlomker/AVFD/plots/BField/Tevolution_5.02TeV_Blifetime_0.2.root");
        TProfile2D *T1local2D;
	TProfile2D *T2local2D;
        TProfile2D *T3local2D;
	TProfile2D *T4local2D;
	Float_t tau1 = 0.600000 + snap1*0.010000;//snap0 = 0.6 -> snapi = 0.6 + i*0.01
	Float_t tau2 = 0.600000 + snap2*0.010000;
	Float_t tau3 = 0.600000 + snap3*0.010000;
	Float_t tau4 = 0.600000 + snap4*0.010000;
	cout<<tau1<<endl;	
	for(Int_t c = 2; c<3; c++){
//cout <<Form("Tlocal2D_tau%f_cent%d0_%d0",tau1,c,c+1)<<endl;//change the read script to produce files with integer ID instead float value ! 
		T1local2D= (TProfile2D*) f->Get(Form("Tlocal2D_tau%d_cent%d0_%d0;1",snap1,c,c+1));
	        T2local2D= (TProfile2D*) f->Get(Form("Tlocal2D_tau%d_cent%d0_%d0;1",snap2,c,c+1));
        	T3local2D= (TProfile2D*) f->Get(Form("Tlocal2D_tau%d_cent%d0_%d0;1",snap3,c,c+1));
        	T4local2D= (TProfile2D*) f->Get(Form("Tlocal2D_tau%d_cent%d0_%d0;1",snap4,c,c+1));
		
		TCanvas *c1 = new TCanvas("c1","Fields");//,200,10,700,500);
		c1->Divide(2,2);
		c1->cd(1);//gStyle->SetPalette(1);//to set colors ... very intense with 1
		T1local2D->SetTitle(Form("local temperature(#tau = %f)",tau1));
		T1local2D->Draw("COLZ");//->SetContour(10);//to group intensities more roughly
		c1->cd(2);
		T2local2D->SetTitle(Form("local temperature(#tau = %f)",tau2));
		T2local2D->Draw("COLZ");
		c1->cd(3);
		T3local2D->SetTitle(Form("local temperature(#tau = %f)",tau3));
		T3local2D->Draw("COLZ");
		c1->cd(4);
		T4local2D->SetTitle(Form("local temperature(#tau = %f)",tau4));
		T4local2D->Draw("COLZ");
		c1->SaveAs(Form("BField/BField0.2/Tevolution_tauID%d-%d_cent%d0_%d0.pdf",snap1,snap4, c, c+1));//double tau..
		T1local2D->Reset();
		T2local2D->Reset();
		T3local2D->Reset();
		T4local2D->Reset();
		delete c1;
	}
}
