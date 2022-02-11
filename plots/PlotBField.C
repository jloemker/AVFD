#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotBField(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	TFile* f = new TFile("/project/alice/users/jlomker/AVFD/plots/BField/B_fields_5.02TeV_lifetime_Panos_BaselineB.root");//B_fields_5.02TeV_lifetime_0.root");
	//TObj* l = (TObj*) f->Get("feBxfield2D_cent20_30");
//f->FindObject("");
//f->Get(
cout<<"here ?"<<endl;
	TProfile2D *feBxfield2D = (TProfile2D*) f->Get("feBxfield2D_cent20_30");
        TProfile2D *feByfield2D = (TProfile2D*) f->Get("feByfield2D[2]");
        TProfile2D *feExfield2D = (TProfile2D*) f->Get("feExfield2D[2]");
        TProfile2D *feEyfield2D = (TProfile2D*) f->Get("feEyfield2D[2]");

	TCanvas *c1 = new TCanvas("c1","Bx",200,10,700,500);
	//gStyle->SetPalette(1);//to set colors ... very intense with 1

	//feBxfield2D->SetContour(10);//to group intensities more roughly
	feBxfield2D->Draw("COLZ");

	//feBxfield2D->Draw("TEXT" "SAME");

	c1->SaveAs("BField/Bx_test.pdf");
}//For one centrality: Put all B and all E fFields in one canvas - then add time evolution - then add multi centrality loop
