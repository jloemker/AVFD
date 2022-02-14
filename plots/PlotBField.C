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
        TProfile2D *feBxfield2D;
        TProfile2D *feByfield2D;
        TProfile2D *feExfield2D;
        TProfile2D *feEyfield2D;
	for(Int_t c = 0; c<7; c++){
		if(c==0){
                feBxfield2D = (TProfile2D*) f->Get(Form("feBxfield2D_cent5_%d0",c+1));
                feByfield2D = (TProfile2D*) f->Get(Form("feByfield2D_cent5_%d0",c+1));
                feExfield2D = (TProfile2D*) f->Get(Form("feExfield2D_cent5_%d0",c+1));
                feEyfield2D = (TProfile2D*) f->Get(Form("feEyfield2D_cent5_%d0",c+1));
		}else{
		feBxfield2D = (TProfile2D*) f->Get(Form("feBxfield2D_cent%d0_%d0",c,c+1));
	        feByfield2D = (TProfile2D*) f->Get(Form("feByfield2D_cent%d0_%d0",c,c+1));
        	feExfield2D = (TProfile2D*) f->Get(Form("feExfield2D_cent%d0_%d0",c,c+1));
        	feEyfield2D = (TProfile2D*) f->Get(Form("feEyfield2D_cent%d0_%d0",c,c+1));
		}
		TCanvas *c1 = new TCanvas("c1","Fields");//,200,10,700,500);
		c1->Divide(2,2);
		c1->cd(1);
		//gStyle->SetPalette(1);//to set colors ... very intense with 1
		//feBxfield2D->SetContour(10);//to group intensities more roughly
		feBxfield2D->Draw("COLZ");
		c1->cd(2);
		feByfield2D->Draw("COLZ");
		c1->cd(3);
		feExfield2D->Draw("COLZ");
		c1->cd(4);
		feEyfield2D->Draw("COLZ");
		c1->SaveAs(Form("BField/Fields_%d0_%d0.pdf", c, c+1));
		feBxfield2D->Reset();
		feByfield2D->Reset();
		feExfield2D->Reset();
		feEyfield2D->Reset();
		delete c1;
	}
}//For one centrality: add time evolution plots/modify the read macro for that or even something in EbE... - then add multi centrality loop
