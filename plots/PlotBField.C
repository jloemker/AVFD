#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "TStyle.h"
#include <math.h>

void PlotBField(){
	gROOT->SetBatch();
	gStyle->SetOptStat(0);
	TFile* f = new TFile("/project/alice/users/jlomker/AVFD/plots/BField/B_fields_tauInit0.4_5.02TeV_lifetime_0.2.root");//B_fields_5.02TeV_lifetime_Panos_BaselineB.root");//B_fields_5.02TeV_lifetime_0.root");
        TProfile2D *feBxfield2D;
        TProfile2D *feByfield2D;
        TProfile2D *feExfield2D;
        TProfile2D *feEyfield2D;
	for(Int_t c = 1; c<7; c++){
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
		TPad *p1 = new TPad("Bx","Bx", 0.,0.5,0.5,1.0);
		p1->SetTopMargin(0.01);
		p1->SetBottomMargin(0.1);
		p1->SetLeftMargin(0.1);
		p1->SetRightMargin(0);	
		p1->Draw();
		p1->cd();
		TLegend *l = new TLegend(0.23,0.73,0.32,0.8);
		l->SetHeader("B_{x}","C");
		l->SetFillStyle(0);
		l->SetBorderSize(0);
		feBxfield2D->SetTitle(" ");
		feBxfield2D->GetXaxis()->SetTitle("x");
		feBxfield2D->GetYaxis()->SetTitle("y");
                feBxfield2D->GetXaxis()->SetTitleSize(0.08);
                feBxfield2D->GetYaxis()->SetTitleSize(0.08);
                feBxfield2D->GetXaxis()->SetTitleOffset(0.5);
                feBxfield2D->GetYaxis()->SetTitleOffset(0.5);
		feBxfield2D->GetXaxis()->SetLabelSize(0.055);
                feBxfield2D->GetYaxis()->SetLabelSize(0.055);
		feBxfield2D->Draw("COL");
		feBxfield2D->GetZaxis()->SetLabelSize(0.001);
		feBxfield2D->GetZaxis()->SetRangeUser(-52,52);
		feBxfield2D->GetZaxis()->SetLabelOffset(2);
		c1->cd();
		l->Draw();
                TPad *p2 = new TPad("By","By", 0.5,0.5,1.0,1.0);
                p2->SetTopMargin(0.01);
		p2->SetBottomMargin(0.1);
		p2->SetLeftMargin(0.);
		p2->SetRightMargin(0.1);
		p2->Draw();
                p2->cd();
                TLegend *l1 = new TLegend(0.68,0.73,0.77,0.8);
                l1->SetHeader("B_{y}","C");
                l1->SetFillStyle(0);
                l1->SetBorderSize(0);
                feByfield2D->SetTitle(" ");
                feByfield2D->GetXaxis()->SetTitle("x");
                feByfield2D->GetYaxis()->SetTitle(" ");
		feByfield2D->GetXaxis()->SetTitleSize(0.08);
                feByfield2D->GetXaxis()->SetTitleOffset(0.5);
                feByfield2D->GetXaxis()->SetLabelSize(0.055);
                feByfield2D->Draw("COLZ");
                feByfield2D->GetZaxis()->SetLabelSize(0.045);
                feByfield2D->GetZaxis()->SetRangeUser(-52,52);
		c1->cd();
		l1->Draw();
		TPad *p3 = new TPad("Ex","Ex", 0.,0.,0.5,0.5);
                p3->SetTopMargin(0.01);
		p3->SetBottomMargin(0.1);
		p3->SetLeftMargin(0.1);
		p3->SetRightMargin(0);
		p3->Draw();
                p3->cd();
                TLegend *l2 = new TLegend(0.23,0.23,0.32,0.3);
                l2->SetHeader("E_{x}","C");
                l2->SetFillStyle(0);
                l2->SetBorderSize(0);
                feExfield2D->SetTitle(" ");
                feExfield2D->GetXaxis()->SetTitle("x");
                feExfield2D->GetYaxis()->SetTitle("y");
                feExfield2D->GetXaxis()->SetTitleSize(0.08);
                feExfield2D->GetYaxis()->SetTitleSize(0.08);
                feExfield2D->GetXaxis()->SetTitleOffset(0.5);
                feExfield2D->GetYaxis()->SetTitleOffset(0.5);
                feExfield2D->GetXaxis()->SetLabelSize(0.055);
                feExfield2D->GetYaxis()->SetLabelSize(0.055);
                feExfield2D->Draw("COL");
                feExfield2D->GetZaxis()->SetLabelSize(0.001);
		feExfield2D->GetZaxis()->SetRangeUser(-52,52);
		c1->cd();
		l2->Draw();
                TPad *p4 = new TPad("Ey","Ey", 0.5,0.,1.0,0.5);
                p4->SetTopMargin(0.01);
                p4->SetBottomMargin(0.1);
                p4->SetLeftMargin(0.);
		p4->SetRightMargin(0.1);
		p4->Draw();
                p4->cd();
                TLegend *l3 = new TLegend(0.68,0.23,0.77,0.3);
                l3->SetHeader("E_{y}","C");
                l3->SetFillStyle(0);
                l3->SetBorderSize(0);
                feEyfield2D->SetTitle(" ");
                feEyfield2D->GetXaxis()->SetTitle("x");
                feEyfield2D->GetYaxis()->SetTitle(" ");
		feEyfield2D->GetXaxis()->SetTitleSize(0.08);
                //feEyfield2D->GetYaxis()->SetTitleSize(0.08);
                feEyfield2D->GetXaxis()->SetTitleOffset(0.5);
                //feEyfield2D->GetYaxis()->SetTitleOffset(0.5);
                feEyfield2D->GetXaxis()->SetLabelSize(0.055);
                //feEyfield2D->GetYaxis()->SetLabelSize(0.055);
                feEyfield2D->Draw("COLZ");
                feEyfield2D->GetZaxis()->SetLabelSize(0.045);
                feEyfield2D->GetZaxis()->SetRangeUser(-52,52);
		c1->cd();
		l3->Draw();
		c1->SaveAs(Form("BField/Initial0.4_BFields0.2_%d0_%d0.png", c, c+1));
		feBxfield2D->Reset();
		feByfield2D->Reset();
		feExfield2D->Reset();
		feEyfield2D->Reset();
		delete c1;
	}
}//For one centrality: add time evolution plots/modify the read macro for that or even something in EbE... - then add multi centrality loop
