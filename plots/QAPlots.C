void QAPlots(Int_t c){
	gStyle -> SetOptStat(0);
	//For B tau = 0 baseline
	//TFile *f = new TFile(Form("// data/alice/jlomker/AVFD/result/dirID-0/NoBField/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",c,c+1));
	TFile *f = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-0/TestEM/tau_init_0.6/BField0.2/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",c,c+1));
	//For B tau = 0.2 baseline
	TFile *b = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-0/TestEM/tau_init_0.1/BField0.2/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",c,c+1));

	TList *l = (TList*) f->Get("QAList;1");
	TList *lb = (TList*) b->Get("QAList;1");
	
	TH1D *N = (TH1D*) l->FindObject("fnEvent");
	TH1D *Multiplicity = (TH1D*) l->FindObject("fMultChargedParticlesDistribution");
	TH1D *Pt = (TH1D*) l ->FindObject("fPtChargedParticlesDistribution");
	TH1D *Eta = (TH1D*) l->FindObject("fEtaChargedParticlesDistribution");
	TH1D *Phi = (TH1D*) l->FindObject("fPhiChargedParticlesDistribution");

	TH1D *Nb = (TH1D*) lb->FindObject("fnEvent");
	TH1D *Multiplicityb = (TH1D*) lb->FindObject("fMultChargedParticlesDistribution");
	TH1D *Ptb = (TH1D*) lb ->FindObject("fPtChargedParticlesDistribution");
	TH1D *Etab = (TH1D*) lb->FindObject("fEtaChargedParticlesDistribution");
	TH1D *Phib = (TH1D*) lb->FindObject("fPhiChargedParticlesDistribution");

	TH1D *newPt = new TH1D("Rebinned Pt", "Rebinned Pt",40,0.2,5);
	TH1D *newPtb = new TH1D("Rebinned PtB", "Rebinned PtB", 40,0.2,5);
	
	for(Int_t iBin = 1; iBin<=Pt->GetNbinsX(); iBin++){
	//cout<<Ptb->GetBinWidth()<<endl;
	Double_t center = Pt->GetBinCenter(iBin);
	Double_t val = Pt->GetBinContent(iBin);
	newPt->Fill(Pt->GetBinCenter(iBin), Pt->GetBinContent(iBin));
        newPtb->Fill(Ptb->GetBinCenter(iBin), Ptb->GetBinContent(iBin));
	}

        auto L1 = new TLegend();
        L1->SetHeader("Pb-Pb, 5.02TeV, pT: 0.2 -5 GeV, #tau B: 0.2","C");
        L1->AddEntry(Multiplicity,Form("#tau 0 = 0.6, Centrality %d0-%d0",c,c+1), "l");
        L1->AddEntry(Multiplicityb,Form("#tau 0 = 0.1, Centrality %d0-%d0",c,c+1), "l");

	TCanvas *c1 = new TCanvas("QA","QAHistograms");
	c1->Divide(2,3);
        c1->cd(1);
	//N->Scale(N->Integral());
        N->SetLineColor(2);
        N->GetYaxis()->SetTitle("Number Of Events");
        N->Draw("h");
        L1->Draw();
        c1->cd(2);
        Nb->GetYaxis()->SetTitle("Number of Events");
        //Nb->Scale(Nb->Integral());
	Nb->Draw("h");
        L1->Draw();
	c1->cd(3);
	Multiplicity->GetXaxis()->SetTitle("Charged Particle Multiplicity");
	Multiplicity->GetXaxis()->SetRangeUser(0,2000);
	Multiplicity->GetYaxis()->SetTitle("Normalized Number Of Events");
	Multiplicity->Scale(1/N->Integral());
	Multiplicity->Rebin(4);
	Multiplicity->SetLineColor(2);
	Multiplicity->Draw("h");
	Multiplicityb->Scale(1/Nb->Integral());
	Multiplicityb->Rebin(4);
	//Multiplicity->SetLineColor(2);
	Multiplicityb->Draw("SAME");
	c1->cd(4);
        newPt->Rebin(2);
        newPtb->Rebin(2);
	newPtb->GetXaxis()->SetTitle("p_{T} [GeV]");
	newPtb->GetXaxis()->SetRangeUser(0,4);
	newPtb->GetYaxis()->SetTitle("Normalized Number Of Events");
	newPtb->Scale(1/Nb->Integral());
	newPtb->Draw("h");
	newPt->SetLineColor(2);
	newPt->Scale(1/N->Integral());
	newPt->Draw("SAME");
	//L1->Draw();
	c1->cd(5);
	Etab->GetXaxis()->SetTitle("#eta");
	Etab->GetXaxis()->SetRangeUser(-0.9,0.9);
	Etab->GetYaxis()->SetTitle("Normalized Number Of Events");
	Etab->Scale(1/Nb->Integral());
	Etab->Draw("h");
	Eta->SetLineColor(2);
	Eta->Scale(1/N->Integral());
	Eta->Draw("SAME");
	c1->cd(6);
	Phib->GetXaxis()->SetTitle("#phi");
	Phib->GetYaxis()->SetTitle("Normalized Number Of Events");
	Phib->Scale(1/Nb->Integral());
	Phib->Draw("h");
	Phi->SetLineColor(2);
	Phi->Scale(1/N->Integral());
	Phi->Draw("SAME");
	c1->SaveAs(Form("controlPlots/QA_0.6_vs_0.1_C%d0_%d0.pdf",c,c+1));
}
