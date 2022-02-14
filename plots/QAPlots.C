void QAPlots(Int_t c){
gStyle -> SetOptStat(0);
TFile *f = new TFile(Form("/data/alice/jlomker/AVFD/result/dirID-0/NoBField/Result_5.02TeV_pT_0_eta_0_Cent%d0_%d0.root",c,c+1));
TList *l = (TList*) f->Get("QAList;1");

//TH1D *N = (TH1D*) l->FindObject("fnEvent");
TH1D *Multiplicity = (TH1D*) l->FindObject("fMultChargedParticlesDistribution");
TH1D *Pt = (TH1D*) l ->FindObject("fPtChargedParticlesDistribution");
TH1D *Eta = (TH1D*) l->FindObject("fEtaChargedParticlesDistribution");
TH1D *Phi = (TH1D*) l->FindObject("fPhiChargedParticlesDistribution");

TCanvas *c1 = new TCanvas("QA","QAHistograms");
c1->Divide(2,2);
c1->cd(1);
//N->GetXaxis()->SetTitle("?");
//N->Draw("h");
Multiplicity->GetXaxis()->SetTitle("Charged Particle Multiplicity");
Multiplicity->Draw("h");
Multiplicity->GetXaxis()->SetRangeUser(0,1500);
c1->cd(2);
Pt->GetXaxis()->SetTitle("p_{T} [GeV]");
Pt->GetXaxis()->SetRangeUser(0,5.5);
Pt->Draw("h");
c1->cd(3);
Eta->GetXaxis()->SetTitle("#eta");
Eta->GetXaxis()->SetRangeUser(-0.9,0.9);
Eta->Draw("h");
c1->cd(4);
Phi->GetXaxis()->SetTitle("#phi");
Phi->Draw("h");
c1->SaveAs(Form("ControlPlots/QA_%d0_%d0.pdf",c,c+1));
}
