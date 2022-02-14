void findSpread(Int_t c, TString obj) {
  TFile *f[10];
  TList *list[10];
  TH1D *fHist[10];
  TH1F *fHistSpread[36];
  for(Int_t i = 0; i < 36; i++)
    fHistSpread[i] = new TH1F(Form("fHistSpreadPtBin%d",i+1),"",1000,0, 0.3);
  for(Int_t iFile = 0; iFile < 10; iFile++) {
    f[iFile] = TFile::Open(Form("/data/alice/jlomker/AVFD/result/dirID-0/split/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    if((!f[iFile])||(!f[iFile]->IsOpen())) {
    cout<<"File "<<iFile<<" not found..."<<endl;
    return;
    }
    
    list[iFile] = dynamic_cast<TList *>(f[iFile]->Get("FlowQCList"));
    if(!list[iFile]) {
      cout<<"Input list of file "<<iFile<<" not found..."<<endl;
      return;
    }
    fHist[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(obj));
    if(!fHist[iFile]) {
      cout<<"Histogram of file "<<iFile<<" not found..."<<endl;
      return;
    }
    fHist[iFile]->GetXaxis()->SetRangeUser(0.2,5);
    fHist[iFile]->SetMarkerColor(iFile+1);
    fHist[iFile]->SetLineColor(iFile+1);
    fHist[iFile]->Draw("ESAME");
    for(Int_t iBin = 1; iBin <= fHist[iFile]->GetNbinsX(); iBin++){ 
	cout<<fHist[iFile]->GetBinContent(iBin)<<endl;
      fHistSpread[iBin-1]->Fill(fHist[iFile]->GetBinContent(iBin));}				
  }//loop over files
  TFile *fOutput = new TFile(Form("spread_"+obj+"_c%d0_%d0.root",c,c+1),"recreate");
  for(Int_t iBin = 1; iBin <= fHist[0]->GetNbinsX(); iBin++){ 
    fHistSpread[iBin-1]->Write();}
    fOutput->Close();
}
