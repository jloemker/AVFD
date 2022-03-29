TH1F *Difference(TH1F *B, TH1F *NoB, Double_t arr[], TString title){
	Double_t NBins = B ->GetNbinsX();
	if(B->GetMean() == 0 || NoB->GetMean() == 0){cout<<" !!! This centrality sample is not yet produced -> These plots are fake news !!!"<<endl;}
	TH1F *difference = new TH1F("difference_"+title,"",NBins, arr);
	for(Int_t i = 1; i <= NBins; i++){
		difference->SetBinContent(i, B->GetBinContent(i) - NoB->GetBinContent(i) );
		difference->SetBinError(i, sqrt(B->GetBinError(i)*B->GetBinError(i) + NoB->GetBinError(i)*NoB->GetBinError(i)) );
	}
	return difference;
}//end difference


void findSpread(Int_t c, Int_t harm,TString spectrum,TString input){
  TFile *f[10];
  TList *list[10];
  TH1D *fHist[10];
  Int_t Nbins = 0;
  if(spectrum == "Pt"){Nbins = 30;}
  else if(spectrum == "Eta"){Nbins = 50;}
  TH1D *fHistPos[10];
  TH1D *fHistNeg[10];
  TH1F *fHistSpreadPos[Nbins];
  TH1F *fHistSpreadNeg[Nbins];
  TH1F *fHistSpreadDiff[Nbins];
 for(Int_t i = 0; i< Nbins; i++){
        fHistSpreadPos[i] = new TH1F(Form("fHistSpreadPos"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
        fHistSpreadNeg[i] = new TH1F(Form("fHistSpreadNeg"+spectrum+"Bin%d",i+1),"",1000,0,0.3);
        fHistSpreadDiff[i] = new TH1F(Form("fHistSpreadDiff"+spectrum+"Bin%d",i+1),"",1000,-0.3,0.3);
  }
  for(Int_t iFile = 0; iFile < 10; iFile++) {
    f[iFile] = TFile::Open(input+Form("/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
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


TH1F *Subsampling(TString spectrum, TString charge, TH1F *result, Int_t c, Int_t harm){
        TH1F *o;
        Int_t Nbins = result->GetNbinsX();
        TFile *output = output-> Open(Form("output_findSpread/spread_"+spectrum+"_c%d0_%d0.root",c,c+1));
        for(Int_t i=1; i<=Nbins; i++){
                Double_t sigma = 0.;
                Double_t mean = 0.;
                TString spreadstring;
                spreadstring = Form("fHistSpread"+charge+spectrum+"Bin%d",i);
                o = (TH1F*) output->Get(spreadstring);
                sigma = o->GetRMS();
                mean = o->GetMean();
                result->SetBinContent(i, mean);
                result->SetBinError(i, sigma);
                o->Reset();
		
        }//end Nbins
	output->Close();
return result;
}//end subsampling

void CompareEffectB(Int_t cent, Int_t cmax, Int_t harm){//add writing to external file to later calcuate the difference for B vs noB
        gROOT->SetBatch();//to avoid opening the plots ak bad wifi struggle
        Double_t fPtDiffNBins = 29;
        Double_t fCRCPtBins[30];
        Double_t PtBins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5.};
        for(Int_t r=0; r<30; r++){fCRCPtBins[r] = PtBins[r];}

        Double_t fEtaDiffNBins = 50;
        Double_t fCRCEtaBins[51]={0};
        Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
        for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}


	for(Int_t c=cent; c<=cmax;c++){
		gROOT->SetBatch();
		gStyle->SetOptStat(0);
		TH1F *DiffPospT, *DiffNegpT, *DiffDelpT;
		TH1F *DiffPosEta, *DiffNegEta, *DiffDelEta;
        	cout<<"Centrality: "<<c<<"0-"<<(c+1)<<"0"<<endl;//actually 5.02
        	TString inputB = "/data/alice/jlomker/AVFD/result/dirID-0/BField0.2/split";
		TString input = "/data/alice/jlomker/AVFD/result/dirID-0/NoBField/split";//something is wrong here too!
                cout<<"Harmonic: "<<(harm+1)<<endl;//apply subsampling to all histograms from pos/neg particles=======================

                findSpread(c,harm,"Pt",input);//finds spread from the 10 splitfiles for pos, neg and delta histogram as input for the subsampling
	        TH1F *PvpT = new TH1F("pvpT", "pvpT", fPtDiffNBins, fCRCPtBins);
		PvpT = Subsampling("Pt", "Pos", PvpT, c, harm);
	        TH1F *NvpT = new TH1F("nvpT","nvpT",fPtDiffNBins,fCRCPtBins);
                NvpT = Subsampling("Pt","Neg",NvpT, c ,harm);
                TH1F *DvpT = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
                DvpT = Subsampling("Pt","Diff", DvpT,c,harm);
		findSpread(c,harm,"Pt",inputB);
                TH1F *BPvpT = new TH1F("BpvpT", "BpvpT", fPtDiffNBins, fCRCPtBins);
                BPvpT = Subsampling("Pt", "Pos", BPvpT, c, harm);
                TH1F *BNvpT = new TH1F("BnvpT","BnvpT",fPtDiffNBins,fCRCPtBins);
                BNvpT = Subsampling("Pt","Neg",BNvpT, c ,harm);
                TH1F *BDvpT = new TH1F("BdvnpT","BdvnpT",fPtDiffNBins,fCRCPtBins);
                BDvpT = Subsampling("Pt","Diff", BDvpT,c,harm);
				
		DiffPospT = Difference(BPvpT,PvpT,fCRCPtBins,"PosPt");
	        DiffNegpT = Difference(BNvpT,NvpT,fCRCPtBins,"NegPt");
		DiffDelpT = Difference(BDvpT,DvpT,fCRCPtBins,"DelPt");
		TCanvas *p1 = new TCanvas("p1","PtDifferencesB",160,200);
		p1->Divide(1,3);
		p1->cd(1);
		DiffPospT->SetTitle(Form("Centrality %d0 - %d0 #Delta [B(#tau = 0.2) - B(#tau = 0)]",c,c+1));
                DiffPospT->GetYaxis()->SetTitle(Form("v_{%d}(+h): B(#tau) - B(0)",harm+1));
                DiffPospT->GetXaxis()->SetTitle("p_{T} [GeV]");
                DiffPospT->GetXaxis()->SetRangeUser(0.,3.4);
		DiffPospT->Draw("h");
		p1->cd(2);
                DiffNegpT->GetYaxis()->SetTitle(Form("v_{%d}(-h): B(#tau) - B(0)",harm+1));
                DiffNegpT->GetXaxis()->SetTitle("p_{T} [GeV]");
		DiffNegpT->GetXaxis()->SetRangeUser(0.,3.4);
		DiffNegpT->Draw("h");
		p1->cd(3);
                DiffDelpT->GetYaxis()->SetTitle(Form("#Delta v_{%d}(+h - -h): B(#tau) - B(0)",harm+1));
                DiffDelpT->GetXaxis()->SetTitle("p_{T} [GeV]");
                DiffDelpT->GetXaxis()->SetRangeUser(0.,3.4);
		DiffDelpT->Draw("h");
		p1->SaveAs(Form("difference/v%d/pT_B0.2_NoB_cen%d.pdf",harm+1,c));
		delete p1;

		DiffPospT->Delete();
		DiffNegpT->Delete();
		DiffDelpT->Delete();

                PvpT->Delete();
                NvpT->Delete();
                DvpT->Delete();
                BPvpT->Delete();
                BNvpT->Delete();
                BDvpT->Delete();

                findSpread(c,harm,"Eta",input);
                TH1F *PvEta = new TH1F("pvEta", "pvEta", fEtaDiffNBins, fCRCEtaBins);
                PvEta = Subsampling("Eta", "Pos", PvEta, c, harm);
                TH1F *NvEta = new TH1F("nvEta", "nvEta", fEtaDiffNBins, fCRCEtaBins);
                NvEta = Subsampling("Eta", "Neg", NvEta, c, harm);
                TH1F *DvEta = new TH1F("deltavneta","deltavneta",fEtaDiffNBins,fCRCEtaBins);
                DvEta = Subsampling( "Eta", "Diff", DvEta, c, harm);

                findSpread(c,harm,"Eta",inputB);
                TH1F *BPvEta = new TH1F("BpvEta", "BpvEta", fEtaDiffNBins, fCRCEtaBins);
                BPvEta = Subsampling("Eta", "Pos", BPvEta, c, harm);
                TH1F *BNvEta = new TH1F("BnvEta", "BnvEta", fEtaDiffNBins, fCRCEtaBins);
                BNvEta = Subsampling("Eta", "Neg", BNvEta, c, harm);
                TH1F *BDvEta = new TH1F("Bdeltavneta","Bdeltavneta",fEtaDiffNBins,fCRCEtaBins);
                BDvEta = Subsampling( "Eta", "Diff", BDvEta, c, harm);

		DiffPosEta = Difference(BPvEta,PvEta,fCRCEtaBins,"PosEta");//difference = result with B - result without B
		DiffNegEta = Difference(BNvEta,NvEta,fCRCEtaBins,"NegEta");
		DiffDelEta = Difference(BDvEta,DvEta,fCRCEtaBins,"DelEta");
		TCanvas *e1 = new TCanvas("e1","EtaDifferencesB",160,200);
		e1->Divide(1,3);
		e1->cd(1);
		DiffPosEta->SetTitle(Form("Centrality %d0 - %d0 #Delta [B(#tau = 0.2) - B(#tau = 0)]",c,c+1));
                DiffPosEta->GetYaxis()->SetTitle(Form("v_{%d}(+h): B(#tau) - B(0)",harm+1));
                DiffPosEta->GetXaxis()->SetTitle("#eta");
		DiffPosEta->GetXaxis()->SetRangeUser(-0.9,0.9);
		DiffPosEta->Draw("h");
		e1->cd(2);
                DiffNegEta->GetYaxis()->SetTitle(Form("v_{%d}(-h): B(#tau) - B(0)",harm+1));
                DiffNegEta->GetXaxis()->SetTitle("#eta");
		DiffNegEta->GetXaxis()->SetRangeUser(-0.9,0.9);
		DiffNegEta->Draw("h");
		e1->cd(3);
		DiffDelEta->GetYaxis()->SetTitle(Form("#Delta v_{%d}(+h - -h): B(#tau) - B(0)",harm+1));
		DiffDelEta->GetXaxis()->SetTitle("#eta");
		DiffDelEta->GetXaxis()->SetRangeUser(-0.9,0.9);
                DiffDelEta->Draw("h");
		e1->SaveAs(Form("difference/v%d/eta_B0.2_NoB_cen%d.pdf",harm+1,c));
		delete e1;

		DiffPosEta->Delete();
	        DiffNegEta->Delete();
		DiffDelEta->Delete();

                PvEta->Delete();
                NvEta->Delete();
                DvEta->Delete();
                BPvEta->Delete();
                BNvEta->Delete();
                BDvEta->Delete();
	}//centrality

}//endvoid

