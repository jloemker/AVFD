TH1F *Difference(TH1F *B, TH1F *NoB, Double_t arr[], TString title){
	Double_t NBins = B ->GetNbinsX();
	if(B->GetMean() == 0 || NoB->GetMean() == 0){cout<<" !!! This is from too poor stats in "+title+" -> This set of plots are worth discussing !!!"<<endl;}
	TH1F *difference = new TH1F("difference_"+title,"",NBins, arr);
	for(Int_t i = 1; i <= NBins; i++){
		difference->SetBinContent(i, B->GetBinContent(i) - NoB->GetBinContent(i) );
		difference->SetBinError(i, sqrt(B->GetBinError(i)*B->GetBinError(i) + NoB->GetBinError(i)*NoB->GetBinError(i)) );
	}
	return difference;
}//end difference


void findSpread(Int_t c, Int_t harm,TString spectrum,TString in){
  TFile *f[20];
  TList *list[20];
  TH1D *fHist[20];
  Int_t Nbins = 0;
  if(spectrum == "Pt"){Nbins = 30;}
  else if(spectrum == "Eta"){Nbins = 50;}
  TH1D *fHistPos[20];
  TH1D *fHistNeg[20];
  TH1F *fHistSpreadPos[Nbins];
  TH1F *fHistSpreadNeg[Nbins];
  TH1F *fHistSpreadDiff[Nbins];
 for(Int_t i = 0; i< Nbins; i++){
        fHistSpreadPos[i] = new TH1F(Form("fHistSpreadPos%d"+spectrum+"Bin%d",harm,i+1),"",1000,0,0.3);
        fHistSpreadNeg[i] = new TH1F(Form("fHistSpreadNeg%d"+spectrum+"Bin%d",harm,i+1),"",1000,0,0.3);
        fHistSpreadDiff[i] = new TH1F(Form("fHistSpreadDiff%d"+spectrum+"Bin%d",harm,i+1),"",1000,-0.3,0.3);
  }
  for(Int_t iFile = 0; iFile < 20; iFile++) {
    f[iFile] = TFile::Open(in+Form("/Results_5.02TeV_pTrange_0_eta_0_Cent%d0_%d0_split_%d.root",c,c+1,iFile));
    if((!f[iFile])||(!f[iFile]->IsOpen())) {
    cout<<"File "<<iFile<<" not found..."<<endl;
   // return;
    }

    list[iFile] = dynamic_cast<TList *>(f[iFile]->Get("FlowQCList"));
    if(!list[iFile]) {
      cout<<"Input list of file "<<iFile<<" not found..."<<endl;
    //  return;
    }
    fHistPos[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[0][%d][%d][0]",c,harm)));
    if(!fHistPos[iFile]) {
      cout<<"Histogram of positive particles from file "<<iFile<<" not found..."<<endl;
    //  return;
    }
    fHistNeg[iFile] = dynamic_cast<TH1D *>(list[iFile]->FindObject(Form("fFlowQCFinal"+spectrum+"DifHist[1][%d][%d][0]",c,harm)));
    if(!fHistNeg[iFile]) {
      cout<<"Histogram of negative particles from file "<<iFile<<" not found..."<<endl;
     // return;
    }

    for(Int_t iBin = 1; iBin <= fHistPos[iFile]->GetNbinsX(); iBin++){
	fHistSpreadPos[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin));}
    for(Int_t iBin = 1; iBin <= fHistNeg[iFile]->GetNbinsX(); iBin++){
      fHistSpreadNeg[iBin-1]->Fill(fHistNeg[iFile]->GetBinContent(iBin));
      if(fabs(fHistPos[iFile]->GetBinContent(iBin))>0 && fabs(fHistNeg[iFile]->GetBinContent(iBin))>0){
      fHistSpreadDiff[iBin-1]->Fill(fHistPos[iFile]->GetBinContent(iBin) - fHistNeg[iFile]->GetBinContent(iBin));}
    }//loop to fill SpreadHists

  }//loop over files
  TString B;
  if(in == "/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.4/split"){B = "B";}
  if(in =="/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.0/split"){B = "NoB";}
  TFile *fOutput = new TFile(Form("output_findSpread/spread%d_"+B+spectrum+"_c%d0_%d0.root",harm,c,c+1),"recreate");
  for(Int_t iBin = 1; iBin <= fHistPos[0]->GetNbinsX(); iBin++){fHistSpreadPos[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadNeg[iBin-1]->Write();}
  for(Int_t iBin = 1; iBin <= fHistNeg[0]->GetNbinsX(); iBin++){fHistSpreadDiff[iBin-1]->Write();}
  fOutput->Close();

 for(Int_t i = 0; i< Nbins; i++){
        delete fHistSpreadPos[i];
        delete fHistSpreadNeg[i];
        delete fHistSpreadDiff[i];
 }
for(Int_t iFile = 0; iFile < 20; iFile++) {
	f[iFile]->Close();
 }
 //for(Int_t i = 0; i<20;i++){
//	delete fHist[i];	
 //}
}//end find spread


TH1F *Subsampling(TString spectrum, TString charge, TH1F *result, Int_t c, Int_t harm, TString B){
	Int_t Nbins = result->GetNbinsX();
	TKey *k;
	TH1F *s;
	Double_t sigma = 0.0;
	Double_t mean = 0.0;
   	TString h = Form("fHistSpread"+charge+"%d",harm);//hmm...something is wrong :( 
    	TFile *of = new TFile(Form("output_findSpread/spread%d_"+B+spectrum+"_c%d0_%d0.root",harm,c,c+1),"read");
	if(!of->IsOpen()){cout<<"closed file"<<endl;}
	for(Int_t i=1; i<=Nbins; i++){
	//	cout<<of<<endl;
		s =dynamic_cast<TH1F *>(of->GetKey(Form(h+spectrum+"Bin%d",i))->ReadObj());	
		TH1F *o = (TH1F *)s->Clone();
		sigma = o->GetRMS();
	//	cout<<"sigma "<<sigma<<endl;
		mean = o->GetMean();//take the file check the values fro 3*sigma        
		result->SetBinContent(i, mean);
		result->SetBinError(i, sigma);
	}//end Nbins
	of->Close();
return result;
}//end subsampling

void CompareEffectB(){//add writing to external file to later calcuate the difference for B vs noB
        gROOT->SetBatch();//to avoid opening the plots ak bad wifi struggle
        Double_t fPtDiffNBins = 29;
        Double_t fCRCPtBins[30];
        Double_t PtBins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5.};
        for(Int_t r=0; r<30; r++){fCRCPtBins[r] = PtBins[r];}

        Double_t fEtaDiffNBins = 50;
        Double_t fCRCEtaBins[51]={0};
        Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
        for(Int_t i=0; i<51; i++){fCRCEtaBins[i] = etabinEdge[i];}

	Int_t c = 2;
//	for(Int_t c=1; c<5;c++){
		gROOT->SetBatch();
		gStyle->SetOptStat(0);
		TH1F *DiffPospT, *DiffNegpT, *DiffDelpT,*DiffPospT2, *DiffNegpT2, *DiffDelpT2, *DiffPospT3, *DiffNegpT3, *DiffDelpT3;
		TH1F *DiffPosEta, *DiffNegEta, *DiffDelEta, *DiffPosEta2, *DiffNegEta2, *DiffDelEta2, *DiffPosEta3, *DiffNegEta3, *DiffDelEta3;
        	cout<<"Centrality: "<<c<<"0-"<<(c+1)<<"0"<<endl;//actually 5.02/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.2

        	TString inputB = "/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.4/split";
		TString input = "/data/alice/jlomker/AVFD/result/dirID-0/tau_init0.4/BField0.0/split";//something is wrong here too!
                //cout<<"Harmonic: "<<(harm+1)<<endl;//apply subsampling to all histograms from pos/neg particles=======================

                findSpread(c,0,"Pt",input);//finds spread from the 10 splitfiles for pos, neg and delta histogram as input for the subsampling
		TH1F *PvpT = new TH1F("pvpT", "pvpT", fPtDiffNBins, fCRCPtBins);
		PvpT = Subsampling("Pt", "Pos", PvpT, c, 0,"NoB");
		TH1F *NvpT = new TH1F("nvpT","nvpT",fPtDiffNBins,fCRCPtBins);
                NvpT = Subsampling("Pt","Neg",NvpT, c ,0,"NoB");
                TH1F *DvpT = new TH1F("dvnpT","dvnpT",fPtDiffNBins,fCRCPtBins);
                DvpT = Subsampling("Pt","Diff", DvpT,c,0,"NoB");
		findSpread(c,0,"Pt",inputB);
                TH1F *BPvpT = new TH1F("BpvpT", "BpvpT", fPtDiffNBins, fCRCPtBins);
                BPvpT = Subsampling("Pt", "Pos", BPvpT, c, 0,"B");
                TH1F *BNvpT = new TH1F("BnvpT","BnvpT",fPtDiffNBins,fCRCPtBins);
                BNvpT = Subsampling("Pt","Neg",BNvpT, c ,0, "B");
                TH1F *BDvpT = new TH1F("BdvnpT","BdvnpT",fPtDiffNBins,fCRCPtBins);
                BDvpT = Subsampling("Pt","Diff", BDvpT,c,0, "B");

		findSpread(c,1,"Pt",input);
                TH1F *PvpT2 = new TH1F("pvpT2", "pvpT2", fPtDiffNBins, fCRCPtBins);
                PvpT2 = Subsampling("Pt", "Pos", PvpT2, c, 1,"NoB");
                TH1F *NvpT2 = new TH1F("nvpT2","nvpT2",fPtDiffNBins,fCRCPtBins);
                NvpT2 = Subsampling("Pt","Neg",NvpT2, c ,1,"NoB");
                TH1F *DvpT2 = new TH1F("dvnpT2","dvnpT2",fPtDiffNBins,fCRCPtBins);
                DvpT2 = Subsampling("Pt","Diff", DvpT2,c,1,"NoB");
                findSpread(c,1,"Pt",inputB);
                TH1F *BPvpT2 = new TH1F("BpvpT2", "BpvpT2", fPtDiffNBins, fCRCPtBins);
                BPvpT2 = Subsampling("Pt", "Pos", BPvpT2, c, 1,"B");
                TH1F *BNvpT2 = new TH1F("BnvpT2","BnvpT2",fPtDiffNBins,fCRCPtBins);
                BNvpT2 = Subsampling("Pt","Neg",BNvpT2, c ,1,"B");
                TH1F *BDvpT2 = new TH1F("BdvnpT2","BdvnpT2",fPtDiffNBins,fCRCPtBins);
                BDvpT2 = Subsampling("Pt","Diff", BDvpT2,c,1,"B");

		findSpread(c,2,"Pt",input);
                TH1F *PvpT3 = new TH1F("pvpT3", "pvpT3", fPtDiffNBins, fCRCPtBins);
                PvpT3 = Subsampling("Pt", "Pos", PvpT3, c, 2,"NoB");
                TH1F *NvpT3 = new TH1F("nvpT3","nvpT3",fPtDiffNBins,fCRCPtBins);
                NvpT3 = Subsampling("Pt","Neg",NvpT3, c ,2,"NoB");
                TH1F *DvpT3 = new TH1F("dvnpT3","dvnpT3",fPtDiffNBins,fCRCPtBins);
                DvpT3 = Subsampling("Pt","Diff", DvpT3,c,2,"NoB");
                findSpread(c,2,"Pt",inputB);
                TH1F *BPvpT3 = new TH1F("BpvpT3", "BpvpT3", fPtDiffNBins, fCRCPtBins);
                BPvpT3 = Subsampling("Pt", "Pos", BPvpT3, c, 2,"B");
                TH1F *BNvpT3 = new TH1F("BnvpT3","BnvpT3",fPtDiffNBins,fCRCPtBins);
                BNvpT3 = Subsampling("Pt","Neg",BNvpT3, c ,2,"B");
                TH1F *BDvpT3 = new TH1F("BdvnpT3","BdvnpT3",fPtDiffNBins,fCRCPtBins);
                BDvpT3 = Subsampling("Pt","Diff", BDvpT3,c,2,"B");

		//make	v1,2,3 for one centrality with low centrality (better statistics) per panel make same color as before but now centrality color -> eta/pT, marker the same, vn color = vn color sheme 	
		DiffPospT = Difference(BPvpT,PvpT,fCRCPtBins,"PosPt");
	        DiffNegpT = Difference(BNvpT,NvpT,fCRCPtBins,"NegPt");
		DiffDelpT = Difference(BDvpT,DvpT,fCRCPtBins,"DelPt");
                DiffPospT2 = Difference(BPvpT2,PvpT2,fCRCPtBins,"PosPt2");
                DiffNegpT2 = Difference(BNvpT2,NvpT2,fCRCPtBins,"NegPt2");
                DiffDelpT2 = Difference(BDvpT2,DvpT2,fCRCPtBins,"DelPt2");
                DiffPospT3 = Difference(BPvpT3,PvpT3,fCRCPtBins,"PosPt3");
                DiffNegpT3 = Difference(BNvpT3,NvpT3,fCRCPtBins,"NegPt3");
                DiffDelpT3 = Difference(BDvpT3,DvpT3,fCRCPtBins,"DelPt3");

                auto L = new TLegend(0.15,0.88,0.9,1);
                L->SetTextSize(0.031);
                L->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 20%-30%, |#eta| #leq 0.8","C");
                L->SetBorderSize(0);
                L->SetFillStyle(0);

                auto L0 = new TLegend(0.32,0.81,0.92,0.91);
                L0->SetTextSize(0.03);
                L0->SetHeader("#bf{#tau_{0} = 0.4 fm/c with #tau_{B} = 0.4 fm/c and #tau_{B'} = 0 fm/c} ","C");
                L0->SetBorderSize(0);
                L0->SetFillStyle(0);


        	auto L1 = new TLegend(0.32,0.84,0.89,0.92);
        	L1->SetTextSize(0.02);//L2->SetTextAlign(11);//or13
        	L1->SetEntrySeparation(0.0006);
        	L1->SetNColumns(3); 
	        L1->AddEntry(DiffPospT,Form("#Delta_{#tau} #nu_{1}(+h)"), "lep");
        	L1->AddEntry(DiffPospT2,Form("#Delta_{#tau} #nu_{2}(+h)"), "lep");
		L1->AddEntry(DiffPospT3,Form("#Delta_{#tau} #nu_{3}(+h)"), "lep");
		L1->SetBorderSize(0);
	        L1->SetFillStyle(0);
                auto L2 = new TLegend(0.32,0.52,0.89,0.67);
                L2->SetTextSize(0.035);//L2->SetTextAlign(11);//or13
                L2->SetEntrySeparation(0.0006);
                L2->SetNColumns(3);
                L2->AddEntry(DiffDelpT,Form("directed"), "lep");
                L2->AddEntry(DiffDelpT2,Form("elliptic"), "lep");
                L2->AddEntry(DiffDelpT3,Form("triangular flow"), "lep");
                L2->SetBorderSize(0);
                L2->SetFillStyle(0);
                auto L3 = new TLegend(0.32,0.28,0.89,0.36);
                L3->SetTextSize(0.02);//L2->SetTextAlign(11);//or13
                L3->SetEntrySeparation(0.0006);
                L3->SetNColumns(3);
                L3->AddEntry(DiffDelpT,Form("directed"), "lep");
                L3->AddEntry(DiffDelpT2,Form("elliptic"), "lep");
                L3->AddEntry(DiffDelpT3,Form("triangular flow"), "lep");
                L3->SetBorderSize(0);
                L3->SetFillStyle(0);

		TCanvas *p1 = new TCanvas("p1","PtDifferencesB",160,200);
		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.64, 1., 0.98);
        	pad1->SetBottomMargin(0.0001); // Upper and lower plot are joined
        	pad1->SetLeftMargin(0.2);
		pad1->SetRightMargin(0.02);
		pad1->SetTopMargin(0.25);
        	pad1->Draw();             // Draw the upper pad: pad1
        	pad1->cd();
		DiffPospT->SetTitle("");
                DiffPospT->GetYaxis()->SetTitle("#Delta_{#tau}#nu_{n}(+h)");
                DiffPospT->GetYaxis()->SetTitleSize(0.13);
                DiffPospT->GetYaxis()->SetTitleOffset(0.6);
                DiffPospT->GetYaxis()->SetLabelSize(0.11);
		DiffPospT->GetYaxis()->SetRangeUser(-0.012,0.042);
                DiffPospT->GetYaxis()->CenterTitle();
		DiffPospT->GetXaxis()->SetLabelOffset(0.4);
                DiffPospT->GetXaxis()->SetRangeUser(0.2,3);
                DiffPospT->GetYaxis()->SetNdivisions(4);
		DiffPospT->SetLineColor(kOrange+2);
		DiffPospT->SetMarkerStyle(24);
		DiffPospT->SetMarkerSize(0.2);
		DiffPospT->SetMarkerColor(kOrange+2);
		DiffPospT->DrawCopy("h");
		DiffPospT2->SetLineColor(kMagenta+2);
		DiffPospT2->SetMarkerStyle(27);
		DiffPospT2->SetMarkerSize(0.2);
		DiffPospT2->SetMarkerColor(kMagenta+2);
		DiffPospT2->DrawCopy("hsame");
		DiffPospT3->SetLineColor(kGreen+2);
		DiffPospT3->SetMarkerStyle(26);
		DiffPospT3->SetMarkerSize(0.2);
		DiffPospT3->SetMarkerColor(kGreen+2);
		DiffPospT3->DrawCopy("hsame");
		p1->cd();	
		//L1->Draw();
		L0->Draw();
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.35, 1., 0.64);
        	pad2->SetTopMargin(0);
		pad2->SetRightMargin(0.02);
		pad2->SetBottomMargin(0.0001); // Upper and lower plot are joined
        	pad2->SetLeftMargin(0.2);
        	pad2->Draw();             // Draw the upper pad: pad1
        	pad2->cd();
                DiffNegpT->GetYaxis()->SetTitle("#Delta_{#tau} #nu_{n}(-h)");
                DiffNegpT->GetYaxis()->SetTitleSize(0.15);
                DiffNegpT->GetYaxis()->SetTitleOffset(0.55);
                DiffNegpT->GetYaxis()->SetLabelSize(0.13);
		DiffNegpT->GetYaxis()->SetRangeUser(-0.012,0.044);
		DiffNegpT->GetYaxis()->CenterTitle();
		DiffNegpT->GetXaxis()->SetRangeUser(0.2,3);
		DiffNegpT->GetXaxis()->SetLabelOffset(0.1);
		DiffNegpT->GetYaxis()->SetNdivisions(4);
                DiffNegpT->SetLineColor(kOrange+2);
		DiffNegpT->SetMarkerStyle(24);
		DiffNegpT->SetMarkerSize(0.2);
		DiffNegpT->SetMarkerColor(kOrange+2);
                DiffNegpT->Draw("h");
                DiffNegpT2->SetLineColor(kMagenta+2);
		DiffNegpT2->SetMarkerStyle(27);
		DiffNegpT2->SetMarkerSize(0.2);
		DiffNegpT2->SetMarkerColor(kMagenta+2);
                DiffNegpT2->Draw("hsame");
                DiffNegpT3->SetLineColor(kGreen+2);
		DiffNegpT3->SetMarkerStyle(26);
		DiffNegpT3->SetMarkerSize(0.2);
		DiffNegpT3->SetMarkerColor(kGreen+2);
		DiffNegpT3->Draw("hsame");
		p1->cd();
		L2->Draw();	
		TPad *pad3 = new TPad("pad3", "pad3", 0, 0.03, 1., 0.35);
		pad3->SetTopMargin(0);
		pad3->SetRightMargin(0.02);
        	pad3->SetBottomMargin(0.25); // Upper and lower plot are joined
        	pad3->SetLeftMargin(0.2);
        	pad3->Draw();             // Draw the upper pad: pad1
        	pad3->cd();
                DiffDelpT->GetYaxis()->SetTitle("#Delta_{#tau} #Delta #nu_{n}");
		DiffDelpT->GetYaxis()->SetRangeUser(-0.022,0.042);
                DiffDelpT->GetYaxis()->SetTitleSize(0.14);
		DiffDelpT->GetYaxis()->SetTitleOffset(0.53);
		DiffDelpT->GetXaxis()->SetTitle("p_{T} [GeV]");
		DiffDelpT->GetYaxis()->CenterTitle();
		DiffDelpT->GetXaxis()->SetTitleOffset(0.9);
		DiffDelpT->GetXaxis()->SetTitleSize(0.12);
                DiffDelpT->GetYaxis()->SetLabelSize(0.11);
		DiffDelpT->GetXaxis()->SetLabelSize(0.11);
                DiffDelpT->GetXaxis()->SetRangeUser(0.2,3);
                DiffDelpT->GetYaxis()->SetNdivisions(6);
                DiffDelpT->SetLineColor(kOrange+2);
		DiffDelpT->SetMarkerStyle(20);
		DiffDelpT->SetMarkerSize(0.2);
		DiffDelpT->SetMarkerColor(kOrange+2);
                DiffDelpT->DrawCopy("h");
                DiffDelpT2->SetLineColor(kMagenta+2);
		DiffDelpT2->SetMarkerStyle(33);
		DiffDelpT2->SetMarkerSize(0.2);
		DiffDelpT2->SetMarkerColor(kMagenta+2);
                DiffDelpT2->DrawCopy("hsame");
                DiffDelpT3->SetLineColor(kGreen+2);
		DiffDelpT3->SetMarkerStyle(22);
		DiffDelpT3->SetMarkerSize(0.2);
		DiffDelpT3->SetMarkerColor(kGreen+2);
		DiffDelpT3->DrawCopy("hsame");
		p1->cd();
		//L3->Draw();
		L->Draw();
		p1->SaveAs(Form("difference/thesis_pT_cen%d.pdf",c));
		
		DiffPospT->Delete();
		DiffNegpT->Delete();
		DiffDelpT->Delete();
                DiffPospT2->Delete();
                DiffNegpT2->Delete();
                DiffDelpT2->Delete();
                DiffPospT3->Delete();
                DiffNegpT3->Delete();
                DiffDelpT3->Delete();
/*                PvpT->Delete();
                NvpT->Delete();
                DvpT->Delete();
                BPvpT->Delete();
                BNvpT->Delete();
                BDvpT->Delete();
                PvpT2->Delete();
                NvpT2->Delete();
                DvpT2->Delete();
                BPvpT2->Delete();
                BNvpT2->Delete();
                BDvpT2->Delete();
                PvpT3->Delete();
                NvpT3->Delete();
                DvpT3->Delete();
                BPvpT3->Delete();
                BNvpT3->Delete();
                BDvpT3->Delete();*/
		delete p1;
                
		findSpread(c,0,"Eta",input);
                TH1F *PvEta = new TH1F("pvEta", "pvEta", fEtaDiffNBins, fCRCEtaBins);
                PvEta = Subsampling("Eta", "Pos", PvEta, c, 0,"NoB");
                TH1F *NvEta = new TH1F("nvEta", "nvEta", fEtaDiffNBins, fCRCEtaBins);
                NvEta = Subsampling("Eta", "Neg", NvEta, c, 0, "NoB");
                TH1F *DvEta = new TH1F("deltavneta","deltavneta",fEtaDiffNBins,fCRCEtaBins);
                DvEta = Subsampling( "Eta", "Diff", DvEta, c, 0,"NoB");
                findSpread(c,0,"Eta",inputB);
                TH1F *BPvEta = new TH1F("BpvEta", "BpvEta", fEtaDiffNBins, fCRCEtaBins);
                BPvEta = Subsampling("Eta", "Pos", BPvEta, c, 0,"B");
                TH1F *BNvEta = new TH1F("BnvEta", "BnvEta", fEtaDiffNBins, fCRCEtaBins);
                BNvEta = Subsampling("Eta", "Neg", BNvEta, c, 0, "B");
                TH1F *BDvEta = new TH1F("Bdeltavneta","Bdeltavneta",fEtaDiffNBins,fCRCEtaBins);
                BDvEta = Subsampling( "Eta", "Diff", BDvEta, c, 0,"B");

                findSpread(c,1,"Eta",input);
                TH1F *PvEta2 = new TH1F("pvEta2", "pvEta2", fEtaDiffNBins, fCRCEtaBins);
                PvEta2 = Subsampling("Eta", "Pos", PvEta2, c,1, "NoB");
                TH1F *NvEta2 = new TH1F("nvEta2", "nvEta2", fEtaDiffNBins, fCRCEtaBins);
                NvEta2 = Subsampling("Eta", "Neg", NvEta2, c, 1, "NoB");
                TH1F *DvEta2 = new TH1F("deltavneta2","deltavneta2",fEtaDiffNBins,fCRCEtaBins);
                DvEta2 = Subsampling( "Eta", "Diff", DvEta2, c, 1, "NoB");
                findSpread(c,1,"Eta",inputB);
                TH1F *BPvEta2 = new TH1F("BpvEta2", "BpvEta2", fEtaDiffNBins, fCRCEtaBins);
                BPvEta2 = Subsampling("Eta", "Pos", BPvEta2, c, 1, "B");
                TH1F *BNvEta2 = new TH1F("BnvEta2", "BnvEta2", fEtaDiffNBins, fCRCEtaBins);
                BNvEta2 = Subsampling("Eta", "Neg", BNvEta2, c, 1, "B");
                TH1F *BDvEta2 = new TH1F("Bdeltavneta2","Bdeltavneta2",fEtaDiffNBins,fCRCEtaBins);
                BDvEta2 = Subsampling( "Eta", "Diff", BDvEta2, c, 1, "B");

                findSpread(c,2,"Eta",input);
                TH1F *PvEta3 = new TH1F("pvEta3", "pvEta3", fEtaDiffNBins, fCRCEtaBins);
                PvEta3 = Subsampling("Eta", "Pos", PvEta3, c, 2, "NoB");
                TH1F *NvEta3 = new TH1F("nvEta3", "nvEta3", fEtaDiffNBins, fCRCEtaBins);
                NvEta3 = Subsampling("Eta", "Neg", NvEta3, c, 2, "NoB");
                TH1F *DvEta3 = new TH1F("deltavneta3","deltavneta3",fEtaDiffNBins,fCRCEtaBins);
                DvEta3 = Subsampling( "Eta", "Diff", DvEta3, c, 2, "NoB");
                findSpread(c,2,"Eta",inputB);
                TH1F *BPvEta3 = new TH1F("BpvEta3", "BpvEta3", fEtaDiffNBins, fCRCEtaBins);
                BPvEta3 = Subsampling("Eta", "Pos", BPvEta3, c, 2, "B");
                TH1F *BNvEta3 = new TH1F("BnvEta3", "BnvEta3", fEtaDiffNBins, fCRCEtaBins);
                BNvEta3 = Subsampling("Eta", "Neg", BNvEta3, c, 2, "B");
                TH1F *BDvEta3 = new TH1F("Bdeltavneta3","Bdeltavneta3",fEtaDiffNBins,fCRCEtaBins);
                BDvEta3 = Subsampling( "Eta", "Diff", BDvEta3, c, 2, "B");

		DiffPosEta = Difference(BPvEta,PvEta,fCRCEtaBins,"PosEta");//difference = result with B - result without B
		DiffNegEta = Difference(BNvEta,NvEta,fCRCEtaBins,"NegEta");
		DiffDelEta = Difference(BDvEta,DvEta,fCRCEtaBins,"DelEta");
                DiffPosEta2 = Difference(BPvEta2,PvEta2,fCRCEtaBins,"PosEta2");//difference = result with B - result without B
                DiffNegEta2 = Difference(BNvEta2,NvEta2,fCRCEtaBins,"NegEta2");
                DiffDelEta2 = Difference(BDvEta2,DvEta2,fCRCEtaBins,"DelEta2");
                DiffPosEta3 = Difference(BPvEta3,PvEta3,fCRCEtaBins,"PosEta3");//difference = result with B - result without B
                DiffNegEta3 = Difference(BNvEta3,NvEta3,fCRCEtaBins,"NegEta3");
                DiffDelEta3 = Difference(BDvEta3,DvEta3,fCRCEtaBins,"DelEta3");

                auto LE = new TLegend(0.15,0.88,0.9,1);
                LE->SetTextSize(0.031);
                LE->SetHeader("#bf{AVFD Pb-Pb} @ #sqrt{s} = 5.02TeV, Cent. 20%-30%, p_{T} #in [0.2,5]GeV","C");
		LE->SetBorderSize(0);
                LE->SetFillStyle(0);

		auto LE1 = new TLegend(0.32,0.79,0.92,0.94);
		LE1->SetTextSize(0.03);
		LE1->SetHeader("#bf{#tau_{0} = 0.4 fm/c with #tau_{B} = 0.4 fm/c and #tau_{B'} = 0 fm/c} ","C");
                LE1->SetBorderSize(0);
                LE1->SetFillStyle(0);

		auto L1e = new TLegend(0.32,0.84,0.89,0.92);
                //L1e->SetTextSize(0.03);
                //L1e->SetHeader("#bf{Pb-Pb} @ #sqrt{s} = 5.02TeV in centrality 20%-30%, #tau = 0.4 fm/c, pT #in {0.2,5} GeV, |#eta| #leq 0.8 ","C");
                L1e->SetTextSize(0.02);//L2->SetTextAlign(11);//or13
                L1e->SetEntrySeparation(0.0006);
                L1e->SetNColumns(3);
                L1e->AddEntry(DiffPosEta,Form("directed"), "l");
                //L1e->AddEntry(DiffPosEta2,Form("#Delta_{#tau} #nu_{2}(+h)"), "lep");
                //L1e->AddEntry(DiffPosEta3,Form("#Delta_{#tau} #nu_{3}(+h)"), "lep");
                L1e->SetBorderSize(0);
                L1e->SetFillStyle(0);
                auto L2e = new TLegend(0.32,0.6,0.89,0.75);
                L2e->SetTextSize(0.035);//L2->SetTextAlign(11);//or13
                L2e->SetEntrySeparation(0.0006);
                L2e->SetNColumns(3);
                L2e->AddEntry(DiffDelEta,Form("directed"), "lep");
                L2e->AddEntry(DiffDelEta2,Form("elliptic"), "lep");
                L2e->AddEntry(DiffDelEta3,Form("triangular flow"), "lep");
                L2e->SetBorderSize(0);
                L2e->SetFillStyle(0);
                auto L3e = new TLegend(0.32,0.28,0.89,0.36);
                L3e->SetTextSize(0.02);//L2->SetTextAlign(11);//or13
                L3e->SetEntrySeparation(0.0006);
                //L3e->SetNColumns(3);
                //L3e->AddEntry(DiffDelEta,Form("#Delta_{#tau} #Delta #nu_{1}"), "lep");
                //L3e->AddEntry(DiffDelEta2,Form("#Delta_{#tau} #Delta #nu_{2}"), "lep");
                L3e->AddEntry(DiffDelEta3,Form("triangular"), "l");
                L3e->SetBorderSize(0);
                L3e->SetFillStyle(0);


                TCanvas *e1 = new TCanvas("e1","PtDifferencesB",160,200);
                TPad *pad1e = new TPad("pad1e", "pad1e", 0, 0.64, 1., 0.98);
                pad1e->SetBottomMargin(0.0001); // Upper and lower plot are joined
                pad1e->SetLeftMargin(0.2);
                pad1e->SetRightMargin(0.02);
                pad1e->SetTopMargin(0.25);
                pad1e->Draw();             // Draw the upper pad: pad1
                pad1e->cd();
                //DiffPosEta->SetTitle(Form("AVFD with #tau_{B} = 0.4 fm/c and #tau_{B'} = 0.0 fm/c"));
                //DiffPosEta->SetTitleSize(0.9);
                DiffPosEta->GetYaxis()->SetTitle("#Delta_{#tau}#nu_{n}(+h)");
                DiffPosEta->GetYaxis()->SetTitleSize(0.13);
                DiffPosEta->GetYaxis()->SetTitleOffset(0.6);
                DiffPosEta->GetYaxis()->SetLabelSize(0.11);
                DiffPosEta->GetYaxis()->SetRangeUser(-0.008,0.026);
                DiffPosEta->GetXaxis()->SetLabelOffset(0.1);
                DiffPosEta->GetXaxis()->SetRangeUser(-0.8,0.8);
                DiffPosEta->GetYaxis()->SetNdivisions(4);
		DiffPosEta->SetLineColor(kOrange+2);
                DiffPosEta->SetMarkerStyle(24);
                DiffPosEta->SetMarkerSize(0.2);
                DiffPosEta->SetMarkerColor(kOrange+2);
                DiffPosEta->DrawCopy("h");
                DiffPosEta2->SetLineColor(kMagenta+2);
                DiffPosEta2->SetMarkerStyle(27);
                DiffPosEta2->SetMarkerSize(0.2);
                DiffPosEta2->SetMarkerColor(kMagenta+2);
                DiffPosEta2->DrawCopy("hsame");
                DiffPosEta3->SetLineColor(kGreen+2);
                DiffPosEta3->SetMarkerStyle(26);
                DiffPosEta3->SetMarkerSize(0.2);
                DiffPosEta3->SetMarkerColor(kGreen+2);
                DiffPosEta3->DrawCopy("hsame");
                e1->cd();
                //L1e->Draw();
		LE1->Draw();
                TPad *pad2e = new TPad("pad2e", "pad2e", 0, 0.35, 1., 0.64);
                pad2e->SetTopMargin(0);
                pad2e->SetRightMargin(0.02);
                pad2e->SetBottomMargin(0.0001); // Upper and lower plot are joined
                pad2e->SetLeftMargin(0.2);
                pad2e->Draw();             // Draw the upper pad: pad1
                pad2e->cd();
                DiffNegEta->GetYaxis()->SetTitle("#Delta_{#tau} #nu_{n}(-h)");
                DiffNegEta->GetYaxis()->SetTitleSize(0.15);
                DiffNegEta->GetYaxis()->SetTitleOffset(0.55);
                DiffNegEta->GetYaxis()->SetLabelSize(0.13);
                DiffNegEta->GetYaxis()->CenterTitle();
		DiffNegEta->GetYaxis()->SetRangeUser(-0.006,0.022);
                DiffNegEta->GetXaxis()->SetRangeUser(-0.8,0.8);
                DiffNegEta->GetXaxis()->SetLabelOffset(0.1);
                DiffNegEta->GetYaxis()->SetNdivisions(4);
                DiffNegEta->SetLineColor(kOrange+2);
                DiffNegEta->SetMarkerStyle(24);
                DiffNegEta->SetMarkerSize(0.2);
                DiffNegEta->SetMarkerColor(kOrange+2);
                DiffNegEta->Draw("h");
                DiffNegEta2->SetLineColor(kMagenta+2);
                DiffNegEta2->SetMarkerStyle(27);
                DiffNegEta2->SetMarkerSize(0.2);
                DiffNegEta2->SetMarkerColor(kMagenta+2);
                DiffNegEta2->Draw("hsame");
                DiffNegEta3->SetLineColor(kGreen+2);
                DiffNegEta3->SetMarkerStyle(26);
                DiffNegEta3->SetMarkerSize(0.2);
                DiffNegEta3->SetMarkerColor(kGreen+2);
                DiffNegEta3->Draw("hsame");
                e1->cd();
                L2e->Draw();
                TPad *pad3e = new TPad("pad3e", "pad3e", 0, 0.03, 1., 0.35);
                pad3e->SetTopMargin(0);
                pad3e->SetRightMargin(0.02);
                pad3e->SetBottomMargin(0.25); // Upper and lower plot are joined
                pad3e->SetLeftMargin(0.2);
                pad3e->Draw();             // Draw the upper pad: pad1
                pad3e->cd();
                DiffDelEta->GetYaxis()->SetTitle("#Delta_{#tau} #Delta #nu_{n}");//#tau_{B}) - #Delta #nu_{n}(#tau_{B'})");
                DiffDelEta->GetYaxis()->SetRangeUser(-0.012,0.022);
                DiffDelEta->GetYaxis()->SetTitleSize(0.14);
		DiffDelEta->GetYaxis()->CenterTitle();
                DiffDelEta->GetYaxis()->SetTitleOffset(0.53);
                DiffDelEta->GetXaxis()->SetTitle("#eta");
                DiffDelEta->GetXaxis()->SetTitleOffset(0.8);
                DiffDelEta->GetXaxis()->SetTitleSize(0.12);
                DiffDelEta->GetYaxis()->SetLabelSize(0.11);
                DiffDelEta->GetXaxis()->SetLabelSize(0.11);
                DiffDelEta->GetXaxis()->SetRangeUser(-0.8,0.8);
		DiffDelEta->GetYaxis()->SetNdivisions(4);
                DiffDelEta->SetLineColor(kOrange+2);
                DiffDelEta->SetMarkerStyle(20);
                DiffDelEta->SetMarkerSize(0.2);
                DiffDelEta->SetMarkerColor(kOrange+2);
                DiffDelEta->DrawCopy("h");
                DiffDelEta2->SetLineColor(kMagenta+2);
                DiffDelEta2->SetMarkerStyle(33);
                DiffDelEta2->SetMarkerSize(0.2);
                DiffDelEta2->SetMarkerColor(kMagenta+2);
                DiffDelEta2->DrawCopy("hsame");
                DiffDelEta3->SetLineColor(kGreen+2);
                DiffDelEta3->SetMarkerStyle(22);
                DiffDelEta3->SetMarkerSize(0.2);
                DiffDelEta3->SetMarkerColor(kGreen+2);
                DiffDelEta3->DrawCopy("hsame");
                e1->cd();
		//LE1->Draw();
                //L3e->Draw();
		LE->Draw();
		e1->SaveAs(Form("difference/thesis_eta_cen%d.pdf",c));

		DiffPosEta->Delete();
	        DiffNegEta->Delete();
		DiffDelEta->Delete();
                DiffPosEta2->Delete();
                DiffNegEta2->Delete();
                DiffDelEta2->Delete();
                DiffPosEta3->Delete();
                DiffNegEta3->Delete();
                DiffDelEta3->Delete();
/*
                PvEta->Delete();
                NvEta->Delete();
                DvEta->Delete();
                BPvEta->Delete();
                BNvEta->Delete();
                BDvEta->Delete();
                PvEta2->Delete();
                NvEta2->Delete();
                DvEta2->Delete();
                BPvEta2->Delete();
                BNvEta2->Delete();
                BDvEta2->Delete();
                PvEta3->Delete();
                NvEta3->Delete();
                DvEta3->Delete();
                BPvEta3->Delete();
                BNvEta3->Delete();
                BDvEta3->Delete();
*/
		delete e1;
//	}//centrality

}//endvoid

