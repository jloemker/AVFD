void drawCorrelatorsFromAVFD() {
  TFile* fFitPb5 = TFile::Open("DataPoints/fitPb5.root", "READ");
  if((!fFitPb5)||(!fFitPb5->IsOpen())) {
    cout<<"The fit parameters for the Pb system are not found..."<<endl;
    return;
  }
  
  TFile* fFitXe = TFile::Open("DataPoints/fitXe.root", "READ");
  if((!fFitXe)||(!fFitXe->IsOpen())) {
    cout<<"The fit parameters for the Xe system are not found..."<<endl;
    return;
  }
 //getting my histograms /data/alice/jlomker/AVFD/result/AnalysisResult_a-0.1_5.44TeV.root
 //For test with reduced sample TFile* File = TFile::Open("/project/alice/users/jlomker/AnalysisResult_a-0.1_5.44TeV.root", "READ");
  TFile* File = TFile::Open("/data/alice/jlomker/AVFD/result/AnalysisResult_a-0.1_5.44TeV.root", "READ");
  if((!File)||(!File->IsOpen())){
  cout<<"The histograms for my comparisons cannot be found..."<<endl;
  return;
  }
  TList *CME = (TList*) File->Get("CMEList;1");
  TH1F* DeltaG112Xe_a_010 = (TH1F*) CME -> FindObject("DeltaG112");
  TH1F* DeltaD11Xe_a_010 = (TH1F*) CME -> FindObject("DeltaD11");
 //declaring arrays for parameter fit
  TF1* fDeltaLCCPb5[7]; // [0]*x + [1]
  TF1* fDeltaCMEPb5[7]; // [0]*x*x + [1]*x + [2]
  TF1* fGammaLCCPb5[7]; // [0]*x + [1]
  TF1* fGammaCMEPb5[7]; // [0]*x*x + [1]*x + [2]
  
  TF1* fDeltaLCCXe[7];
  TF1* fDeltaCMEXe[7];
  TF1* fGammaLCCXe[7];
  TF1* fGammaCMEXe[7];
  
  Double_t fDeltaLCCPb5Par[7][2];
  Double_t fDeltaLCCPb5ParErr[7][2];
  Double_t fDeltaCMEPb5Par[7][3];
  Double_t fDeltaCMEPb5ParErr[7][3];
  Double_t fGammaLCCPb5Par[7][2];
  Double_t fGammaLCCPb5ParErr[7][2];
  Double_t fGammaCMEPb5Par[7][3];
  Double_t fGammaCMEPb5ParErr[7][3];
  
  Double_t fDeltaLCCXePar[7][2];
  Double_t fDeltaLCCXeParErr[7][2];
  Double_t fDeltaCMEXePar[7][3];
  Double_t fDeltaCMEXeParErr[7][3];
  Double_t fGammaLCCXePar[7][2];
  Double_t fGammaLCCXeParErr[7][2];
  Double_t fGammaCMEXePar[7][3];
  Double_t fGammaCMEXeParErr[7][3];
  
  Double_t fDeltaD11Pb5_LCC15[7];
  Double_t fDeltaG112Pb5_LCC15[7];
  Double_t fDeltaD11Pb5_LCC33[7];
  Double_t fDeltaG112Pb5_LCC33[7];
  Double_t fDeltaD11Pb5_LCC50[7];
  Double_t fDeltaG112Pb5_LCC50[7];
  
  Double_t fDeltaD11Pb5_LCC15Err[7];
  Double_t fDeltaG112Pb5_LCC15Err[7];
  Double_t fDeltaD11Pb5_LCC33Err[7];
  Double_t fDeltaG112Pb5_LCC33Err[7];
  Double_t fDeltaD11Pb5_LCC50Err[7];
  Double_t fDeltaG112Pb5_LCC50Err[7];
  
  Double_t fDeltaD11Xe_LCC15[7];
  Double_t fDeltaG112Xe_LCC15[7];
  Double_t fDeltaD11Xe_LCC33[7];
  Double_t fDeltaG112Xe_LCC33[7];
  Double_t fDeltaD11Xe_LCC50[7];
  Double_t fDeltaG112Xe_LCC50[7];
  
  Double_t fDeltaD11Xe_LCC15Err[7];
  Double_t fDeltaG112Xe_LCC15Err[7];
  Double_t fDeltaD11Xe_LCC33Err[7];
  Double_t fDeltaG112Xe_LCC33Err[7];
  Double_t fDeltaD11Xe_LCC50Err[7];
  Double_t fDeltaG112Xe_LCC50Err[7];
  
  for(Int_t i = 0; i < 7; i++) { // 0 to 70% centrality
    fDeltaLCCPb5[i] = (TF1*) fFitPb5->FindObjectAny(Form("fDeltaLCC%dTo%d", i*10, (i+1)*10));
    fDeltaCMEPb5[i] = (TF1*) fFitPb5->FindObjectAny(Form("fDeltaCME%dTo%d", i*10, (i+1)*10));
    fGammaLCCPb5[i] = (TF1*) fFitPb5->FindObjectAny(Form("fGammaLCC%dTo%d", i*10, (i+1)*10));
    fGammaCMEPb5[i] = (TF1*) fFitPb5->FindObjectAny(Form("fGammaCME%dTo%d", i*10, (i+1)*10));
    
    fDeltaLCCXe[i] = (TF1*) fFitXe->FindObjectAny(Form("fDeltaLCC%dTo%d", i*10, (i+1)*10));
    fDeltaCMEXe[i] = (TF1*) fFitXe->FindObjectAny(Form("fDeltaCME%dTo%d", i*10, (i+1)*10));
    fGammaLCCXe[i] = (TF1*) fFitXe->FindObjectAny(Form("fGammaLCC%dTo%d", i*10, (i+1)*10));
    fGammaCMEXe[i] = (TF1*) fFitXe->FindObjectAny(Form("fGammaCME%dTo%d", i*10, (i+1)*10));
    
    for(Int_t j = 0; j < 2; j++) {
      fDeltaLCCPb5Par[i][j] = fDeltaLCCPb5[i]->GetParameter(j);
      fDeltaLCCPb5ParErr[i][j] = fDeltaLCCPb5[i]->GetParError(j);
      fGammaLCCPb5Par[i][j] = fGammaLCCPb5[i]->GetParameter(j);
      fGammaLCCPb5ParErr[i][j] = fGammaLCCPb5[i]->GetParError(j);
      
      fDeltaLCCXePar[i][j] = fDeltaLCCXe[i]->GetParameter(j);
      fDeltaLCCXeParErr[i][j] = fDeltaLCCXe[i]->GetParError(j);
      fGammaLCCXePar[i][j] = fGammaLCCXe[i]->GetParameter(j);
      fGammaLCCXeParErr[i][j] = fGammaLCCXe[i]->GetParError(j);
    }
    
    for(Int_t j = 0; j < 3; j++) {
      fDeltaCMEPb5Par[i][j] = fDeltaCMEPb5[i]->GetParameter(j);
      fDeltaCMEPb5ParErr[i][j] = fDeltaCMEPb5[i]->GetParError(j);
      fGammaCMEPb5Par[i][j] = fGammaCMEPb5[i]->GetParameter(j);
      fGammaCMEPb5ParErr[i][j] = fGammaCMEPb5[i]->GetParError(j);
      
      fDeltaCMEXePar[i][j] = fDeltaCMEXe[i]->GetParameter(j);
      fDeltaCMEXeParErr[i][j] = fDeltaCMEXe[i]->GetParError(j);
      fGammaCMEXePar[i][j] = fGammaCMEXe[i]->GetParameter(j);
      fGammaCMEXeParErr[i][j] = fGammaCMEXe[i]->GetParError(j);
    }
    
    if (i>=2) { // loaded for i>=2 for Pb5
      //cout<<"==> fDeltaLCCPb5Par[i][0]*0.15+fDeltaLCCPb5Par[i][1] = "<<fDeltaLCCPb5Par[i][0]*0.15<<"+"<<fDeltaLCCPb5Par[i][1]<<endl;
      fDeltaD11Pb5_LCC15[i] = fDeltaLCCPb5Par[i][0]*15+fDeltaLCCPb5Par[i][1];
      fDeltaG112Pb5_LCC15[i] = fGammaLCCPb5Par[i][0]*15+fGammaLCCPb5Par[i][1];
      fDeltaD11Pb5_LCC33[i] = fDeltaLCCPb5Par[i][0]*33+fDeltaLCCPb5Par[i][1];
      fDeltaG112Pb5_LCC33[i] = fGammaLCCPb5Par[i][0]*33+fGammaLCCPb5Par[i][1];
      fDeltaD11Pb5_LCC50[i] = fDeltaLCCPb5Par[i][0]*50+fDeltaLCCPb5Par[i][1];
      fDeltaG112Pb5_LCC50[i] = fGammaLCCPb5Par[i][0]*50+fGammaLCCPb5Par[i][1];
      
      fDeltaD11Pb5_LCC15Err[i] = fDeltaLCCPb5ParErr[i][0]*15+fDeltaLCCPb5ParErr[i][1];
      fDeltaG112Pb5_LCC15Err[i] = fGammaLCCPb5ParErr[i][0]*15+fGammaLCCPb5ParErr[i][1];
      fDeltaD11Pb5_LCC33Err[i] = fDeltaLCCPb5ParErr[i][0]*33+fDeltaLCCPb5ParErr[i][1];
      fDeltaG112Pb5_LCC33Err[i] = fGammaLCCPb5ParErr[i][0]*33+fGammaLCCPb5ParErr[i][1];
      fDeltaD11Pb5_LCC50Err[i] = fDeltaLCCPb5ParErr[i][0]*50+fDeltaLCCPb5ParErr[i][1];
      fDeltaG112Pb5_LCC50Err[i] = fGammaLCCPb5ParErr[i][0]*50+fGammaLCCPb5ParErr[i][1];
    }
    else {
      fDeltaD11Pb5_LCC15[i] = 0;
      fDeltaG112Pb5_LCC15[i] = 0;
      fDeltaD11Pb5_LCC33[i] = 0;
      fDeltaG112Pb5_LCC33[i] = 0;
      fDeltaD11Pb5_LCC50[i] = 0;
      fDeltaG112Pb5_LCC50[i] = 0;
      
      fDeltaD11Pb5_LCC15Err[i] = 0;
      fDeltaG112Pb5_LCC15Err[i] = 0;
      fDeltaD11Pb5_LCC33Err[i] = 0;
      fDeltaG112Pb5_LCC33Err[i] = 0;
      fDeltaD11Pb5_LCC50Err[i] = 0;
      fDeltaG112Pb5_LCC50Err[i] = 0;
    }
    
    if (i>=1) { // loaded for i>=1 for Xe
      fDeltaD11Xe_LCC15[i] = fDeltaLCCXePar[i][0]*15+fDeltaLCCXePar[i][1];
      fDeltaG112Xe_LCC15[i] = fGammaLCCXePar[i][0]*15+fGammaLCCXePar[i][1];
      fDeltaD11Xe_LCC33[i] = fDeltaLCCXePar[i][0]*33+fDeltaLCCXePar[i][1];
      fDeltaG112Xe_LCC33[i] = fGammaLCCXePar[i][0]*33+fGammaLCCXePar[i][1];
      fDeltaD11Xe_LCC50[i] = fDeltaLCCXePar[i][0]*50+fDeltaLCCXePar[i][1];
      fDeltaG112Xe_LCC50[i] = fGammaLCCXePar[i][0]*50+fGammaLCCXePar[i][1];
      
      fDeltaD11Xe_LCC15Err[i] = fDeltaLCCXeParErr[i][0]*15+fDeltaLCCXeParErr[i][1];
      fDeltaG112Xe_LCC15Err[i] = fGammaLCCXeParErr[i][0]*15+fGammaLCCXeParErr[i][1];
      fDeltaD11Xe_LCC33Err[i] = fDeltaLCCXeParErr[i][0]*33+fDeltaLCCXeParErr[i][1];
      fDeltaG112Xe_LCC33Err[i] = fGammaLCCXeParErr[i][0]*33+fGammaLCCXeParErr[i][1];
      fDeltaD11Xe_LCC50Err[i] = fDeltaLCCXeParErr[i][0]*50+fDeltaLCCXeParErr[i][1];
      fDeltaG112Xe_LCC50Err[i] = fGammaLCCXeParErr[i][0]*50+fGammaLCCXeParErr[i][1];
    }
    else {
      fDeltaD11Xe_LCC15[i] = 0;
      fDeltaG112Xe_LCC15[i] = 0;
      fDeltaD11Xe_LCC33[i] = 0;
      fDeltaG112Xe_LCC33[i] = 0;
      fDeltaD11Xe_LCC50[i] = 0;
      fDeltaG112Xe_LCC50[i] = 0;
      
      fDeltaD11Xe_LCC15Err[i] = 0;
      fDeltaG112Xe_LCC15Err[i] = 0;
      fDeltaD11Xe_LCC33Err[i] = 0;
      fDeltaG112Xe_LCC33Err[i] = 0;
      fDeltaD11Xe_LCC50Err[i] = 0;
      fDeltaG112Xe_LCC50Err[i] = 0;
    }
  }
  
  // load data point for others
  // =====> DeltaD11_Baseline_544TeV
  Double_t fDeltaD11Xe_Baseline[7] = {-0.000331522, -0.000126885, -2.62533e-05, 0.000155665, 0.000363263, 0.000652082, 0.00130889};
  Double_t fDeltaD11Xe_BaselineErr[7] = {1.0032e-05, 6.45053e-06, 8.81636e-06, 1.29124e-05, 2.0362e-05, 3.45409e-05, 0.000108742};
  // =====> DeltaG112_Baseline_544TeV
  Double_t fDeltaG112Xe_Baseline[7] = {-3.97637e-05, -7.25526e-05, -9.3799e-05, -4.42574e-05, 1.88029e-06, 3.46793e-05, -0.000110355};
  Double_t fDeltaG112Xe_BaselineErr[7] = {9.63893e-06, 6.29834e-06, 8.70659e-06, 1.27697e-05, 2.01926e-05, 3.44832e-05, 0.000107899};
  // =====> DeltaD11_alpha_0.05_544TeV
  Double_t fDeltaD11Xe_alpha005[7] = {0, -0.00021634, -0.000121582, 2.34761e-05, 0.000103792, 0.000171474, 0.000565419};
  Double_t fDeltaD11Xe_alpha005Err[7] = {0, 9.48569e-06, 1.26704e-05, 1.8274e-05, 2.90105e-05, 4.93592e-05, 9.32429e-05};
  // =====> DeltaG112_alpha_0.05_544TeV
  Double_t fDeltaG112Xe_alpha005[7] = {0, -1.14352e-05, 1.81063e-05, 8.95028e-05, 0.000220084, 0.000440227, 0.000726613};
  Double_t fDeltaG112Xe_alpha005Err[7] = {0, 9.04085e-06, 1.25813e-05, 1.81599e-05, 2.89246e-05, 4.92971e-05, 9.31471e-05};
  // =====> DeltaD11_alpha_0.07_544TeV
  Double_t fDeltaD11Xe_alpha007[7] = {0, -0.000164271, -0.000190888, -0.000163383, -0.000124993, -0.000174209, -0.00018937};
  Double_t fDeltaD11Xe_alpha007Err[7] = {0, 9.37509e-06, 1.26357e-05, 1.86275e-05, 2.95154e-05, 5.01936e-05, 9.41364e-05};
  // =====> DeltaG112_alpha_0.07_544TeV
  Double_t fDeltaG112Xe_alpha007[7] = {0, -5.13023e-05, 0.000123957, 0.000190992, 0.000433107, 0.000678548, 0.00143401};
  Double_t fDeltaG112Xe_alpha007Err[7] = {0, 9.24852e-06, 1.2494e-05, 1.84723e-05, 2.92967e-05, 4.97479e-05, 9.35928e-05};
  // =====> DeltaD11_alpha_0.10_544TeV
  Double_t fDeltaD11Xe_alpha010[7] = {0, -0.000337974, -0.000429772, -0.000398192, -0.000630153, -0.00102805, -0.00155161}; 
  Double_t fDeltaD11Xe_alpha010Err[7] = {0, 9.38455e-06, 1.31214e-05, 1.89359e-05, 2.99713e-05, 5.16764e-05, 9.7365e-05};
  // =====> DeltaG112_alpha_0.10_544TeV
  Double_t fDeltaG112Xe_alpha010[7] = {0, 0.00010246, 0.000297959, 0.000487409, 0.00085607, 0.00147142, 0.00243475};
  Double_t fDeltaG112Xe_alpha010Err[7] = {0, 9.20617e-06, 1.30119e-05, 1.88362e-05, 2.98305e-05, 5.14506e-05, 9.68323e-05};
  // =====> DeltaD11_Baseline_502TeV
  Double_t fDeltaD11Pb5_Baseline[7] = {-0.000109385, -5.34092e-05, 1.27933e-05, 7.62681e-05, 0.000306204, 0.000530114, 0.000950383};
  Double_t fDeltaD11Pb5_BaselineErr[7] = {7.58608e-06, 4.74629e-06, 6.86963e-06, 9.7541e-06, 1.51568e-05, 2.55744e-05, 5.04002e-05};
  // =====> DeltaG112_Baseline_502TeV
  Double_t fDeltaG112Pb5_Baseline[7] = {-6.31493e-05, -6.30012e-05, -5.61324e-05, -3.51231e-05, -2.44125e-05, 2.00313e-05, 0.000109326};
  Double_t fDeltaG112Pb5_BaselineErr[7] = {7.38216e-06, 4.63792e-06, 6.7758e-06, 9.65316e-06, 1.51612e-05, 2.54402e-05, 5.02408e-05};
  // =====> DeltaD11_alpha_0.05_502TeV
  Double_t fDeltaD11Pb5_alpha005[7] = {0, 0, -3.91463e-05, -3.24798e-05, 0.000117713, 0.00021277, 0.000424306};
  Double_t fDeltaD11Pb5_alpha005Err[7] = {0, 0, 6.51407e-06, 9.70705e-06, 1.47769e-05, 2.58273e-05, 5.10005e-05};
  // =====> DeltaG112_alpha_005_502TeV
  Double_t fDeltaG112Pb5_alpha005[7] = {0, 0, 4.38116e-06, 1.96707e-05, 0.000155515, 0.000292368, 0.000588499}; 
  Double_t fDeltaG112Pb5_alpha005Err[7] = {0, 0, 6.42795e-06, 9.64005e-06, 1.47959e-05, 2.5706e-05, 5.084e-05};
  // =====> DeltaD11_alpha_007_502TeV
  Double_t fDeltaD11Pb5_alpha007[7] = {0, 0, -0.000120007, -0.000113831, -5.1251e-05, -9.87425e-05, -0.000119845};
  Double_t fDeltaD11Pb5_alpha007Err[7] = {0, 0, 6.65681e-06, 9.73736e-06, 1.5245e-05, 2.65579e-05, 5.05958e-05};
  // =====> DeltaG112_alpha_007_502TeV
  Double_t fDeltaG112Pb5_alpha007[7] = {0, 0, 5.92242e-05, 0.000173973, 0.000288812, 0.000556749, 0.000948951};
  Double_t fDeltaG112Pb5_alpha007Err[7] = {0, 0, 6.61801e-06, 9.69286e-06, 1.52306e-05, 2.65984e-05, 5.03822e-05};
  // =====> DeltaD11_alpha_010_502TeV
  Double_t fDeltaD11Pb5_alpha010[7] = {0, 0, -0.000267316, -0.000385626, -0.000443574, -0.00064701, -0.00117045};
  Double_t fDeltaD11Pb5_alpha010Err[7] = {0, 0, 6.87774e-06, 1.02838e-05, 1.57763e-05, 2.69816e-05, 5.20458e-05};
  // =====> DeltaG112_alpha_010_502TeV
  Double_t fDeltaG112Pb5_alpha010[7] = {0, 0, 0.000170681, 0.00034874, 0.000573241, 0.00113717, 0.00213};
  Double_t fDeltaG112Pb5_alpha010Err[7] = {0, 0, 6.82479e-06, 1.01803e-05, 1.5729e-05, 2.68995e-05, 5.17827e-05};
 
//here I try to get bin contents and errors as doubles to pass them to my TGraphErrors in the next step !!! I have only 5 centrality bins !!! 
  Double_t JDeltaG112Xe_a_010[7];
  Double_t JDeltaG112Xe_a_010Err[7];
  Double_t JDeltaD11Xe_a_010[7];
  Double_t JDeltaD11Xe_a_010Err[7];
  for(int i = 0; i < 8; i++){
	 JDeltaG112Xe_a_010[i] = DeltaG112Xe_a_010 -> GetBinContent(i);
 	 JDeltaG112Xe_a_010Err[i] = DeltaG112Xe_a_010 -> GetBinError(i);
 	 JDeltaD11Xe_a_010[i] = DeltaD11Xe_a_010 -> GetBinContent(i);
 	 JDeltaD11Xe_a_010Err[i] = DeltaD11Xe_a_010 -> GetBinError(i);
  }
//here I caluclate 1) the differences J - f and the errors
  Double_t DiffD11Xe_a_010[7];
  Double_t DiffD11Xe_a_010Err[7];
  Double_t DiffG112Xe_a_010[7];
  Double_t DiffG112Xe_a_010Err[7];
  for(int i = 0; i < 2; i++){
  DiffD11Xe_a_010[i] = 0;
  DiffD11Xe_a_010Err[i] = 0;
  DiffG112Xe_a_010[i] = 0;
  DiffG112Xe_a_010Err[i] = 0;
  }
  for(int i = 2; i < 8; i++){
	 DiffD11Xe_a_010[i] = JDeltaD11Xe_a_010[i] - fDeltaD11Xe_alpha010[i];
         DiffD11Xe_a_010Err[i] = TMath::Sqrt(TMath::Power(JDeltaG112Xe_a_010Err[i],2) + TMath::Power(fDeltaD11Xe_alpha010Err[i],2));
         DiffG112Xe_a_010[i] = JDeltaG112Xe_a_010[i] - fDeltaG112Xe_alpha010[i];
         DiffG112Xe_a_010Err[i] = TMath::Sqrt(TMath::Power(JDeltaG112Xe_a_010Err[i],2) + TMath::Power(fDeltaG112Xe_alpha010Err[i],2));  
  cout<<"DiffD11: "<<DiffD11Xe_a_010[i]<<endl;
  cout<<"J: "<<JDeltaG112Xe_a_010[i]<<" PSJ: "<<fDeltaG112Xe_alpha010[i]<<endl;
  cout<<"DiffG112: "<<DiffG112Xe_a_010[i]<<endl;
  }
//here I calculate 2) the ratios of J/f and the errors
  Double_t RatioD11Xe_a_010[7];
  Double_t RatioD11Xe_a_010Err[7];
  Double_t RatioG112Xe_a_010[7];
  Double_t RatioG112Xe_a_010Err[7];
  for(int i = 0; i < 2; i++){
  RatioD11Xe_a_010[i] = 0;
  RatioD11Xe_a_010Err[i] = 0;
  RatioG112Xe_a_010[i] = 0;
  RatioG112Xe_a_010Err[i] = 0;
  }
  for(int i = 2; i < 8; i++){
	 RatioD11Xe_a_010[i] = (1/JDeltaD11Xe_a_010[i])*fDeltaD11Xe_alpha010[i];
	 RatioD11Xe_a_010Err[i] = RatioD11Xe_a_010[i]*TMath::Sqrt(TMath::Power(JDeltaD11Xe_a_010Err[i]/JDeltaD11Xe_a_010[i],2)+ TMath::Power(fDeltaD11Xe_alpha010Err[i]/fDeltaD11Xe_alpha010[i],2));
	 RatioG112Xe_a_010[i] = (1/JDeltaG112Xe_a_010[i])*fDeltaG112Xe_alpha010[i];
	 RatioG112Xe_a_010Err[i] = RatioG112Xe_a_010[i]*TMath::Sqrt( TMath::Power(JDeltaG112Xe_a_010Err[i]/JDeltaG112Xe_a_010[i],2)+ TMath::Power(fDeltaG112Xe_alpha010Err[i]/fDeltaG112Xe_alpha010[i],2));
  cout<<"RatioD11: "<<RatioD11Xe_a_010[i]<<"J: "<< JDeltaD11Xe_a_010[i]<<" PSJ: "<<fDeltaD11Xe_alpha010[i]<<endl;
  cout<<"RatioG112: "<<RatioG112Xe_a_010[i]<<"J: "<< JDeltaG112Xe_a_010[i]<<" PSJ: "<<fDeltaG112Xe_alpha010[i]<<endl;
  }

  // loaded for i>=2 for Pb5
  // loaded for i>=1 for Xe
  Double_t gCentralityAVFD[7] = { 5.,15.,25.,35.,45.,55.,65.};
  Double_t gCentralityAVFDError[7] = {0.,0.,0.,0.,0.,0.,0.};
  
  Double_t gDelta1AVFD0[7] = {3.2981916e-07,3.3153323e-05,0.00015727837,0.00021617493,0.00037009935,0.00062122336,0.0012135616};
  Double_t gDelta1AVFD0Error[7] = {7.4361248e-06,6.7075673e-06,8.2698811e-06,1.1559705e-05,1.7802286e-05,2.9499850e-05,5.5649497e-05};
  
  TGraphErrors *grDelta1AVFD0 = new TGraphErrors(7,gCentralityAVFD,gDelta1AVFD0,gCentralityAVFDError,gDelta1AVFD0Error);
  grDelta1AVFD0->SetLineColor(kGreen+2);
  grDelta1AVFD0->SetLineWidth(1);
  grDelta1AVFD0->SetFillColor(kGreen+2);
  grDelta1AVFD0->SetFillStyle(1);

//here I have to pass the extracted values from above to the TGraphErrors -  have only 5 centrality bins!
  TGraphErrors* JGraphDiffDeltaD11Xe_a_010 = new TGraphErrors(7,gCentralityAVFD, DiffD11Xe_a_010, gCentralityAVFDError, DiffD11Xe_a_010Err);
  TGraphErrors* JGraphDiffDeltaG112Xe_a_010 = new TGraphErrors(7, gCentralityAVFD, DiffG112Xe_a_010, gCentralityAVFDError, DiffG112Xe_a_010Err);
  TGraphErrors* JGraphRatioDeltaD11Xe_a_010 = new TGraphErrors(7, gCentralityAVFD, RatioD11Xe_a_010, gCentralityAVFDError, RatioD11Xe_a_010Err); 
  TGraphErrors* JGraphRatioDeltaG112Xe_a_010 = new TGraphErrors(7, gCentralityAVFD, RatioG112Xe_a_010, gCentralityAVFDError, RatioG112Xe_a_010Err);

  TGraphErrors* JGraphErrorsDeltaD11Xe_a_010 = new TGraphErrors(7, gCentralityAVFD, JDeltaD11Xe_a_010, gCentralityAVFDError, JDeltaD11Xe_a_010Err );
  TGraphErrors* JGraphErrorsDeltaG112Xe_a_010 = new TGraphErrors(7, gCentralityAVFD, JDeltaG112Xe_a_010, gCentralityAVFDError, JDeltaG112Xe_a_010Err );

  TGraphErrors* fGraphErrorsDeltaD11Pb5_Baseline = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_Baseline, gCentralityAVFDError, fDeltaD11Pb5_BaselineErr);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_Baseline = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_Baseline, gCentralityAVFDError, fDeltaG112Pb5_BaselineErr);
  
  TGraphErrors* fGraphErrorsDeltaD11Pb5_LCC15 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_LCC15, gCentralityAVFDError, fDeltaD11Pb5_LCC15Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_LCC15 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_LCC15, gCentralityAVFDError, fDeltaG112Pb5_LCC15Err);
  TGraphErrors* fGraphErrorsDeltaD11Pb5_LCC33 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_LCC33, gCentralityAVFDError, fDeltaD11Pb5_LCC33Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_LCC33 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_LCC33, gCentralityAVFDError, fDeltaG112Pb5_LCC33Err);
  TGraphErrors* fGraphErrorsDeltaD11Pb5_LCC50 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_LCC50, gCentralityAVFDError, fDeltaD11Pb5_LCC50Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_LCC50 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_LCC50, gCentralityAVFDError, fDeltaG112Pb5_LCC50Err);
  
  TGraphErrors* fGraphErrorsDeltaD11Pb5_alpha005 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_alpha005, gCentralityAVFDError, fDeltaD11Pb5_alpha005Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_alpha005 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_alpha005, gCentralityAVFDError, fDeltaG112Pb5_alpha005Err);
  TGraphErrors* fGraphErrorsDeltaD11Pb5_alpha007 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_alpha007, gCentralityAVFDError, fDeltaD11Pb5_alpha007Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_alpha007 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_alpha007, gCentralityAVFDError, fDeltaG112Pb5_alpha007Err);
  TGraphErrors* fGraphErrorsDeltaD11Pb5_alpha010 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Pb5_alpha010, gCentralityAVFDError, fDeltaD11Pb5_alpha010Err);
  TGraphErrors* fGraphErrorsDeltaG112Pb5_alpha010 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Pb5_alpha010, gCentralityAVFDError, fDeltaG112Pb5_alpha010Err);
  
  TGraphErrors* fGraphErrorsDeltaD11Xe_Baseline = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_Baseline, gCentralityAVFDError, fDeltaD11Xe_BaselineErr);
  TGraphErrors* fGraphErrorsDeltaG112Xe_Baseline = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_Baseline, gCentralityAVFDError, fDeltaG112Xe_BaselineErr);
  
  TGraphErrors* fGraphErrorsDeltaD11Xe_LCC15 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_LCC15, gCentralityAVFDError, fDeltaD11Xe_LCC15Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_LCC15 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_LCC15, gCentralityAVFDError, fDeltaG112Xe_LCC15Err);
  TGraphErrors* fGraphErrorsDeltaD11Xe_LCC33 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_LCC33, gCentralityAVFDError, fDeltaD11Xe_LCC33Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_LCC33 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_LCC33, gCentralityAVFDError, fDeltaG112Xe_LCC33Err);
  TGraphErrors* fGraphErrorsDeltaD11Xe_LCC50 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_LCC50, gCentralityAVFDError, fDeltaD11Xe_LCC50Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_LCC50 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_LCC50, gCentralityAVFDError, fDeltaG112Xe_LCC50Err);
  
  TGraphErrors* fGraphErrorsDeltaD11Xe_alpha005 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_alpha005, gCentralityAVFDError, fDeltaD11Xe_alpha005Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_alpha005 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_alpha005, gCentralityAVFDError, fDeltaG112Xe_alpha005Err);
  TGraphErrors* fGraphErrorsDeltaD11Xe_alpha007 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_alpha007, gCentralityAVFDError, fDeltaD11Xe_alpha007Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_alpha007 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_alpha007, gCentralityAVFDError, fDeltaG112Xe_alpha007Err);
  TGraphErrors* fGraphErrorsDeltaD11Xe_alpha010 = new TGraphErrors(7, gCentralityAVFD, fDeltaD11Xe_alpha010, gCentralityAVFDError, fDeltaD11Xe_alpha010Err);
  TGraphErrors* fGraphErrorsDeltaG112Xe_alpha010 = new TGraphErrors(7, gCentralityAVFD, fDeltaG112Xe_alpha010, gCentralityAVFDError, fDeltaG112Xe_alpha010Err);
  ///////////////////////////////////////////////////
  fGraphErrorsDeltaD11Pb5_Baseline->SetMarkerStyle(20);
  fGraphErrorsDeltaD11Pb5_Baseline->SetMarkerColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_Baseline->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_Baseline->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_Baseline->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_Baseline->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_Baseline->SetFillStyle(3003);
  fGraphErrorsDeltaG112Pb5_Baseline->SetMarkerStyle(20);
  fGraphErrorsDeltaG112Pb5_Baseline->SetMarkerColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_Baseline->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_Baseline->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_Baseline->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_Baseline->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_Baseline->SetFillStyle(3003);
  
  fGraphErrorsDeltaD11Pb5_LCC15->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC15->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_LCC15->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_LCC15->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC15->SetFillStyle(1000);
  fGraphErrorsDeltaG112Pb5_LCC15->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC15->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_LCC15->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC15->SetFillStyle(1000);
  fGraphErrorsDeltaD11Pb5_LCC33->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC33->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_LCC33->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_LCC33->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC33->SetFillStyle(3001);
  fGraphErrorsDeltaG112Pb5_LCC33->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC33->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_LCC33->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC33->SetFillStyle(3001);
  fGraphErrorsDeltaD11Pb5_LCC50->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC50->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_LCC50->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_LCC50->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_LCC50->SetFillStyle(3008);
  fGraphErrorsDeltaG112Pb5_LCC50->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC50->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_LCC50->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_LCC50->SetFillStyle(3008);
  
  fGraphErrorsDeltaD11Xe_Baseline->SetMarkerStyle(24);
  fGraphErrorsDeltaD11Xe_Baseline->SetMarkerColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_Baseline->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_Baseline->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_Baseline->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_Baseline->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_Baseline->SetFillStyle(3144);
  fGraphErrorsDeltaG112Xe_Baseline->SetMarkerStyle(24);
  fGraphErrorsDeltaG112Xe_Baseline->SetMarkerColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_Baseline->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_Baseline->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_Baseline->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_Baseline->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_Baseline->SetFillStyle(3144);
  
  fGraphErrorsDeltaD11Xe_LCC15->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC15->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_LCC15->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_LCC15->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC15->SetFillStyle(1000);
  fGraphErrorsDeltaG112Xe_LCC15->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC15->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_LCC15->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC15->SetFillStyle(1000);
  fGraphErrorsDeltaD11Xe_LCC33->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC33->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_LCC33->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_LCC33->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC33->SetFillStyle(3001);
  fGraphErrorsDeltaG112Xe_LCC33->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC33->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_LCC33->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC33->SetFillStyle(3001);
  fGraphErrorsDeltaD11Xe_LCC50->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC50->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_LCC50->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_LCC50->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_LCC50->SetFillStyle(3008);
  fGraphErrorsDeltaG112Xe_LCC50->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC50->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_LCC50->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_LCC50->SetFillStyle(3008);
  
  // ====================================================
  fGraphErrorsDeltaD11Pb5_alpha005->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha005->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_alpha005->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_alpha005->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha005->SetFillStyle(1000);
  fGraphErrorsDeltaG112Pb5_alpha005->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha005->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_alpha005->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha005->SetFillStyle(1000);
  fGraphErrorsDeltaD11Pb5_alpha007->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha007->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_alpha007->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_alpha007->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha007->SetFillStyle(3008);
  fGraphErrorsDeltaG112Pb5_alpha007->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha007->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_alpha007->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha007->SetFillStyle(3008);
  fGraphErrorsDeltaD11Pb5_alpha010->SetLineColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha010->SetLineWidth(1);
  fGraphErrorsDeltaD11Pb5_alpha010->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Pb5_alpha010->SetFillColor(kBlue+2);
  fGraphErrorsDeltaD11Pb5_alpha010->SetFillStyle(1);
  fGraphErrorsDeltaG112Pb5_alpha010->SetLineColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha010->SetLineWidth(1);
  fGraphErrorsDeltaG112Pb5_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Pb5_alpha010->SetFillColor(kBlue+2);
  fGraphErrorsDeltaG112Pb5_alpha010->SetFillStyle(1);
//here I add mine:
  JGraphErrorsDeltaG112Xe_a_010->SetLineColor(kRed+2);
  JGraphErrorsDeltaG112Xe_a_010->SetLineWidth(1);
  JGraphErrorsDeltaG112Xe_a_010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  JGraphErrorsDeltaG112Xe_a_010->SetFillColor(kRed+2);  
  JGraphErrorsDeltaG112Xe_a_010->SetFillStyle(1);
  
  JGraphErrorsDeltaD11Xe_a_010->SetLineColor(kRed+2);   
  JGraphErrorsDeltaD11Xe_a_010->SetLineWidth(1);
  JGraphErrorsDeltaD11Xe_a_010->GetYaxis()->SetTitle("#Delta#delta_{11}");                    
  JGraphErrorsDeltaD11Xe_a_010->SetFillColor(kRed+2);
  JGraphErrorsDeltaD11Xe_a_010->SetFillStyle(100);


  fGraphErrorsDeltaD11Xe_alpha005->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha005->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_alpha005->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_alpha005->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha005->SetFillStyle(1000);
  fGraphErrorsDeltaG112Xe_alpha005->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha005->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_alpha005->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha005->SetFillStyle(1000);
  fGraphErrorsDeltaD11Xe_alpha007->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha007->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_alpha007->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_alpha007->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha007->SetFillStyle(3008);
  fGraphErrorsDeltaG112Xe_alpha007->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha007->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_alpha007->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha007->SetFillStyle(3008);
  fGraphErrorsDeltaD11Xe_alpha010->SetLineColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha010->SetLineWidth(1);
  fGraphErrorsDeltaD11Xe_alpha010->GetYaxis()->SetTitle("#Delta#delta_{11}");
  fGraphErrorsDeltaD11Xe_alpha010->SetFillColor(kGreen+2);
  fGraphErrorsDeltaD11Xe_alpha010->SetFillStyle(1);
  fGraphErrorsDeltaG112Xe_alpha010->SetLineColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha010->SetLineWidth(1);
  fGraphErrorsDeltaG112Xe_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphErrorsDeltaG112Xe_alpha010->SetFillColor(kGreen+2);
  fGraphErrorsDeltaG112Xe_alpha010->SetFillStyle(1);
//here I add the differences and ratios
  JGraphDiffDeltaG112Xe_a_010->SetLineColor(kRed+2);
  JGraphDiffDeltaG112Xe_a_010->SetLineWidth(1);
  JGraphDiffDeltaG112Xe_a_010->GetYaxis()->SetTitle("#Delta#gamma_{112}(J) - #Delta#gamma_{112}");
  JGraphDiffDeltaG112Xe_a_010->SetFillColor(kRed+2);
  JGraphDiffDeltaG112Xe_a_010->SetFillStyle(1);
  JGraphDiffDeltaD11Xe_a_010->SetLineColor(kRed+4);
  JGraphDiffDeltaD11Xe_a_010->SetLineWidth(1);
  JGraphDiffDeltaD11Xe_a_010->GetYaxis()->SetTitle("#Delta#delta_{11}(J) - #Delta#delta_{11}");
  JGraphDiffDeltaD11Xe_a_010->SetFillColor(kRed+4);
  JGraphDiffDeltaD11Xe_a_010->SetFillStyle(1000);

  JGraphRatioDeltaG112Xe_a_010->SetLineColor(kRed+2);
  JGraphRatioDeltaG112Xe_a_010->SetLineWidth(1);
  JGraphRatioDeltaG112Xe_a_010->GetYaxis()->SetTitle("#Delta#gamma_{112}(PSJ/J)");
  JGraphRatioDeltaG112Xe_a_010->SetFillColor(kRed+2);
  JGraphRatioDeltaG112Xe_a_010->SetFillStyle(1);
  JGraphRatioDeltaD11Xe_a_010->SetLineColor(kRed+4);
  JGraphRatioDeltaD11Xe_a_010->SetLineWidth(1);
  JGraphRatioDeltaD11Xe_a_010->GetYaxis()->SetTitle("#Delta#delta_{11}(PSJ/J)");
  JGraphRatioDeltaD11Xe_a_010->SetFillColor(kRed+4);
  JGraphRatioDeltaD11Xe_a_010->SetFillStyle(1000);


  //====================vs dN/deta===========================//
  //Pb-Pb
  Double_t gAlice502[7] = {0.25*(2035.43+1850.15+1666.32+1505.11), 1180.48, 786.45, 511.66, 318.15, 183.33, 96.27};
  Double_t gAlice502Error[7] = {0.25*(TMath::Sqrt(TMath::Power(1.08,2) + TMath::Power(52.43,2)) + TMath::Sqrt(TMath::Power(1.19,2) + TMath::Power(55.23,2)) + TMath::Sqrt(TMath::Power(1.12,2) + TMath::Power(47.89,2)) + TMath::Sqrt(TMath::Power(1.07,2) + TMath::Power(43.71,2))),
    TMath::Sqrt(TMath::Power(0.46,2) + TMath::Power(31.18,2)),
    TMath::Sqrt(TMath::Power(0.36,2) + TMath::Power(20.44,2)),
    TMath::Sqrt(TMath::Power(0.28,2) + TMath::Power(15.2,2)),
    TMath::Sqrt(TMath::Power(0.22,2) + TMath::Power(12.27,2)),
    TMath::Sqrt(TMath::Power(0.16,2) + TMath::Power(8.07,2)),
    TMath::Sqrt(TMath::Power(0.12,2) + TMath::Power(5.81,2))};
  
  TGraphErrors* fGraphvsMultDeltaD11Pb5_Baseline = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_Baseline, gAlice502Error, fDeltaD11Pb5_BaselineErr);
  fGraphvsMultDeltaD11Pb5_Baseline->SetMarkerStyle(20);
  fGraphvsMultDeltaD11Pb5_Baseline->SetMarkerColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_Baseline->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_Baseline->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_Baseline->GetYaxis()->SetTitle("#Delta#delta_{11}");
  
  TGraphErrors* fGraphvsMultDeltaG112Pb5_Baseline = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_Baseline, gAlice502Error, fDeltaG112Pb5_BaselineErr);
  fGraphvsMultDeltaG112Pb5_Baseline->SetMarkerStyle(20);
  fGraphvsMultDeltaG112Pb5_Baseline->SetMarkerColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_Baseline->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_Baseline->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_Baseline->GetYaxis()->SetTitle("#Delta#gamma_{112}");

  //LCC-delta
  TGraphErrors* fGraphvsMultDeltaD11Pb5_LCC15 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_LCC15, gAlice502Error, fDeltaD11Pb5_LCC15Err);
  fGraphvsMultDeltaD11Pb5_LCC15->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC15->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_LCC15->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC15->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaD11Pb5_LCC33 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_LCC33, gAlice502Error, fDeltaD11Pb5_LCC33Err);
  fGraphvsMultDeltaD11Pb5_LCC33->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC33->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_LCC33->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC33->SetFillStyle(3001);

  TGraphErrors* fGraphvsMultDeltaD11Pb5_LCC50 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_LCC50, gAlice502Error, fDeltaD11Pb5_LCC50Err);
  fGraphvsMultDeltaD11Pb5_LCC50->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC50->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_LCC50->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_LCC50->SetFillStyle(3008);

  //n5/s-delta
  TGraphErrors* fGraphvsMultDeltaD11Pb5_alpha005 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_alpha005, gAlice502Error, fDeltaD11Pb5_alpha005Err);
  fGraphvsMultDeltaD11Pb5_alpha005->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha005->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_alpha005->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha005->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaD11Pb5_alpha007 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_alpha007, gAlice502Error, fDeltaD11Pb5_alpha007Err);
  fGraphvsMultDeltaD11Pb5_alpha007->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha007->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_alpha007->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha007->SetFillStyle(3008);

  TGraphErrors* fGraphvsMultDeltaD11Pb5_alpha010 = new TGraphErrors(7, gAlice502, fDeltaD11Pb5_alpha010, gAlice502Error, fDeltaD11Pb5_alpha010Err);
  fGraphvsMultDeltaD11Pb5_alpha010->SetLineColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha010->SetLineWidth(1);
  fGraphvsMultDeltaD11Pb5_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Pb5_alpha010->SetFillColor(kBlue+2);
  fGraphvsMultDeltaD11Pb5_alpha010->SetFillStyle(1);

  //LCC-gamma
  TGraphErrors* fGraphvsMultDeltaG112Pb5_LCC15 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_LCC15, gAlice502Error, fDeltaG112Pb5_LCC15Err);
  fGraphvsMultDeltaG112Pb5_LCC15->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC15->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_LCC15->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC15->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaG112Pb5_LCC33 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_LCC33, gAlice502Error, fDeltaG112Pb5_LCC33Err);
  fGraphvsMultDeltaG112Pb5_LCC33->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC33->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_LCC33->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC33->SetFillStyle(3001);

  TGraphErrors* fGraphvsMultDeltaG112Pb5_LCC50 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_LCC50, gAlice502Error, fDeltaG112Pb5_LCC50Err);
  fGraphvsMultDeltaG112Pb5_LCC50->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC50->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_LCC50->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_LCC50->SetFillStyle(3008);

  //n5/s-gamma
  TGraphErrors* fGraphvsMultDeltaG112Pb5_alpha005 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_alpha005, gAlice502Error, fDeltaG112Pb5_alpha005Err);
  fGraphvsMultDeltaG112Pb5_alpha005->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha005->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_alpha005->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha005->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaG112Pb5_alpha007 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_alpha007, gAlice502Error, fDeltaG112Pb5_alpha007Err);
  fGraphvsMultDeltaG112Pb5_alpha007->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha007->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_alpha007->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha007->SetFillStyle(3008);

  TGraphErrors* fGraphvsMultDeltaG112Pb5_alpha010 = new TGraphErrors(7, gAlice502, fDeltaG112Pb5_alpha010, gAlice502Error, fDeltaG112Pb5_alpha010Err);
  fGraphvsMultDeltaG112Pb5_alpha010->SetLineColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha010->SetLineWidth(1);
  fGraphvsMultDeltaG112Pb5_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Pb5_alpha010->SetFillColor(kBlue+2);
  fGraphvsMultDeltaG112Pb5_alpha010->SetFillStyle(1);

  //Xe-Xe
  Double_t gAlice544[7] = {0.5*(1167+939),706,478,315,198,118,64.7};
  Double_t gAlice544Error[7] = {0.5*TMath::Sqrt(26.*26+24.*24),17,11,8,5,3,2};

  TGraphErrors* fGraphvsMultDeltaD11Xe_Baseline = new TGraphErrors(7, gAlice502, fDeltaD11Xe_Baseline, gAlice502Error, fDeltaD11Xe_BaselineErr);
  fGraphvsMultDeltaD11Xe_Baseline->SetMarkerStyle(24);
  fGraphvsMultDeltaD11Xe_Baseline->SetMarkerColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_Baseline->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_Baseline->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_Baseline->GetYaxis()->SetTitle("#Delta#delta_{11}");
  
  TGraphErrors* fGraphvsMultDeltaG112Xe_Baseline = new TGraphErrors(7, gAlice502, fDeltaG112Xe_Baseline, gAlice502Error, fDeltaG112Xe_BaselineErr);
  fGraphvsMultDeltaG112Xe_Baseline->SetMarkerStyle(24);
  fGraphvsMultDeltaG112Xe_Baseline->SetMarkerColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_Baseline->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_Baseline->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_Baseline->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  
  //LCC-delta
  TGraphErrors* fGraphvsMultDeltaD11Xe_LCC15 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_LCC15, gAlice544Error, fDeltaD11Xe_LCC15Err);
  fGraphvsMultDeltaD11Xe_LCC15->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC15->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_LCC15->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC15->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaD11Xe_LCC33 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_LCC33, gAlice544Error, fDeltaD11Xe_LCC33Err);
  fGraphvsMultDeltaD11Xe_LCC33->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC33->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_LCC33->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC33->SetFillStyle(3001);

  TGraphErrors* fGraphvsMultDeltaD11Xe_LCC50 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_LCC50, gAlice544Error, fDeltaD11Xe_LCC50Err);
  fGraphvsMultDeltaD11Xe_LCC50->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC50->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_LCC50->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_LCC50->SetFillStyle(3008);
  
  //n5/s-delta
  TGraphErrors* fGraphvsMultDeltaD11Xe_alpha005 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_alpha005, gAlice544Error, fDeltaD11Xe_alpha005Err);
  fGraphvsMultDeltaD11Xe_alpha005->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha005->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_alpha005->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha005->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaD11Xe_alpha007 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_alpha007, gAlice544Error, fDeltaD11Xe_alpha007Err);
  fGraphvsMultDeltaD11Xe_alpha007->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha007->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_alpha007->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha007->SetFillStyle(3008);

  TGraphErrors* fGraphvsMultDeltaD11Xe_alpha010 = new TGraphErrors(7, gAlice544, fDeltaD11Xe_alpha010, gAlice544Error, fDeltaD11Xe_alpha010Err);
  fGraphvsMultDeltaD11Xe_alpha010->SetLineColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha010->SetLineWidth(1);
  fGraphvsMultDeltaD11Xe_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaD11Xe_alpha010->SetFillColor(kGreen+2);
  fGraphvsMultDeltaD11Xe_alpha010->SetFillStyle(1);

  //LCC-gamma
  TGraphErrors* fGraphvsMultDeltaG112Xe_LCC15 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_LCC15, gAlice544Error, fDeltaG112Xe_LCC15Err);
  fGraphvsMultDeltaG112Xe_LCC15->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC15->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_LCC15->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_LCC15->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC15->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaG112Xe_LCC33 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_LCC33, gAlice544Error, fDeltaG112Xe_LCC33Err);
  fGraphvsMultDeltaG112Xe_LCC33->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC33->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_LCC33->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_LCC33->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC33->SetFillStyle(3001);

  TGraphErrors* fGraphvsMultDeltaG112Xe_LCC50 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_LCC50, gAlice544Error, fDeltaG112Xe_LCC50Err);
  fGraphvsMultDeltaG112Xe_LCC50->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC50->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_LCC50->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_LCC50->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_LCC50->SetFillStyle(3008);
  
  //n5/s-gamma
  TGraphErrors* fGraphvsMultDeltaG112Xe_alpha005 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_alpha005, gAlice544Error, fDeltaG112Xe_alpha005Err);
  fGraphvsMultDeltaG112Xe_alpha005->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha005->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_alpha005->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_alpha005->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha005->SetFillStyle(1000);

  TGraphErrors* fGraphvsMultDeltaG112Xe_alpha007 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_alpha007, gAlice544Error, fDeltaG112Xe_alpha007Err);
  fGraphvsMultDeltaG112Xe_alpha007->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha007->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_alpha007->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_alpha007->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha007->SetFillStyle(3008);

  TGraphErrors* fGraphvsMultDeltaG112Xe_alpha010 = new TGraphErrors(7, gAlice544, fDeltaG112Xe_alpha010, gAlice544Error, fDeltaG112Xe_alpha010Err);
  fGraphvsMultDeltaG112Xe_alpha010->SetLineColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha010->SetLineWidth(1);
  fGraphvsMultDeltaG112Xe_alpha010->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  fGraphvsMultDeltaG112Xe_alpha010->SetFillColor(kGreen+2);
  fGraphvsMultDeltaG112Xe_alpha010->SetFillStyle(1);

  //=============Draw the results============================//
  const Int_t nCentralityBins = 7;
  TH2F *gEmpty0 = new TH2F("gEmpty0",";dN/d#eta;", 1000,0,2500,1000,-1.5e-03,0.015);
  gEmpty0->SetStats(kFALSE);
  gEmpty0->GetYaxis()->SetTitleSize(0.07);
  gEmpty0->GetYaxis()->SetTitleOffset(1.2);
  gEmpty0->GetYaxis()->SetNdivisions(10);
  gEmpty0->GetXaxis()->SetNdivisions(5);
  gEmpty0->GetYaxis()->CenterTitle();
  gEmpty0->GetXaxis()->CenterTitle();

  TH2F *gEmpty = new TH2F("gEmpty",";centrality, %;", nCentralityBins,0,72,1000,-1.5e-03,0.015);
  gEmpty->SetStats(kFALSE);
  gEmpty->GetYaxis()->SetTitleSize(0.07);
  gEmpty->GetYaxis()->SetTitleOffset(0.95);
  gEmpty->GetYaxis()->SetNdivisions(10);
  gEmpty->GetXaxis()->SetNdivisions(10);
  gEmpty->GetYaxis()->CenterTitle();
  gEmpty->GetXaxis()->CenterTitle();
  
  
  TLatex *text = new TLatex();
  text->SetTextFont(42);
  text->SetTextSize(0.07);

  TLegend *legend1 =new TLegend(0.20,0.4,0.45,0.865);
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(0.058);
  legend1->SetHeader("AVFD");
  legend1->AddEntry(fGraphErrorsDeltaG112Pb5_Baseline,"Pb-Pb Baseline","P");
  legend1->AddEntry(fGraphErrorsDeltaG112Xe_Baseline,"Xe-Xe Baseline","P");
  legend1->AddEntry(fGraphErrorsDeltaG112Pb5_LCC15,"Pb-Pb n_{5}/s=0.0-LCC=15%","F");
  legend1->AddEntry(fGraphErrorsDeltaG112Xe_LCC15,"Xe-Xe n_{5}/s=0.0-LCC=15%","F");
  //legend1->AddEntry(fGraphErrorsDeltaG112Pb5_LCC33,"Pb-Pb n_{5}/s=0.0-LCC=33%","F");
  //legend1->AddEntry(fGraphErrorsDeltaG112Xe_LCC33,"Xe-Xe n_{5}/s=0.0-LCC=33%","F");
  legend1->AddEntry(fGraphErrorsDeltaG112Pb5_LCC50,"Pb-Pb n_{5}/s=0.0-LCC=50%","F");
  legend1->AddEntry(fGraphErrorsDeltaG112Xe_LCC50,"Xe-Xe n_{5}/s=0.0-LCC=50%","F");

  TLegend *legend2 =new TLegend(0.20,0.66,0.45,0.97);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.058);
  legend2->AddEntry(fGraphErrorsDeltaG112Pb5_alpha005,"Pb-Pb n_{5}/s=0.05-LCC=0%","F");
  legend2->AddEntry(fGraphErrorsDeltaG112Xe_alpha005,"Xe-Xe n_{5}/s=0.05-LCC=0%","F");
  legend2->AddEntry(fGraphErrorsDeltaG112Pb5_alpha007,"Pb-Pb n_{5}/s=0.07-LCC=0%","F");
  legend2->AddEntry(fGraphErrorsDeltaG112Xe_alpha007,"Xe-Xe n_{5}/s=0.07-LCC=0%","F");
  legend2->AddEntry(fGraphErrorsDeltaG112Pb5_alpha010,"Pb-Pb n_{5}/s=0.10-LCC=0%","F");
  legend2->AddEntry(fGraphErrorsDeltaG112Xe_alpha010,"Xe-Xe n_{5}/s=0.10-LCC=0%","F");
//here mine:
  legend2->AddEntry(JGraphErrorsDeltaG112Xe_a_010, "J: Xe-Xe n_{5}/s=0.10-LCC=0%","F");

  TLegend *legend3 =new TLegend(0.55,0.4,0.85,0.865);
  legend3->SetBorderSize(0);
  legend3->SetFillColor(0);
  legend3->SetTextFont(42);
  legend3->SetTextSize(0.058);
  legend3->SetHeader("AVFD");
  legend3->AddEntry(fGraphErrorsDeltaG112Pb5_Baseline,"Pb-Pb Baseline","P");
  legend3->AddEntry(fGraphErrorsDeltaG112Xe_Baseline,"Xe-Xe Baseline","P");
  legend3->AddEntry(fGraphErrorsDeltaG112Pb5_LCC15,"Pb-Pb n_{5}/s=0.0-LCC=15%","F");
  legend3->AddEntry(fGraphErrorsDeltaG112Xe_LCC15,"Xe-Xe n_{5}/s=0.0-LCC=15%","F");
  //legend3->AddEntry(fGraphErrorsDeltaG112Pb5_LCC33,"Pb-Pb n_{5}/s=0.0-LCC=33%","F");
  //legend3->AddEntry(fGraphErrorsDeltaG112Xe_LCC33,"Xe-Xe n_{5}/s=0.0-LCC=33%","F");
  legend3->AddEntry(fGraphErrorsDeltaG112Pb5_LCC50,"Pb-Pb n_{5}/s=0.0-LCC=50%","F");
  legend3->AddEntry(fGraphErrorsDeltaG112Xe_LCC50,"Xe-Xe n_{5}/s=0.0-LCC=50%","F");
  
  TLegend *legend4 =new TLegend(0.55,0.66,0.85,0.97);
  legend4->SetBorderSize(0);
  legend4->SetFillColor(0);
  legend4->SetTextFont(42);
  legend4->SetTextSize(0.045);
  legend4->AddEntry(fGraphErrorsDeltaG112Pb5_alpha005,"Pb-Pb n_{5}/s=0.05-LCC=0%","F");
  legend4->AddEntry(fGraphErrorsDeltaG112Xe_alpha005,"Xe-Xe n_{5}/s=0.05-LCC=0%","F");
  legend4->AddEntry(fGraphErrorsDeltaG112Pb5_alpha007,"Pb-Pb n_{5}/s=0.07-LCC=0%","F");
  legend4->AddEntry(fGraphErrorsDeltaG112Xe_alpha007,"Xe-Xe n_{5}/s=0.07-LCC=0%","F");
  legend4->AddEntry(fGraphErrorsDeltaG112Pb5_alpha010,"Pb-Pb n_{5}/s=0.10-LCC=0%","F");
  legend4->AddEntry(fGraphErrorsDeltaG112Xe_alpha010,"Xe-Xe n_{5}/s=0.10-LCC=0%","F");
  
//here mine:      
  TLegend *legend5 =new TLegend(0.20,0.56,0.45,0.7);
  legend5->SetBorderSize(0);
  legend5->SetFillColor(0);
  legend5->SetTextFont(42);
  legend5->SetTextSize(0.035);//from 0.045
  legend5->AddEntry(JGraphDiffDeltaD11Xe_a_010,"Xe-Xe n_{5}/s=0.1-LCC=0%","F");

  TLegend *legend6 =new TLegend(0.20,0.66,0.45,0.8);
  legend6->SetBorderSize(0);
  legend6->SetFillColor(0);
  legend6->SetTextFont(42);
  legend6->SetTextSize(0.035);
  legend6->AddEntry(JGraphDiffDeltaG112Xe_a_010,"Xe-Xe n_{5}/s=0.1-LCC=0%","F");

  TLegend *legend7 =new TLegend(0.20,0.66,0.45,0.8);
  legend7->SetBorderSize(0);
  legend7->SetFillColor(0);
  legend7->SetTextFont(42);
  legend7->SetTextSize(0.035);
  legend7->AddEntry(JGraphRatioDeltaD11Xe_a_010,"Xe-Xe n_{5}/s=0.1-LCC=0%","F");

  TLegend *legend8 =new TLegend(0.20,0.66,0.45,0.8);
  legend8->SetBorderSize(0);
  legend8->SetFillColor(0);
  legend8->SetTextFont(42);
  legend8->SetTextSize(0.035);
  legend8->AddEntry(JGraphRatioDeltaG112Xe_a_010,"Xe-Xe n_{5}/s=0.1-LCC=0%","F");

  TF1 *f1 = new TF1("f1","0",0,1000);
  f1->SetLineColor(1); 
  f1->SetLineStyle(1); 
  f1->SetLineWidth(1);
  
 //____________________Differences_______________________//
  TCanvas *d1 = new TCanvas("d1","Comparison: J - PSJ",0,0,700,800);
  d1->SetFillColor(10);
  d1->SetHighLightColor(10);
  d1->Divide(1,2,0.99,0.0,10);
 
 //===============Delta delta=============//
  d1->cd(1)->SetLeftMargin(0.19);
  d1->cd(1)->SetRightMargin(0.01);
  d1->cd(1)->SetTopMargin(0.083);
  gEmpty->GetYaxis()->SetLabelSize(0.075);
  gEmpty->GetXaxis()->SetLabelSize(0.075);
  gEmpty->GetYaxis()->SetTitleSize(0.095);
  gEmpty->GetYaxis()->SetTitleOffset(0.99);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(4);
  gEmpty->GetXaxis()->SetNdivisions(0);
  gEmpty->GetYaxis()->SetTitle("#Delta#delta_{1}(J - PSJ)");
  gEmpty->GetYaxis()->SetRangeUser(-0.00053,0.00053); // delta D11
  //gEmpty->GetYaxis()->SetRangeUser(-0.0004,0.001); // delta G112
  gEmpty->DrawCopy();
  f1->Draw("same");
  JGraphDiffDeltaD11Xe_a_010->Draw("P");
  legend5->Draw();
  //============Delta gamma===============//
  d1->cd(2)->SetLeftMargin(0.19);
  d1->cd(2)->SetRightMargin(0.01);
  d1->cd(2)->SetBottomMargin(0.183);
  gEmpty->GetYaxis()->SetTitleOffset(0.99);
  gEmpty->GetYaxis()->SetLabelSize(0.07);
  gEmpty->GetXaxis()->SetLabelSize(0.07);
  gEmpty->GetYaxis()->SetTitleSize(0.09);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(3);
  gEmpty->GetXaxis()->SetNdivisions(10);
  gEmpty->GetYaxis()->SetRangeUser(-0.002,0.0001); // delta G112
  gEmpty->GetYaxis()->SetTitle("#Delta#gamma_{1,1}(J - PSJ)");
  gEmpty->DrawCopy();
  f1->Draw("same");
  JGraphDiffDeltaG112Xe_a_010->Draw("P");
  legend6->Draw();
  //d1->SaveAs("graphs/deltaJ_PSJ_AVFD.eps");
  //d1->SaveAs("graphs/deltaJ_PSJ_AVFD.png");

  //_____________________Ratios___________________________//
  TCanvas *d2 = new TCanvas("d2","Ratios: PSJ/J",0,0,700,800);
  d2->SetFillColor(10);
  d2->SetHighLightColor(10);
  d2->Divide(1,2,0.99,0.0,10);

  //===============Delta delta=============//
  d2->cd(1)->SetLeftMargin(0.19);
  d2->cd(1)->SetRightMargin(0.01);
  d2->cd(1)->SetTopMargin(0.083);
  gEmpty->GetYaxis()->SetLabelSize(0.075);
  gEmpty->GetXaxis()->SetLabelSize(0.075);
  gEmpty->GetYaxis()->SetTitleSize(0.095);
  gEmpty->GetYaxis()->SetTitleOffset(0.9);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(4);
  gEmpty->GetXaxis()->SetNdivisions(0);
  gEmpty->GetYaxis()->SetRangeUser(0.5,2); // delta D11 doesn't work... why ?
  //gEmpty->GetYaxis()->SetMinimum(0.5);
  //gEmpty->GetYaxis()->SetMaximum(1.8);
  gEmpty->GetYaxis()->SetTitle("#Delta#delta_{1}(PSJ/J)");
  gEmpty->DrawCopy();
  f1->Draw("same");
  JGraphRatioDeltaD11Xe_a_010->Draw("P");
  legend7->Draw();
  //============Delta gamma===============//
  d2->cd(2)->SetLeftMargin(0.19);
  d2->cd(2)->SetRightMargin(0.01);
  d2->cd(2)->SetBottomMargin(0.183);
  gEmpty->GetYaxis()->SetTitleOffset(0.85);
  gEmpty->GetYaxis()->SetLabelSize(0.07);
  gEmpty->GetXaxis()->SetLabelSize(0.07);
  gEmpty->GetYaxis()->SetTitleSize(0.09);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(5);
  gEmpty->GetXaxis()->SetNdivisions(10);
  gEmpty->GetYaxis()->SetRangeUser(1.5,3);
  //gEmpty->GetYaxis()->SetMinimum(1.5); // delta G112
  //gEmpty->GetYaxis()->SetMaximum(3);// No clue why the SetRangeUser doesn't work here -.-
  gEmpty->GetYaxis()->SetTitle("#Delta#gamma_{1,1}(PSJ/J)");
  gEmpty->GetYaxis()->SetRangeUser(1.5,3);
  gEmpty->DrawCopy();
  f1->Draw("same");
  JGraphRatioDeltaG112Xe_a_010->GetYaxis()->SetRangeUser(1.5,3);
  JGraphRatioDeltaG112Xe_a_010->Draw("P");
  legend8->Draw();

  //d2->SaveAs("graphs/ratioJ_PSJ_AVFD.eps");
  //d2->SaveAs("graphs/ratioJ_PSJ_AVFD.png");

  //____________________Delta gamma_______________________//
  TCanvas *c1 = new TCanvas("c1","Centrality dependence: Delta gamma",0,0,700,800);
  c1->SetFillColor(10); 	
  c1->SetHighLightColor(10);
  c1->Divide(1,2,0.99,0.0,10);

  //===============gamma vs LCC=============//
  c1->cd(1)->SetLeftMargin(0.19); 
  c1->cd(1)->SetRightMargin(0.01);
  c1->cd(1)->SetTopMargin(0.083);
  gEmpty->GetYaxis()->SetLabelSize(0.075);
  gEmpty->GetXaxis()->SetLabelSize(0.075);
  gEmpty->GetYaxis()->SetTitleSize(0.095);
  gEmpty->GetYaxis()->SetTitleOffset(0.9);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(4);
  gEmpty->GetXaxis()->SetNdivisions(0);
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha}-#varphi_{#beta}) #GT");
  gEmpty->GetYaxis()->SetTitle("#Delta#gamma_{1,1}");
  //gEmpty->GetYaxis()->SetRangeUser(-0.0006,0.0039); // delta D11
  gEmpty->GetYaxis()->SetRangeUser(-0.0004,0.001); // delta G112
  gEmpty->DrawCopy();

  f1->Draw("same");
  fGraphErrorsDeltaG112Xe_Baseline->Draw("P");
  fGraphErrorsDeltaG112Pb5_Baseline->Draw("P");
  fGraphErrorsDeltaG112Xe_LCC15->Draw("P3");
  fGraphErrorsDeltaG112Pb5_LCC15->Draw("P3");
  //fGraphErrorsDeltaG112Pb5_LCC33->Draw("P3");
  //fGraphErrorsDeltaG112Xe_LCC33->Draw("P3");
  fGraphErrorsDeltaG112Xe_LCC50->Draw("P3");
  fGraphErrorsDeltaG112Pb5_LCC50->Draw("P3");
	
  legend1->Draw();
	
  //============gamma vs n5/s===============//
  c1->cd(2)->SetLeftMargin(0.19); 	
  c1->cd(2)->SetRightMargin(0.01);
  c1->cd(2)->SetBottomMargin(0.183);
  gEmpty->GetYaxis()->SetTitleOffset(0.85);
  gEmpty->GetYaxis()->SetLabelSize(0.07);
  gEmpty->GetXaxis()->SetLabelSize(0.07);
  gEmpty->GetYaxis()->SetTitleSize(0.09);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(5);
  gEmpty->GetXaxis()->SetNdivisions(10);
  //gEmpty->GetYaxis()->SetRangeUser(-0.0048,0.0013); // delta D11
  gEmpty->GetYaxis()->SetRangeUser(-0.0005,0.0025); // delta G112
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha} + #varphi_{#beta} - 2#Psi_{RP}) #GT");
  gEmpty->GetYaxis()->SetTitle("#Delta#gamma_{1,1}");
  gEmpty->DrawCopy();
  f1->Draw("same");
//here mine:
  JGraphErrorsDeltaG112Xe_a_010->Draw("P3");
  fGraphErrorsDeltaG112Pb5_alpha005->Draw("P3");
  fGraphErrorsDeltaG112Xe_alpha005->Draw("P3");
  fGraphErrorsDeltaG112Pb5_alpha007->Draw("P3");
  fGraphErrorsDeltaG112Xe_alpha007->Draw("P3");
  fGraphErrorsDeltaG112Pb5_alpha010->Draw("P3");
  fGraphErrorsDeltaG112Xe_alpha010->Draw("P3");

  legend2->Draw();

  //c1->SaveAs("graphs/deltaGammaAVFD.eps");
  //c1->SaveAs("graphs/deltaGammaAVFD.png");
 
  //____________________Delta delta_______________________//
  TCanvas *c2 = new TCanvas("c2","Centrality dependence: Delta delta",100,100,700,800);
  c2->SetFillColor(10); 	
  c2->SetHighLightColor(10);
  c2->Divide(1,2,0.99,0.0,10);

  //===============gamma vs LCC=============//
  c2->cd(1)->SetLeftMargin(0.19); 
  c2->cd(1)->SetRightMargin(0.01);
  c2->cd(1)->SetTopMargin(0.083);
  gEmpty->GetYaxis()->SetLabelSize(0.075);
  gEmpty->GetXaxis()->SetLabelSize(0.075);
  gEmpty->GetYaxis()->SetTitleSize(0.095);
  gEmpty->GetYaxis()->SetTitleOffset(0.9);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(4);
  gEmpty->GetXaxis()->SetNdivisions(0);
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha}-#varphi_{#beta}) #GT");
  gEmpty->GetYaxis()->SetTitle("#Delta#delta_{1}");
  gEmpty->GetYaxis()->SetRangeUser(-0.0006,0.0039); // delta D11
  //gEmpty->GetYaxis()->SetRangeUser(-0.0004,0.001); // delta G112
  gEmpty->DrawCopy();

  f1->Draw("same");
  fGraphErrorsDeltaD11Xe_Baseline->Draw("P");
  fGraphErrorsDeltaD11Pb5_Baseline->Draw("P");
  fGraphErrorsDeltaD11Xe_LCC15->Draw("P3");
  fGraphErrorsDeltaD11Pb5_LCC15->Draw("P3");
  //fGraphErrorsDeltaD11Pb5_LCC33->Draw("P3");
  //fGraphErrorsDeltaD11Xe_LCC33->Draw("P3");
  fGraphErrorsDeltaD11Xe_LCC50->Draw("P3");
  fGraphErrorsDeltaD11Pb5_LCC50->Draw("P3");
	
  legend1->Draw();
	
  //============gamma vs n5/s===============//
  c2->cd(2)->SetLeftMargin(0.19); 	
  c2->cd(2)->SetRightMargin(0.01);
  c2->cd(2)->SetBottomMargin(0.183);
  gEmpty->GetYaxis()->SetTitleOffset(0.9);
  gEmpty->GetYaxis()->SetLabelSize(0.07);
  gEmpty->GetXaxis()->SetLabelSize(0.07);
  gEmpty->GetYaxis()->SetTitleSize(0.09);
  gEmpty->GetXaxis()->SetTitleSize(0.075);
  gEmpty->GetYaxis()->SetNdivisions(5);
  gEmpty->GetXaxis()->SetNdivisions(10);
  gEmpty->GetYaxis()->SetRangeUser(-0.0048,0.0013); // delta D11
  //gEmpty->GetYaxis()->SetRangeUser(-0.0005,0.0025); // delta G112
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha} + #varphi_{#beta} - 2#Psi_{RP}) #GT");
  gEmpty->GetYaxis()->SetTitle("#Delta#delta_{1}");
  gEmpty->DrawCopy();
  f1->Draw("same");
//here mine:
  JGraphErrorsDeltaD11Xe_a_010->Draw("P3");
  fGraphErrorsDeltaD11Pb5_alpha005->Draw("P3");
  fGraphErrorsDeltaD11Xe_alpha005->Draw("P3");
  fGraphErrorsDeltaD11Pb5_alpha007->Draw("P3");
  fGraphErrorsDeltaD11Xe_alpha007->Draw("P3");
  fGraphErrorsDeltaD11Pb5_alpha010->Draw("P3");
  fGraphErrorsDeltaD11Xe_alpha010->Draw("P3");

  legend2->Draw();

  //c2->SaveAs("graphs/deltaDeltaAVFD.eps");
  //c2->SaveAs("graphs/deltaDeltaAVFD.png");

  //____________________Delta gamma_______________________//
  TCanvas *c3 = new TCanvas("c3","Centrality dependence: Delta gamma",200,200,700,800);
  c3->SetFillColor(10); 	
  c3->SetHighLightColor(10);
  c3->Divide(1,2,0.99,0.0,10);

  //===============gamma vs LCC=============//
  c3->cd(1)->SetLeftMargin(0.19); 
  c3->cd(1)->SetRightMargin(0.01);
  c3->cd(1)->SetTopMargin(0.083);
  gEmpty0->GetYaxis()->SetLabelSize(0.075);
  gEmpty0->GetXaxis()->SetLabelSize(0.075);
  gEmpty0->GetYaxis()->SetTitleSize(0.095);
  gEmpty0->GetYaxis()->SetTitleOffset(0.9);
  gEmpty0->GetXaxis()->SetTitleSize(0.075);
  gEmpty0->GetYaxis()->SetNdivisions(4);
  gEmpty0->GetXaxis()->SetNdivisions(0);
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha}-#varphi_{#beta}) #GT");
  gEmpty0->GetYaxis()->SetTitle("#Delta#gamma_{1,1}");
  gEmpty0->GetXaxis()->SetRangeUser(13,990); 
  //gEmpty0->GetYaxis()->SetRangeUser(-0.0006,0.0039); // delta D11
  gEmpty0->GetYaxis()->SetRangeUser(-0.0004,0.001); // delta G112
  gEmpty0->DrawCopy();

  fGraphvsMultDeltaG112Pb5_Baseline->Draw("P");
  fGraphvsMultDeltaG112Xe_Baseline->Draw("P");
  fGraphvsMultDeltaG112Pb5_LCC15->Draw("P3");
  fGraphvsMultDeltaG112Xe_LCC15->Draw("P3");
  fGraphvsMultDeltaG112Pb5_LCC33->Draw("P3");
  fGraphvsMultDeltaG112Xe_LCC33->Draw("P3");
  fGraphvsMultDeltaG112Pb5_LCC50->Draw("P3");
  fGraphvsMultDeltaG112Xe_LCC50->Draw("P3");

  legend3->Draw();
   
  //============gamma vs n5/s===============//
  c3->cd(2)->SetLeftMargin(0.19); 	
  c3->cd(2)->SetRightMargin(0.01);
  c3->cd(2)->SetBottomMargin(0.183);
  gEmpty0->GetYaxis()->SetTitleOffset(0.9);
  gEmpty0->GetYaxis()->SetLabelSize(0.07);
  gEmpty0->GetXaxis()->SetLabelSize(0.07);
  gEmpty0->GetYaxis()->SetTitleSize(0.09);
  gEmpty0->GetXaxis()->SetTitleSize(0.075);
  gEmpty0->GetYaxis()->SetNdivisions(5);
  gEmpty0->GetXaxis()->SetNdivisions(5);
  gEmpty0->GetXaxis()->SetRangeUser(13,990); 
  gEmpty0->GetYaxis()->SetRangeUser(-0.0005,0.0025); // delta G112
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha} + #varphi_{#beta} - 2#Psi_{RP}) #GT");
  gEmpty0->GetYaxis()->SetTitle("#Delta#gamma_{1,1}");
  gEmpty0->DrawCopy();
  f1->Draw("same");
  //here better not mine :P
  //JGraphErrorsDeltaG112Xe_a_010->Draw("P3");
  fGraphvsMultDeltaG112Pb5_alpha005->Draw("P3");
  fGraphvsMultDeltaG112Xe_alpha005->Draw("P3");
  fGraphvsMultDeltaG112Pb5_alpha007->Draw("P3");
  fGraphvsMultDeltaG112Xe_alpha007->Draw("P3");
  fGraphvsMultDeltaG112Pb5_alpha010->Draw("P3");
  fGraphvsMultDeltaG112Xe_alpha010->Draw("P3");

  legend4->Draw();

  //c3->SaveAs("/graphs/deltaGammaAVFDVsNch.eps");
  //c3->SaveAs("/graphs/deltaGammaAVFDVsNch.png");

   //____________________Delta gamma_______________________//
  TCanvas *c4 = new TCanvas("c4","Centrality dependence: Delta delta",300,300,700,800);
  c4->SetFillColor(10); 	
  c4->SetHighLightColor(10);
  c4->Divide(1,2,0.99,0.0,10);

  //===============gamma vs LCC=============//
  c4->cd(1)->SetLeftMargin(0.19); 
  c4->cd(1)->SetRightMargin(0.01);
  c4->cd(1)->SetTopMargin(0.083);
  gEmpty0->GetYaxis()->SetLabelSize(0.075);
  gEmpty0->GetXaxis()->SetLabelSize(0.075);
  gEmpty0->GetYaxis()->SetTitleSize(0.095);
  gEmpty0->GetYaxis()->SetTitleOffset(0.9);
  gEmpty0->GetXaxis()->SetTitleSize(0.075);
  gEmpty0->GetYaxis()->SetNdivisions(4);
  gEmpty0->GetXaxis()->SetNdivisions(0);
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha}-#varphi_{#beta}) #GT");
  gEmpty0->GetYaxis()->SetTitle("#Delta#delta_{1}");
  gEmpty0->GetXaxis()->SetRangeUser(13,990); 
  gEmpty0->GetYaxis()->SetRangeUser(-0.0006,0.0039); // delta D11
  //gEmpty0->GetYaxis()->SetRangeUser(-0.0004,0.001); // delta G112
  gEmpty0->DrawCopy();

  fGraphvsMultDeltaD11Pb5_Baseline->Draw("P");
  fGraphvsMultDeltaD11Xe_Baseline->Draw("P");
  fGraphvsMultDeltaD11Pb5_LCC15->Draw("P3");
  fGraphvsMultDeltaD11Xe_LCC15->Draw("P3");
  fGraphvsMultDeltaD11Pb5_LCC33->Draw("P3");
  fGraphvsMultDeltaD11Xe_LCC33->Draw("P3");
  fGraphvsMultDeltaD11Pb5_LCC50->Draw("P3");
  fGraphvsMultDeltaD11Xe_LCC50->Draw("P3");

  legend3->Draw();
   
  //============gamma vs n5/s===============//
  c4->cd(2)->SetLeftMargin(0.19); 	
  c4->cd(2)->SetRightMargin(0.01);
  c4->cd(2)->SetBottomMargin(0.183);
  gEmpty0->GetYaxis()->SetTitleOffset(0.9);
  gEmpty0->GetYaxis()->SetLabelSize(0.07);
  gEmpty0->GetXaxis()->SetLabelSize(0.07);
  gEmpty0->GetYaxis()->SetTitleSize(0.09);
  gEmpty0->GetXaxis()->SetTitleSize(0.075);
  gEmpty0->GetYaxis()->SetNdivisions(5);
  gEmpty0->GetXaxis()->SetNdivisions(5);
  gEmpty0->GetXaxis()->SetRangeUser(13,990); 
  gEmpty0->GetYaxis()->SetRangeUser(-0.0048,0.0013); // delta D11
  //gEmpty->GetYaxis()->SetTitle("#LT cos(#varphi_{#alpha} + #varphi_{#beta} - 2#Psi_{RP}) #GT");
  gEmpty0->GetYaxis()->SetTitle("#Delta#delta_{1}");
  gEmpty0->DrawCopy();
  f1->Draw("same");
  //here better not mine :P
  //JGraphErrorsDeltaD11Xe_a_010->Draw("P3");  
  fGraphvsMultDeltaD11Pb5_alpha005->Draw("P3");
  fGraphvsMultDeltaD11Xe_alpha005->Draw("P3");
  fGraphvsMultDeltaD11Pb5_alpha007->Draw("P3");
  fGraphvsMultDeltaD11Xe_alpha007->Draw("P3");
  fGraphvsMultDeltaD11Pb5_alpha010->Draw("P3");
  fGraphvsMultDeltaD11Xe_alpha010->Draw("P3");

  legend4->Draw();

  //c4->SaveAs("/graphs/deltaDeltaAVFDVsNch.eps");
  //c4->SaveAs("/graphs/deltaDeltaAVFDVsNch.png");

}
