#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "Event.h"
#include "Particle.h"
#include "CalculateFlowCME.h"

using namespace std;

CalculateFlowCME::CalculateFlowCME () 
{
	std::cout<<"default CalculateFlowCME constructor"<<'\n';

}

//=====================================================================================================

CalculateFlowCME::CalculateFlowCME(const char* name):fQAList(NULL),fFlowQCList(NULL),fFlowGFList(NULL),fFlowFromBWList(NULL),fCMEList(NULL),fCMWList(NULL),fCMWQAList(NULL),fFlowQCCenBin(10),fReQGF(NULL),fImQGF(NULL)
{
	std::cout<<"CalculateFlowCME constructor"<<'\n';

	fFlowQCList = new TList();
    fFlowQCList->SetName("fFlowQCList");
    fFlowQCList->SetOwner(kTRUE);
    
    fFlowGFList = new TList();
    fFlowGFList->SetName("fFlowGFList");
    fFlowGFList->SetOwner(kTRUE);
    
    fFlowFromBWList = new TList();
    fFlowFromBWList->SetName("fFlowFromBWList");
    fFlowFromBWList->SetOwner(kTRUE);
    
    fCMEList = new TList();
    fCMEList->SetName("fCMEList");
    fCMEList->SetOwner(kTRUE);
    
    fCMWList = new TList();
    fCMWList->SetName("fCMWList");
    fCMWList->SetOwner(kTRUE);
    
    fCMWQAList = new TList();
    fCMWQAList->SetName("fCMWQAList");
    fCMWQAList->SetOwner(kTRUE);
	
	fQAList = new TList();
    fQAList->SetName("fQAList");
    fQAList->SetOwner(kTRUE);
	
	InitializeArraysForFlowQC();
	InitializeArraysForFlowGF();
	InitializeArraysForCMW();
	InitializeArraysForQA();
	InitializeArraysForFlowFromBW();
	InitializeArraysForCME();
	
	for(Int_t i=0; i<fkGFPtB; i++) {
		fReQGFPt[i] = NULL;
		fImQGFPt[i] = NULL;
	}
	
}

void CalculateFlowCME::InitializeArraysForQA()
{
	fnEvent = NULL;
	fMultChargedParticlesDistribution = NULL;
	fPtChargedParticlesDistribution = NULL;
	fEtaChargedParticlesDistribution = NULL;
	fPhiChargedParticlesDistribution = NULL;
	fEtaChargedParticlesDistributionPaperbinning = NULL;
	fPionsPtSpectra = NULL;
	fPionsEtaSpectra = NULL;
	fPionsPhiSpectra = NULL;
	fPosPionsPtSpectra = NULL;
	fPosPionsEtaSpectra = NULL;
	fPosPionsPhiSpectra = NULL;
	fAntiPionsPtSpectra = NULL;
	fAntiPionsEtaSpectra = NULL;
	fAntiPionsPhiSpectra = NULL;
	fKaonsPtSpectra = NULL;
	fKaonsEtaSpectra = NULL;
	fKaonsPhiSpectra = NULL;
	fPosKaonsPtSpectra = NULL;
	fPosKaonsEtaSpectra = NULL;
	fPosKaonsPhiSpectra = NULL;
	fAntiKaonsPtSpectra = NULL;
	fAntiKaonsEtaSpectra = NULL;
	fAntiKaonsPhiSpectra = NULL;
	fProtonsPtSpectra = NULL;
	fProtonsEtaSpectra = NULL;
	fProtonsPhiSpectra = NULL;
	fPosProtonsPtSpectra = NULL;
	fPosProtonsEtaSpectra = NULL;
	fPosProtonsPhiSpectra = NULL;
	fAntiProtonsPtSpectra = NULL;
	fAntiProtonsEtaSpectra = NULL;
	fAntiProtonsPhiSpectra = NULL;
}
	
void CalculateFlowCME::InitializeArraysForFlowQC()
{
 for(Int_t ch=0; ch < charge; ch++){
	for (Int_t c=0;c<fQVecPower;c++) {
		for (Int_t h=0;h<fFlowNHarmMax;h++) {
			fPOIPtDiffQRe[ch][c][h] = NULL;
			fPOIPtDiffQIm[ch][c][h] = NULL;
			fPOIPtDiffMul[ch][c][h] = NULL;
		}
	}
	
	for(Int_t i=0; i<fFlowNHarm; i++) {
		for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
			fFlowQCIntCorPro[ch][i][j] = NULL;
			fFlowQCIntCorHist[ch][i][j] = NULL;
			fFlowQCIntCumHist[ch][i][j] = NULL;
		}
	}
	
	// reference flow
	for(Int_t i=0; i<fFlowNHarm; i++) {
		for(Int_t j=0; j<fFlowQCNRef; j++) {
			fFlowQCRefCorPro[ch][i][j] = NULL;
			fFlowQCRefCorHist[ch][i][j] = NULL;
		}
		for(Int_t j=0; j<4; j++) {
			fFlowQCRefCorFinal[ch][i][j] = NULL;
		}
	}	
  
    // differential flow
	for (Int_t h=0; h<fCRCMaxnCen; h++) {
		for(Int_t i=0; i<fFlowNHarm; i++) {
			for(Int_t j=0; j<fFlowQCNPro; j++) {
				fFlowQCCorPro[ch][h][i][j] = NULL;
				fFlowQCCorHist[ch][h][i][j] = NULL;
			}
			for(Int_t k=0; k<fFlowQCNCov; k++) {
				fFlowQCCorCovPro[ch][h][i][k] = NULL;
				fFlowQCCorCovHist[ch][h][i][k] = NULL;
				fFlowQCFinalPtDifHist[ch][h][i][k] = NULL;
			}
		}
	}
 }
	//final hist for differential \Delta v_1
 for(Int_t h = 0; h<fCRCMaxnCen; h++){
 	for(Int_t j=0; j<fFlowNHarm; j++){
 		for(Int_t k=0; k<fFlowQCNCov; k++) {
			fFlowQCFinalPtDifDeltaHist[h][j][k]; 
		}
 	} 
 }
}
//=====================================================================================================
void CalculateFlowCME::InitializeArraysForFlowGF()
{
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNOrde; i++) {
      fFlowGFIntCorPro[h][i] = NULL;
      fFlowGFIntCorHist[h][i] = NULL;
      fFlowGFIntCumHist[h][i] = NULL;
      fFlowGFIntFinalHist[h][i] = NULL;
      for(Int_t k=0; k<fkFlowGFNOrde; k++) {
        fFlowGFIntCovPro[h][i][k] = NULL;
        fFlowGFIntCovHist[h][i][k] = NULL;
      }
      for(Int_t s=0; s<fkGFPtB; s++) {
        fFlowGFIntCorProPtB[s][h][i] = NULL;
        fFlowGFIntCorHistPtB[s][h][i] = NULL; 
       for(Int_t k=0; k<fkFlowGFNOrde; k++) {
          fFlowGFIntCovProPtB[s][h][i][k] = NULL;
          fFlowGFIntCovHistPtB[s][h][i][k] = NULL;
        }
      }
    }
  }
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNHarm; i++) {
      fFlowGFMixedCorPro[h][i] = NULL;
      fFlowGFMixedCorHist[h][i] = NULL;
      fFlowGFMixedFinalHist[h][i] = NULL;
    }
  }
 hMpn = NULL;

} // end of CalculateFlowCME::InitializeArraysForFlowGF()

void CalculateFlowCME::InitializeArraysForCMW()
{
  for(int i=0;i<2;i++){
    for(int j=0;j<9;j++){
      fHistv2AchChrgPosEtaNeg[i][j] = NULL;
      fHistv2AchChrgPosEtaPos[i][j] = NULL;
      //fHistv2AchChrgNegEtaNeg[i][j] = NULL;
      fHistv2AchChrgNegEtaPos[i][j] = NULL;
      
      fHistv2AchChrgPosChrgNeg[i][j] = NULL;
      fHistv2AchChrgNegChrgPos[i][j] = NULL;
    }
  }

  for(int i=0; i<9; i++){
    fHistEPResolutionAch[i] = NULL;
  }

  fHistAChrgVsCent = NULL;
  fHistEPResolution = NULL;

  for(int i=0; i<9; i++){
    fHistv2cumAchChrgAll[i] = NULL;
  }
  
  fHistAveV2PosAch = NULL;
  fHistAveV2Pos = NULL;
  fHistAveV2NegAch = NULL;
  fHistAveV2Neg = NULL;
 
  //fHistAveV1Pos = NULL;
  //fHistAveV1Neg = NULL;  
  // <v2^{+}A>-<A><v2> = <<v2^{+}A>>/sqrt(<<v2^{+}>>)-sqrt(<<v2^{+}>>)<A>
  fCMWThreeParticleCorrelator[0] = NULL;
  fCMWThreeParticleCorrelator[1] = NULL;
  
  // CMW QA hists
  fEvPlTPC = NULL;
}
//=====================================================================================================

void CalculateFlowCME::InitializeArraysForFlowFromBW() 
{
  for(int h=0; h<nHar; h++) {
    V2IntPro[h] = NULL;
    V2IntProQC[h] = NULL;
  for(int i=0; i<10; i++) {
      V2PtDiffPro[h][i] = NULL;
	 }
  }
  
  for(int h=0; h<nHar; h++) {
    hRepn[h] = NULL;
    hImpn[h] = NULL;
  } 
}

//=====================================================================================================

void CalculateFlowCME::InitializeArraysForCME()
{
  //--------- individual Q-vector terms -----------
  for(Int_t i=0; i<nQVector; i++) {
    fReQVectorPro[i]  = NULL;
    fImQVectorPro[i]  = NULL;
  }
  
  //------- 2particle correlator -------------
  for(Int_t i=0; i<nDnncor; i++) {
    fDnnPro[i] = NULL;
    fDnnHist[i] = NULL;
  }
  
  //--------- 3 particle correlator -----------
  for(Int_t i=0; i<nCMEcor; i++) {
    fCMEPro[i] = NULL;
    fCMEHist[i] = NULL;
  }
  
  DeltaD11 = NULL;
  DeltaG112 = NULL;
  DeltaG132 = NULL;
  DeltaG123 = NULL;
  DeltaG224 = NULL;
}

//=====================================================================================================

void CalculateFlowCME::UserCreateOutputObjects() {
	fReQGF = new TMatrixD(21,9);
	fImQGF = new TMatrixD(21,9);
	for(Int_t i=0; i<fkGFPtB; i++) {
		fReQGFPt[i] = new TMatrixD(21,9);
		fImQGFPt[i] = new TMatrixD(21,9);
	}
 
    Double_t pTbinEdge[59] = {0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20};
    
    
    Double_t pTPionbinEdge[42] = {0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3};
    
	Double_t pTKaonbinEdge[37] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3};
	
	Double_t pTProtonbinEdge[43] = {0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6};

	Double_t etabinEdge[51] = {-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.25,2.5,2.75,3,3.25,3.5};
	
    // QA Hists
    fnEvent = new TH1D("fnEvent", "fnEvent", 1, 0, 1);
    fMultChargedParticlesDistribution = new TH1D("fMultChargedParticlesDistribution", "fMultChargedParticlesDistribution", 200, 0, 4000);
	//fPtChargedParticlesDistribution = new TH1D("fPtChargedParticlesDistribution", "fPtChargedParticlesDistribution", 500, 0, 20);
	fPtChargedParticlesDistribution = new TH1D("fPtChargedParticlesDistribution", "fPtChargedParticlesDistribution", 58, pTbinEdge);
	fEtaChargedParticlesDistribution = new TH1D("fEtaChargedParticlesDistribution", "fEtaChargedParticlesDistribution", 200, -5, 5);
	fPhiChargedParticlesDistribution = new TH1D("fPhiChargedParticlesDistribution", "fPhiChargedParticlesDistribution", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	
	fEtaChargedParticlesDistributionPaperbinning = new TH1D("fEtaChargedParticlesDistributionPaperbinning", "fEtaChargedParticlesDistributionPaperbinning", 50, etabinEdge);
	
	//fPionsPtSpectra = new TH1D("fPionsPtSpectra", "fPionsPtSpectra", 200, 0, 2);
	fPionsPtSpectra = new TH1D("fPionsPtSpectra", "fPionsPtSpectra", 41, pTPionbinEdge);
	fPionsPtSpectra->Sumw2();
	fPionsEtaSpectra = new TH1D("fPionsEtaSpectra", "fPionsEtaSpectra", 200, -5, 5);
	fPionsEtaSpectra->Sumw2();
	fPionsPhiSpectra = new TH1D("fPionsPhiSpectra", "fPionsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fPionsPhiSpectra->Sumw2();
	
	fPosPionsPtSpectra = new TH1D("fPosPionsPtSpectra", "fPosPionsPtSpectra", 41, pTPionbinEdge);
	fPosPionsPtSpectra->Sumw2();
	fPosPionsEtaSpectra = new TH1D("fPosPionsEtaSpectra", "fPosPionsEtaSpectra", 200, -5, 5);
	fPosPionsEtaSpectra->Sumw2();
	fPosPionsPhiSpectra = new TH1D("fPosPionsPhiSpectra", "fPosPionsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fPosPionsPhiSpectra->Sumw2();
	
	fAntiPionsPtSpectra = new TH1D("fAntiPionsPtSpectra", "fAntiPionsPtSpectra", 41, pTPionbinEdge);
	fAntiPionsPtSpectra->Sumw2();
	fAntiPionsEtaSpectra = new TH1D("fAntiPionsEtaSpectra", "fAntiPionsEtaSpectra", 200, -5, 5);
	fAntiPionsEtaSpectra->Sumw2();
	fAntiPionsPhiSpectra = new TH1D("fAntiPionsPhiSpectra", "fAntiPionsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fAntiPionsPhiSpectra->Sumw2();

	//fKaonsPtSpectra = new TH1D("fKaonsPtSpectra", "fKaonsPtSpectra", 300, 0, 3);
	fKaonsPtSpectra = new TH1D("fKaonsPtSpectra", "fKaonsPtSpectra", 36, pTKaonbinEdge);
	fKaonsPtSpectra->Sumw2();
	fKaonsEtaSpectra = new TH1D("fKaonsEtaSpectra", "fKaonsEtaSpectra", 200, -5, 5);
	fKaonsEtaSpectra->Sumw2();
	fKaonsPhiSpectra = new TH1D("fKaonsPhiSpectra", "fKaonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fKaonsPhiSpectra->Sumw2();
	
	fPosKaonsPtSpectra = new TH1D("fPosKaonsPtSpectra", "fPosKaonsPtSpectra", 36, pTKaonbinEdge);
	fPosKaonsPtSpectra->Sumw2();
	fPosKaonsEtaSpectra = new TH1D("fPosKaonsEtaSpectra", "fPosKaonsEtaSpectra", 200, -5, 5);
	fPosKaonsEtaSpectra->Sumw2();
	fPosKaonsPhiSpectra = new TH1D("fPosKaonsPhiSpectra", "fPosKaonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fPosKaonsPhiSpectra->Sumw2();
	
	fAntiKaonsPtSpectra = new TH1D("fAntiKaonsPtSpectra", "fAntiKaonsPtSpectra", 36, pTKaonbinEdge);
	fAntiKaonsPtSpectra->Sumw2();
	fAntiKaonsEtaSpectra = new TH1D("fAntiKaonsEtaSpectra", "fAntiKaonsEtaSpectra", 200, -5, 5);
	fAntiKaonsEtaSpectra->Sumw2();
	fAntiKaonsPhiSpectra = new TH1D("fAntiKaonsPhiSpectra", "fAntiKaonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fAntiKaonsPhiSpectra->Sumw2();

	//fProtonsPtSpectra = new TH1D("fProtonsPtSpectra", "fProtonsPtSpectra", 400, 0, 4);
	fProtonsPtSpectra = new TH1D("fProtonsPtSpectra", "fProtonsPtSpectra", 42, pTProtonbinEdge);
	fProtonsPtSpectra->Sumw2();
	fProtonsEtaSpectra = new TH1D("fProtonsEtaSpectra", "fProtonsEtaSpectra", 200, -5, 5);
	fProtonsEtaSpectra->Sumw2();
	fProtonsPhiSpectra = new TH1D("fProtonsPhiSpectra", "fProtonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fProtonsPhiSpectra->Sumw2();
	
	fPosProtonsPtSpectra = new TH1D("fPosProtonsPtSpectra", "fPosProtonsPtSpectra", 42, pTProtonbinEdge);
	fPosProtonsPtSpectra->Sumw2();
	fPosProtonsEtaSpectra = new TH1D("fPosProtonsEtaSpectra", "fPosProtonsEtaSpectra", 200, -5, 5);
	fPosProtonsEtaSpectra->Sumw2();
	fPosProtonsPhiSpectra = new TH1D("fPosProtonsPhiSpectra", "fPosProtonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fPosProtonsPhiSpectra->Sumw2();
	
	fAntiProtonsPtSpectra = new TH1D("fAntiProtonsPtSpectra", "fAntiProtonsPtSpectra", 42, pTProtonbinEdge);
	fAntiProtonsPtSpectra->Sumw2();
	fAntiProtonsEtaSpectra = new TH1D("fAntiProtonsEtaSpectra", "fAntiProtonsEtaSpectra", 200, -5, 5);
	fAntiProtonsEtaSpectra->Sumw2();
	fAntiProtonsPhiSpectra = new TH1D("fAntiProtonsPhiSpectra", "fAntiProtonsPhiSpectra", 200, -TMath::Pi()-0.1, TMath::Pi()+0.1);
	fAntiProtonsPhiSpectra->Sumw2();

	fQAList->Add(fnEvent);
	fQAList->Add(fMultChargedParticlesDistribution);
	fQAList->Add(fPtChargedParticlesDistribution);
	fQAList->Add(fEtaChargedParticlesDistribution);
	fQAList->Add(fPhiChargedParticlesDistribution);
	fQAList->Add(fEtaChargedParticlesDistributionPaperbinning);
	fQAList->Add(fPionsPtSpectra);
	fQAList->Add(fPionsEtaSpectra);
	fQAList->Add(fPionsPhiSpectra);
	fQAList->Add(fPosPionsPtSpectra);
	fQAList->Add(fPosPionsEtaSpectra);
	fQAList->Add(fPosPionsPhiSpectra);
	fQAList->Add(fAntiPionsPtSpectra);
	fQAList->Add(fAntiPionsEtaSpectra);
	fQAList->Add(fAntiPionsPhiSpectra);
	fQAList->Add(fKaonsPtSpectra);
	fQAList->Add(fKaonsEtaSpectra);
	fQAList->Add(fKaonsPhiSpectra);
	fQAList->Add(fPosKaonsPtSpectra);
	fQAList->Add(fPosKaonsEtaSpectra);
	fQAList->Add(fPosKaonsPhiSpectra);
	fQAList->Add(fAntiKaonsPtSpectra);
	fQAList->Add(fAntiKaonsEtaSpectra);
	fQAList->Add(fAntiKaonsPhiSpectra);
	fQAList->Add(fProtonsPtSpectra);
	fQAList->Add(fProtonsEtaSpectra);
	fQAList->Add(fProtonsPhiSpectra);
	fQAList->Add(fPosProtonsPtSpectra);
	fQAList->Add(fPosProtonsEtaSpectra);
	fQAList->Add(fPosProtonsPhiSpectra);
	fQAList->Add(fAntiProtonsPtSpectra);
	fQAList->Add(fAntiProtonsEtaSpectra);
	fQAList->Add(fAntiProtonsPhiSpectra);
	
	// Calculate FlowQC Hists
	fPtDiffNBins = 36;
	fCRCPtBins = new Double_t[37];
	Double_t PtBins[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,17.,20.,25.,30.,40.,50.};
	for(Int_t r=0; r<37; r++) {
		fCRCPtBins[r] = PtBins[r];
	}
	//To get eta for the loop
	fEtaDiffNBins = 50;
        fCRCEtaBins = new Double_t[51];
	for(Int_t i=0; i<51; i++){
		fCRCEtaBins[i] = etabinEdge[i];
//cout << "fill fCRCEtaBins: " << fCRCEtaBins[i] << endl;
	}

  for(Int_t ch = 0; ch<charge; ch++){
	for (Int_t c=0;c<fQVecPower;c++) {
		for (Int_t h=0;h<fFlowNHarmMax;h++) {
		    fPOIPtDiffQRe[ch][c][h] = new TProfile(Form("fPOIPtDiffQRe[%d][%d][%d]",ch,c,h),Form("fPOIPtDiffQRe[%d][%d][%d]",ch,c,h),fPtDiffNBins,fCRCPtBins,"s");
		    fFlowQCList->Add(fPOIPtDiffQRe[ch][c][h]);
		    fPOIPtDiffQIm[ch][c][h] = new TProfile(Form("fPOIPtDiffQIm[%d][%d][%d]",ch,c,h),Form("fPOIPtDiffQIm[%d][%d][%d]",ch,c,h),fPtDiffNBins,fCRCPtBins,"s");
		    fFlowQCList->Add(fPOIPtDiffQIm[ch][c][h]);
		    fPOIPtDiffMul[ch][c][h] = new TProfile(Form("fPOIPtDiffMul[%d][%d][%d]",ch,c,h),Form("fPOIPtDiffMul[%d][%d][%d]",ch,c,h),fPtDiffNBins,fCRCPtBins, "s");
		    fFlowQCList->Add(fPOIPtDiffMul[ch][c][h]);
		}
	}

	for(Int_t i=0; i<fFlowNHarm; i++) {
		for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
			fFlowQCIntCorPro[ch][i][j] = new TProfile(Form("fFlowQCIntCorPro[%d][%d][%d]",ch,i,j),Form("fFlowQCIntCorPro[%d][%d][%d]",ch,i,j),fFlowQCCenBin,0.,100.,"s");
			fFlowQCIntCorPro[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCIntCorPro[ch][i][j]);
			fFlowQCIntCorHist[ch][i][j] = new TH1D(Form("fFlowQCIntCorHist[%d][%d][%d]",ch,i,j),Form("fFlowQCIntCorHist[%d][%d][%d]",ch,i,j),fFlowQCCenBin,0.,100.);
			fFlowQCIntCorHist[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCIntCorHist[ch][i][j]);
			fFlowQCIntCumHist[ch][i][j] = new TH1D(Form("fFlowQCIntCumHist[%d][%d][%d]",ch,i,j),Form("fFlowQCIntCumHist[%d][%d][%d]",ch,i,j),fFlowQCCenBin,0.,100.);
			fFlowQCIntCumHist[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCIntCumHist[ch][i][j]);
		}
	}
	
	
	// reference flow
	Double_t xbins[12] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
	for(Int_t i=0; i<fFlowNHarm; i++) {
		for(Int_t j=0; j<fFlowQCNRef; j++) {
			fFlowQCRefCorPro[ch][i][j] = new TProfile(Form("fFlowQCRefCorPro[[%d]%d][%d]",ch,i,j),Form("fFlowQCRefCorPro[%d][%d][%d]",ch,i,j),11,xbins,"s");
			fFlowQCRefCorPro[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCRefCorPro[ch][i][j]);
			fFlowQCRefCorHist[ch][i][j] = new TH1D(Form("fFlowQCRefCorHist[%d][%d][%d]",ch,i,j),Form("fFlowQCRefCorHist[%d][%d][%d]",ch,i,j),11,xbins);
			fFlowQCRefCorHist[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCRefCorHist[ch][i][j]);
		}
		for(Int_t j=0; j<4; j++) {
			fFlowQCRefCorFinal[ch][i][j] = new TH1D(Form("fFlowQCRefCorFinal[%d][%d][%d]",ch,i,j),Form("fFlowQCRefCorFinal[%d][%d][%d]",ch,i,j),11,xbins);
			fFlowQCRefCorFinal[ch][i][j]->Sumw2();
			fFlowQCList->Add(fFlowQCRefCorFinal[ch][i][j]);
		}
	}

    // differential flow
	for (Int_t h=0; h<fCRCnCen; h++) {
		for(Int_t i=0; i<fFlowNHarm; i++) {
			for(Int_t j=0; j<fFlowQCNPro; j++) {
				fFlowQCCorPro[ch][h][i][j] = new TProfile(Form("fFlowQCCorPro[%d][%d][%d][%d]",ch,h,i,j),Form("fFlowQCCorPro[%d][%d][%d][%d]",ch,h,i,j),fPtDiffNBins,fCRCPtBins,"s");
				fFlowQCCorPro[ch][h][i][j]->Sumw2();
				fFlowQCList->Add(fFlowQCCorPro[ch][h][i][j]);
				fFlowQCCorHist[ch][h][i][j] = new TH1D(Form("fFlowQCCorHist[%d][%d][%d][%d]",ch,h,i,j),Form("fFlowQCCorHist[%d][%d][%d][%d]",ch,h,i,j),fPtDiffNBins,fCRCPtBins);
				fFlowQCCorHist[ch][h][i][j]->Sumw2();
				fFlowQCList->Add(fFlowQCCorHist[ch][h][i][j]);
			}
			for(Int_t k=0; k<fFlowQCNCov; k++) {
				fFlowQCCorCovPro[ch][h][i][k] = new TProfile(Form("fFlowQCCorCovPro[%d][%d][%d][%d]",ch,h,i,k),Form("fFlowQCCorCovPro[%d][%d][%d][%d]",ch,h,i,k),fPtDiffNBins,fCRCPtBins,"s");
				fFlowQCCorCovPro[ch][h][i][k]->Sumw2();
				fFlowQCList->Add(fFlowQCCorCovPro[ch][h][i][k]);
				fFlowQCCorCovHist[ch][h][i][k] = new TH1D(Form("fFlowQCCorCovHist[%d][%d][%d][%d]",ch,h,i,k),Form("fFlowQCCorCovHist[%d][%d][%d][%d]",ch,h,i,k),fPtDiffNBins,fCRCPtBins);
				fFlowQCCorCovHist[ch][h][i][k]->Sumw2();
				fFlowQCList->Add(fFlowQCCorCovHist[ch][h][i][k]);
				fFlowQCFinalPtDifHist[ch][h][i][k] = new TH1D(Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d]",ch,h,i,k),Form("fFlowQCFinalPtDifHist[%d][%d][%d][%d]",ch,h,i,k),fPtDiffNBins,fCRCPtBins);
				fFlowQCFinalPtDifHist[ch][h][i][k]->Sumw2();
				fFlowQCList->Add(fFlowQCFinalPtDifHist[ch][h][i][k]);
			}
		}
	}
  }
	//here the delta v1 hist
 for(Int_t h=0; h<fCRCnCen; h++) {
          for(Int_t i=0; i<fFlowNHarm; i++) {
                for(Int_t k=0; k<fFlowQCNCov; k++) {
			fFlowQCFinalPtDifDeltaHist[h][i][k] = new TH1D(Form("fFlowQCFinalPtDifDeltaHist[%d][%d][%d]",h,i,k),Form("fFlowQCFinalPtDifDeltaHist[%d][%d][%d]",h,i,k),fPtDiffNBins,fCRCPtBins);
                        fFlowQCFinalPtDifDeltaHist[h][i][k]->Sumw2();
                        fFlowQCList->Add(fFlowQCFinalPtDifDeltaHist[h][i][k]);
			}
		}
	}	

  // for CalculateFlowGF
  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNOrde; i++) {
      fFlowGFIntCorPro[h][i] = new TProfile(Form("fFlowGFIntCorPro[%d][%d]",h,i),Form("fFlowGFIntCorPro[%d][%d]",h,i),fFlowGFCenBin,0.,100.,"s");
      fFlowGFIntCorPro[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFIntCorPro[h][i]);
      fFlowGFIntCorHist[h][i] = new TH1D(Form("fFlowGFIntCorHist[%d][%d]",h,i),Form("fFlowGFIntCorHist[%d][%d]",h,i),fFlowGFCenBin,0.,100.);
      fFlowGFIntCorHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFIntCorHist[h][i]);
      fFlowGFIntCumHist[h][i] = new TH1D(Form("fFlowGFIntCumHist[%d][%d]",h,i),Form("fFlowGFIntCumHist[%d][%d]",h,i),fFlowGFCenBin,0.,100.);
      fFlowGFIntCumHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFIntCumHist[h][i]);
      fFlowGFIntFinalHist[h][i] = new TH1D(Form("fFlowGFIntFinalHist[%d][%d]",h,i),Form("fFlowGFIntFinalHist[%d][%d]",h,i),fFlowGFCenBin,0.,100.);
      fFlowGFIntFinalHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFIntFinalHist[h][i]);

      for(Int_t k=0; k<fkFlowGFNOrde; k++) {
        fFlowGFIntCovPro[h][i][k] = new TProfile(Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovPro[%d][%d][%d]",h,i,k),fFlowGFCenBin,0.,100.,"s");
        fFlowGFIntCovPro[h][i][k]->Sumw2();
        fFlowGFList->Add(fFlowGFIntCovPro[h][i][k]);
        fFlowGFIntCovHist[h][i][k] = new TH1D(Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),Form("fFlowGFIntCovHist[%d][%d][%d]",h,i,k),fFlowGFCenBin,0.,100.);
        fFlowGFIntCovHist[h][i][k]->Sumw2();
        fFlowGFList->Add(fFlowGFIntCovHist[h][i][k]);
//for Pt
      }
      for(Int_t s=0; s<fkGFPtB; s++) {
        fFlowGFIntCorProPtB[s][h][i] = new TProfile(Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorProPtB[%d][%d][%d]",s,h,i),fFlowGFCenBin,0.,100.,"s");
        fFlowGFIntCorProPtB[s][h][i]->Sumw2();
        fFlowGFList->Add(fFlowGFIntCorProPtB[s][h][i]);
        fFlowGFIntCorHistPtB[s][h][i] = new TH1D(Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),Form("fFlowGFIntCorHistPtB[%d][%d][%d]",s,h,i),fFlowGFCenBin,0.,100.);
        fFlowGFIntCorHistPtB[s][h][i]->Sumw2();
        fFlowGFList->Add(fFlowGFIntCorHistPtB[s][h][i]);
        for(Int_t k=0; k<fkFlowGFNOrde; k++) {
          fFlowGFIntCovProPtB[s][h][i][k] = new TProfile(Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovProPtB[%d][%d][%d][%d]",s,h,i,k),fFlowGFCenBin,0.,100.,"s");
          fFlowGFIntCovProPtB[s][h][i][k]->Sumw2();
          fFlowGFList->Add(fFlowGFIntCovProPtB[s][h][i][k]);
          fFlowGFIntCovHistPtB[s][h][i][k] = new TH1D(Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),Form("fFlowGFIntCovHistPtB[%d][%d][%d][%d]",s,h,i,k),fFlowGFCenBin,0.,100.);
          fFlowGFIntCovHistPtB[s][h][i][k]->Sumw2();
          fFlowGFList->Add(fFlowGFIntCovHistPtB[s][h][i][k]);
        }
      }
    }
  }

  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNHarm; i++) {
      fFlowGFMixedCorPro[h][i] = new TProfile(Form("fFlowGFMixedCorPro[%d][%d]",h,i),Form("fFlowGFMixedCorPro[%d][%d]",h,i),fFlowGFCenBin,0.,100.,"s");
      fFlowGFMixedCorPro[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFMixedCorPro[h][i]);
      fFlowGFMixedCorHist[h][i] = new TH1D(Form("fFlowGFMixedCorHist[%d][%d]",h,i),Form("fFlowGFMixedCorHist[%d][%d]",h,i),fFlowGFCenBin,0.,100.);
      fFlowGFMixedCorHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFMixedCorHist[h][i]);
      fFlowGFMixedFinalHist[h][i] = new TH1D(Form("fFlowGFMixedFinalHist[%d][%d]",h,i),Form("fFlowGFMixedFinalHist[%d][%d]",h,i),fFlowGFCenBin,0.,100.);
      fFlowGFMixedFinalHist[h][i]->Sumw2();
      fFlowGFList->Add(fFlowGFMixedFinalHist[h][i]);
    }
  }
  hMpn = new TH1F("hMpn","M(p,n)", fPtDiffNBins,fCRCPtBins);                                          
  fFlowGFList->Add(hMpn);                                                                         
  
  // CMW histograms
  Char_t name[100];
  Char_t title[100];
 
  Double_t centRange[10] = {0,10,20,30,40,50,60,70,80,90};

  fHistAChrgVsCent = new TProfile("fHistAChrgVsCent","Ach vs Cent",10,0,100);
  fCMWList->Add(fHistAChrgVsCent);

  fHistAveV2PosAch = new TProfile("fHistAveV2PosAch","<<v_{2}^{+}A>>",10,0,100);
  fCMWList->Add(fHistAveV2PosAch);
  
  fHistAveV2Pos = new TProfile("fHistAveV2Pos","<<v_{2}^{+}>>",10,0,100);
  fCMWList->Add(fHistAveV2Pos);
          
  fHistAveV2NegAch = new TProfile("fHistAveV2NegAch"," <<v_{2}^{-}A>>",10,0,100);
  fCMWList->Add(fHistAveV2NegAch);
  
  fHistAveV2Neg = new TProfile("fHistAveV2Neg","<<v_{2}^{-}>>",10,0,100);
  fCMWList->Add(fHistAveV2Neg);
  
  fv2plusminus = new TProfile("fv2plusminus", "<<v_{2}^{-}>> and <<v_{2}^{+}>>", 2, 0, 2);
  fCMWList->Add(fv2plusminus);

  // <v2^{+}A>-<A><v2> = <<v2^{+}A>>/sqrt(<<v2^{+}>>)-sqrt(<<v2^{+}>>)<A>
  
  fCMWThreeParticleCorrelator[0] = new TProfile("fCMWThreeParticleCorrelator[0]","<v_{2}^{+}A>-<A><v_{2}^{+}>",10,0,100);
  fCMWList->Add(fCMWThreeParticleCorrelator[0]);
  
  fCMWThreeParticleCorrelator[1] = new TProfile("fCMWThreeParticleCorrelator[1]","<v_{2}^{-}A>-<A><v_{2}^{-}>",10,0,100);
  fCMWList->Add(fCMWThreeParticleCorrelator[1]);

  // v2 vs Ach
  for(int i=0;i<2;i++){
    for(int j=0;j<9;j++){
      ////Charge:
      sprintf(name,"fHistv2AchChrgPosEtaNeg_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgPosEtaNeg[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
      fHistv2AchChrgPosEtaNeg[i][j]->Sumw2();
      fCMWList->Add(fHistv2AchChrgPosEtaNeg[i][j]);
      //sprintf(name,"fHistv2AchChrgNeg_Method%d_Cent%d",i,j);
      //sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
      //fHistv2AchChrgNegEtaNeg[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
      //fHistv2AchChrgNegEtaNeg[i][j]->Sumw2();
      //fCMWList->Add(fHistv2AchChrgNegEtaNeg[i][j]);
      sprintf(name,"fHistv2AchChrgPosEtaPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgPosEtaPos[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
      fHistv2AchChrgPosEtaPos[i][j]->Sumw2();
      fCMWList->Add(fHistv2AchChrgPosEtaPos[i][j]);
      sprintf(name,"fHistv2AchChrgNegEtaPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgNegEtaPos[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
      fHistv2AchChrgNegEtaPos[i][j]->Sumw2();
      fCMWList->Add(fHistv2AchChrgNegEtaPos[i][j]);
      
      sprintf(name,"fHistv2AchChrgPosChrgNeg_Method%d_Cent%d",i,j);
	  sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
	  fHistv2AchChrgPosChrgNeg[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
	  fHistv2AchChrgPosChrgNeg[i][j]->Sumw2();
	  fCMWList->Add(fHistv2AchChrgPosChrgNeg[i][j]);
	  sprintf(name,"fHistv2AchChrgNegChrgPos_Method%d_Cent%d",i,j);
	  sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[j],centRange[j+1]);
	  fHistv2AchChrgNegChrgPos[i][j] = new TProfile(name,title,10,-0.1,0.1,"");
	  fHistv2AchChrgNegChrgPos[i][j]->Sumw2();
	  fCMWList->Add(fHistv2AchChrgNegChrgPos[i][j]);
    }
  }

  for(int i=0; i<9; i++){
    ////Charge:
    sprintf(name,"fHistResolutionvsAch_Cent%d",i);
    sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; Resolution",centRange[i],centRange[i+1]);
    fHistEPResolutionAch[i] = new TProfile(name,title,10,-0.1,0.1,"");
    fHistEPResolutionAch[i]->Sumw2();
    fCMWList->Add(fHistEPResolutionAch[i]);
  }

  sprintf(name,"fHistEPResolution_Method");
  sprintf(title,"EP Resolution; centrality; Reso");
  fHistEPResolution = new TProfile(name,title,9,centRange,"");
  fHistEPResolution->Sumw2();
  fCMWList->Add(fHistEPResolution);

  for(int i=0; i<9; i++){
    ////Charge:
    sprintf(name,"fHistv2cumAchChrgAllQcumCent%d",i);
    sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; v_{2}",centRange[i],centRange[i+1]);
    fHistv2cumAchChrgAll[i] = new TProfile(name,title,10,-0.1,0.1,"");
    fHistv2cumAchChrgAll[i]->Sumw2();
    fCMWList->Add(fHistv2cumAchChrgAll[i]);
  }

  //CMW QA hists
  fEvPlTPC = new TH1D("fEvPlTPC", "EvPlTPC distribution", 200, -3, 3);
  fCMWQAList->Add(fEvPlTPC);

  //Flow from BW hists 
  for(int h=0; h<nHar; h++) {
    V2IntPro[h] = new TProfile(Form("V%dIntPro",h+1),Form("V%dIntPro",h+1),10,0,100,"s");
    V2IntPro[h]->Sumw2();
    fFlowFromBWList->Add(V2IntPro[h]);
    V2IntProQC[h] = new TProfile(Form("V%dIntProQC",h+1),Form("V%dIntProQC",h+1),10,0,100,"s");
    V2IntProQC[h]->Sumw2();
    fFlowFromBWList->Add(V2IntProQC[h]);
    for(int i=0; i<10; i++) { // i index is the centrality
      V2PtDiffPro[h][i] = new TProfile(Form("V%dPtDiffPro_Cent%d",h+1,i),Form("V%dPtDiffPro_Cent%d",h+1,i),fPtDiffNBins,fCRCPtBins,"s");
      V2PtDiffPro[h][i]->Sumw2();
      fFlowFromBWList->Add(V2PtDiffPro[h][i]);
    }
  }
  hMpn = new TH1F("hMpn","M(p,n)", fPtDiffNBins,fCRCPtBins);
  fFlowFromBWList->Add(hMpn);
  
  for(int h=0; h<nHar; h++) {
    hRepn[h] = new TH1F(Form("hRep%d",h),Form("Rep%d",h),fPtDiffNBins,fCRCPtBins);
    fFlowFromBWList->Add(hRepn[h]);
    hImpn[h] = new TH1F(Form("hImp%d",h),Form("Imp%d",h),fPtDiffNBins,fCRCPtBins);
    fFlowFromBWList->Add(hImpn[h]);
  }
  //CME hists
  //--------- individual Q-vector terms -----------
  for(Int_t i=0; i<nQVector/2; i++) {
    fReQVectorPro[2*i]  = new TProfile(Form("fReQ%dPPro",i+1),Form("fReQ%dPPro",i+1),10,0,100,"s");
    fReQVectorPro[2*i]->Sumw2();
    fCMEList->Add(fReQVectorPro[2*i]);
    fReQVectorPro[2*i+1]  = new TProfile(Form("fReQ%dNPro",i+1),Form("fReQ%dNPro",i+1),10,0,100,"s");
    fReQVectorPro[2*i+1]->Sumw2();
    fCMEList->Add(fReQVectorPro[2*i+1]);

    fImQVectorPro[2*i]  = new TProfile(Form("fImQ%dPPro",i+1),Form("fImQ%dPPro",i+1),10,0,100,"s");
    fImQVectorPro[2*i]->Sumw2();
    fCMEList->Add(fImQVectorPro[2*i]);
    fImQVectorPro[2*i+1]  = new TProfile(Form("fImQ%dNPro",i+1),Form("fImQ%dNPro",i+1),10,0,100,"s");
    fImQVectorPro[2*i+1]->Sumw2();
    fCMEList->Add(fImQVectorPro[2*i+1]);
  }
 
 
  //------- 2particle correlator -------------
  
  fDnnPro[0] = new TProfile("fDnnPro_0_D11_SS","fDnnPro_0_D11_SS",10,0,100,"s");
  fDnnPro[0]->Sumw2();
  fDnnPro[1] = new TProfile("fDnnPro_1_D11_OS","fDnnPro_1_D11_OS",10,0,100,"s");
  fDnnPro[1]->Sumw2();
  
  fDnnHist[0] = new TH1F("fDnnHist_0_D11_SS","fDnnHist_0_D11_SS",10,0,100);
  fDnnHist[0]->Sumw2();
  fDnnHist[1] = new TH1F("fDnnHist_1_D11_OS","fDnnHist_1_D11_OS",10,0,100);
  fDnnHist[1]->Sumw2();
    
  for(Int_t i=0; i<nDnncor; i++) {
    fCMEList->Add(fDnnPro[i]);
    fCMEList->Add(fDnnHist[i]);
  }
   
  //--------- 3 particle correlator -----------
  fCMEPro[0]  = new TProfile("fCMEPro_0_C112_SS","fCMEPro_0_C112_SS",10,0,100,"s");
  fCMEPro[0]->Sumw2();
  fCMEPro[1]  = new TProfile("fCMEPro_1_C112_OS","fCMEPro_1_C112_OS",10,0,100,"s");
  fCMEPro[1]->Sumw2();
  fCMEPro[2]  = new TProfile("fCMEPro_2_C132_SS","fCMEPro_2_C132 SS",10,0,100,"s");
  fCMEPro[2]->Sumw2();
  fCMEPro[3]  = new TProfile("fCMEPro_3_C132_OS","fCMEPro_3_C132_OS",10,0,100,"s");
  fCMEPro[3]->Sumw2();
  fCMEPro[4]  = new TProfile("fCMEPro_4_C123_SS","fCMEPro_4_C123_SS",10,0,100,"s");
  fCMEPro[4]->Sumw2();
  fCMEPro[5]  = new TProfile("fCMEPro_5_C123_OS","fCMEPro_5_C123_OS",10,0,100,"s");
  fCMEPro[5]->Sumw2();
  fCMEPro[6]  = new TProfile("fCMEPro_6_C224_SS","fCMEPro_6_C224_SS",10,0,100,"s");
  fCMEPro[6]->Sumw2();
  fCMEPro[7]  = new TProfile("fCMEPro_7_C224_OS","fCMEPro_7_C224_OS",10,0,100,"s");
  fCMEPro[7]->Sumw2();
  
  fCMEHist[0]  = new TH1F("fCMEHist_0_C112_SS","fCMEHist_0_C112_SS",10,0,100);
  fCMEHist[0]->Sumw2();
  fCMEHist[1]  = new TH1F("fCMEHist_1_C112_OS","fCMEHist_1_C112_OS",10,0,100);
  fCMEHist[1]->Sumw2();
  fCMEHist[2]  = new TH1F("fCMEHist_2_C132_SS","fCMEHist_2_C132_SS",10,0,100);
  fCMEHist[2]->Sumw2();
  fCMEHist[3]  = new TH1F("fCMEHist_3_C132_OS","fCMEHist_3_C132_OS",10,0,100);
  fCMEHist[3]->Sumw2();
  fCMEHist[4]  = new TH1F("fCMEHist_4_C123_SS","fCMEHist_4_C123_SS",10,0,100);
  fCMEHist[4]->Sumw2();
  fCMEHist[5]  = new TH1F("fCMEHist_5_C123_OS","fCMEHist_5_C123_OS",10,0,100);
  fCMEHist[5]->Sumw2();
  fCMEHist[6]  = new TH1F("fCMEHist_6_C224_SS","fCMEHist_6_C224_SS",10,0,100);
  fCMEHist[6]->Sumw2();
  fCMEHist[7]  = new TH1F("fCMEHist_7_C224_OS","fCMEHist_7_C224_OS",10,0,100);
  fCMEHist[7]->Sumw2();
  
  for(Int_t i=0; i<nCMEcor; i++) {
    fCMEList->Add(fCMEPro[i]);
    fCMEList->Add(fCMEHist[i]);
  }
  
  
  
}

//=====================================================================================================

void CalculateFlowCME::UserExec() 
{
	Make(fEvent);
}

//=====================================================================================================

void CalculateFlowCME::Make(Event* anEvent) {
	Int_t nPrim = 0;
	Int_t Pid = -99999;
	Double_t dPhi = 0.; // azimuthal angle in the laboratory frame
	Double_t dPt  = 0.; // transverse momentum
	Double_t dEta = 0.; // pseudorapidity
	Int_t dCharge = -2; // charge
	fCenBin = GetCRCCenBin(fCentralityEBE);
	nPrim = anEvent->getNtrack();
	Int_t cw = 0; 
	Int_t nChargedParticles = 0;
	// multiplicity for charged particles 
	
	Double_t fNumOfPos = 0;
	Double_t fNumOfNeg = 0;
	
	if (nPrim < minNtracks) return;

	fnEvent->Fill(0.5);
	
	for(Int_t i=0;i<nPrim;i++) {
		Pid = anEvent->getParticle(i).getPid();
		dPhi = anEvent->getParticle(i).getPhi();
		dPt = anEvent->getParticle(i).getPt();
		dEta = anEvent->getParticle(i).getRapidity();  // in MC, the rapidity is stored rather than pesudorapidity (eta)
		dCharge = anEvent->getParticle(i).getCharge();

		if (dCharge == 0) continue;
		cw = (dCharge > 0. ? 0 : 1);
		
		if (dPt > maxPtCut) continue;
		if (dPt < minPtCut) continue;
		if (TMath::Abs(dEta) > maxEtaCut) continue;
		
		nChargedParticles++;
		// ====== Plot some QA ============
		// spectra for pions, kaons and protons
		if (doQA) {
			fPtChargedParticlesDistribution->Fill(dPt);
			fEtaChargedParticlesDistribution->Fill(dEta);
			fPhiChargedParticlesDistribution->Fill(dPhi);
			fEtaChargedParticlesDistributionPaperbinning->Fill(dEta);
			// Pions
			if (Pid == 211 || Pid == -211) {
				fPionsPtSpectra->Fill(dPt);
				fPionsEtaSpectra->Fill(dEta);
				fPionsPhiSpectra->Fill(dPhi);
				if (Pid == 211) {
					fPosPionsPtSpectra->Fill(dPt);
					fPosPionsEtaSpectra->Fill(dEta);
					fPosPionsPhiSpectra->Fill(dPhi);
				}
				if (Pid == -211) {
					fAntiPionsPtSpectra->Fill(dPt);
					fAntiPionsEtaSpectra->Fill(dEta);
					fAntiPionsPhiSpectra->Fill(dPhi);
				}
			}
			// Kaons
			if (Pid == 321 || Pid == -321) {
				fKaonsPtSpectra->Fill(dPt);
				fKaonsEtaSpectra->Fill(dEta);
				fKaonsPhiSpectra->Fill(dPhi);
				if (Pid == 321) {
					fPosKaonsPtSpectra->Fill(dPt);
					fPosKaonsEtaSpectra->Fill(dEta);
					fPosKaonsPhiSpectra->Fill(dPhi);
				}
				if (Pid == -321) {
					fAntiKaonsPtSpectra->Fill(dPt);
					fAntiKaonsEtaSpectra->Fill(dEta);
					fAntiKaonsPhiSpectra->Fill(dPhi);
				}
			}
			// Protons
			if (Pid == 2212 || Pid == -2212) {
				fProtonsPtSpectra->Fill(dPt);
				fProtonsEtaSpectra->Fill(dEta);
				fProtonsPhiSpectra->Fill(dPhi);
				if (Pid == 2212) {
					fPosProtonsPtSpectra->Fill(dPt);
					fPosProtonsEtaSpectra->Fill(dEta);
					fPosProtonsPhiSpectra->Fill(dPhi);
				}
				if (Pid == -2212) {
					fAntiProtonsPtSpectra->Fill(dPt);
					fAntiProtonsEtaSpectra->Fill(dEta);
					fAntiProtonsPhiSpectra->Fill(dPhi);
				}
			}
		}
		
		
		// ====== for calculateFlowGF (generic framework) =====
		// Generic Framework: Calculate Re[Q_{m*n,k}] and Im[Q_{m*n,k}] for this event (m = 1,2,...,12, k = 0,1,...,8):
        for(Int_t m=0;m<21;m++) {
			for(Int_t k=0;k<9;k++) {
				(*fReQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
				(*fImQGF)(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);	
			 }
        }
        
		for(Int_t ptb=0; ptb<fkGFPtB; ptb++) {
			if(ptb==0 && dPt>0.5) continue;
			if(ptb==1 && (dPt<0.5 || dPt>1.)) continue;
			if(ptb==2 && (dPt<1. || dPt>2.)) continue;
			if(ptb==3 && dPt<2.) continue;
			if(ptb==4 && (dPt<2. || dPt>2.5)) continue;
			if(ptb==5 && dPt<2.5) continue;
			if(ptb==6 && (dPt<2.5 || dPt>3.)) continue;
			if(ptb==7 && dPt<3.) continue;
			for(Int_t m=0;m<21;m++) {
				for(Int_t k=0;k<9;k++) {
					(*fReQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Cos(m*dPhi);
					(*fImQGFPt[ptb])(m,k) += pow(wPhiEta*wPhi*wPt*wEta*wTrack,k)*TMath::Sin(m*dPhi);
				}
			}
        }
		// ====== for calculateFlowQC ========= // extra dimension for charge !
		Int_t ch = 0;//to store all, and +- separately
		if(dCharge>0){ch=1.;}//this is not optmal 0 will never be fillew
		if(dCharge<0){ch=2.;}
		for (Int_t k=0; k<fQVecPower; k++) {
			for (Int_t h=0;h<fFlowNHarmMax;h++) {//h=0 for n=1
                                fPOIPtDiffQRe[0][k][h]->Fill(dPt,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));//for powers of the weight from 0th to xth *cos(n*dPhi) 
                                fPOIPtDiffQIm[0][k][h]->Fill(dPt,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
                                fPOIPtDiffMul[0][k][h]->Fill(dPt,pow(wPhiEta,k));

				fPOIPtDiffQRe[ch][k][h]->Fill(dPt,pow(wPhiEta,k)*TMath::Cos((h+1.)*dPhi));//for powers of the weight from 0th to xth *cos(n*dPhi) 
				fPOIPtDiffQIm[ch][k][h]->Fill(dPt,pow(wPhiEta,k)*TMath::Sin((h+1.)*dPhi));
				fPOIPtDiffMul[ch][k][h]->Fill(dPt,pow(wPhiEta,k));
			}
		}
		// ====== for calculateFlowfromBW ======= 
		if(dCharge > 0) {//here are concrete numbers
          ReQ1P += cos(dPhi);
          ImQ1P += sin(dPhi);
          ReQ2P += cos(2*dPhi);
          ImQ2P += sin(2*dPhi);
          ReQ3P += cos(3*dPhi);
          ImQ3P += sin(3*dPhi);
          ReQ4P += cos(4*dPhi);
          ImQ4P += sin(4*dPhi);
	      MQP++;
        } else if(dCharge < 0) {
          ReQ1N += cos(dPhi);
          ImQ1N += sin(dPhi);
          ReQ2N += cos(2*dPhi);
          ImQ2N += sin(2*dPhi);
          ReQ3N += cos(3*dPhi);
          ImQ3N += sin(3*dPhi);
          ReQ4N += cos(4*dPhi);
          ImQ4N += sin(4*dPhi);

          MQN++;
        }
	//cout <<"after fill ReQ1n"<<ReQ1N<<endl;
        for(int h=0; h<nHar; h++) {
          hRepn[h]->Fill(dPt,cos((h+1.)*dPhi)); // -> change to h+1 to get the v1 
      	  hImpn[h]->Fill(dPt,sin((h+1.)*dPhi));   
   
          ReQn[h] += cos((h+1.)*dPhi);//for integrated vn
          ImQn[h] += sin((h+1.)*dPhi);	
	}
		hMpn->Fill(dPt);
		MQ++;
		// ====== for calculateCMW ========		//fTPCQn2xWhole->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi)); // pt vs. cos(2phi)
		//fTPCQn2yWhole->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi)); // pt vs. sin(2phi)
		//fWgtWhole->Fill(dPt, trkWgt);
		//fWgt2Whole->Fill(dPt, pow(trkWgt,2));
				
		if (dEta < fEtaGapNeg) {
			//fTPCQn2xEtaNegWhole->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
			//fTPCQn2yEtaNegWhole->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
			//fWgtEtaNegWhole->Fill(dPt, trkWgt);
			//fWgt2EtaNegWhole->Fill(dPt, pow(trkWgt,2));
			
			fSumTPCQn2xEtaNegWhole += trkWgt*TMath::Cos(gPsiN*dPhi);
			fSumTPCQn2yEtaNegWhole += trkWgt*TMath::Sin(gPsiN*dPhi);
			fSumWgtEtaNegWhole  += trkWgt;
		} else if (dEta > fEtaGapPos) {
			//fTPCQn2xEtaPosWhole->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
			//fTPCQn2yEtaPosWhole->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
			//fWgtEtaPosWhole->Fill(dPt, trkWgt);
			//fWgt2EtaPosWhole->Fill(dPt, pow(trkWgt,2));
			
			fSumTPCQn2xEtaPosWhole += trkWgt*TMath::Cos(gPsiN*dPhi);
			fSumTPCQn2yEtaPosWhole += trkWgt*TMath::Sin(gPsiN*dPhi);
			fSumWgtEtaPosWhole  += trkWgt;
		}

		if (dEta < fEtaGapNeg) {
			if (dCharge>0) {
				//fTPCQn2xEtaNegChPos->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
				//fTPCQn2yEtaNegChPos->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
				//fWgtEtaNegChPos->Fill(dPt, trkWgt);
				//fWgt2EtaNegChPos->Fill(dPt, pow(trkWgt,2));
				
				fSumTPCQn2xEtaNegChPos += trkWgt*TMath::Cos(gPsiN*dPhi);
				fSumTPCQn2yEtaNegChPos += trkWgt*TMath::Sin(gPsiN*dPhi);
				fSumWgtEtaNegChPos += trkWgt;
			} else if (dCharge<0) {
				//fTPCQn2xEtaNegChNeg->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
				//fTPCQn2yEtaNegChNeg->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
				//fWgtEtaNegChNeg->Fill(dPt, trkWgt);
				//fWgt2EtaNegChNeg->Fill(dPt, pow(trkWgt,2));
				
				fSumTPCQn2xEtaNegChNeg += trkWgt*TMath::Cos(gPsiN*dPhi);
				fSumTPCQn2yEtaNegChNeg += trkWgt*TMath::Sin(gPsiN*dPhi);
				fSumWgtEtaNegChNeg  += trkWgt;
			}
		} else if (dEta > fEtaGapPos) {
			if (dCharge>0) {
				//fTPCQn2xEtaPosChPos->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
				//fTPCQn2yEtaPosChPos->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
				//fWgtEtaPosChPos->Fill(dPt, trkWgt);
				//fWgt2EtaPosChPos->Fill(dPt, pow(trkWgt,2));
				
				fSumTPCQn2xEtaPosChPos += trkWgt*TMath::Cos(gPsiN*dPhi);
				fSumTPCQn2yEtaPosChPos += trkWgt*TMath::Sin(gPsiN*dPhi);
				fSumWgtEtaPosChPos  += trkWgt;
			} else if (dCharge<0) {
				//fTPCQn2xEtaPosChNeg->Fill(dPt, trkWgt*TMath::Cos(gPsiN*dPhi));
				//fTPCQn2yEtaPosChNeg->Fill(dPt, trkWgt*TMath::Sin(gPsiN*dPhi));
				//fWgtEtaPosChNeg->Fill(dPt, trkWgt);
				//fWgt2EtaPosChNeg->Fill(dPt, pow(trkWgt,2));
				
				fSumTPCQn2xEtaPosChNeg += trkWgt*TMath::Cos(gPsiN*dPhi);
				fSumTPCQn2yEtaPosChNeg += trkWgt*TMath::Sin(gPsiN*dPhi);
				fSumWgtEtaPosChNeg  += trkWgt;
			}
		}
		

		if (dCharge > 0) {	  
			fNumOfPos += trkWgt;
		} else if (dCharge < 0) {
			fNumOfNeg += trkWgt;
		}
	}
	/*
	/////////////////////////////////////////////////////////////////////////
	if ((fNumOfPos + fNumOfNeg)!=0) fAchrgNet = (fNumOfPos - fNumOfNeg)/(fNumOfPos + fNumOfNeg); // A
	Float_t fWgtEvent = 1.0;
	fHistAChrgVsCent->Fill(fCentralityEBE, fAchrgNet, fWgtEvent); // <A>

	Double_t ResolutionEP = fSumTPCQn2xEtaNegChPos*fSumTPCQn2xEtaPosChPos + fSumTPCQn2yEtaNegChPos*fSumTPCQn2yEtaPosChPos;
	if ((fSumWgtEtaPosChPos*fSumWgtEtaNegChPos)!=0) ResolutionEP  = ResolutionEP/(fSumWgtEtaPosChPos*fSumWgtEtaNegChPos);

	Double_t ResWgt = TMath::Sqrt(fSumWgtEtaPosChPos*fSumWgtEtaNegChPos);

	/// Resolution vs cent:
	fHistEPResolution->Fill(fCentralityEBE,ResolutionEP,ResWgt*fWgtEvent);

	/// Resolution vs Ach for each cent:
	fHistEPResolutionAch[fCenBin]->Fill(fAchrgNet,ResolutionEP,ResWgt*fWgtEvent);
	*/
	/*
	////////-----> Starting 2nd track Loop -----------
		
	for (Int_t i=0;i<nPrim;i++) {
	///----------------------------------------------------------------	
		
		uqRe = TMath::Cos(gPsiN*dPhi); // one track
		uqIm = TMath::Sin(gPsiN*dPhi);
		
		//----> Remove Track from EP calculation ------
		sumQxTPC = fSumTPCQn2xEtaNegWhole+fSumTPCQn2xEtaPosWhole; // Cos(2*dPhi)_{eta-} + Cos(2*dPhi)_{eta+} 
		sumQyTPC = fSumTPCQn2yEtaNegWhole+fSumTPCQn2yEtaPosWhole; // Sin(2*dPhi)_{eta-} + Sin(2*dPhi)_{eta+} 
		sumWgtTPC = fSumWgtEtaNegWhole+fSumWgtEtaPosWhole; // # of tracks {eta-} + # of tracks {eta+}
		
		sumQxTPCEtaNegChPos = fSumTPCQn2xEtaNegChPos;   
		sumQyTPCEtaNegChPos = fSumTPCQn2yEtaNegChPos;
		sumQxTPCEtaPosChPos = fSumTPCQn2xEtaPosChPos;
		sumQyTPCEtaPosChPos = fSumTPCQn2yEtaPosChPos;
		
		sumWgtEtaNegChPos = fSumWgtEtaNegChPos;
		sumWgtEtaPosChPos = fSumWgtEtaPosChPos;

		sumQxTPC -= trkWgt*uqRe; // [∑cos(2phi)] - cos(2phi_i)
		sumQyTPC -= trkWgt*uqIm; // [∑sin(2phi)] - sin(2phi_i)
		sumWgtTPC -= trkWgt; // [∑tracks] - 1
		//// remove AutoCorrelation:
		if (dEta < fEtaGapNeg) {
		  sumQxTPCEtaNegChPos -= trkWgt*uqRe;
		  sumQyTPCEtaNegChPos -= trkWgt*uqIm; 
		  sumWgtEtaNegChPos   -= trkWgt;
		} else if (dEta > fEtaGapPos) {
		  sumQxTPCEtaPosChPos -= trkWgt*uqRe;
		  sumQyTPCEtaPosChPos -= trkWgt*uqIm; 
		  sumWgtEtaPosChPos   -= trkWgt;
		}
		//---------------------------------------------

		///---------- Prepare components for three-particle correlator -------------
		
		if (dCharge > 0) {
		
		  fHistv2AchChrgPosEtaNeg[0][fCenBin]->Fill(fAchrgNet, (uqRe*sumQxTPCEtaNegChPos + uqIm*sumQyTPCEtaNegChPos)/sumWgtEtaNegChPos, trkWgt); //This Trk weigth is for uQ	  
          fHistv2AchChrgPosEtaPos[0][fCenBin]->Fill(fAchrgNet, (uqRe*sumQxTPCEtaPosChPos + uqIm*sumQyTPCEtaPosChPos)/sumWgtEtaPosChPos, trkWgt);
          
          // calculation of <<v2^{+}A>> and <<v2^{+}>>: uqRe = cos(2phi_i), sumQxTPC = cos(2Phi') = [∑cos(2phi)] - cos(2phi_i), uqIm = sin(2phi_i), sumQyTPC = sin(2Phi') = [∑sin(2phi)] - sin(2phi_i)
          fHistAveV2PosAch->Fill(fCentralityEBE, (uqRe*sumQxTPC + uqIm*sumQyTPC)/sumWgtTPC*fAchrgNet, trkWgt); // <<v2^{+}A>> = <cos(2phi_i)*cos(2Phi')+sin(2phi_1)*sin(2Phi')> * A
          fHistAveV2Pos->Fill(fCentralityEBE, (uqRe*sumQxTPC + uqIm*sumQyTPC)/sumWgtTPC); // <<v2^{+}>> 
          
		  if (dEta > fEtaGapPos) {
			sumQ2xChrgPosEtaPos += trkWgt*uqRe;
			sumQ2yChrgPosEtaPos += trkWgt*uqIm;
			NumOfChrgPosEtaPos  += trkWgt;
		  } else if (dEta < fEtaGapNeg) {
			sumQ2xChrgPosEtaNeg += trkWgt*uqRe;
			sumQ2yChrgPosEtaNeg += trkWgt*uqIm;
			NumOfChrgPosEtaNeg  += trkWgt;
		  }
		  
		  ///+ve Ch done	
		} else if (dCharge < 0) {  //-Ve charge
		  
		  fHistv2AchChrgNegEtaNeg[0][fCenBin]->Fill(fAchrgNet, (uqRe*sumQxTPCEtaNegChPos + uqIm*sumQyTPCEtaNegChPos)/sumWgtEtaNegChPos, trkWgt);
		  fHistv2AchChrgNegEtaPos[0][fCenBin]->Fill(fAchrgNet, (uqRe*sumQxTPCEtaPosChPos + uqIm*sumQyTPCEtaPosChPos)/sumWgtEtaPosChPos, trkWgt);

          fHistAveV2NegAch->Fill(fCentralityEBE, (uqRe*sumQxTPC + uqIm*sumQyTPC)/sumWgtTPC*fAchrgNet, trkWgt); // <<v2^{-}A>>
          fHistAveV2Neg->Fill(fCentralityEBE, (uqRe*sumQxTPC + uqIm*sumQyTPC)/sumWgtTPC); // <<v2^{-}>>

		  if (dEta > fEtaGapPos) {
			sumQ2xChrgNegEtaPos += trkWgt*uqRe;
			sumQ2yChrgNegEtaPos += trkWgt*uqIm;
			NumOfChrgNegEtaPos  += trkWgt;
		  } else if (dEta < fEtaGapNeg) {
			sumQ2xChrgNegEtaNeg += trkWgt*uqRe;
			sumQ2yChrgNegEtaNeg += trkWgt*uqIm;
			NumOfChrgNegEtaNeg  += trkWgt;
		  }

		}/// if -ve Particle
		//----------- v2 vs Ach filled ---------
    
	}///------> 2nd Track loop Ends here.<--------
		/////////////////////////////////////////////////////////////////////////
		
	*/

	
	if (doQA) {
		fMultChargedParticlesDistribution->Fill(nChargedParticles); // x-axis no. charge particles, y-axis no. of events
	}
	CalculateFlowQC();
	CalculateFlowGF();
	CalculateFlowFromBW();
	//CalculateCMW();
	CalculateCME();
	ResetEventByEventQuantities();
	
}                               

//=====================================================================================================

void CalculateFlowCME::ResetEventByEventQuantities() 
{
	//FlowGF
    fReQGF->Zero();
    fImQGF->Zero();
    for(Int_t i=0; i<fkGFPtB; i++) {
		fReQGFPt[i]->Zero();
		fImQGFPt[i]->Zero();
	}

    //FlowQC
    for(Int_t ch=0; ch<charge;ch++){
    for(Int_t c=0;c<fQVecPower;c++) {
      for (Int_t h=0;h<fFlowNHarmMax;h++) {
        if(fPOIPtDiffQRe[ch][c][h]) fPOIPtDiffQRe[ch][c][h]->Reset();
        if(fPOIPtDiffQIm[ch][c][h]) fPOIPtDiffQIm[ch][c][h]->Reset();
        if(fPOIPtDiffMul[ch][c][h]) fPOIPtDiffMul[ch][c][h]->Reset();
      }
    }
    }
    //FlowFromBW
    // reset Q-vectors
    ReQ1P = 0; ReQ1N = 0;
    ImQ1P = 0; ImQ1N = 0;
    ReQ2P = 0; ImQ2P = 0;
    ReQ2N = 0; ImQ2N = 0;
    ReQ3P = 0; ImQ3P = 0;
    ReQ3N = 0; ImQ3N = 0;
    ReQ4P = 0; ImQ4P = 0;
    ReQ4N = 0; ImQ4N = 0;
    MQP = 0; MQN = 0;
    
    for(int h=0; h<nHar; h++) {
      ReQn[h] = 0;
      ImQn[h] = 0;
    }
    MQ = 0;
    PMQ = 0;
    NMQ = 0;
    
    //CMW
    // first loop
    fSumTPCQn2xEtaNegWhole = 0;
	fSumTPCQn2yEtaNegWhole = 0;
	fSumWgtEtaNegWhole = 0;
	fSumTPCQn2xEtaPosWhole = 0;
	fSumTPCQn2yEtaPosWhole = 0;
	fSumWgtEtaPosWhole = 0;
	fSumTPCQn2xEtaNegChPos = 0;
	fSumTPCQn2yEtaNegChPos = 0;
	fSumWgtEtaNegChPos = 0;
	fSumTPCQn2xEtaNegChNeg = 0;
	fSumTPCQn2yEtaNegChNeg = 0;
	fSumWgtEtaNegChNeg = 0;
	fSumTPCQn2xEtaPosChPos = 0;
	fSumTPCQn2yEtaPosChPos = 0;
	fSumWgtEtaPosChPos = 0;
	fSumTPCQn2xEtaPosChNeg = 0;
	fSumTPCQn2yEtaPosChNeg = 0;
	fSumWgtEtaPosChNeg = 0;
	// second loop
    uqRe=0; 
	uqIm=0;
	sumQxTPC=0;
	sumQyTPC=0;
	sumWgtTPC=0;
	sumQxTPCEtaNegChPos=0; 
	sumQyTPCEtaNegChPos=0;
	sumQxTPCEtaPosChPos=0;
	sumQyTPCEtaPosChPos=0;
	sumWgtEtaNegChPos=0;
	sumWgtEtaPosChPos=0;
	NumOfChrgPosEtaPos=0;
	NumOfChrgPosEtaNeg=0;
	sumQxTPCEtaNegChPos=0;
	sumQyTPCEtaNegChPos=0;
	sumWgtEtaNegChPos=0;
	sumQxTPCEtaPosChPos=0;
	sumQyTPCEtaPosChPos=0;
	sumWgtEtaPosChPos=0;
	sumQ2xChrgPosEtaPos=0;
	sumQ2yChrgPosEtaPos=0;
	NumOfChrgPosEtaPos=0;
	sumQ2xChrgPosEtaNeg=0;
	sumQ2yChrgPosEtaNeg=0;
	NumOfChrgPosEtaNeg=0;
	sumQ2xChrgNegEtaPos=0;
	sumQ2yChrgNegEtaPos=0;
	NumOfChrgNegEtaPos=0;
	sumQ2xChrgNegEtaNeg=0;
	sumQ2yChrgNegEtaNeg=0;
	NumOfChrgNegEtaNeg=0;
	
}
//=====================================================================================================

void CalculateFlowCME::Terminate()
{
	//FinalizeFlowQC(); //only for the integrated results relevant ... diifpt from calculate function
	FinalizeFlowGF();
	//FinalizeCMW();
	FinalizeCME();
}

//=====================================================================================================

void CalculateFlowCME::CalculateFlowFromBW()
{
	Double_t EvPlTPC = TMath::ATan2(fSumTPCQn2yEtaNegWhole+fSumTPCQn2yEtaPosWhole,fSumTPCQn2xEtaNegWhole+fSumTPCQn2xEtaPosWhole)/2.;
	Double_t pRe = 0, pIm = 0, mp = 0, posmp = 0, negmp = 0, ep = 0, posep = 0, negep = 0, pt = 0, eta=0;
	Double_t v1plus = 0, v1minus = 0, v1plusErr = 0, v1minusErr = 0, delta = 0, deltaErr = 0;
	for(int h=1; h<nHar; h++) {
//for pt diff
        for(Int_t k=1; k<=hRepn[h]->GetNbinsX(); k++) {//another loop for pos and neg with already existing hists
        pRe = hRepn[h]->GetBinContent(k);
        pIm = hImpn[h]->GetBinContent(k);
        mp = hMpn->GetBinContent(k);
        pt = hMpn->GetBinCenter(k);
        if(mp>0) {
          V2PtDiffPro[h][fCenBin]->Fill(pt,(pRe*cos((h+1.)*EvPlTPC)+pIm*sin((h+1.)*EvPlTPC))/mp);
        }
       
      if(MQ>0) {//WHAT IS MQ ?? -> here for VInt & CO

        V2IntPro[h]->Fill(fCentralityEBE,(ReQn[h]*cos((h+1.)*EvPlTPC)+ImQn[h]*sin((h+1.)*EvPlTPC))/MQ);
        V2IntProQC[h]->Fill(fCentralityEBE,(ReQn[h]*ReQn[h]+ImQn[h]*ImQn[h]-MQ)/(MQ*(MQ-1.)),MQ*(MQ-1.));
       	 }
        }
       }	
}
//=====================================================================================================

void CalculateFlowCME::CalculateCME()
{
	Double_t Psi = 0;
	Double_t CMEres = 0;
    /////////  delta_NN calculations ///////////////
    //Delta11:
    CMEres = (ReQ1P*ReQ1P+ImQ1P*ImQ1P-MQP)/(MQP*MQP-MQP);
    fDnnPro[0]->Fill(fCentralityEBE,CMEres,MQP*MQP-MQP);
    CMEres = (ReQ1P*ReQ1N+ImQ1P*ImQ1N)/(MQP*MQN);
    fDnnPro[1]->Fill(fCentralityEBE,CMEres,MQP*MQN);
    
    //////////// Q-vector components ////////////
    fReQVectorPro[0]->Fill(fCentralityEBE,ReQ1P);
    fReQVectorPro[1]->Fill(fCentralityEBE,ReQ1N);
    fReQVectorPro[2]->Fill(fCentralityEBE,ReQ2P);
    fReQVectorPro[3]->Fill(fCentralityEBE,ReQ2N);
    fReQVectorPro[4]->Fill(fCentralityEBE,ReQ3P);
    fReQVectorPro[5]->Fill(fCentralityEBE,ReQ3N);
    fReQVectorPro[6]->Fill(fCentralityEBE,ReQ4P);
    fReQVectorPro[7]->Fill(fCentralityEBE,ReQ4N);
    
    fImQVectorPro[0]->Fill(fCentralityEBE,ImQ1P);
    fImQVectorPro[1]->Fill(fCentralityEBE,ImQ1N);
    fImQVectorPro[2]->Fill(fCentralityEBE,ImQ2P);
    fImQVectorPro[3]->Fill(fCentralityEBE,ImQ2N);
    fImQVectorPro[4]->Fill(fCentralityEBE,ImQ3P);
    fImQVectorPro[5]->Fill(fCentralityEBE,ImQ3N);
    fImQVectorPro[6]->Fill(fCentralityEBE,ImQ4P);
    fImQVectorPro[7]->Fill(fCentralityEBE,ImQ4N);
    
    //////////// gamma_m,n calculations ////////////
    //C112
    CMEres = ((ReQ1P*ReQ1P-ImQ1P*ImQ1P-ReQ2P)*cos(2.*Psi)+(2.*ReQ1P*ImQ1P-ImQ2P)*sin(2.*Psi))/(MQP*MQP-MQP);
    fCMEPro[0]->Fill(fCentralityEBE,CMEres,MQP*MQP-MQP); //Same-Sign
    CMEres = ((ReQ1P*ReQ1N-ImQ1P*ImQ1N)*cos(2.*Psi)+(ReQ1P*ImQ1N+ImQ1P*ReQ1N)*sin(2.*Psi))/(MQP*MQN);
    fCMEPro[1]->Fill(fCentralityEBE,CMEres,MQP*MQN);     //Oppo-Sign
    
    //C132
    CMEres = ((ReQ1P*ReQ3P+ImQ1P*ImQ3P-ReQ2P)*cos(2.*Psi)-(ReQ3P*ImQ1P-ImQ3P*ReQ1P+ImQ2P)*sin(2.*Psi))/(MQP*(MQP-1.));
    fCMEPro[2]->Fill(fCentralityEBE,CMEres,MQP*(MQP-1.)); //Same-Sign
    CMEres = ((ReQ1P*ReQ3N+ImQ1P*ImQ3N)*cos(2.*Psi)-(ImQ1P*ReQ3N-ReQ1P*ImQ3N)*sin(2.*Psi))/(MQP*MQN);
    fCMEPro[3]->Fill(fCentralityEBE,CMEres,MQP*MQN);     //Oppo-Sign
    
    //C123 added for Sergei
    CMEres = ((ReQ1P*ReQ2P-ImQ1P*ImQ2P-ReQ3P)*cos(3.*Psi)+(ReQ1P*ImQ2P+ReQ2P*ImQ1P-ImQ3P)*sin(3.*Psi))/(MQP*(MQP-1.));
    fCMEPro[4]->Fill(fCentralityEBE,CMEres,MQP*(MQP-1.)); //Same-Sign
    CMEres = ((ReQ1P*ReQ2N-ImQ1P*ImQ2N)*cos(3.*Psi)+(ReQ1P*ImQ2N+ImQ1P*ReQ2N)*sin(3.*Psi))/(MQP*MQN);
    fCMEPro[5]->Fill(fCentralityEBE,CMEres,MQP*MQN);     //Oppo-Sign
    
    //C224
    CMEres = ((ReQ2P*ReQ2P-ImQ2P*ImQ2P-ReQ4P)*cos(4.*Psi)+(2.*ReQ2P*ImQ2P-ImQ4P)*sin(4.*Psi))/(MQP*(MQP-1.));
    fCMEPro[6]->Fill(fCentralityEBE,CMEres,MQP*(MQP-1.)); //Same-Sign
    CMEres = ((ReQ2P*ReQ2N-ImQ2P*ImQ2N)*cos(4.*Psi)+(ReQ2P*ImQ2N+ImQ2P*ReQ2N)*sin(4.*Psi))/(MQP*MQN);
    fCMEPro[7]->Fill(fCentralityEBE,CMEres,MQP*MQN);     //Oppo-Sign
    
    // use actual reaction plane
	
}

//=====================================================================================================

void CalculateFlowCME::CalculateCMW()
{			
  Double_t EvPlTPC = TMath::ATan2(fSumTPCQn2yEtaNegWhole+fSumTPCQn2yEtaPosWhole,fSumTPCQn2xEtaNegWhole+fSumTPCQn2xEtaPosWhole)/2.; // atan2(sin(2phi), cos(2phi))/2
  fEvPlTPC->Fill(EvPlTPC);
  
  Double_t QRePosCh = fSumTPCQn2xEtaNegChPos + fSumTPCQn2xEtaPosChPos;
  Double_t QImPosCh = fSumTPCQn2yEtaNegChPos + fSumTPCQn2yEtaPosChPos;
  Double_t QMuPosCh = fSumWgtEtaNegChPos + fSumWgtEtaPosChPos;
  
  Double_t QReNegCh = fSumTPCQn2xEtaNegChNeg + fSumTPCQn2xEtaPosChNeg;
  Double_t QImNegCh = fSumTPCQn2yEtaNegChNeg + fSumTPCQn2yEtaPosChNeg;
  Double_t QMuNegCh = fSumWgtEtaNegChNeg + fSumWgtEtaPosChNeg;
  
  Double_t v2plus = (QRePosCh * TMath::Cos(EvPlTPC) + QImPosCh * TMath::Sin(EvPlTPC))/QMuPosCh;
  Double_t v2minus = (QReNegCh * TMath::Cos(EvPlTPC) + QImNegCh * TMath::Sin(EvPlTPC))/QMuNegCh;
//either here or in CME  

  fv2plusminus->Fill(0.5, v2plus);
  fv2plusminus->Fill(1.5, v2minus);
  
  
  //////////////////////////////////////////////////////////////////////////////////
	  /// For cumulant method:

  /// Charge All(+-):


  Double_t c2WeightChrg     = fSumWgtEtaPosWhole*fSumWgtEtaNegWhole;
  if (c2WeightChrg!=0.0) {
    Double_t c2cumulantChrg   = (fSumTPCQn2xEtaPosWhole*fSumTPCQn2xEtaNegWhole + fSumTPCQn2yEtaPosWhole*fSumTPCQn2yEtaNegWhole)/c2WeightChrg;
    fHistv2cumAchChrgAll[fCenBin]->Fill(fAchrgNet, c2cumulantChrg, c2WeightChrg);   /// for denominator
  }

  Double_t c2WeightChrgPos   =  NumOfChrgPosEtaPos*fSumWgtEtaNegChPos;
  Double_t c2WeightChrgNeg   =  NumOfChrgNegEtaPos*fSumWgtEtaNegChNeg;
  
  Double_t c2WeightChrgPosChrgNeg   =  NumOfChrgPosEtaPos*fSumWgtEtaNegChNeg; //positive charge correlation with negative charge in opp subevent
  Double_t c2WeightChrgNegChrgPos   =  NumOfChrgNegEtaPos*fSumWgtEtaNegChPos;
  
  
  if((NumOfChrgPosEtaPos*fSumWgtEtaNegChPos)!=0.0){
    ///Chrg+:  
    //Double_t c2WeightChrgPos   =  NumOfChrgPosEtaPos*fSumWgtEtaNeg;
    Double_t c2cumulantChrgPos =  (sumQ2xChrgPosEtaPos*fSumTPCQn2xEtaNegChPos + sumQ2yChrgPosEtaPos*fSumTPCQn2yEtaNegChPos)/c2WeightChrgPos;
    fHistv2AchChrgPosEtaNeg[1][fCenBin]->Fill(fAchrgNet, c2cumulantChrgPos, c2WeightChrgPos);   /// for denominator
    
  }
  
  if((NumOfChrgNegEtaPos*fSumWgtEtaNegChNeg)!=0.0){
    ///Chrg-:  
    //Double_t c2WeightChrgNeg   =  NumOfChrgNegEtaPos*fSumWgtEtaNegChNeg;
    Double_t c2cumulantChrgNeg =  (sumQ2xChrgNegEtaPos*fSumTPCQn2xEtaNegChNeg + sumQ2yChrgNegEtaPos*fSumTPCQn2yEtaNegChNeg)/c2WeightChrgNeg;
    //fHistv2AchChrgNegEtaNeg[1][fCenBin]->Fill(fAchrgNet, c2cumulantChrgNeg, c2WeightChrgNeg);   /// for denominator
  }
  


  if((c2WeightChrgPosChrgNeg)!=0.0){
    //Chrg+: opp correlation
    Double_t c2cumulantChrgPosChrgNeg =  (sumQ2xChrgPosEtaPos*fSumTPCQn2xEtaNegChNeg + sumQ2yChrgPosEtaPos*fSumTPCQn2yEtaNegChNeg)/c2WeightChrgPosChrgNeg;
    fHistv2AchChrgPosChrgNeg[1][fCenBin]->Fill(fAchrgNet, c2cumulantChrgPosChrgNeg, c2WeightChrgPosChrgNeg);   /// for denominator
  }


  if((c2WeightChrgNegChrgPos)!=0.0){
    ///Chrg-:  
    Double_t c2cumulantChrgNegChrgPos =  (sumQ2xChrgNegEtaPos*fSumTPCQn2xEtaNegChPos + sumQ2yChrgNegEtaPos*fSumTPCQn2yEtaNegChPos)/c2WeightChrgNegChrgPos;
    fHistv2AchChrgNegChrgPos[1][fCenBin]->Fill(fAchrgNet, c2cumulantChrgNegChrgPos, c2WeightChrgNegChrgPos);   /// for denominator
  }
  
	
} // end of CalculateCMW()
//=====================================================================================================

void CalculateFlowCME::CalculateFlowGF()//here is the relevant filling
{
	Int_t order[fkFlowGFNOrde] = {0};
	// I ignore the mixed and wide Pt bins for now and trz to get the differential v1 pt
	for(Int_t k=0; k<fkFlowGFNOrde; k++) {
		order[k] = (k+1);
	}

	Double_t dMult = (*fReQGF)(0,0);//implement for pos and neg!

	for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {

		Double_t CorrOrd[fkFlowGFNOrde] = {0.};
		Double_t WeigOrd[fkFlowGFNOrde] = {0.};

		for(Int_t no=0; no<fkFlowGFNOrde; no++) {

			if(dMult<order[no]) continue;

			TArrayI harmonics(order[no]);
			Int_t halforder = (Int_t)order[no]/2;
			for(Int_t k=0; k<order[no]; k++) {//rethink this part too
				if(k<halforder) harmonics[k] = -(hr+2);//+2
				else            harmonics[k] = hr+2;//+2
			}
			TArrayI emptyness(order[no]);
			for(Int_t k=0; k<order[no]; k++) emptyness[k] = 0;

			std::complex<double> N = this->ucN(order[no], harmonics, -1);
			std::complex<double> D = this->ucN(order[no], emptyness, -1);

			if(D.real()>0.) {
				fFlowGFIntCorPro[hr][no]->Fill(fCentralityEBE,N.real()/D.real(),D.real()*fCenWeightEbE);
				CorrOrd[no] = N.real()/D.real();                     
				WeigOrd[no] = D.real();
			}

		} // end of for(Int_t no=0; no<fkFlowGFNOrde; no++)

		for(Int_t no=0; no<fkFlowGFNOrde; no++) {
			for(Int_t no2=0; no2<fkFlowGFNOrde; no2++) {
				if(WeigOrd[no]>0. && WeigOrd[no2]>0.) {
					fFlowGFIntCovPro[hr][no][no2]->Fill(fCentralityEBE,CorrOrd[no]*CorrOrd[no2],WeigOrd[no]*WeigOrd[no2]*fCenWeightEbE*fCenWeightEbE);
				}
			}
		}

	} // end of for(Int_t hr=0; hr<fkFlowGFNHarm; hr++)

	for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {
		for(Int_t hr2=0; hr2<fkFlowGFNHarm; hr2++) {

			if(dMult<4) continue;

			TArrayI harmonics(4);// rethink this part !!
			harmonics[0] = hr+2;//2
			harmonics[1] = hr2+2;//2
			harmonics[2] = -(hr+2);//2
			harmonics[3] = -(hr2+2);//2

			TArrayI emptyness(4);
			for(Int_t k=0; k<4; k++) emptyness[k] = 0;

			std::complex<double> N = this->ucN(4, harmonics, -1);
			std::complex<double> D = this->ucN(4, emptyness, -1);

			if(D.real()>0.) {
				fFlowGFMixedCorPro[hr][hr2]->Fill(fCentralityEBE,N.real()/D.real(),D.real()*fCenWeightEbE);
			}
		}
	}

	// in wide pt bins
	for(Int_t i=0; i<fkGFPtB; i++) {
		Double_t dMult = (*fReQGFPt[i])(0,0);

		for(Int_t hr=0; hr<fkFlowGFNHarm; hr++) {

			Double_t CorrOrd[fkFlowGFNOrde] = {0.};
			Double_t WeigOrd[fkFlowGFNOrde] = {0.};

			for(Int_t no=0; no<fkFlowGFNOrde; no++) {

				if(dMult<order[no]) continue;

				TArrayI harmonics(order[no]);
				Int_t halforder = (Int_t)order[no]/2;
				for(Int_t k=0; k<order[no]; k++) {//maybe this should be +2
					if(k<halforder) harmonics[k] = -(hr+2);
					else            harmonics[k] = hr+2;
				}
				TArrayI emptyness(order[no]);
				for(Int_t k=0; k<order[no]; k++) emptyness[k] = 0;
			
                        	 std::complex<double> N = this->ucN(order[no], harmonics, -1);
                        	 std::complex<double> D = this->ucN(order[no], emptyness, -1);
				if(D.real()>0.) {
					fFlowGFIntCorProPtB[i][hr][no]->Fill(fCentralityEBE,N.real()/D.real(),D.real()*fCenWeightEbE);
					CorrOrd[no] = N.real()/D.real();
					WeigOrd[no] = D.real();
				}

			} // end of for(Int_t no=0; no<fkFlowGFNOrde; no++)

			for(Int_t no=0; no<fkFlowGFNOrde; no++) {
				for(Int_t no2=0; no2<fkFlowGFNOrde; no2++) {
					if(WeigOrd[no]>0. && WeigOrd[no2]>0.) {
						fFlowGFIntCovProPtB[i][hr][no][no2]->Fill(fCentralityEBE,CorrOrd[no]*CorrOrd[no2],WeigOrd[no]*WeigOrd[no2]*fCenWeightEbE*fCenWeightEbE);
					}
				}
			}

		} // end of for(Int_t hr=0; hr<fkFlowGFNHarm; hr++)
	} // end of for(Int_t i=0; i<fkGFPtB; i++)

} // end of CalculateFlowCME::CalculateFlowGF()


void CalculateFlowCME::CalculateFlowQC() 
{//add the charge loop then think about the matrices and arrays - should the fPOIPTDiffQ rather not be a TProfile ?...hmmm..!
	Double_t FillPtBin = 0.;
	Double_t IQC2[fFlowNHarm] = {0.};
	Double_t IQC4[fFlowNHarm] = {0.};
	Double_t IQM2=0., IQM4=0.;
	Double_t WQM2=0., WQM4=0.;
	Double_t dQC2=0., dQC4=0.;
	Double_t dQM2=0., dQM4=0.;
	Double_t WdQM2=0., WdQM4=0., WdQM2EG=0., WdQM2EGB=0.;
	Double_t QRe=0., QIm=0., Q2Re2=0., Q2Im2=0., QRe3=0., QIm3=0., QM0=0., QM=0., QM2=0., QM3=0., QM4=0.;
	Double_t qpRe0=0., qpIm0=0., qpRe2=0., qpIm2=0., qp2Re=0., qp2Im=0., qpM0=0., qpM=0., qpM2=0., qpM3=0.;
    Double_t WqpM0=0., WqpAM=0.;
	Bool_t Q2f=kFALSE, Q4f=kFALSE, dQ2f=kFALSE, dQ4f=kFALSE;
	Bool_t WeigMul = kTRUE; //(fCorrWeightTPC==kMultiplicity ? kTRUE : kFALSE);
 for(Int_t ch = 0; ch<charge; ch++){
	for(Int_t hr=0; hr<fFlowNHarm; hr++) {
		// ********************************************************************
		// pT-integrated: {2}, {4} ********************************************
		// ********************************************************************

		// store reference flow (2 and 4p) ***********************************
		QRe=0.; QIm=0.; Q2Re2=0.; Q2Im2=0.; QRe3=0.; QIm3=0.;
		QM0=0.; QM=0.; QM2=0.; QM3=0.; QM4=0.;
		Q2f=kFALSE; Q4f=kFALSE;

		for(Int_t pt=0; pt<fPtDiffNBins; pt++) {//now 2nd [] for the order of the weight-> the weight is set to 0 for now; hr = 0 for n=1 harmonic
			QRe += fPOIPtDiffQRe[ch][1][hr]->GetBinContent(pt+1); // Cos((hr+1.)*dPhi)
			QIm += fPOIPtDiffQIm[ch][1][hr]->GetBinContent(pt+1); // Sin((hr+1.)*dPhi)
			//the next Q vectors are only for the cn4 calculation - 1. ask why the confusion with the [2*hr+3] vs note in the comment ... 
			Q2Re2 += fPOIPtDiffQRe[ch][2][2*hr+3]->GetBinContent(pt+1); // w^2*Cos((hr+2.)*dPhi) // WTF is this vector doing ?? double check with shi tomorrow !! 
			Q2Im2 += fPOIPtDiffQIm[ch][2][2*hr+3]->GetBinContent(pt+1); // w^2*Sin((hr+2.)*dPhi)
			QRe3 += fPOIPtDiffQRe[ch][3][hr+1]->GetBinContent(pt+1); // w^3*Cos((hr+2.)*dPhi)
			QIm3 += fPOIPtDiffQIm[ch][3][hr+1]->GetBinContent(pt+1); // w^3*Sin((hr+2.)*dPhi)

			QM0 += fPOIPtDiffMul[ch][0][0]->GetBinContent(pt+1); // w^0
			QM  += fPOIPtDiffMul[ch][1][0]->GetBinContent(pt+1); // w^1
			QM2 += fPOIPtDiffMul[ch][2][0]->GetBinContent(pt+1); // w^2
			QM3 += fPOIPtDiffMul[ch][3][0]->GetBinContent(pt+1); // w^3
			QM4 += fPOIPtDiffMul[ch][4][0]->GetBinContent(pt+1); // w^4
		}

		IQM2 = QM*QM-QM2;
		WQM2 = (WeigMul? IQM2 : 1.);
		if(QM0>1) {
			IQC2[hr] = (QRe*QRe+QIm*QIm-QM2)/IQM2; // <2> = |Q2,1|^2-S1,2/(S2,1-S1,2)
			fFlowQCIntCorPro[ch][hr][0]->Fill(fCentralityEBE,IQC2[hr],WQM2*fCenWeightEbE); // cent vs. <2>
			fFlowQCRefCorPro[ch][hr][0]->Fill(fCentralityEBE,IQC2[hr],WQM2*fCenWeightEbE); // cent vs. <2>
			Q2f = kTRUE;
		}

		IQM4 = QM*QM*QM*QM - 6.*QM2*QM*QM + 8.*QM3*QM + 3.*QM2*QM2 - 6.*QM4;
		WQM4 = (WeigMul? IQM4 : 1.);
		if(QM0>3) {
			IQC4[hr] = ((QRe*QRe+QIm*QIm)*(QRe*QRe+QIm*QIm)
						- 2.*(QRe*QRe*Q2Re2+2.*QRe*QIm*Q2Im2-QIm*QIm*Q2Re2)
						+ 8.*(QRe3*QRe+QIm3*QIm)
						+ (Q2Re2*Q2Re2+Q2Im2*Q2Im2)
						- 4.*QM2*(QRe*QRe+QIm*QIm)
						- 6.*QM4 + 2.*QM2*QM2) / IQM4;

			fFlowQCIntCorPro[ch][hr][1]->Fill(fCentralityEBE,IQC4[hr],WQM4*fCenWeightEbE);
			fFlowQCRefCorPro[ch][hr][1]->Fill(fCentralityEBE,IQC4[hr],WQM4*fCenWeightEbE);
			Q4f = kTRUE;
		}

		// product of correlations or covariances
		if(Q2f && Q4f) {
		  fFlowQCIntCorPro[ch][hr][2]->Fill(fCentralityEBE,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
		  fFlowQCRefCorPro[ch][hr][13]->Fill(fCentralityEBE,IQC2[hr]*IQC4[hr],WQM2*WQM4*fCenWeightEbE);
		}


		// ********************************************************************
		// pT-differential: {2}, {4} ******************************************
		// ********************************************************************

		// store pt-differential flow ****************************************
		for(Int_t pt=0; pt<fPtDiffNBins; pt++) {

		  FillPtBin = fPOIPtDiffQRe[ch][0][0]->GetBinCenter(pt+1);
		  qpRe0=0.; qpIm0=0.; qpRe2=0.; qpIm2=0.; qp2Re=0.; qp2Im=0.; qpM0=0.; qpM=0.; qpM2=0.; qpM3=0.;

		  qpRe0 = fPOIPtDiffQRe[ch][0][hr]->GetBinContent(pt+1);//for consistency also changed from [hr+1] to hr to get the 0th weight*cos(n=1*dPhi)
		  qpIm0 = fPOIPtDiffQIm[ch][0][hr]->GetBinContent(pt+1);
		//this is only relevant if I include c4 and needs some more understanding of the theory ... upsi wupsi
		  qpRe2 = fPOIPtDiffQRe[ch][2][hr+1]->GetBinContent(pt+1);//and this should be then [hr+1] -> hr
		  qpIm2 = fPOIPtDiffQIm[ch][2][hr+1]->GetBinContent(pt+1);
		  qp2Re = fPOIPtDiffQRe[ch][1][2*hr+3]->GetBinContent(pt+1);// would be next order previously: hr=0 -> n=2 and and second n=3
		  qp2Im = fPOIPtDiffQIm[ch][1][2*hr+3]->GetBinContent(pt+1);// now hr=0 -> n=1 and second n=2 ->should be hr+1 

		  qpM0 = fPOIPtDiffMul[ch][0][0]->GetBinContent(pt+1);
		  qpM  = fPOIPtDiffMul[ch][1][0]->GetBinContent(pt+1);
		  qpM2 = fPOIPtDiffMul[ch][2][0]->GetBinContent(pt+1);
		  qpM3 = fPOIPtDiffMul[ch][3][0]->GetBinContent(pt+1);

		  //if(hr==0) {//////////////////////////////////////////////////////////////////////////////////
		  // fFlowQCSpectra->Fill(fCentralityEBE,FillPtBin,qpM*fCenWeightEbE);
		  //}

		  dQM2 = qpM0*QM-qpM;
		  WdQM2 = (WeigMul? dQM2 : 1.);

		  if(qpM0>0 && QM0>0) {
			dQC2 = (qpRe0*QRe+qpIm0*QIm-qpM)/dQM2;
			fFlowQCCorPro[ch][fCenBin][hr][1]->Fill(FillPtBin,dQC2,WdQM2*fCenWeightEbE);
			dQ2f = kTRUE;
		  }

		  dQM4 = qpM0*(QM*QM*QM-3.*QM*QM2+2.*QM3)-3.*(qpM*(QM*QM-QM2)+2.*(qpM3-qpM2*QM));
		  WdQM4 = (WeigMul? dQM4 : 1.);
		  //@Shi dQM4 has to be nonzero
		  if(qpM0>0 && QM0>3 && dQM4!=0) { 
			dQC4 = ((pow(QRe,2.)+pow(QIm,2.))*(qpRe0*QRe+qpIm0*QIm)
							 - qp2Re*(pow(QRe,2.)-pow(QIm,2.))
							 - 2.*qp2Im*QRe*QIm
							 - qpRe0*(QRe*Q2Re2+QIm*Q2Im2)
							 + qpIm0*(QIm*Q2Re2-QRe*Q2Im2)
							 - 2.*QM2*(qpRe0*QRe+qpIm0*QIm)
							 - 2.*(pow(QRe,2.)+pow(QIm,2.))*qpM
							 + 6.*(qpRe2*QRe+qpIm2*QIm)
							 + 1.*(qp2Re*Q2Re2+qp2Im*Q2Im2)
							 + 2.*(qpRe0*QRe3+qpIm0*QIm3)
							 + 2.*qpM*QM2
							 - 6.*qpM3) / dQM4;
			fFlowQCCorPro[ch][fCenBin][hr][2]->Fill(FillPtBin,dQC4,WdQM4*fCenWeightEbE);
			dQ4f = kTRUE;
		  }

		  // product of correlations or covariances
		  if(Q2f && dQ2f) fFlowQCCorCovPro[ch][fCenBin][hr][0]->Fill(FillPtBin,IQC2[hr]*dQC2,WQM2*WdQM2*fCenWeightEbE);
		  if(Q4f && dQ2f) fFlowQCCorCovPro[ch][fCenBin][hr][1]->Fill(FillPtBin,IQC4[hr]*dQC2,WQM4*WdQM2*fCenWeightEbE);
		  if(Q2f && dQ4f) fFlowQCCorCovPro[ch][fCenBin][hr][2]->Fill(FillPtBin,IQC2[hr]*dQC4,WQM2*WdQM4*fCenWeightEbE);
		  if(dQ2f && dQ4f) fFlowQCCorCovPro[ch][fCenBin][hr][3]->Fill(FillPtBin,dQC2*dQC4,WdQM2*WdQM4*fCenWeightEbE);
		  if(Q4f && dQ4f) fFlowQCCorCovPro[ch][fCenBin][hr][4]->Fill(FillPtBin,IQC4[hr]*dQC4,WQM4*WdQM4*fCenWeightEbE);

		} // end of for(Int_t pt=0; pt<fCRCnPtBin; pt++)

	} // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
    // pt diff get fFlowQCCorHist ==============================


    // Pt-DIFFERENTIAL
  for(Int_t hr=0; hr<fFlowNHarm; hr++) {
    for(Int_t j=0; j<fFlowQCNRef; j++) {
      for(Int_t pt=1;pt<=fFlowQCRefCorPro[ch][hr][j]->GetNbinsX();pt++) {
        Double_t stats[6]={0.};
        fFlowQCRefCorPro[ch][hr][j]->GetXaxis()->SetRange(pt,pt);
        fFlowQCRefCorPro[ch][hr][j]->GetStats(stats);
        Double_t SumWeig   = stats[0];
        Double_t SumWeigSq  = stats[1];
        Double_t SumTwo  = stats[4];
        Double_t SumTwoSq = stats[5];

        if(SumWeig>0.) {
          Double_t Corr = SumTwo/SumWeig;
          Double_t SqCorr = SumTwoSq/SumWeig;
          Double_t Weig = SumWeig;
          Double_t SqWeig = SumWeigSq;
          Double_t spread=0., termA=0., termB=0.;
          if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
          if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
          if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
          Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
          if(CorrErr) {
            fFlowQCRefCorHist[ch][hr][j]->SetBinContent(pt,Corr);
            fFlowQCRefCorHist[ch][hr][j]->SetBinError(pt,CorrErr);
          }
        }
      } // end of for(Int_t pt=1;pt<=100;pt++)
      fFlowQCRefCorPro[ch][hr][j]->GetXaxis()->SetRange(1,fFlowQCRefCorPro[ch][hr][j]->GetNbinsX());
    } // end of for(Int_t j=0; j<5; j++)
    
    
    for (Int_t h=0; h<fCRCnCen; h++) {

      // STORE IN HISTOGRAMS

      for(Int_t j=0; j<fFlowQCNPro; j++) {
        for(Int_t pt=1;pt<=fPtDiffNBins;pt++) {

          Double_t stats[6]={0.};
          fFlowQCCorPro[ch][h][hr][j]->GetXaxis()->SetRange(pt,pt);
          fFlowQCCorPro[ch][h][hr][j]->GetStats(stats);
          Double_t SumWeig   = stats[0];
          Double_t SumWeigSq  = stats[1];
          Double_t SumTwo  = stats[4];
          Double_t SumTwoSq = stats[5];

          if(SumWeig>0.) {
            Double_t Corr = SumTwo/SumWeig;
            Double_t SqCorr = SumTwoSq/SumWeig;
            Double_t Weig = SumWeig;
            Double_t SqWeig = SumWeigSq;
            Double_t spread=0., termA=0., termB=0.;
            if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
            if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
            if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
            Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
            if(CorrErr) {
              fFlowQCCorHist[ch][h][hr][j]->SetBinContent(pt,Corr);
              fFlowQCCorHist[ch][h][hr][j]->SetBinError(pt,CorrErr);
            }
          }

        } // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)
        fFlowQCCorPro[ch][h][hr][j]->GetXaxis()->SetRange(1,fPtDiffNBins);
      }
      
      // reference flow
      // 2- and 4-particle cumulants
      Double_t QC2    = fFlowQCRefCorHist[ch][hr][0]->GetBinContent(h+1);
      Double_t QC2E = fFlowQCRefCorHist[ch][hr][0]->GetBinError(h+1);
      Double_t QC4    = fFlowQCRefCorHist[ch][hr][1]->GetBinContent(h+1);
      Double_t QC4E = fFlowQCRefCorHist[ch][hr][1]->GetBinError(h+1);
      Double_t Cn2 = QC2;
      Double_t Cn2E = QC2E;
      Double_t wCov24 = fFlowQCRefCorHist[ch][hr][13]->GetBinContent(h+1);
      Double_t Cn4 = QC4-2.*QC2*QC2;
      Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;

      Double_t Cos1 = fFlowQCRefCorHist[ch][hr][3]->GetBinContent(h+1); // <<cos(n*phi1)>>
      Double_t Sin1 = fFlowQCRefCorHist[ch][hr][4]->GetBinContent(h+1); // <<sin(n*phi1)>>
      Double_t Sin1P2 = fFlowQCRefCorHist[ch][hr][5]->GetBinContent(h+1);
      Double_t Cos1P2 = fFlowQCRefCorHist[ch][hr][6]->GetBinContent(h+1);
      Double_t Sin1M2M3 = fFlowQCRefCorHist[ch][hr][7]->GetBinContent(h+1);
      Double_t Cos1M2M3 = fFlowQCRefCorHist[ch][hr][8]->GetBinContent(h+1);
      // change vocabulary, to be changed
      Double_t cosP1nPhi = fFlowQCRefCorHist[ch][hr][3]->GetBinContent(h+1); // <<cos(n*phi1)>>
      Double_t sinP1nPhi = fFlowQCRefCorHist[ch][hr][4]->GetBinContent(h+1); // <<sin(n*phi1)>>
      Double_t sinP1nPhi1P1nPhi2 = fFlowQCRefCorHist[ch][hr][5]->GetBinContent(h+1); //sin(n*(phi1+phi2))
      Double_t cosP1nPhi1P1nPhi2 = fFlowQCRefCorHist[ch][hr][6]->GetBinContent(h+1);  //cos(n*(phi1+phi2))
      Double_t sinP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[ch][hr][7]->GetBinContent(h+1);  //sin(n*(phi1-phi2-phi3))
      Double_t cosP1nPhi1M1nPhi2M1nPhi3 = fFlowQCRefCorHist[ch][hr][8]->GetBinContent(h+1); //cos(n*(phi1-phi2-phi3))
      
      fFlowQCRefCorFinal[ch][hr][0]->SetBinContent(h+1,Cn2);
      fFlowQCRefCorFinal[ch][hr][0]->SetBinError(h+1,Cn2E);

      if(Cn4Esq>0.) {
        Double_t Cn4E = pow(Cn4Esq,0.5);
        fFlowQCRefCorFinal[ch][hr][3]->SetBinContent(h+1,Cn4);
        fFlowQCRefCorFinal[ch][hr][3]->SetBinError(h+1,Cn4E);
        if(Cn4<0.) {
          Double_t Flow4 = pow(fabs(Cn4),0.25);
          Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
          fFlowQCRefCorFinal[ch][hr][1]->SetBinContent(h+1,Flow4);
          fFlowQCRefCorFinal[ch][hr][1]->SetBinError(h+1,Flow4E);
        }
      }
      
      // pt-differential
      for(Int_t pt=1; pt<=fPtDiffNBins; pt++) {

        Double_t qp2    = fFlowQCCorHist[ch][h][hr][1]->GetBinContent(pt);
        Double_t qp2E = fFlowQCCorHist[ch][h][hr][1]->GetBinError(pt);
        Double_t qp4    = fFlowQCCorHist[ch][h][hr][2]->GetBinContent(pt);
        Double_t qp4E = fFlowQCCorHist[ch][h][hr][2]->GetBinError(pt);
        Double_t Dn2 = qp2;
        Double_t Dn2E = qp2E;
        Double_t Dn4 = qp4-2.*qp2*QC2;
        Double_t wCovTwoFourReduced = fFlowQCCorCovHist[ch][h][hr][1]->GetBinContent(pt);
        Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[ch][h][hr][4]->GetBinContent(pt);
        Double_t Dn4Esq = 4.*pow(QC2,2.)*pow(qp2E,2) + 4.*pow(qp2,2.)*pow(QC2E,2) + pow(qp4E,2.) - 4.*qp2*wCovTwoFourReduced - 4.*QC2*wCovTwoReducedFourReduced;


       fFlowQCFinalPtDifHist[ch][h][hr][5]->SetBinContent(pt,Dn2);
       fFlowQCFinalPtDifHist[ch][h][hr][5]->SetBinError(pt,Dn2E);

        if(abs(Cn2)>0) {
          Double_t Flow2 = Dn2/sqrt(fabs(Cn2));
          Double_t Flow2E = 0.;
          // change vocabulary, to be changed
          Double_t two = QC2;
          Double_t twoError = QC2E;
          Double_t twoReduced = qp2;
          Double_t twoReducedError = qp2E;
          Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[ch][h][hr][0]->GetBinContent(pt);
          Double_t v2PrimeErrorSquared = abs((1./4.)*pow(two,-3.)*(pow(twoReduced,2.)*pow(twoError,2.)
                                + 4.*pow(two,2.)*pow(twoReducedError,2.)
                                - 4.*two*twoReduced*wCovTwoTwoReduced));//absolute value 
          if(v2PrimeErrorSquared>0.){Flow2E = pow(v2PrimeErrorSquared,0.5);}

          if(Flow2E>0.) {
	    cout<<"pt "<<pt<<" FLow2 "<<Flow2<<endl;
            fFlowQCFinalPtDifHist[ch][h][hr][0]->SetBinContent(pt,Flow2);
            fFlowQCFinalPtDifHist[ch][h][hr][0]->SetBinError(pt,Flow2E);
          }
        }

        if(Dn4Esq>0.) {
          Double_t Dn4E = pow(Dn4Esq,0.5);
         fFlowQCFinalPtDifHist[ch][h][hr][6]->SetBinContent(pt,Dn4);
         fFlowQCFinalPtDifHist[ch][h][hr][6]->SetBinError(pt,Dn4E);
        }

        if(Cn4Esq>0.) {
          Double_t Flow4 = - Dn4/pow(fabs(Cn4),0.75);
          Double_t Flow4E = 0.;
          // change vocabulary, to be changed
          Double_t two = QC2;
          Double_t twoError = QC2E;
          Double_t twoReduced = qp2;
          Double_t twoReducedError = qp2E;
          Double_t four = QC4;
          Double_t fourError = QC4E;
          Double_t fourReduced = qp4;
          Double_t fourReducedError = qp4E;
          Double_t wCovTwoTwoReduced = fFlowQCCorCovHist[ch][h][hr][0]->GetBinContent(pt);
          Double_t wCovTwoFourReduced = fFlowQCCorCovHist[ch][h][hr][1]->GetBinContent(pt);
          Double_t wCovFourTwoReduced = fFlowQCCorCovHist[ch][h][hr][2]->GetBinContent(pt);
          Double_t wCovFourFourReduced = fFlowQCCorCovHist[ch][h][hr][3]->GetBinContent(pt);
          Double_t wCovTwoReducedFourReduced = fFlowQCCorCovHist[ch][h][hr][4]->GetBinContent(pt);
          Double_t wCovTwoFour = fFlowQCRefCorHist[ch][hr][13]->GetBinContent(h+1);

          Double_t v4PrimeErrorSquared = 0.;
          if(2.*pow(two,2.)-four>0.) {
            v4PrimeErrorSquared = pow(2.*pow(two,2.)-four,-7./2.)
            * (pow(2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced,2.)*pow(twoError,2.)
            + (9./16.)*pow(2.*two*twoReduced-fourReduced,2.)*pow(fourError,2.)
            + 4.*pow(two,2.)*pow(2.*pow(two,2.)-four,2.)*pow(twoReducedError,2.)
            + pow(2.*pow(two,2.)-four,2.)*pow(fourReducedError,2.)
            - (3./2.)*(2.*two*twoReduced-fourReduced)
            * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFour
            - 4.*two*(2.*pow(two,2.)-four)
            * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoTwoReduced
            + 2.*(2.*pow(two,2.)-four)
            * (2.*pow(two,2.)*twoReduced-3.*two*fourReduced+2.*four*twoReduced)*wCovTwoFourReduced
            + 3.*two*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourTwoReduced
            - (3./2.)*(2.*pow(two,2.)-four)*(2.*two*twoReduced-fourReduced)*wCovFourFourReduced
            - 4.*two*pow(2.*pow(two,2.)-four,2.)*wCovTwoReducedFourReduced);
          }
          if(v4PrimeErrorSquared>0.){Flow4E = pow(v4PrimeErrorSquared,0.5);}

          if(Flow4E>0.) {
            fFlowQCFinalPtDifHist[ch][h][hr][1]->SetBinContent(pt,Flow4);
            fFlowQCFinalPtDifHist[ch][h][hr][1]->SetBinError(pt,Flow4E);
          }
        }
      } // end of for(Int_t pt=1; pt<=fPtDiffNBins; pt++) {
    } // end of for (Int_t h=0; h<fCRCnCen; h++) {
  } // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
}//end of charge loop

for(Int_t pt=1; pt<=fPtDiffNBins; pt++){
	for(Int_t h = 0; h<fCRCnCen; h++){
		for(Int_t hr=0; hr<fFlowNHarm; hr++){
		//variables for calculation
		Double_t vpos = 0.;
		Double_t vneg = 0.;
		Double_t deltav = 0.;
		//statistical errors
		Double_t vposError = 0.;
		Double_t vnegError = 0.;
		Double_t deltavError = 0.;
 		//getting values for the calculation - focus on v1
		vpos = fFlowQCFinalPtDifHist[1][h][hr][0]->GetBinContent(pt);//[1] for pos charge [h] for centralities [h] for harmonic [hr=0] for differential flow v1 from QC 2 -> this means v2 is in here as well ! For hr=1 ...	
		vneg = fFlowQCFinalPtDifHist[2][h][hr][0]->GetBinContent(pt);//for neg charge the last is actually from cov loop ! No clue what that means for me, but 0 stores the flow from cn2 and not cn4 
		//getting errors for propagation as fully correlated variables
		vposError = fFlowQCFinalPtDifHist[1][h][hr][0]->GetBinError(pt);
		vnegError = fFlowQCFinalPtDifHist[2][h][hr][0]->GetBinError(pt);
		//calculate delta v1
		deltav = vpos - vneg;
		cout<<"pt "<<pt<<" vpos "<<vpos<<" vneg "<<vneg<<" deltav "<<deltav<<endl;
		//propagete error on v1
		deltavError = sqrt(abs(pow(vposError,2)-pow(vnegError,2)));
		//setting bin content
		fFlowQCFinalPtDifDeltaHist[h][hr][0]->SetBinContent(pt,deltav);//centrality, harmonic, cov
		fFlowQCFinalPtDifDeltaHist[h][hr][0]->SetBinError(pt, deltavError);
		}//end of hr<fkFlowNHarm
	}//end of h<fCRCnCen	
}//end of pt<=fPtDiffNBins
	
//use diffrence of final pt hist + error propagation to get the \Delta v1 hist filled
}

//=======================================================================================================================

void CalculateFlowCME::FinalizeCME()
{
  // finalise CME corr
  for(Int_t k=0;k<nCMEcor;k++) {
    for(Int_t c=1;c<=fCMEPro[k]->GetNbinsX();c++) { //centrality loop
      
      Double_t stats[6]={0.};
      
      fCMEPro[k]->GetXaxis()->SetRange(c,c);
      fCMEPro[k]->GetStats(stats);
      
      Double_t SumWeig    = stats[0];
      Double_t SumWeigSq  = stats[1];
      Double_t SumTwo     = stats[4];
      Double_t SumTwoSq   = stats[5];
      
      if(!SumWeig) continue;
      
      Double_t Corr   = SumTwo/SumWeig;
      Double_t SqCorr = SumTwoSq/SumWeig;
      Double_t Weig   = SumWeig;
      Double_t SqWeig = SumWeigSq;
      Double_t spread=0., termA=0., termB=0.;
      if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
      if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
      if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
      Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)

      fCMEHist[k]->SetBinContent(c,Corr);
      fCMEHist[k]->SetBinError(c,CorrErr);
    } // end of for(Int_t c=1;c<=fCMEPro[k]->GetNbinsX();c++)
  } // end of for(Int_t k=0;k<4;k++)



  ///// Get (OS - SS) for \gamma correlator:
  
  DeltaG112 = (TH1F*)(fCMEHist[1]->Clone("DeltaG112"));
  DeltaG112->Add(fCMEHist[0],-1);
  DeltaG132 = (TH1F*)(fCMEHist[3]->Clone("DeltaG132"));
  DeltaG132->Add(fCMEHist[2],-1);
  DeltaG123 = (TH1F*)(fCMEHist[5]->Clone("DeltaG123"));
  DeltaG123->Add(fCMEHist[4],-1);
  DeltaG224 = (TH1F*)(fCMEHist[6]->Clone("DeltaG224"));
  DeltaG224->Add(fCMEHist[6],-1);


  // finalise Dnn corr
  for(Int_t k=0;k<nDnncor;k++) {
    for(Int_t c=1;c<=fDnnPro[k]->GetNbinsX();c++) { //centrality loop

      Double_t stats[6]={0.};

      fDnnPro[k]->GetXaxis()->SetRange(c,c);
      fDnnPro[k]->GetStats(stats);

      Double_t SumWeig    = stats[0];
      Double_t SumWeigSq  = stats[1];
      Double_t SumTwo     = stats[4];
      Double_t SumTwoSq   = stats[5];

      if(!SumWeig) continue;

      Double_t Corr   = SumTwo/SumWeig;
      Double_t SqCorr = SumTwoSq/SumWeig;
      Double_t Weig   = SumWeig;
      Double_t SqWeig = SumWeigSq;
      Double_t spread=0., termA=0., termB=0.;
      if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
      if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
      if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
      Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)

      fDnnHist[k]->SetBinContent(c,Corr);
      fDnnHist[k]->SetBinError(c,CorrErr);
    } // end of for(Int_t c=1;c<=fCMEPro[k]->GetNbinsX();c++)
  } // end of for(Int_t k=0;k<4;k++)



  
  //------ OS-SS for \Delta_NN ----------
  DeltaD11 = (TH1F*)(fDnnHist[1]->Clone("DeltaD11"));
  DeltaD11->Add(fDnnHist[0],-1);
  
  //-------
  fCMEList->Add(DeltaD11);
  fCMEList->Add(DeltaG112);
  fCMEList->Add(DeltaG132);
  fCMEList->Add(DeltaG123);
  fCMEList->Add(DeltaG224);
}
//=======================================================================================================================

void CalculateFlowCME::FinalizeFlowGF()
{
  cout << "*************************************" << "\n";
  cout << "\n";
  cout << "Finalizing Flow with Generic Framework";
  cout << "\n";
  cout << "\n";

  for (Int_t h=0; h<fkFlowGFNHarm; h++){
    for(Int_t i=0; i<fkFlowGFNOrde; i++) {
      for(Int_t pt=1; pt<=fFlowGFIntCorPro[h][i]->GetNbinsX(); pt++) {
	 Double_t stats[6]={0.};
        fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(pt,pt);
        fFlowGFIntCorPro[h][i]->GetStats(stats);
        Double_t SumWeig   = stats[0];
        Double_t SumWeigSq  = stats[1];
        Double_t SumTwo  = stats[4];
        Double_t SumTwoSq = stats[5];
        if(SumWeig>0.) {
          Double_t Corr = SumTwo/SumWeig;
          Double_t SqCorr = SumTwoSq/SumWeig;
          Double_t Weig = SumWeig;
          Double_t SqWeig = SumWeigSq;
          Double_t spread=0., termA=0., termB=0.;
          if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
          if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
          if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
          Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
          if(CorrErr) {
            fFlowGFIntCorHist[h][i]->SetBinContent(pt,Corr);
            fFlowGFIntCorHist[h][i]->SetBinError(pt,CorrErr);
          }
        }
      } // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)
      fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][i]->GetNbinsX());
    } // end of for(Int_t i=0; i<fkFlowGFNOrde; i++)
  }

  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNOrde; i++) {
      for(Int_t k=0; k<fkFlowGFNOrde; k++) {
        for(Int_t pt=1; pt<=fFlowGFIntCovPro[h][i][k]->GetNbinsX(); pt++) {	// correlations:
          Double_t A = fFlowGFIntCorHist[h][i]->GetBinContent(pt); // <<A>>
          Double_t B = fFlowGFIntCorHist[h][k]->GetBinContent(pt); // <<B>>	 // sum of weights for correlation:
          Double_t sumOfWeightsForA = GetSumPro(fFlowGFIntCorPro[h][i],pt); // sum_{i=1}^{N} w_{<A>}
          Double_t sumOfWeightsForB = GetSumPro(fFlowGFIntCorPro[h][k],pt); // sum_{i=1}^{N} w_{<B>} // products for correlations:
          Double_t AB = fFlowGFIntCovPro[h][i][k]->GetBinContent(pt); // <<A><B>>  // sum of weights for product of correlation:
          Double_t productOfWeightsForAB = GetSumPro(fFlowGFIntCovPro[h][i][k],pt); // sum_{i=1}^{N} w_{<A>}w_{<B>} 	// <A>,<B>:
          Double_t term1 = productOfWeightsForAB;
          Double_t term2 = sumOfWeightsForA;
          Double_t term3 = sumOfWeightsForB;
          if(term2*term3>0.)
          {
            Double_t denominator = 1.-term1/(term2*term3);
            Double_t prefactor = term1/(term2*term3);
            if(TMath::Abs(denominator)>1.e-6)
            {
              Double_t covAB = (AB-A*B)/denominator;
              Double_t wCovAB = covAB*prefactor;
              fFlowGFIntCovHist[h][i][k]->SetBinContent(pt,wCovAB);
            }
          }
        } // end of for(Int_t pt=1; pt<=fFlowGFIntCovPro[h][i][k]->GetNbinsX(); pt++)
        fFlowGFIntCovPro[h][i][k]->GetXaxis()->SetRange(1,fFlowGFIntCovPro[h][i][k]->GetNbinsX());
        fFlowGFIntCorPro[h][k]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][k]->GetNbinsX());
      } // end of for(Int_t k=0; k<fkFlowGFNOrde; k++)
      fFlowGFIntCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorPro[h][i]->GetNbinsX());
    }
  }

  for(Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t pt=1; pt<=fFlowGFIntCorHist[h][0]->GetNbinsX(); pt++) { // Correlations: 
    
      Double_t two = fFlowGFIntCorHist[h][0]->GetBinContent(pt); // <<2>>
      Double_t four = fFlowGFIntCorHist[h][1]->GetBinContent(pt); // <<4>>
      Double_t six = fFlowGFIntCorHist[h][2]->GetBinContent(pt); // <<6>>
      Double_t eight = fFlowGFIntCorHist[h][3]->GetBinContent(pt); // <<8>> // Statistical errors of average 2-, 4-, 6- and 8-particle azimuthal correlations:
      Double_t twoError = fFlowGFIntCorHist[h][0]->GetBinError(pt); // statistical error of <2>
      Double_t fourError = fFlowGFIntCorHist[h][1]->GetBinError(pt); // statistical error of <4>
      Double_t sixError = fFlowGFIntCorHist[h][2]->GetBinError(pt); // statistical error of <6>
      Double_t eightError = fFlowGFIntCorHist[h][3]->GetBinError(pt); // statistical error of <8> // Q-cumulants:
      Double_t qc1 = 0.;
      Double_t qc2 = 0.; // QC{2}
      Double_t qc4 = 0.; // QC{4}
      Double_t qc6 = 0.; // QC{6}
      Double_t qc8 = 0.; // QC{8}
      if(TMath::Abs(two) > 0.){qc2 = two;}
      if(TMath::Abs(four) > 0.){qc4 = four-2.*pow(two,2.);}
      if(TMath::Abs(six) > 0.){qc6 = six-9.*two*four+12.*pow(two,3.);}
      if(TMath::Abs(eight) > 0.){qc8 = eight-16.*two*six-18.*pow(four,2.)+144.*pow(two,2.)*four-144.*pow(two,4.);} // Statistical errors of Q-cumulants:
      Double_t qc2Error = 0.;
      Double_t qc4Error = 0.;
      Double_t qc6Error = 0.;
      Double_t qc8Error = 0.; // Squared statistical errors of Q-cumulants:
      Double_t qc2ErrorSquared = 0.;
      Double_t qc4ErrorSquared = 0.;
      Double_t qc6ErrorSquared = 0.;
      Double_t qc8ErrorSquared = 0.; // covariances:
      Double_t wCov24 = fFlowGFIntCovHist[h][0][1]->GetBinContent(pt);
      Double_t wCov26 = fFlowGFIntCovHist[h][0][2]->GetBinContent(pt);
      Double_t wCov28 = fFlowGFIntCovHist[h][0][3]->GetBinContent(pt);
      Double_t wCov46 = fFlowGFIntCovHist[h][1][2]->GetBinContent(pt);
      Double_t wCov48 = fFlowGFIntCovHist[h][1][3]->GetBinContent(pt);
      Double_t wCov68 = fFlowGFIntCovHist[h][2][3]->GetBinContent(pt); // Statistical error of QC{2}:
      qc2Error = twoError; // Statistical error of QC{4}:
      qc4ErrorSquared = 16.*pow(two,2.)*pow(twoError,2.)+pow(fourError,2.)
      - 8.*two*wCov24;
      if(qc4ErrorSquared>0.) {
        qc4Error = pow(qc4ErrorSquared,0.5);
      } // Statistical error of QC{6}:
      qc6ErrorSquared = 81.*pow(4.*pow(two,2.)-four,2.)*pow(twoError,2.)
      + 81.*pow(two,2.)*pow(fourError,2.)
      + pow(sixError,2.)
      - 162.*two*(4.*pow(two,2.)-four)*wCov24
      + 18.*(4.*pow(two,2.)-four)*wCov26
      - 18.*two*wCov46;
      if(qc6ErrorSquared>0.) {
        qc6Error = pow(qc6ErrorSquared,0.5);
      } // Statistical error of QC{8}:
      qc8ErrorSquared = 256.*pow(36.*pow(two,3.)-18.*four*two+six,2.)*pow(twoError,2.)
      + 1296.*pow(4.*pow(two,2.)-four,2.)*pow(fourError,2.)
      + 256.*pow(two,2.)*pow(sixError,2.)
      + pow(eightError,2.)
      - 1152.*(36.*pow(two,3.)-18.*four*two+six)*(4.*pow(two,2.)-four)*wCov24
      + 512.*two*(36.*pow(two,3.)-18.*four*two+six)*wCov26
      - 32.*(36.*pow(two,3.)-1/8.*four*two+six)*wCov28
      - 1152.*two*(4.*pow(two,2.)-four)*wCov46
      + 72.*(4.*pow(two,2.)-four)*wCov48
      - 32.*two*wCov68;
      if(qc8ErrorSquared>0.) {
        qc8Error = pow(qc8ErrorSquared,0.5);
      }      // Store the cumulants:
      fFlowGFIntCumHist[h][0]->SetBinContent(pt,qc2);
      fFlowGFIntCumHist[h][0]->SetBinError(pt,qc2Error);
      fFlowGFIntCumHist[h][1]->SetBinContent(pt,qc4);
      fFlowGFIntCumHist[h][1]->SetBinError(pt,qc4Error);
      fFlowGFIntCumHist[h][2]->SetBinContent(pt,qc6);
      fFlowGFIntCumHist[h][2]->SetBinError(pt,qc6Error);
      fFlowGFIntCumHist[h][3]->SetBinContent(pt,qc8);
      fFlowGFIntCumHist[h][3]->SetBinError(pt,qc8Error); // Reference flow estimates:
      Double_t v2 = 0.; // v{2,QC}
      Double_t v4 = 0.; // v{4,QC}
      Double_t v6 = 0.; // v{6,QC}
      Double_t v8 = 0.; // v{8,QC} // Reference flow statistical errors:
      Double_t v2Error = 0.; // v{2,QC} stat. error
      Double_t v4Error = 0.; // v{4,QC} stat. error
      Double_t v6Error = 0.; // v{6,QC} stat. error
      Double_t v8Error = 0.; // v{8,QC} stat. error // calculate flow
      if(qc2>=0.){v2 = pow(qc2,0.5);}
      if(qc4<=0.){v4 = pow(-1.*qc4,1./4.);}
      if(qc6>=0.){v6 = pow((1./4.)*qc6,1./6.);}
      if(qc8<=0.){v8 = pow((-1./33.)*qc8,1./8.);} // Calculate stat. error for reference flow estimates from stat. error of Q-cumulants:
      if(qc2>0.){v2Error = (1./2.)*pow(qc2,-0.5)*qc2Error;}
      if(qc4<0.){v4Error = (1./4.)*pow(-qc4,-3./4.)*qc4Error;}
      if(qc6>0.){v6Error = (1./6.)*pow(2.,-1./3.)*pow(qc6,-5./6.)*qc6Error;}
      if(qc8<0.){v8Error = (1./8.)*pow(33.,-1./8.)*pow(-qc8,-7./8.)*qc8Error;}	// Store the results:
      if(qc2>0.) {
        fFlowGFIntFinalHist[h][0]->SetBinContent(pt,v2);
        fFlowGFIntFinalHist[h][0]->SetBinError(pt,v2Error);
	 }
      if(qc4<0.) {
        fFlowGFIntFinalHist[h][1]->SetBinContent(pt,v4);
        fFlowGFIntFinalHist[h][1]->SetBinError(pt,v4Error);
	 }
      if(qc6>0.) {
        fFlowGFIntFinalHist[h][2]->SetBinContent(pt,v6);
        fFlowGFIntFinalHist[h][2]->SetBinError(pt,v6Error);
	}
      if(qc8<0.) {
        fFlowGFIntFinalHist[h][3]->SetBinContent(pt,v8);
        fFlowGFIntFinalHist[h][3]->SetBinError(pt,v8Error);
	 }
     }//here ends the pt loop
  }

  // MIXED HARMONICS ***********************************************************

  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNHarm; i++) {

      for(Int_t pt=1; pt<=fFlowGFMixedCorPro[h][i]->GetNbinsX(); pt++) {

        Double_t stats[6]={0.};
        fFlowGFMixedCorPro[h][i]->GetXaxis()->SetRange(pt,pt);
        fFlowGFMixedCorPro[h][i]->GetStats(stats);
        Double_t SumWeig   = stats[0];
        Double_t SumWeigSq  = stats[1];
        Double_t SumTwo  = stats[4];
        Double_t SumTwoSq = stats[5];

        if(SumWeig>0.) {
          Double_t Corr = SumTwo/SumWeig;
          Double_t SqCorr = SumTwoSq/SumWeig;
          Double_t Weig = SumWeig;
          Double_t SqWeig = SumWeigSq;
          Double_t spread=0., termA=0., termB=0.;
          if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
          if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
          if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
          Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
          if(CorrErr) {
            fFlowGFMixedCorHist[h][i]->SetBinContent(pt,Corr);
            fFlowGFMixedCorHist[h][i]->SetBinError(pt,CorrErr);
          }
        }

      } // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)
      fFlowGFMixedCorPro[h][i]->GetXaxis()->SetRange(1,fFlowGFMixedCorPro[h][i]->GetNbinsX());
    }
  }

  for (Int_t h=0; h<fkFlowGFNHarm; h++) {
    for(Int_t i=0; i<fkFlowGFNHarm; i++) {

      for(Int_t pt=1; pt<=fFlowGFIntCorHist[h][0]->GetNbinsX(); pt++) {
        // Correlations:
        Double_t twoA = fFlowGFIntCumHist[h][0]->GetBinContent(pt); // <<2A>>
        Double_t twoB = fFlowGFIntCumHist[i][0]->GetBinContent(pt); // <<2B>>
        Double_t four = fFlowGFMixedCorHist[h][i]->GetBinContent(pt); // <<4>>
        // Statistical errors:
        Double_t twoAError = fFlowGFIntCumHist[h][0]->GetBinError(pt); // statistical error of <2A>
        Double_t twoBError = fFlowGFIntCumHist[i][0]->GetBinError(pt); // statistical error of <2B>
        Double_t fourError = fFlowGFMixedCorHist[h][i]->GetBinError(pt); // statistical error of <4>
        // Symmetric Cumulants:
        Double_t SC = four - twoA*twoB;
        Double_t SCErrorSquared = pow(fourError,2.) + pow(twoA*twoBError,2.) + pow(twoB*twoAError,2.); // TBI
        // Store the results:
        if(SCErrorSquared>0.) {
          fFlowGFMixedFinalHist[h][i]->SetBinContent(pt,SC);
          fFlowGFMixedFinalHist[h][i]->SetBinError(pt,pow(SCErrorSquared,0.5));
        }
      }  // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)

    }
  }

  // in wide pt bins
  for(Int_t s=0; s<fkGFPtB; s++) {

    if(!fFlowGFIntCorProPtB[0][0][0]) continue;

    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
      for(Int_t i=0; i<fkFlowGFNOrde; i++) {
        for(Int_t pt=1; pt<=fFlowGFIntCorProPtB[s][h][i]->GetNbinsX(); pt++) {
          Double_t stats[6]={0.};
          fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(pt,pt);
          fFlowGFIntCorProPtB[s][h][i]->GetStats(stats);
          Double_t SumWeig   = stats[0];
          Double_t SumWeigSq  = stats[1];
          Double_t SumTwo  = stats[4];
          Double_t SumTwoSq = stats[5];
          if(SumWeig>0.) {
            Double_t Corr = SumTwo/SumWeig;
            Double_t SqCorr = SumTwoSq/SumWeig;
            Double_t Weig = SumWeig;
            Double_t SqWeig = SumWeigSq;
            Double_t spread=0., termA=0., termB=0.;
            if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
            if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
            if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
            Double_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
            if(CorrErr) {
              fFlowGFIntCorHistPtB[s][h][i]->SetBinContent(pt,Corr);
              fFlowGFIntCorHistPtB[s][h][i]->SetBinError(pt,CorrErr);
            }
          }
        } // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)
        fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][i]->GetNbinsX());
      } // end of for(Int_t i=0; i<fkFlowGFNOrde; i++)
    }

    for (Int_t h=0; h<fkFlowGFNHarm; h++) {
      for(Int_t i=0; i<fkFlowGFNOrde; i++) {
        for(Int_t k=0; k<fkFlowGFNOrde; k++) {
          for(Int_t pt=1; pt<=fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX(); pt++) {
            // correlations:
            Double_t A = fFlowGFIntCorHistPtB[s][h][i]->GetBinContent(pt); // <<A>>
            Double_t B = fFlowGFIntCorHistPtB[s][h][k]->GetBinContent(pt); // <<B>>
            // sum of weights for correlation:
            Double_t sumOfWeightsForA = GetSumPro(fFlowGFIntCorProPtB[s][h][i],pt); // sum_{i=1}^{N} w_{<A>}
            Double_t sumOfWeightsForB = GetSumPro(fFlowGFIntCorProPtB[s][h][k],pt); // sum_{i=1}^{N} w_{<B>}
            // products for correlations:
            Double_t AB = fFlowGFIntCovProPtB[s][h][i][k]->GetBinContent(pt); // <<A><B>>
            // sum of weights for product of correlation:
            Double_t productOfWeightsForAB = GetSumPro(fFlowGFIntCovProPtB[s][h][i][k],pt); // sum_{i=1}^{N} w_{<A>}w_{<B>}
            // <A>,<B>:
            Double_t term1 = productOfWeightsForAB;
            Double_t term2 = sumOfWeightsForA;
            Double_t term3 = sumOfWeightsForB;
            if(term2*term3>0.)
            {
              Double_t denominator = 1.-term1/(term2*term3);
              Double_t prefactor = term1/(term2*term3);
              if(TMath::Abs(denominator)>1.e-6)
              {
                Double_t covAB = (AB-A*B)/denominator;
                Double_t wCovAB = covAB*prefactor;
                fFlowGFIntCovHistPtB[s][h][i][k]->SetBinContent(pt,wCovAB);
              }
            }
          } // end of for(Int_t pt=1; pt<=fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX(); pt++)
          fFlowGFIntCovProPtB[s][h][i][k]->GetXaxis()->SetRange(1,fFlowGFIntCovProPtB[s][h][i][k]->GetNbinsX());
          fFlowGFIntCorProPtB[s][h][k]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][k]->GetNbinsX());
        } // end of for(Int_t k=0; k<fkFlowGFNOrde; k++)
        fFlowGFIntCorProPtB[s][h][i]->GetXaxis()->SetRange(1,fFlowGFIntCorProPtB[s][h][i]->GetNbinsX());
      }
    }
	
/*    for(Int_t h = 0; h<fkFlowGFNHarm; h++){
	//for(Int_t i=0; i<fkFlowGFNOrde; i++){//for third argument now 0
	    for(Int_t c = 1; c<=fFlowGFIntCorHistPtB[s][h][0]->GetNbinsX(); c++){//c for centbin  
	//correlation
		Double_t two = fFlowGFIntCorHistPtB[s][h][0]->GetBinContent(c);//<<2>>
		Double_t four = fFlowGFIntCorHistPtB[s][h][1]->GetBinContent(c);//<<4>>
		Double_t six = fFlowGFIntCorHistPtB[s][h][2]->GetBinContent(c);//<<6>>
		Double_t eight = fFlowGFIntCorHistPtB[s][h][3]->GetBinContent(c);//<<8>>
//from here reduce calculation to the two particle correlator and ask tomorrow how many I should include - are there only even ones that I can/should use + my change of k to include n = 1 may mess the calculation for higher correlatons up ... mayb adjust that part accordingly ??
//Different  order  cumulants  provide  independent  estimates  for  the  same  reference  harmonic this means with k-> k+1 we keep harmonic 1,2 and 4 -> maybe make [0] separate and start the loop frm one ->thus harmonic 1,2,4,6 would be present - I assume v3 would follow... so changing this s appropriate ... but what about the hr+2 part ?  
	//statistical errors of average 2-
		Double_t twoError = fFlowGFIntCorHistPtB[s][h][0]->GetBinError(c);
	//Q-cumulants:
	Double_t qc2 = 0.; //QC{2}
	if(TMath::Abs(two)>0.){qc2 = two;}
	//declaration statistical error of Q-cumulants
	Double_t qc2Error = 0.;
	//declaration squared statistical error of Q-cumulants not for qc2
	
	//declaration covariances - not for qc2
	
	qc2Error = twoError;	
	//declaration reference flow estimates
	Double_t v2 = 0.;// v{2}_n = sqrt(qc{2}_n) -> integrated
// aaaahhh... we use [0] this is where the harmonic n is fixed ...uff that took ages
	//declaration reference flow stat. errors
	Double_t v2Error = 0.;//v{2} stat. error
	//calculate flow:
	if(qc2>=0.){v2 = pow(qc2, 0.5);}
	//calculate stat. error for reference flow estimate from stat. error of Q-cumulants 
	if(qc2>=0.){v2Error = (1./2.)*pow(qc2,-0.5)*qc2Error;}
	//store the results
	Double_t pt = hPt[s];
//	cout<<"pt "<<pt<<"c "<<c<<endl;
	if(qc2>0.){
	//final hist to get the centrality dependent one with fill in (c,v2)
	//here in pt dependence
	cout<<"hPt[s] "<<pt<<" v2 "<<v2<<" c "<<c<<" h "<<h<<" s "<<s<<endl;
	//v1pt[h][0]->SetBinContent(pt,v2);
	v1pt[h][0]->SetBinError(pt, v2Error);
	}
	
	}//end of c<=fFlowGFIntCorHistPtB[s][h][0]->GetNBinsX()
	//}//end of i<fkFlowGFNOrde	
    }//end h<FlowGFHarm  
*/
  }//end s<fkGFPtB loop
  cout << "*************************************" << "\n";
  cout << "\n";

} // end of void CalculateFlowCME::FinalizeFlowGF()

//=====================================================================================================

void CalculateFlowCME::FinalizeFlowQC()
{
 Int_t ch = 0.;
	std::cout << "Finalizing Flow QC"<<'\n'; 
	for(Int_t hr=0; hr<fFlowNHarm; hr++) {
		// Pt-INTEGRATED
		// STORE IN HISTOGRAMS
		// 2- and 4-particle cumulants
		for(Int_t j=0; j<fkFlowQCnIntCorPro; j++) {
			for(Int_t pt=1;pt<=fFlowQCIntCorPro[ch][hr][j]->GetNbinsX();pt++) {
				Double_t stats[6]={0.};
				fFlowQCIntCorPro[ch][hr][j]->GetXaxis()->SetRange(pt,pt);
				fFlowQCIntCorPro[ch][hr][j]->GetStats(stats);	
				LongDouble_t SumWeig   = stats[0];
				LongDouble_t SumWeigSq  = stats[1];
				LongDouble_t SumTwo  = stats[4];
				LongDouble_t SumTwoSq = stats[5];
				if(SumWeig>0.){
				  LongDouble_t Corr = SumTwo/SumWeig;
				  LongDouble_t SqCorr = SumTwoSq/SumWeig;
				  LongDouble_t Weig = SumWeig;
				  LongDouble_t SqWeig = SumWeigSq;
				  LongDouble_t spread=0., termA=0., termB=0.;
				  if(SqCorr-pow(Corr,2.)>=0.) { spread = pow(SqCorr-pow(Corr,2.),0.5); }
				  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
				  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
				  LongDouble_t CorrErr = termA*spread*termB; // final error (unbiased estimator for standard deviation)
				  if(CorrErr) {
					fFlowQCIntCorHist[ch][hr][j]->SetBinContent(pt,Corr);
					fFlowQCIntCorHist[ch][hr][j]->SetBinError(pt,CorrErr);
				  }
				}
			} // end of for(Int_t pt=1;pt<=100;pt++)
			fFlowQCIntCorPro[ch][hr][j]->GetXaxis()->SetRange(1,fFlowQCIntCorPro[ch][hr][j]->GetNbinsX());
		} // end of for(Int_t j=0; j<5; j++)

	} // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
	
	
	
	// FINALISE (calculate flow)

	for(Int_t hr=0; hr<fFlowNHarm; hr++) {
		// calculate covariance
		for(Int_t pt=1; pt<=fFlowQCIntCorHist[ch][hr][0]->GetNbinsX(); pt++) {
			// average reduced correlations:
			Double_t two = fFlowQCIntCorHist[ch][hr][0]->GetBinContent(pt); // <<2>>
			Double_t four = fFlowQCIntCorHist[ch][hr][1]->GetBinContent(pt); // <<4>>
			// sum of weights for reduced correlation:
			Double_t sumOfWeightsForTwo = GetSumPro(fFlowQCIntCorPro[ch][hr][0],pt); // sum_{i=1}^{N} w_{<2>}
			Double_t sumOfWeightsForFour = GetSumPro(fFlowQCIntCorPro[ch][hr][1],pt); // sum_{i=1}^{N} w_{<4>}
			// product of weights for reduced correlation:
			Double_t productOfWeightsForTwoFour = GetSumPro(fFlowQCIntCorPro[ch][hr][2],pt); // sum_{i=1}^{N} w_{<2>}w_{<4>}
			// products for differential flow:
			Double_t twoFour = fFlowQCIntCorHist[ch][hr][2]->GetBinContent(pt); // <<2><4>>

			// <2>,<4>:
			Double_t term1 = productOfWeightsForTwoFour;
			Double_t term2 = sumOfWeightsForTwo;
			Double_t term3 = sumOfWeightsForFour;
			if(term2*term3>0.)
			{
				Double_t denominator = 1.-term1/(term2*term3);
				Double_t prefactor = term1/(term2*term3);
				if(TMath::Abs(denominator)>1.e-6)
				{
					Double_t covTwoFour = (twoFour-two*four)/denominator;
					Double_t wCovTwoFour = covTwoFour*prefactor;
					fFlowQCIntCorHist[ch][hr][2]->SetBinContent(pt,wCovTwoFour);
				}
			}
		} // end of for(Int_t pt=1;pt<=fPtDiffNBins;pt++)

		// 2- and 4-particle cumulants
		for(Int_t pt=1; pt<=fFlowQCIntCorHist[ch][hr][0]->GetNbinsX(); pt++) {
			Double_t QC2    = fFlowQCIntCorHist[ch][hr][0]->GetBinContent(pt);
			Double_t QC2E   = fFlowQCIntCorHist[ch][hr][0]->GetBinError(pt);
			Double_t QC4    = fFlowQCIntCorHist[ch][hr][1]->GetBinContent(pt);
			Double_t QC4E   = fFlowQCIntCorHist[ch][hr][1]->GetBinError(pt);
			Double_t wCov24 = fFlowQCIntCorHist[ch][hr][2]->GetBinContent(pt);
			Double_t Cn2 = QC2;
			Double_t Cn2E = QC2E;
			Double_t Cn4 = QC4-2.*QC2*QC2;
			Double_t Cn4Esq = 16.*pow(QC2,2.)*pow(QC2E,2) + pow(QC4E,2.) - 8.*QC2*wCov24;

			fFlowQCIntCumHist[ch][hr][0]->SetBinContent(pt,Cn2);
			fFlowQCIntCumHist[ch][hr][0]->SetBinError(pt,Cn2E);

			if(Cn4Esq>0.) {
				Double_t Cn4E = pow(Cn4Esq,0.5);
				fFlowQCIntCumHist[ch][hr][1]->SetBinContent(pt,Cn4);
				fFlowQCIntCumHist[ch][hr][1]->SetBinError(pt,Cn4E);
				if (Cn4<0.) {
					Double_t Flow4 = pow(fabs(Cn4),0.25);
					Double_t Flow4E = fabs(Flow4/(4.*Cn4))*Cn4E;
					fFlowQCIntCorHist[ch][hr][2]->SetBinContent(pt,Flow4);
					fFlowQCIntCorHist[ch][hr][2]->SetBinError(pt,Flow4E);
				} else {
					fFlowQCIntCorHist[ch][hr][2]->SetBinContent(pt,0.);
					fFlowQCIntCorHist[ch][hr][2]->SetBinError(pt,0.);
				}
			}
		}
	} // end of for(Int_t hr=0; hr<fFlowNHarm; hr++)
}

void CalculateFlowCME::FinalizeCMW()
{
  // fHistAveV2PosAch //<<v2^{+}A>>
  for(Int_t i = 1; i<=fHistAveV2PosAch->GetNbinsX(); i++) {
	  cout<<"fHistAveV2PosAch->GetBinContent("<<i<<") = "<<fHistAveV2PosAch->GetBinContent(i)<<endl;
  }
  cout<<"fHistAveV2PosAch->FindBin("<<fCentralityEBE<<") == "<<fHistAveV2PosAch->FindBin(fCentralityEBE)<<endl;
  Double_t statsAveV2PosAch[6]={0.};
  fHistAveV2PosAch->GetXaxis()->SetRange(fHistAveV2PosAch->FindBin(fCentralityEBE),fHistAveV2PosAch->FindBin(fCentralityEBE));
  fHistAveV2PosAch->GetStats(statsAveV2PosAch);
  Double_t SumWeigAveV2PosAch   = statsAveV2PosAch[0]; // sumw
  Double_t SumWeigSqAveV2PosAch  = statsAveV2PosAch[1]; // sumw2
  Double_t SumTwoAveV2PosAch  = statsAveV2PosAch[4]; // sumwy
  Double_t SumTwoSqAveV2PosAch = statsAveV2PosAch[5]; // sumwy2
  
  Double_t CorrAveV2PosAch = 0;
  Double_t CorrErrAveV2PosAch = 0;
  if(SumWeigAveV2PosAch>0.) {
	  CorrAveV2PosAch = SumTwoAveV2PosAch/SumWeigAveV2PosAch; // sumwy/sumw
	  Double_t SqCorr = SumTwoSqAveV2PosAch/SumWeigAveV2PosAch;
	  Double_t Weig = SumWeigAveV2PosAch;
	  Double_t SqWeig = SumWeigSqAveV2PosAch;
	  Double_t spread=0., termA=0., termB=0.;
	  if(SqCorr-pow(CorrAveV2PosAch,2.)>=0.) { spread = pow(SqCorr-pow(CorrAveV2PosAch,2.),0.5); }
	  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
	  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
	  CorrErrAveV2PosAch = termA*spread*termB; // final error (unbiased estimator for standard deviation)
  }

  // fHistAveV2Pos // <<v2^{+}>>
  Double_t statsAveV2Pos[6]={0.};
  fHistAveV2Pos->GetXaxis()->SetRange(fHistAveV2Pos->FindBin(fCentralityEBE),fHistAveV2Pos->FindBin(fCentralityEBE));
  fHistAveV2Pos->GetStats(statsAveV2Pos);
  Double_t SumWeigAveV2Pos   = statsAveV2Pos[0];
  Double_t SumWeigSqAveV2Pos  = statsAveV2Pos[1];
  Double_t SumTwoAveV2Pos  = statsAveV2Pos[4]; // sumwy
  Double_t SumTwoSqAveV2Pos = statsAveV2Pos[5]; // sumwy2
  
  Double_t CorrAveV2Pos = 0;
  Double_t CorrErrAveV2Pos = 0;
  if(SumWeigAveV2Pos>0.) {
	  CorrAveV2Pos = SumTwoAveV2Pos/SumWeigAveV2Pos; // sumwy/sumw
	  Double_t SqCorr = SumTwoSqAveV2Pos/SumWeigAveV2Pos;
	  Double_t Weig = SumWeigAveV2Pos;
	  Double_t SqWeig = SumWeigSqAveV2Pos;
	  Double_t spread=0., termA=0., termB=0.;
	  if(SqCorr-pow(CorrAveV2Pos,2.)>=0.) { spread = pow(SqCorr-pow(CorrAveV2Pos,2.),0.5); }
	  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
	  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
	  CorrErrAveV2Pos = termA*spread*termB; // final error (unbiased estimator for standard deviation)
  }

  // fHistAveV2NegAch // <<v2^{-}A>>
  Double_t statsAveV2NegAch[6]={0.};
  fHistAveV2NegAch->GetXaxis()->SetRange(fHistAveV2NegAch->FindBin(fCentralityEBE),fHistAveV2NegAch->FindBin(fCentralityEBE));
  fHistAveV2NegAch->GetStats(statsAveV2NegAch);
  Double_t SumWeigAveV2NegAch   = statsAveV2NegAch[0];
  Double_t SumWeigSqAveV2NegAch  = statsAveV2NegAch[1];
  Double_t SumTwoAveV2NegAch  = statsAveV2NegAch[4]; // sumwy
  Double_t SumTwoSqAveV2NegAch = statsAveV2NegAch[5]; // sumwy2
  
  Double_t CorrAveV2NegAch = 0;
  Double_t CorrErrAveV2NegAch = 0;
  if(SumWeigAveV2NegAch>0.) {
	  CorrAveV2NegAch = SumTwoAveV2NegAch/SumWeigAveV2NegAch; // sumwy/sumw
	  Double_t SqCorr = SumTwoSqAveV2NegAch/SumWeigAveV2NegAch;
	  Double_t Weig = SumWeigAveV2NegAch;
	  Double_t SqWeig = SumWeigSqAveV2NegAch;
	  Double_t spread=0., termA=0., termB=0.;
	  if(SqCorr-pow(CorrAveV2NegAch,2.)>=0.) { spread = pow(SqCorr-pow(CorrAveV2NegAch,2.),0.5); }
	  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
	  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
	  CorrErrAveV2NegAch = termA*spread*termB; // final error (unbiased estimator for standard deviation)
  }

  // fHistAveV2Neg // <<v2^{-}>>
  Double_t statsAveV2Neg[6]={0.};
  fHistAveV2Neg->GetXaxis()->SetRange(fHistAveV2Neg->FindBin(fCentralityEBE),fHistAveV2Neg->FindBin(fCentralityEBE));
  fHistAveV2Neg->GetStats(statsAveV2Neg);
  Double_t SumWeigAveV2Neg   = statsAveV2Neg[0];
  Double_t SumWeigSqAveV2Neg  = statsAveV2Neg[1];
  Double_t SumTwoAveV2Neg  = statsAveV2Neg[4]; // sumwy
  Double_t SumTwoSqAveV2Neg = statsAveV2Neg[5]; // sumwy2
  
  Double_t CorrAveV2Neg = 0;
  Double_t CorrErrAveV2Neg = 0;
  if(SumWeigAveV2Neg>0.) {
	  CorrAveV2Neg = SumTwoAveV2Neg/SumWeigAveV2Neg; // sumwy/sumw
	  Double_t SqCorr = SumTwoSqAveV2Neg/SumWeigAveV2Neg;
	  Double_t Weig = SumWeigAveV2Neg;
	  Double_t SqWeig = SumWeigSqAveV2Neg;
	  Double_t spread=0., termA=0., termB=0.;
	  if(SqCorr-pow(CorrAveV2Neg,2.)>=0.) { spread = pow(SqCorr-pow(CorrAveV2Neg,2.),0.5); }
	  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
	  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
	  CorrErrAveV2Neg = termA*spread*termB; // final error (unbiased estimator for standard deviation)
  }

  // fHistAChrgVsCent // <A>
  Double_t statsAChrgVsCent[6]={0.};
  fHistAChrgVsCent->GetXaxis()->SetRange(fHistAChrgVsCent->FindBin(fCentralityEBE),fHistAChrgVsCent->FindBin(fCentralityEBE));
  fHistAChrgVsCent->GetStats(statsAChrgVsCent);
  Double_t SumWeigAChrgVsCent   = statsAChrgVsCent[0];
  Double_t SumWeigSqAChrgVsCent  = statsAChrgVsCent[1];
  Double_t SumTwoAChrgVsCent  = statsAChrgVsCent[4]; // sumwy
  Double_t SumTwoSqAChrgVsCent = statsAChrgVsCent[5]; // sumwy2
  
  Double_t CorrAChrgVsCent = 0;
  Double_t CorrErrAChrgVsCent = 0;
  if(SumWeigAChrgVsCent>0.) {
	  CorrAChrgVsCent = SumTwoAChrgVsCent/SumWeigAChrgVsCent; // sumwy/sumw
	  Double_t SqCorr = SumTwoSqAChrgVsCent/SumWeigAChrgVsCent;
	  Double_t Weig = SumWeigAChrgVsCent;
	  Double_t SqWeig = SumWeigSqAChrgVsCent;
	  Double_t spread=0., termA=0., termB=0.;
	  if(SqCorr-pow(CorrAChrgVsCent,2.)>=0.) { spread = pow(SqCorr-pow(CorrAChrgVsCent,2.),0.5); }
	  if(TMath::Abs(Weig)>0.) { termA = (pow(SqWeig,0.5)/Weig); }
	  if(1.-pow(termA,2.)>0.) { termB = 1./pow(1.-pow(termA,2.),0.5); }
	  CorrErrAChrgVsCent = termA*spread*termB; // final error (unbiased estimator for standard deviation)
  }

cout<<"==> CorrAveV2PosAch = "<<CorrAveV2PosAch<<endl;
cout<<"==> CorrAveV2Pos    = "<<CorrAveV2Pos<<endl;
cout<<"==> CorrAveV2NegAch = "<<CorrAveV2NegAch<<endl;
cout<<"==> CorrAveV2Neg    = "<<CorrAveV2Neg<<endl;
cout<<"==> CorrAChrgVsCent = "<<CorrAChrgVsCent<<endl;

  ///------- Calculation of three-particls correlator: <<cos[2(phi1-Phi2)]A>>/sqrt(<<cos[2(phi1-Phi2)]>>)-sqrt(<<cos[2(phi1-Phi2)]>>)<<A>>


  //// <v2^{+}A>-<A><v2> = <<v2^{+}A>>/sqrt(<<v2^{+}>>)-sqrt(<<v2^{+}>>)<A>
  //Double_t threePartCorrelatorPos = CorrAveV2PosAch/sqrt(CorrAveV2Pos) - sqrt(CorrAveV2Pos)*CorrAChrgVsCent;
  //Double_t threePartCorrelatorPosErr = sqrt(pow(CorrErrAveV2PosAch/sqrt(CorrAveV2Pos),2)+pow(CorrAveV2PosAch/(2*pow(CorrAveV2Pos,1.5))*CorrErrAveV2Pos,2)+pow(CorrAChrgVsCent/(2*sqrt(CorrAveV2Pos))*CorrErrAveV2Pos,2)+pow(sqrt(CorrAveV2Pos)*CorrErrAChrgVsCent,2));
  
  //// <v2^{-}A>-<A><v2> = <<v2^{-}A>>/sqrt(<<v2^{-}>>)-sqrt(<<v2^{-}>>)<A>
  //Double_t threePartCorrelatorNeg = CorrAveV2NegAch/sqrt(CorrAveV2Neg) - sqrt(CorrAveV2Neg)*CorrAChrgVsCent;
  //Double_t threePartCorrelatorNegErr = sqrt(pow(CorrErrAveV2NegAch/sqrt(CorrAveV2Neg),2)+pow(CorrAveV2NegAch/(2*pow(CorrAveV2Neg,1.5))*CorrErrAveV2Neg,2)+pow(CorrAChrgVsCent/(2*sqrt(CorrAveV2Neg))*CorrErrAveV2Neg,2)+pow(sqrt(CorrAveV2Neg)*CorrErrAChrgVsCent,2));
  
  Double_t threePartCorrelatorPos = 1;
  Double_t threePartCorrelatorPosErr = 1;
  Double_t threePartCorrelatorNeg = 1;
  Double_t threePartCorrelatorNegErr = 1;
  
  cout<<"threePartCorrelatorPos="<<threePartCorrelatorPos<<endl;
  cout<<"threePartCorrelatorPosErr="<<threePartCorrelatorPosErr<<endl;
  cout<<"threePartCorrelatorNeg="<<threePartCorrelatorNeg<<endl;
  cout<<"threePartCorrelatorNegErr="<<threePartCorrelatorNegErr<<endl;

  fCMWThreeParticleCorrelator[0]->SetBinContent(fCenBin+1, threePartCorrelatorPos);
  fCMWThreeParticleCorrelator[0]->SetBinError(fCenBin+1, threePartCorrelatorPosErr);
  
  fCMWThreeParticleCorrelator[1]->SetBinContent(fCenBin+1, threePartCorrelatorNeg);
  fCMWThreeParticleCorrelator[1]->SetBinError(fCenBin+1, threePartCorrelatorNegErr);
}
//=====================================================================================================

Int_t CalculateFlowCME::GetCRCCenBin(Double_t Centrality)
{
	Int_t CenBin=-1;
	if (Centrality>0. && Centrality<10.) CenBin=0;
	if (Centrality>10. && Centrality<20.) CenBin=1;
	if (Centrality>20. && Centrality<30.) CenBin=2;
	if (Centrality>30. && Centrality<40.) CenBin=3;
	if (Centrality>40. && Centrality<50.) CenBin=4;
	if (Centrality>50. && Centrality<60.) CenBin=5;
	if (Centrality>60. && Centrality<70.) CenBin=6;
	if (Centrality>70. && Centrality<80.) CenBin=7;
	if (Centrality>80. && Centrality<90.) CenBin=8;
	if (Centrality>90. && Centrality<100.) CenBin=9;
	if (CenBin>=fCRCnCen) CenBin=-1;
	if (fCRCnCen==1) CenBin=0;
	return CenBin;
} // end of CalculateFlowCME::GetCRCCenBin(Double_t Centrality)

//=====================================================================================================

Double_t CalculateFlowCME::GetSumPro(TProfile *pro, Int_t bin)
{
  Double_t stats[6]={0.};
  pro->GetXaxis()->SetRange(bin,bin);
  pro->GetStats(stats);
  pro->GetXaxis()->SetRange(1,pro->GetNbinsX());
  return stats[0];
}

//=====================================================================================================

std::complex<double> CalculateFlowCME::ucN(const Int_t n, const TArrayI& h, Int_t ptb=-1)
{
  TArrayI cnt(n);
  for (Int_t i = 0; i < n; i++) {
    cnt[i] = 1;
  }
  TArrayI hh(h);

  return ucN2(n, hh, cnt, ptb);
}

std::complex<double> CalculateFlowCME::ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, Int_t ptb=-1)//implement for the pos and neg!
{
  Int_t j = n-1;
  std::complex<double> c;
  if(ptb<0) {
    if(h[j] >= 0) {
      c = std::complex<double>((*fReQGF)(h[j],cnt[j]),(*fImQGF)(h[j],cnt[j]));
    } else {
      c = std::complex<double>((*fReQGF)(-h[j],cnt[j]),-(*fImQGF)(-h[j],cnt[j]));
    }
  } else {
    if(h[j] >= 0) {
      c = std::complex<double>((*fReQGFPt[ptb])(h[j],cnt[j]),(*fImQGFPt[ptb])(h[j],cnt[j]));
    } else {
      c = std::complex<double>((*fReQGFPt[ptb])(-h[j],cnt[j]),-(*fImQGFPt[ptb])(-h[j],cnt[j]));
    }
  }

  if (n == 1) return c;

  c *= ucN2(j, h, cnt, ptb);

  if (cnt[j] > 1) return c;

  for (Int_t i = 0; i < (n-1); i++) {
    h[i] += h[j];
    cnt[i] = cnt[i] + 1;
    Double_t factor = 1.*(cnt[i]-1);
    c -= factor * ucN2(j, h, cnt, ptb);
    cnt[i]--;
    h[i] -= h[j];
  }
  return c;
}
