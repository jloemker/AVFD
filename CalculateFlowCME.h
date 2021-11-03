#ifndef CalculateFlowCME_H
#define CalculateFlowCME_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <complex>
#include <cmath>
#include "TList.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectoryFile.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TMath.h"
#include "Event.h"


class CalculateFlowCME 
{
	public:
	
	CalculateFlowCME();
	CalculateFlowCME(const char* name);
	~CalculateFlowCME() = default;
	
	virtual void UserCreateOutputObjects();
    virtual void UserExec();
    virtual void Terminate();
    
    virtual void Make(Event* fEvent);
    void InitializeArraysForFlowQC();
    void InitializeArraysForFlowGF();
    void InitializeArraysForQA();
    void InitializeArraysForCMW();
    void InitializeArraysForFlowFromBW();
    void InitializeArraysForCME();
    virtual void ResetEventByEventQuantities();
    
    virtual void CalculateFlowQC();
    virtual void FinalizeFlowQC();
    virtual void CalculateFlowGF();
    virtual void FinalizeFlowGF();
    virtual void CalculateFlowFromBW();
    virtual void CalculateCMW();
    virtual void FinalizeCMW();
    virtual void CalculateCME();
    virtual void FinalizeCME();
    
    Int_t GetCRCCenBin(Double_t Centrality);
    Double_t GetSumPro(TProfile *pro, Int_t bin);
    
    void SetCentralityEBE(Double_t const c) {this->fCentralityEBE = c;};
	Double_t GetCentralityEBE() const {return this->fCentralityEBE;};
    void SetEvent(Event* e) {this->fEvent = e;};
    TList* GetFlowQCList() {return this->fFlowQCList;};
    TList* GetFlowGFList() {return this->fFlowGFList;};
    TList* GetFlowFromBWList() {return this->fFlowFromBWList;};
    TList* GetCMEList() {return this->fCMEList;};
    TList* GetQAList() {return this->fQAList;};
    TList* GetCMWList() {return this->fCMWList;};
    TList* GetCMWQAList() {return this->fCMWQAList;};
    Event* GetEvent() {return this->fEvent;};
    virtual std::complex<double> ucN(const Int_t n, const TArrayI& h, Int_t ptb);
    virtual std::complex<double> ucN2(const Int_t n, TArrayI& h, TArrayI& cnt, Int_t ptb);
    
    void SetmaxPtCut(Double_t maxPt) {this->maxPtCut = maxPt;};
    void SetminPtCut(Double_t minPt) {this->minPtCut = minPt;};
    void SetminNtrackCut(Double_t minNtrk) {this->minNtracks = minNtrk;};
    void SetmaxEtaCut(Double_t maxEta) {this->maxEtaCut = maxEta;};
    void SetdoQA(Bool_t bflag) {this->doQA = bflag;};
    
    //CMW 
    void SetEtaGapNeg(Double_t negEtaGap) {this->fEtaGapNeg = negEtaGap;};
    void SetEtaGapPos(Double_t posEtaGap) {this->fEtaGapPos = posEtaGap;};
    
	private:
	Event* fEvent;
	
	Int_t fCenBin = -1;
	Int_t fCRCnCen = 10;
	Double_t fCentralityEBE;
	Double_t fCenWeightEbE = 1; // In MC, set to 1 for now.
	Double_t wPhiEta = 1; // In MC, set to 1 for now.
	Double_t wPhi = 1; // In MC, set to 1 for now.
	Double_t wPt = 1; // In MC, set to 1 for now.
	Double_t wEta = 1; // In MC, set to 1 for now.
	Double_t wTrack = 1; // In MC, set to 1 for now.
	Double_t *fCRCPtBins;
	Double_t *PtB;
	//for weight in v1pt and v1eta differential
	Double_t fPtWeight = 1;
	Double_t fEtaWeight = 1;
//for v1(eta)
	Double_t *fCRCEtaBins;
	Bool_t doQA = kFALSE;
	Double_t trkWgt = 1;
	
	// QA Histograms
	TList *fQAList;
	TH1D *fnEvent;
	TH1D *fMultChargedParticlesDistribution;
	TH1D *fPtChargedParticlesDistribution;
	TH1D *fEtaChargedParticlesDistribution;
	TH1D *fPhiChargedParticlesDistribution;
	TH1D *fEtaChargedParticlesDistributionPaperbinning;
	TH1D *fPionsPtSpectra;
	TH1D *fPionsEtaSpectra;
	TH1D *fPionsPhiSpectra;
	TH1D *fPosPionsPtSpectra;
	TH1D *fPosPionsEtaSpectra;
	TH1D *fPosPionsPhiSpectra;
	TH1D *fAntiPionsPtSpectra;
	TH1D *fAntiPionsEtaSpectra;
	TH1D *fAntiPionsPhiSpectra;
	TH1D *fKaonsPtSpectra;
	TH1D *fKaonsEtaSpectra;
	TH1D *fKaonsPhiSpectra;
	TH1D *fPosKaonsPtSpectra;
	TH1D *fPosKaonsEtaSpectra;
	TH1D *fPosKaonsPhiSpectra;
	TH1D *fAntiKaonsPtSpectra;
	TH1D *fAntiKaonsEtaSpectra;
	TH1D *fAntiKaonsPhiSpectra;
	TH1D *fProtonsPtSpectra;
	TH1D *fProtonsEtaSpectra;
	TH1D *fProtonsPhiSpectra;
	TH1D *fPosProtonsPtSpectra;
	TH1D *fPosProtonsEtaSpectra;
	TH1D *fPosProtonsPhiSpectra;
	TH1D *fAntiProtonsPtSpectra;
	TH1D *fAntiProtonsEtaSpectra;
	TH1D *fAntiProtonsPhiSpectra;
	// Flow GF part
	TList *fFlowGFList;
	const static Int_t fkFlowGFNHarm = 4;//Increased from 4 to 5 !
	const static Int_t fkFlowGFNOrde = 4;
	const static Int_t fFlowGFCenBin = 10;
	const static Int_t fkGFPtB = 8;
	TMatrixD *fReQGF; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
	TMatrixD *fImQGF; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
	TMatrixD *fReQGFPt[fkGFPtB]; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
	TMatrixD *fImQGFPt[fkGFPtB]; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
	//for pos/neg and thus \Delta v1
	TMatrixD *fReQGFPtPos[fkGFPtB]; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
        TMatrixD *fImQGFPtPos[fkGFPtB]; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
        TMatrixD *fReQGFPtNeg[fkGFPtB]; // fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
        TMatrixD *fImQGFPtNeg[fkGFPtB]; // fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})

	TProfile *fFlowGFIntCorPro[fkFlowGFNHarm][fkFlowGFNOrde]; //initialized with centrality bins
	TH1D *fFlowGFIntCorHist[fkFlowGFNHarm][fkFlowGFNOrde]; //
	TH1D *fFlowGFIntCumHist[fkFlowGFNHarm][fkFlowGFNOrde]; //
	TH1D *fFlowGFIntFinalHist[fkFlowGFNHarm][fkFlowGFNOrde]; // from pt integrated vn      
	//for differential pt_____________________ 
	TProfile *fFlowGFPtDifCorPro[fkFlowGFNHarm][fkFlowGFNOrde];//initialized with ptdiff bins	
	TH1D *fFlowGFPtDifCorHist[fkFlowGFNHarm][fkFlowGFNOrde];
	TH1D *fFlowGFPtDifCumHist[fkFlowGFNHarm][fkFlowGFNOrde];
	TH1D *fFlowGFPtDifFinalHist[fkFlowGFNHarm][fkFlowGFNOrde];

	TProfile *fFlowGFPtDifCovPro[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde];
	TH1D *fFlowGFPtDifCovHist[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde];
	
	Double_t hPt[fkGFPtB];
	TH1F *hPosPt;
	TH1F *hNegPt;
	//for diffeential eta____________________
	TProfile *fFlowGFEtaDifCorPro[fkFlowGFNHarm][fkFlowGFNOrde];//initialized with Etadiff bins       
        TH1D *fFlowGFEtaDifCorHist[fkFlowGFNHarm][fkFlowGFNOrde];
        TH1D *fFlowGFEtaDifCumHist[fkFlowGFNHarm][fkFlowGFNOrde];
        TH1D *fFlowGFEtaDifFinalHist[fkFlowGFNHarm][fkFlowGFNOrde];

        TProfile *fFlowGFEtaDifCovPro[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde];
        TH1D *fFlowGFEtaDifCovHist[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde];
	
	TH1F *hEta;
	TH1F *hPosEta;
	TH1F *hNegEta;
	//________________________________________
	TProfile *fFlowGFIntCovPro[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
	TH1D *fFlowGFIntCovHist[fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
  
	TProfile *fFlowGFMixedCorPro[fkFlowGFNHarm][fkFlowGFNHarm]; //
	TH1D *fFlowGFMixedCorHist[fkFlowGFNHarm][fkFlowGFNHarm]; //
	TH1D *fFlowGFMixedFinalHist[fkFlowGFNHarm][fkFlowGFNHarm]; //
	
	TProfile *fFlowGFIntCorProPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde]; //
	TH1D *fFlowGFIntCorHistPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde]; //
	TProfile *fFlowGFIntCovProPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
	TH1D *fFlowGFIntCovHistPtB[fkGFPtB][fkFlowGFNHarm][fkFlowGFNOrde][fkFlowGFNOrde]; //
	//TH1D *fFlowGFIntFinalHistPtB[fkFlowGFNHarm][fkFlowGFNOrde];
	//here is the v1pt Histogram 
	TH1D *v1pt[fkFlowGFNHarm][fkFlowGFNOrde];






	// Flow QC part
	TList *fFlowQCList;
	const static Int_t fFlowNHarm = 6;
	const static Int_t fFlowNHarmMax = 14; // WARNING: MIN (2*fFlowNHarm+2)
	const static Int_t fQVecPower = 5;
	Int_t fPtDiffNBins; //
	//Double_t PtBins[37];
//for v1(eta)	
	Int_t fEtaDiffNBins;
	TH1D *fPOIPtDiffQRe[fQVecPower][fFlowNHarmMax]; // real part
	TH1D *fPOIPtDiffQIm[fQVecPower][fFlowNHarmMax]; // imaginary part
	TH1D *fPOIPtDiffMul[fQVecPower][fFlowNHarmMax]; // imaginary part

	const static Int_t fkFlowQCnIntCorPro = 5;
	TProfile *fFlowQCIntCorPro[fFlowNHarm][fkFlowQCnIntCorPro]; //
	TH1D *fFlowQCIntCorHist[fFlowNHarm][fkFlowQCnIntCorPro]; //
	TH1D *fFlowQCIntCumHist[fFlowNHarm][fkFlowQCnIntCorPro];

	const static Int_t fFlowQCNPro = 4;
	const static Int_t fCRCMaxnCen = 10;
	const static Int_t fFlowQCNCov = 8;
	TProfile *fFlowQCCorPro[fCRCMaxnCen][fFlowNHarm][fFlowQCNPro];
	TProfile *fFlowQCCorCovPro[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov];
	TH1D *fFlowQCCorHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNPro]; // <<2'>>, [CRCBin][eg]
	TH1D *fFlowQCCorCovHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov]; // histo for covariances
	TH1D *fFlowQCFinalPtDifHist[fCRCMaxnCen][fFlowNHarm][fFlowQCNCov]; //

	Int_t fFlowQCCenBin;
	
	const static Int_t fFlowQCNRef = 14;
	TProfile *fFlowQCRefCorPro[fFlowNHarm][fFlowQCNRef]; //
	TH1D *fFlowQCRefCorHist[fFlowNHarm][fFlowQCNRef]; //
	TH1D *fFlowQCRefCorFinal[fFlowNHarm][4]; //
	
	// Flow from BW
	TList *fFlowFromBWList;
	const static Int_t fFlowFromBWCenBin = 10;
	Double_t ReQn[fFlowFromBWCenBin] = {0.}, ImQn[fFlowFromBWCenBin] = {0.};
        //for integrated deltav1
        Double_t PosReQn[fFlowFromBWCenBin] = {0.}, PosImQn[fFlowFromBWCenBin] = {0.};
	Double_t NegReQn[fFlowFromBWCenBin] = {0.}, NegImQn[fFlowFromBWCenBin] = {0.};   
    const static Int_t nHar = 4;//was 3, but to inclde v1 and to still account for the others - I think this messes up things ... - still doesn't matter for v1
	TProfile *V2IntPro[nHar];
	TProfile *V2IntDelta[nHar];
    TProfile *V2IntProQC[nHar];
    TProfile *V2PtDiffPro[nHar][10];
	TProfile *V2PtDiffPos[nHar][10];
	TProfile *V2PtDiffNeg[nHar][10];
	TProfile *V2PtDiffDelta[nHar][10];
//and for eta
	TProfile *V2EtaDiffPro[nHar][10];
        TProfile *V2EtaDiffPos[nHar][10];       
        TProfile *V2EtaDiffNeg[nHar][10];                          
        TProfile *V2EtaDiffDelta[nHar][10];
    TH1F *hRepn[nHar], *hImpn[nHar];//declarations for v1 caluclation(s)
    TH1F *hReEtan[nHar], *hImEtan[nHar];
    TH1F *hPosRepn[nHar], *hNegRepn[nHar];
    TH1F *hPosReEtan[nHar], *hNegReEtan[nHar];
    TH1F *hPosImpn[nHar], *hNegImpn[nHar];
    TH1F *hPosImEtan[nHar], *hNegImEtan[nHar];
//here are ins for pT an eta stored in MakeEvent---old
	TH1F *hMpn;
	TH1F *hPosMpn;
	TH1F *hNegMpn;

	TH1F *hMen;
	TH1F *hPosMen;
	TH1F *hNegMen;
	//v1 differential (vs pT/Eta) hists - maybe add the Integrated one for \Delta v1(Cent)
 	TH1F *Posv1pt;
	TH1F *Negv1pt;
	TH1F *Deltav1pt;

	TH1F *V1eta;
	TH1F *Posv1eta;
	TH1F *Negv1eta;
	TH1F *Deltav1eta;    

	// CMW part
	TList *fCMWList;
	Double_t fEtaGapNeg = -0.1;
	Double_t fEtaGapPos = 0.1;
	Int_t gPsiN = 2;
	Double_t fAchrgNet = -1;
	
	TProfile *fHistAChrgVsCent;
	TProfile *fHistEPResolution;
	TProfile *fHistEPResolutionAch[9];
	TProfile *fHistv2AchChrgPosEtaNeg[2][9];
	TProfile *fHistv2AchChrgPosEtaPos[2][9];
	TProfile *fHist2AchChrgNegEtaNeg[2][9];
	TProfile *fHistv2AchChrgNegEtaPos[2][9];
	TProfile *fHistv2AchChrgPosChrgNeg[2][9];
	TProfile *fHistv2AchChrgNegChrgPos[2][9];
	TProfile *fHistv2cumAchChrgAll[9];  // Charge inclusive
	
	TProfile *fHistAveV2PosAch;
    TProfile *fHistAveV2Pos;
          
	TProfile *fHistAveV2NegAch;
    TProfile *fHistAveV2Neg;
    
    TProfile *fv2plusminus;	
	TH1D *fCMWThreeParticleCorrelator[2];
	  
	Double_t fSumTPCQn2xEtaNegWhole = 0;
	Double_t fSumTPCQn2yEtaNegWhole = 0;
	Double_t fSumWgtEtaNegWhole = 0;
	Double_t fSumTPCQn2xEtaPosWhole = 0;
	Double_t fSumTPCQn2yEtaPosWhole = 0;
	Double_t fSumWgtEtaPosWhole = 0;
	Double_t fSumTPCQn2xEtaNegChPos = 0;
	Double_t fSumTPCQn2yEtaNegChPos = 0;
	Double_t fSumWgtEtaNegChPos = 0;
	Double_t fSumTPCQn2xEtaNegChNeg = 0;
	Double_t fSumTPCQn2yEtaNegChNeg = 0;
	Double_t fSumWgtEtaNegChNeg = 0;
	Double_t fSumTPCQn2xEtaPosChPos = 0;
	Double_t fSumTPCQn2yEtaPosChPos = 0;
	Double_t fSumWgtEtaPosChPos = 0;
	Double_t fSumTPCQn2xEtaPosChNeg = 0;
	Double_t fSumTPCQn2yEtaPosChNeg = 0;
	Double_t fSumWgtEtaPosChNeg = 0;
	
	Double_t uqRe=0, uqIm=0;
	Double_t sumQxTPC=0, sumQyTPC=0, sumWgtTPC=0;
	Double_t sumQxTPCEtaNegChPos=0, sumQyTPCEtaNegChPos=0, sumQxTPCEtaPosChPos=0, sumQyTPCEtaPosChPos=0, sumWgtEtaNegChPos=0, sumWgtEtaPosChPos=0;
	Double_t sumQ2xChrgPosEtaPos=0, sumQ2yChrgPosEtaPos=0, NumOfChrgPosEtaPos=0, sumQ2xChrgPosEtaNeg=0, sumQ2yChrgPosEtaNeg=0, NumOfChrgPosEtaNeg=0, sumQ2xChrgNegEtaPos=0, sumQ2yChrgNegEtaPos=0, NumOfChrgNegEtaPos=0, sumQ2xChrgNegEtaNeg=0, sumQ2yChrgNegEtaNeg=0, NumOfChrgNegEtaNeg=0;
	
	//CMW QA hist
	TList *fCMWQAList;
	TH1D *fEvPlTPC;
	
	//CME part
	TList *fCMEList;
	Double_t ReQ1P=0, ReQ1N=0, ImQ1P=0, ImQ1N=0;
    Double_t ReQ2P=0, ImQ2P=0, ReQ2N=0, ImQ2N=0;
    Double_t ReQ3P=0, ImQ3P=0, ReQ3N=0, ImQ3N=0;
    Double_t ReQ4P=0, ImQ4P=0, ReQ4N=0, ImQ4N=0;
    Double_t MQP=0, MQN=0, MQ=0, PMQ = 0, NMQ = 0;

	//--------- individual Q-vector terms -----------
	const static Int_t nQVector = 8;
	TProfile *fReQVectorPro[nQVector];
	TProfile *fImQVectorPro[nQVector];
	
    //------- 2particle correlator -------------
    const static Int_t nDnncor = 2;
    TProfile *fDnnPro[nDnncor];
    TH1F     *fDnnHist[nDnncor];
  
	//--------- 3 particle correlator -----------
    const static Int_t nCMEcor = 8;
    TProfile *fCMEPro[nCMEcor];
    TH1F     *fCMEHist[nCMEcor];
    
    //--------- final hists --------------------
    TH1F* DeltaD11;
    TH1F* DeltaG112;
    TH1F* DeltaG132;
    TH1F* DeltaG123;
    TH1F* DeltaG224;
    
	// Cuts:
	Double_t maxPtCut = 9999; // the data is in GeV? Double check.
	Double_t minPtCut = 0;
	Double_t maxEtaCut = 99; //|eta|<maxEtaCut
	Int_t minNtracks = 0;
};
#endif



//////////////////////////////////////////


	//TH1D *fTPCQn2xWhole;
	//TH1D *fTPCQn2yWhole;
	//TH1D *fWgtWhole;
	//TH1D *fWgt2Whole;
	//TH1D *fTPCQn2xEtaNegWhole;
	//TH1D *fTPCQn2yEtaNegWhole;
	//TH1D *fWgtEtaNegWhole;
	//TH1D *fWgt2EtaNegWhole;
	//TH1D *fTPCQn2xEtaPosWhole;
	//TH1D *fTPCQn2yEtaPosWhole;
	//TH1D *fWgtEtaPosWhole;
	//TH1D *fWgt2EtaPosWhole;
	//TH1D *fTPCQn2xEtaNegChPos;
	//TH1D *fTPCQn2yEtaNegChPos;
	//TH1D *fWgtEtaNegChPos;
	//TH1D *fWgt2EtaNegChPos;
	//TH1D *fTPCQn2xEtaNegChNeg;
	//TH1D *fTPCQn2yEtaNegChNeg;
	//TH1D *fWgtEtaNegChNeg;
	//TH1D *fWgt2EtaNegChNeg;
	//TH1D *fTPCQn2xEtaPosChPos;
	//TH1D *fTPCQn2yEtaPosChPos;
	//TH1D *fWgtEtaPosChPos;
	//TH1D *fWgt2EtaPosChPos;
	//TH1D *fTPCQn2xEtaPosChNeg;
	//TH1D *fTPCQn2yEtaPosChNeg;
	//TH1D *fWgtEtaPosChNeg;
	//TH1D *fWgt2EtaPosChNeg;
