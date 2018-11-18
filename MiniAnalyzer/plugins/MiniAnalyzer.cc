// -*- C++ -*-
//
// Package:    PythiaAnalysis/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc PythiaAnalysis/MiniAnalyzer/plugins/MiniAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Meng Xiao
//         Created:  Tue, 24 Jul 2018 20:00:32 GMT
//
//


// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>

#include "DataFormats/JetReco/interface/GenJet.h"
#include <MelaAnalytics/GenericMEComputer/interface/GMECHelperFunctions.h>


#include "PythiaAnalysis/MiniAnalyzer/interface/EwkCorrections.h"
#include "PythiaAnalysis/MiniAnalyzer/interface/LHEHandler.h"
#include "PythiaAnalysis/MiniAnalyzer/src/kFactors.C"
#include "PythiaAnalysis/MiniAnalyzer/interface/DaughterDataHelpers.h"
#include "PythiaAnalysis/MiniAnalyzer/interface/HZZ4lNtupleFactory.h"
#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>
#include <MelaAnalytics/GenericMEComputer/interface/GMECHelperFunctions.h>

#include <PythiaAnalysis/MiniAnalyzer/interface/MCHistoryTools.h>
#include <PythiaAnalysis/MiniAnalyzer/interface/FinalStates.h>



using namespace BranchHelpers;

enum FailedTreeLevel {
	//only used if skipEmptyEvents == false
	noFailedTree = 0,      //events with no candidate are skipped entirely
	minimalFailedTree = 1, //events with no candidate are written to a separate tree with minimal information: Gen Higgs momentum, lepton flavors, and LHE level ME's
	LHEFailedTree = 2,     //same as (1) + full LHE level information
	fullFailedTree = 3,    //same as (2) + full pythia level information + reco information that is available without a candidate, e.g. jets
};
namespace {

	//List of variables with default values
	Int_t RunNumber  = 0;
	Long64_t EventNumber  = 0;
	Int_t LumiNumber  = 0;
	Short_t NRecoMu  = 0;
	Short_t NRecoEle  = 0;
	Short_t Nvtx  = 0;
	Short_t NObsInt  = 0;
	Float_t NTrueInt  = 0;
	Float_t PUWeight  = 0;

	Float_t KFactor_QCD_ggZZ_Nominal = 0;
	Float_t KFactor_EW_qqZZ = 0;
	Float_t KFactor_QCD_qqZZ_dPhi = 0;
	Float_t KFactor_QCD_qqZZ_M = 0;
	Float_t KFactor_QCD_qqZZ_Pt = 0;

	Float_t PFMET  =  -99;
	Short_t nCleanedJets  =  0;
	Short_t nCleanedJetsPt30  = 0;
	Short_t trigWord  = 0;
	Float_t ZZMass  = 0;
	Short_t ZZsel  = 0;
	Float_t ZZPt  = 0;
	Float_t ZZEta  = 0;
	Float_t ZZPhi  = 0;
	Int_t CRflag  = 0;
	Float_t Z1Mass  = 0;
	Float_t Z1Pt  = 0;
	Short_t Z1Flav  = 0;
	Float_t Z2Mass  = 0;
	Float_t Z2Pt  = 0;
	Short_t Z2Flav  = 0;
	Float_t costhetastar  = 0;
	Float_t helphi  = 0;
	Float_t helcosthetaZ1  = 0;
	Float_t helcosthetaZ2  = 0;
	Float_t phistarZ1  = 0;
	Float_t phistarZ2  = 0;
	Float_t xi  = 0;
	Float_t xistar  = 0;
	Short_t evtPassMETTrigger = 0;

	std::vector<float> LepPt;
	std::vector<float> LepEta;
	std::vector<float> LepPhi;
	std::vector<short> LepLepId;



	std::vector<float> JetPt ;
	std::vector<float> JetEta ;
	std::vector<float> JetPhi ;
	std::vector<float> JetMass ;


	Float_t DiJetMass  = -99;
	Float_t DiJetDEta  = -99;
	Short_t nExtraLep  = 0;
	Short_t nExtraZ  = 0;
	std::vector<float> ExtraLepPt;
	std::vector<float> ExtraLepEta;
	std::vector<float> ExtraLepPhi ;
	std::vector<short> ExtraLepLepId;
	Short_t genFinalState  = 0;
	Int_t genProcessId  = 0;
	Float_t genHEPMCweight  = 0;
	Float_t genHEPMCweight_NNLO  = 0;
	Float_t genHEPMCweight_POWHEGonly = 0;


	std::vector<float> LHEMotherPz;
	std::vector<float> LHEMotherE;
	std::vector<short> LHEMotherId;
	std::vector<float> LHEDaughterPt;
	std::vector<float> LHEDaughterEta;
	std::vector<float> LHEDaughterPhi;
	std::vector<float> LHEDaughterMass;
	std::vector<short> LHEDaughterId;
	std::vector<float> LHEAssociatedParticlePt;
	std::vector<float> LHEAssociatedParticleEta;
	std::vector<float> LHEAssociatedParticlePhi;
	std::vector<float> LHEAssociatedParticleMass;
	std::vector<short> LHEAssociatedParticleId;

	Float_t LHEPDFScale = 0;


	Short_t genExtInfo  = 0;
	Float_t xsection  = 0;
	Float_t genxsection = 0;
	Float_t genbranchingratio = 0;
	Float_t dataMCWeight  = 0;
	Float_t trigEffWeight  = 0;
	Float_t HqTMCweight  = 0;
	Float_t ZXFakeweight  = 0;
	Float_t overallEventWeight  = 0;
	Float_t GenHMass  = 0;
	Float_t GenHPt  = 0;
	Float_t GenHRapidity  = 0;
	Float_t GenZ1Mass  = 0;
	Float_t GenZ1Eta  = 0;
	Float_t GenZ1Pt  = 0;
	Float_t GenZ1Phi  = 0;
	Float_t GenZ1Flav  = 0;
	Float_t GenZ2Mass  = 0;
	Float_t GenZ2Eta  = 0;
	Float_t GenZ2Pt  = 0;
	Float_t GenZ2Phi  = 0;
	Float_t GenZ2Flav  = 0;
	Float_t GenLep1Pt  = 0;
	Float_t GenLep1Eta  = 0;
	Float_t GenLep1Phi  = 0;
	Short_t GenLep1Id  = 0;
	Float_t GenLep2Pt  = 0;
	Float_t GenLep2Eta  = 0;
	Float_t GenLep2Phi  = 0;
	Short_t GenLep2Id  = 0;
	Float_t GenLep3Pt  = 0;
	Float_t GenLep3Eta  = 0;
	Float_t GenLep3Phi  = 0;
	Short_t GenLep3Id  = 0;
	Float_t GenLep4Pt  = 0;
	Float_t GenLep4Eta  = 0;
	Float_t GenLep4Phi  = 0;
	Short_t GenLep4Id  = 0;
	Float_t GenAssocLep1Pt  = 0;
	Float_t GenAssocLep1Eta  = 0;
	Float_t GenAssocLep1Phi  = 0;
	Short_t GenAssocLep1Id  = 0;
	Float_t GenAssocLep2Pt  = 0;
	Float_t GenAssocLep2Eta  = 0;
	Float_t GenAssocLep2Phi  = 0;
	Short_t GenAssocLep2Id  = 0;

	//Int_t   htxsNJets = 0;
	//Float_t htxsHPt = 0;
	//Float_t ggH_NNLOPS_weight = 0;
	std::vector<float> qcd_ggF_uncertSF;


	//FIXME: temporary fix to the mismatch of charge() and sign(pdgId()) for muons with BTT=4
	int getPdgId(const reco::Candidate* p) {
		int id = p->pdgId();
		if (id!=22 && //for TLEs
				signbit(id) && p->charge()<0) id*=-1; // negative pdgId must be positive charge
		return id;
	}

}
void set_bit( int& mask, unsigned int iBit ) { mask |= (1<<iBit); }



// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit MiniAnalyzer(const edm::ParameterSet&);
		~MiniAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		static void addweight(float &weight, float weighttoadd);

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
		void FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
				const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
		void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
				const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
		void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);
		virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
		void BookAllBranches();
		float EvalSpline(TSpline3* const& sp, float xval);
		void FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps);
		FailedTreeLevel failedTreeLevel;
		//counters
		Float_t Nevt_Gen;
		Float_t Nevt_Gen_lumiBlock;

		Float_t gen_ZZ4mu;
		Float_t gen_ZZ4e;
		Float_t gen_ZZ2mu2e;
		Float_t gen_ZZ2l2tau;
		Float_t gen_ZZ2emu2tau;
		Float_t gen_ZZ4tau;
		Float_t gen_ZZ4mu_EtaAcceptance;
		Float_t gen_ZZ4mu_LeptonAcceptance;
		Float_t gen_ZZ4e_EtaAcceptance;
		Float_t gen_ZZ4e_LeptonAcceptance;
		Float_t gen_ZZ2mu2e_EtaAcceptance;
		Float_t gen_ZZ2mu2e_LeptonAcceptance;
		Float_t gen_BUGGY;
		Float_t gen_Unknown;

		Float_t gen_sumPUWeight;
		Float_t gen_sumGenMCWeight;
		Float_t gen_sumWeights;

		edm::EDGetTokenT<edm::View<reco::Candidate> > prunedGenToken_;
		edm::Handle<edm::View<reco::Candidate>> pruned;
		edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
		edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;

		LHEHandler* lheHandler;

		int apply_K_NNLOQCD_ZZGG; // 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
		bool apply_K_NNLOQCD_ZZQQB;
		bool apply_K_NLOEW_ZZQQB;

		void buildMELABranches();
		void computeMELABranches(MELACandidate* cand);
		void updateMELAClusters_Common(const string clustertype);
		void updateMELAClusters_NoInitialQ(const string clustertype);
		void updateMELAClusters_NoInitialG(const string clustertype);
		void updateMELAClusters_NoAssociatedG(const string clustertype);
		void updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype);
		void updateMELAClusters_BestLOAssociatedZ(const string clustertype);
		void updateMELAClusters_BestLOAssociatedW(const string clustertype);
		void updateMELAClusters_BestLOAssociatedVBF(const string clustertype);
		void updateMELAClusters_BestNLOVHApproximation(const string clustertype);
		void updateMELAClusters_BestNLOVBFApproximation(const string clustertype);
		void pushRecoMELABranches(const pat::CompositeCandidate& cand);
		void pushLHEMELABranches();
		void clearMELABranches();
		virtual void FillLHECandidate();
		Float_t getHqTWeight(double mH, double genPt) const;


		virtual void FillJet(const reco::GenJet& jet);


		Mela mela;
		std::vector<std::string> recoMElist;
		std::vector<MELAOptionParser*> recome_originalopts;
		std::vector<MELAOptionParser*> recome_copyopts;
		std::vector<std::string> lheMElist;
		//std::vector<MELAOptionParser*> lheme_originalopts;
		std::vector<MELAOptionParser*> lheme_copyopts;
		std::vector<MELAHypothesis*> lheme_units;
		std::vector<MELAHypothesis*> lheme_aliased_units;
		std::vector<MELAComputation*> lheme_computers;
		std::vector<MELACluster*> lheme_clusters;
		string sampleName;

		std::vector<std::vector<float> > ewkTable;

		TSpline3* spkfactor_ggzz_nnlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp
		TSpline3* spkfactor_ggzz_nlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp

		HZZ4lNtupleFactory *myTree;
		bool printedLHEweightwarning;

		Float_t xsec;
		Float_t genxsec;
		Float_t genbr;
		std::string theCandLabel;
		TH1F *hCounter;
		edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
		TString theFileName;
		  int year;

		bool skipEmptyEvents;




		// ----------member data ---------------------------
};

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
	//prunedGenToken_(consumes<edm::View<reco::GenParticle> >( edm::InputTag("prunedGenParticles"))),
	failedTreeLevel(FailedTreeLevel(iConfig.getParameter<int>("failedTreeLevel"))),
	prunedGenToken_(consumes<edm::View<reco::Candidate> >( iConfig.getParameter<edm::InputTag>("pruned"))),
	genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjet"))),
	apply_K_NNLOQCD_ZZGG(iConfig.getParameter<int>("Apply_K_NNLOQCD_ZZGG")),
	apply_K_NNLOQCD_ZZQQB(iConfig.getParameter<bool>("Apply_K_NNLOQCD_ZZQQB")),
	apply_K_NLOEW_ZZQQB(iConfig.getParameter<bool>("Apply_K_NLOEW_ZZQQB")),
	mela(13., 125., TVar::ERROR),
	recoMElist(iConfig.getParameter<std::vector<std::string>>("recoProbabilities")),
	lheMElist(iConfig.getParameter<std::vector<std::string>>("lheProbabilities")),
	sampleName(iConfig.getParameter<string>("sampleName")),
	myTree(nullptr),
	xsec(iConfig.getParameter<double>("xsec")),
	genxsec(iConfig.getParameter<double>("GenXSEC")),
	genbr(iConfig.getParameter<double>("GenBR")),
	theCandLabel(iConfig.getUntrackedParameter<string>("CandCollection")), // Name of input ZZ collection

	theFileName(iConfig.getUntrackedParameter<string>("fileName")),
	 year(iConfig.getParameter<int>("setup")),

	skipEmptyEvents(iConfig.getParameter<bool>("skipEmptyEvents"))




{
	//now do what ever initialization is needed

	lheHandler = new LHEHandler(iConfig.getParameter<int>("VVMode"), iConfig.getParameter<int>("VVDecayMode"), true, year);
	genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
	consumesMany<LHEEventProduct>();
	candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));

	Nevt_Gen = 0;
	Nevt_Gen_lumiBlock = 0;

	//For Efficiency studies
	gen_ZZ4mu = 0;
	gen_ZZ4e = 0;
	gen_ZZ2mu2e = 0;
	gen_ZZ2l2tau = 0;
	gen_ZZ2emu2tau = 0;
	gen_ZZ4tau = 0;
	gen_ZZ4mu_EtaAcceptance = 0;
	gen_ZZ4mu_LeptonAcceptance = 0;
	gen_ZZ4e_EtaAcceptance = 0;
	gen_ZZ4e_LeptonAcceptance = 0;
	gen_ZZ2mu2e_EtaAcceptance = 0;
	gen_ZZ2mu2e_LeptonAcceptance = 0;
	gen_BUGGY = 0;
	gen_Unknown = 0;

	gen_sumPUWeight = 0.f;
	gen_sumGenMCWeight = 0.f;
	gen_sumWeights =0.f;

	std::string fipPath;

	// Read EWK K-factor table from file
	edm::FileInPath ewkFIP("PythiaAnalysis/MiniAnalyzer/data/kfactors/ZZ_EwkCorrections.dat");
	fipPath=ewkFIP.fullPath();
	ewkTable = EwkCorrections::readFile_and_loadEwkTable(fipPath.data());

	// Read the ggZZ k-factor shape from file
	TString strZZGGKFVar[9]={
	  "Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"
	};
	edm::FileInPath ggzzFIP_NNLO("PythiaAnalysis/MiniAnalyzer/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
	fipPath=ggzzFIP_NNLO.fullPath();
	TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
	for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nnlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NNLO", strZZGGKFVar[ikf].Data()));
	ggZZKFactorFile->Close();
	edm::FileInPath ggzzFIP_NLO("PythiaAnalysis/MiniAnalyzer/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root");
	fipPath=ggzzFIP_NLO.fullPath();
	ggZZKFactorFile = TFile::Open(fipPath.data());
	for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NLO", strZZGGKFVar[ikf].Data()));
	ggZZKFactorFile->Close();

	usesResource("TFileService");

}


MiniAnalyzer::~MiniAnalyzer()
{

	clearMELABranches(); // Cleans LHE branches
	delete lheHandler;

}
//
// constants, enums and typedefs
//

//
// static data member definitions
void MiniAnalyzer::computeMELABranches(MELACandidate* cand){
	mela.setCurrentCandidate(cand);
	// Sequantial computation
	updateMELAClusters_Common("Common");
	updateMELAClusters_NoInitialQ("NoInitialQ");
	updateMELAClusters_NoInitialG("NoInitialG");
	updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZ");
	updateMELAClusters_BestLOAssociatedW("BestLOAssociatedW");
	updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBF");
	updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximation");
	updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximation");
	updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximation");
	updateMELAClusters_NoAssociatedG("NoAssociatedG");
	updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedG");
	// Reverse sequence
	updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedGLast");
	updateMELAClusters_NoAssociatedG("NoAssociatedGLast");
	updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximationLast");
	updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximationLast");
	updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximationLast");
	updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBFLast");
	updateMELAClusters_BestLOAssociatedW("BestLOAssociatedWLast");
	updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZLast");
	updateMELAClusters_NoInitialG("NoInitialGLast");
	updateMELAClusters_NoInitialQ("NoInitialQLast");
	updateMELAClusters_Common("CommonLast");
	// Reset mela
	mela.resetInputEvent();
}
// Common ME computations that do not manipulate the LHE candidate
void MiniAnalyzer::updateMELAClusters_Common(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}
}
// ME computations that require no quark initial state
void MiniAnalyzer::updateMELAClusters_NoInitialQ(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of quark mothers
	std::vector<int> motherIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAQuark(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluon initial state
void MiniAnalyzer::updateMELAClusters_NoInitialG(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> motherIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluons as associated particles
void MiniAnalyzer::updateMELAClusters_NoAssociatedG(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> ajetIds;
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
		ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
		if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require no gluon initial state and no gluons as associated particles
void MiniAnalyzer::updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> motherIds;
	std::vector<int> ajetIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
		ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
		if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require best Z, W or VBF topology at LO (no gluons)
void MiniAnalyzer::updateMELAClusters_BestLOAssociatedZ(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> motherIds;
	std::vector<int> ajetIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
		ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
		if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
	}
	// Give precedence to leptonic V decays
	bool hasALepV=false;
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			if (
					PDGHelpers::isALepton(Vtmp->getDaughter(0)->id)
					||
					PDGHelpers::isANeutrino(Vtmp->getDaughter(0)->id)
			   ){
				hasALepV=true;
			}
		}
	}
	int bestVbyMass=-1;
	float bestVMassDiff=1e5;
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			if (
					PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
					&& hasALepV
			   ){
				for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(false);
			}
			else if (fabs(Vtmp->m()-PDGHelpers::Zmass)<bestVMassDiff){
				bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Zmass);
				bestVbyMass = iv;
			}
		}
	}
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected((iv==bestVbyMass));
		}
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(true);
		}
	}
}
void MiniAnalyzer::updateMELAClusters_BestLOAssociatedW(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> motherIds;
	std::vector<int> ajetIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
		ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
		if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
	}
	// Give precedence to leptonic V decays
	bool hasALepV=false;
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			if (
					PDGHelpers::isALepton(Vtmp->getDaughter(0)->id)
					||
					PDGHelpers::isANeutrino(Vtmp->getDaughter(0)->id)
			   ){
				hasALepV=true;
			}
		}
	}
	int bestVbyMass=-1;
	float bestVMassDiff=1e5;
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			if (
					PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
					&& hasALepV
			   ){
				for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(false);
			}
			else if (fabs(Vtmp->m()-PDGHelpers::Wmass)<bestVMassDiff){
				bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Wmass);
				bestVbyMass = iv;
			}
		}
	}
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected((iv==bestVbyMass));
		}
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
	for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
		MELAParticle* Vtmp = melaCand->getSortedV(iv);
		if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
			for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(true);
		}
	}
}
void MiniAnalyzer::updateMELAClusters_BestLOAssociatedVBF(const string clustertype){
	// Same as updateMELAClusters_NoInitialGNoAssociatedG, but keep a separate function for future studies
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Manipulate the candidate
	// Assign 0 to the id of gluon mothers
	std::vector<int> motherIds;
	std::vector<int> ajetIds;
	for (int imot=0; imot<melaCand->getNMothers(); imot++){
		motherIds.push_back(melaCand->getMother(imot)->id);
		if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
	}
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
		ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
		if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
	}

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	// Restore the candidate properties
	for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
	for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that can approximate the NLO QCD (-/+ MiNLO extra jet) phase space to LO QCD in signal VBF or VH
// Use these for POWHEG samples
// MELACandidateRecaster has very specific use cases, so do not use these functions for other cases.
void MiniAnalyzer::updateMELAClusters_BestNLOVHApproximation(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Check if any clusters request this computation
	bool clustersRequest=false;
	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			clustersRequest=true;
			break;
		}
	}
	if (!clustersRequest) return;

	// Need one recaster for each of ZH and WH, so distinguish by the cluster name
	TVar::Production candScheme;
	if (clustertype.find("BestNLOZHApproximation")!=string::npos) candScheme = TVar::Had_ZH;
	else if (clustertype.find("BestNLOWHApproximation")!=string::npos) candScheme = TVar::Had_WH;
	else return;

	MELACandidateRecaster recaster(candScheme);
	MELACandidate* candModified=nullptr;
	MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(melaCand, candScheme);
	if (bestAV){
		recaster.copyCandidate(melaCand, candModified);
		recaster.deduceLOVHTopology(candModified);
		mela.setCurrentCandidate(candModified);
	}
	else return; // No associated Vs found. The algorithm won't work.

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	delete candModified;
	mela.setCurrentCandidate(melaCand); // Go back to the original candidate
}
void MiniAnalyzer::updateMELAClusters_BestNLOVBFApproximation(const string clustertype){
	MELACandidate* melaCand = mela.getCurrentCandidate();
	if (melaCand==0) return;

	// Check if any clusters request this computation
	bool clustersRequest=false;
	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			clustersRequest=true;
			break;
		}
	}
	if (!clustersRequest) return;

	// Need one recaster for VBF
	TVar::Production candScheme;
	if (clustertype.find("BestNLOVBFApproximation")!=string::npos) candScheme = TVar::JJVBF;
	else return;

	MELACandidateRecaster recaster(candScheme);
	MELACandidate* candModified=nullptr;
	recaster.copyCandidate(melaCand, candModified);
	recaster.reduceJJtoQuarks(candModified);
	mela.setCurrentCandidate(candModified);

	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
		MELACluster* theCluster = lheme_clusters.at(ic);
		if (theCluster->getName()==clustertype){
			// Re-compute all related hypotheses first...
			theCluster->computeAll();
			// ...then update the cluster
			theCluster->update();
		}
	}

	delete candModified;
	mela.setCurrentCandidate(melaCand); // Go back to the original candidate
}


void MiniAnalyzer::pushRecoMELABranches(const pat::CompositeCandidate& cand){
	std::vector<MELABranch*>* recome_branches = myTree->getRecoMELABranches();
	// Pull + push...
	for (unsigned int ib=0; ib<recome_branches->size(); ib++){
		std::string branchname = recome_branches->at(ib)->bname.Data();
//		std::cout <<"pushing branch "<< branchname<< "\t"<<cand.userFloat(branchname)<<std::endl;
		if (cand.hasUserFloat(branchname)) recome_branches->at(ib)->setValue((Float_t)cand.userFloat(branchname));
		else cerr << "MiniAnalyzer::pushRecoMELABranches: Candidate does not contain the reco ME " << branchname << " it should have calculated!" << endl;
	}
}
void MiniAnalyzer::pushLHEMELABranches(){
	std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
	// Pull + push...
	for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->setVal();
	// ...then reset
	for (unsigned int ic=0; ic<lheme_clusters.size(); ic++) lheme_clusters.at(ic)->reset();
}
void MiniAnalyzer::clearMELABranches(){
	for (unsigned int it=0; it<lheme_clusters.size(); it++) delete lheme_clusters.at(it);
	for (unsigned int it=0; it<lheme_computers.size(); it++) delete lheme_computers.at(it);
	for (unsigned int it=0; it<lheme_copyopts.size(); it++) delete lheme_copyopts.at(it);
	//for (unsigned int it=0; it<lheme_aliased_units.size(); it++) delete lheme_aliased_units.at(it); // DO NOT DELETE THIS, WILL BE DELETED WITH lheme_units!
	for (unsigned int it=0; it<lheme_units.size(); it++) delete lheme_units.at(it);

	for (unsigned int it=0; it<recome_copyopts.size(); it++) delete recome_copyopts.at(it);
	for (unsigned int it=0; it<recome_originalopts.size(); it++) delete recome_originalopts.at(it);
}
//

//
// constructors and destructor


void MiniAnalyzer::FillLHECandidate(){
	LHEMotherPz.clear();
	LHEMotherE.clear();
	LHEMotherId.clear();
	LHEDaughterPt.clear();
	LHEDaughterEta.clear();
	LHEDaughterPhi.clear();
	LHEDaughterMass.clear();
	LHEDaughterId.clear();
	LHEAssociatedParticlePt.clear();
	LHEAssociatedParticleEta.clear();
	LHEAssociatedParticlePhi.clear();
	LHEAssociatedParticleMass.clear();
	LHEAssociatedParticleId.clear();

	LHEPDFScale = 0;
	genHEPMCweight_POWHEGonly = 0;

	MELACandidate* cand = lheHandler->getBestCandidate();
	if (cand!=0){
		for (int imot=0; imot<cand->getNMothers(); imot++){
			MELAParticle* apart = cand->getMother(imot);
			if (apart==0){ LHEMotherPz.clear(); LHEMotherE.clear(); LHEMotherId.clear(); break; } // Something went wrong
			LHEMotherPz.push_back(apart->z());
			LHEMotherE.push_back(apart->t());
			LHEMotherId.push_back((short)apart->id);
		}

		for (int iV=0; iV<min(2, cand->getNSortedVs()); iV++){
			MELAParticle* Vi = cand->getSortedV(iV);
			if (Vi!=0){
				for (int iVj=0; iVj<Vi->getNDaughters(); iVj++){
					MELAParticle* Vij = Vi->getDaughter(iVj);
					if (Vij!=0){
						LHEDaughterPt.push_back(Vij->pt());
						if (abs(Vij->pt() / Vij->z()) > 2e-8) {
							LHEDaughterEta.push_back(Vij->eta());
						} else {
							edm::LogWarning("ZeroPt") << "pt = 0!  Using eta = +/-1e10\n" << Vij->id << " " << Vij->x() << " " << Vij->y() << " " << Vij->z() << " " << Vij->t();
							LHEDaughterEta.push_back(copysign(1e10, Vij->z()));
						}
						LHEDaughterPhi.push_back(Vij->phi());
						LHEDaughterMass.push_back(Vij->m());
						LHEDaughterId.push_back((short)Vij->id);
					}
				}
			}
		}

		vector<MELAParticle*> AssociatedParticle;
		vector<MELAParticle*> tmpAssociatedParticle;
		for (int aa=0; aa<cand->getNAssociatedJets(); aa++){
			MELAParticle* apart = cand->getAssociatedJet(aa);
			tmpAssociatedParticle.push_back(apart);
		}
		for (int aa=0; aa<cand->getNAssociatedLeptons(); aa++){
			MELAParticle* apart = cand->getAssociatedLepton(aa);
			if (!PDGHelpers::isANeutrino(apart->id)) tmpAssociatedParticle.push_back(apart);
		}
		for (int aa=0; aa<cand->getNAssociatedNeutrinos(); aa++){
			MELAParticle* apart = cand->getAssociatedNeutrino(aa);
			tmpAssociatedParticle.push_back(apart);
		}
		while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
			MELAParticle* tmpPart=0;
			int pos=0;
			for (unsigned int el=0; el<tmpAssociatedParticle.size(); el++){
				if (tmpPart==0){ tmpPart = tmpAssociatedParticle.at(el); pos=el; }
				else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){ tmpPart = tmpAssociatedParticle.at(el); pos=el; } // Safer to do in two steps
			}
			AssociatedParticle.push_back(tmpPart);
			tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
		}
		for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
			MELAParticle* apart = AssociatedParticle.at(aa);
			if (apart!=0){
				LHEAssociatedParticlePt.push_back(apart->pt());
				if (abs(apart->pt() / apart->z()) > 2e-8) {
					LHEAssociatedParticleEta.push_back(apart->eta());
				} else {
					edm::LogWarning("ZeroPt") << "pt = 0!  Using eta = +/-1e10\n" << apart->id << " " << apart->x() << " " << apart->y() << " " << apart->z() << " " << apart->t();
					LHEAssociatedParticleEta.push_back(copysign(1e10, apart->z()));
				}
				LHEAssociatedParticlePhi.push_back(apart->phi());
				LHEAssociatedParticleMass.push_back(apart->m());
				LHEAssociatedParticleId.push_back((short)apart->id);
			}
		}

		/*
		   cout << "NEW EVENT:" << endl;
		   cout << "Mothers:" << endl;
		   for (unsigned int ipart=0; ipart<LHEMotherId.size(); ipart++) cout << "\t Mot" << ipart << " (pz, E, id) = " << LHEMotherPz.at(ipart) << " " << LHEMotherE.at(ipart) << " " << LHEMotherId.at(ipart) << endl;
		   cout << "Daughters:" << endl;
		   for (unsigned int ipart=0; ipart<LHEDaughterId.size(); ipart++) cout << "\t Dau" << ipart << " (pt, eta, phi, m, id) = " << LHEDaughterPt.at(ipart) << " " << LHEDaughterEta.at(ipart) << " " << LHEDaughterPhi.at(ipart) << " " << LHEDaughterMass.at(ipart) << " " << LHEDaughterId.at(ipart) << endl;
		   cout << "Associated:" << endl;
		   for (unsigned int ipart=0; ipart<LHEAssociatedParticleId.size(); ipart++) cout << "\t APart" << ipart << " (pt, eta, phi, m, id) = " << LHEAssociatedParticlePt.at(ipart) << " " << LHEAssociatedParticleEta.at(ipart) << " " << LHEAssociatedParticlePhi.at(ipart) << " " << LHEAssociatedParticleMass.at(ipart) << " " << LHEAssociatedParticleId.at(ipart) << endl;
		   cout << endl;
		   */

		computeMELABranches(cand);
		pushLHEMELABranches();
	}

	LHEPDFScale = lheHandler->getPDFScale();
	if (genHEPMCweight==1.) {
		genHEPMCweight_NNLO = genHEPMCweight = lheHandler->getLHEOriginalWeight();
		if (!printedLHEweightwarning && genHEPMCweight!=1) {
			printedLHEweightwarning = true;
			edm::LogWarning("InconsistentWeights") << "Gen weight is 1, LHE weight is " << genHEPMCweight;
		}
	}
		if (year == 2017) {
	genHEPMCweight *= lheHandler->reweightNNLOtoNLO();
		}

	genHEPMCweight_POWHEGonly = lheHandler->getPowhegOriginalWeight();
}



//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	using namespace edm;

	using namespace reco;
	using namespace pat;

	iEvent.getByToken(prunedGenToken_,pruned);

	Handle<reco::GenJetCollection> cleanedJets;
	iEvent.getByToken(genJetToken_,cleanedJets);

	edm::Handle<GenEventInfoProduct> genInfo;
	iEvent.getByToken(genInfoToken, genInfo);

//	std::cout<<"event pruned"<<std::endl;
	SimpleParticleCollection_t Leptons;
	SimpleParticleCollection_t theAssociatedV;
	TLorentzVector Jets;

	const reco::Candidate * genH = 0;
	std::vector<const reco::Candidate *> genZLeps;
	std::vector<const reco::Candidate *> genAssocLeps;

	edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
	iEvent.getByToken(candToken, candHandle);
	if(candHandle.failedToGet()) {
		edm::LogError("") << "ZZ collection not found in non-loose electron flow. This should never happen";
	}
	const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

	myTree->InitializeVariables();

	MCHistoryTools mch(iEvent, sampleName, pruned, genInfo);
	genFinalState = mch.genFinalState();
	genProcessId = mch.getProcessID();
	genHEPMCweight_NNLO = genHEPMCweight = mch.gethepMCweight(); // Overridden by LHEHandler if genHEPMCweight==1.
	// For 2017 MC, genHEPMCweight is reweighted later from NNLO to NLO


	genExtInfo = mch.genAssociatedFS();

	genH = mch.genH();
	genZLeps     = mch.sortedGenZZLeps();
	genAssocLeps = mch.genAssociatedLeps();
	//genFSR       = mch.genFSR();

	if(genH != 0){
		FillHGenInfo(genH->p4(),1.);
	}
	else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
		FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
	}

	if (genFinalState!=BUGGY) {

		if (genZLeps.size()==4) {

			// "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
			FillZGenInfo(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId(),
					genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId(),
					genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
					genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

			// Gen leptons
			FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
					genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4());

		}

		if (genZLeps.size()==3) {
			FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
					genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector));
		}
		if (genZLeps.size()==2) {
			FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
					genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
		}

		addweight(gen_sumPUWeight, PUWeight);
		addweight(gen_sumGenMCWeight, genHEPMCweight);
		addweight(gen_sumWeights, PUWeight*genHEPMCweight);
	}

	addweight(Nevt_Gen_lumiBlock, 1); // Needs to be outside the if-block

	edm::Handle<LHEEventProduct> lhe_evt;
	vector<edm::Handle<LHEEventProduct> > lhe_handles;
	iEvent.getManyByType(lhe_handles);
	if (lhe_handles.size()>0){
		lhe_evt = lhe_handles.front();
		lheHandler->setHandle(&lhe_evt);
		lheHandler->extract();
		FillLHECandidate(); // Also writes weights
		lheHandler->clear();
	}
	bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
	bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

	mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);

	if (genFinalState == EEEE) {
		addweight(gen_ZZ4e, 1);
		if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4e_EtaAcceptance, 1);
		if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4e_LeptonAcceptance, 1);
	} else if (genFinalState == MMMM) {
		addweight(gen_ZZ4mu, 1);
		if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4mu_EtaAcceptance, 1);
		if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4mu_LeptonAcceptance, 1);
	} else if (genFinalState == EEMM) {
		addweight(gen_ZZ2mu2e, 1);
		if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ2mu2e_EtaAcceptance, 1);
		if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ2mu2e_LeptonAcceptance, 1);
	} else if (genFinalState == llTT){
		addweight(gen_ZZ2emu2tau, 1);
		addweight(gen_ZZ2l2tau, 1);
	} else if (genFinalState == TTTT){
		addweight(gen_ZZ4tau, 1);
		addweight(gen_ZZ2l2tau, 1);
	} else if (genFinalState == BUGGY){  //handle MCFM ZZ->4tau mZ<2mtau bug
		addweight(gen_BUGGY, 1);
		return; // BUGGY events are skipped
	} else {
		addweight(gen_Unknown, 1);
	}

	//for( View<Candidate>::const_iterator par = pruned->begin(); par != pruned->end(); ++ par ) {
	//	//		for (unsigned int i=0;i<pruned.size();i++){
	//	//					reco::Candidate par = pruned.at(i);
	//	int pdgid = par->pdgId();
	//	int status= par->status();
	//	//if(status==1 && (abs(pdgid)==11 || abs(pdgid==13) ))
	//	if( (abs(pdgid)==11 || abs(pdgid)==13) && status==1){
	//		int hasmother=(par->mother()!=0);
	//		int motherid =hasmother? par->mother()->pdgId():-1;
	//		if(motherid == 23 || motherid==pdgid || motherid==25){
	//			std::cout<<pdgid<<"\t"<<status<<"\t"<< par->pt()<<"\t"<<motherid<<std::endl;
	//			//Leptons.push_back(&*par);
	//		}
	//	}
	//	else if(  abs(pdgid)==15 && status==2){
	//		int hasmother=(par->mother()!=0);
	//		int motherid =hasmother? par->mother()->pdgId():-1;
	//		std::cout<<pdgid<<"\t"<<status<<"\t"<< par->pt()<<"\t"<<motherid<<std::endl;
	//		//Leptons.push_back(&*par);
	//	}
	//	else if ( (abs(pdgid)<6|| pdgid==21) && status==23){
	//		std::cout<<pdgid<<"\t"<<status<<"\t"<< par->pt()<<"\t"<<std::endl;
	//	}
	//	//		else if (pdgid==24 || pdgid==23 ) {
	//	//			if (par->mother()!=0 && par->mother()->pdgId()!=pdgid) {
	//	//				int pid = getParentCode((const GenParticle*)&par);
	//	//				if (pid!=25) theAssociatedV.push_back(par);
	//	//			}
	//	//		}

	//}

	FillKFactors(genInfo, genZLeps);

	RunNumber=iEvent.id().run();
	LumiNumber=iEvent.luminosityBlock();
	EventNumber=iEvent.id().event();
	xsection=xsec;
	genxsection=genxsec;
	genbranchingratio=genbr;

    for (std::vector<reco::GenJet>::const_iterator jet = cleanedJets->begin(); jet != cleanedJets->end(); ++jet){
        reco::GenJet jetc = (reco::GenJet)(*jet);
		++nCleanedJets;
		float pt_nominal = jetc.pt();

		if(pt_nominal>30){
			++nCleanedJetsPt30;
		}

		FillJet(jetc); // No additional pT cut (for JEC studies)
	}
	int nFilled=0;
	vector<Int_t> CRFLAG(cands->size());

	for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
		size_t icand= cand-cands->begin();

		if (!(  (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)

		//For the SR, also fold information about acceptance in CRflag.
		if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG[icand],28);
		if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG[icand],29);
		FillCandidate(*cand, true, iEvent, CRFLAG[icand]);

		// Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
		// there could be several candidates per event.
		myTree->FillCurrentTree(true);
		++nFilled;
	}
	if (nFilled==0) {
		if (skipEmptyEvents==false)
			myTree->FillCurrentTree(true);
		else
			myTree->FillCurrentTree(false); //puts it in the failed tree if there is one
	}


}




// ------------ method called once each job just before starting event loop  ------------
	void 
MiniAnalyzer::beginJob()
{
	edm::Service<TFileService> fs;
	TTree *candTree = fs->make<TTree>(theFileName,"Event Summary");
	TTree *candTree_failed = 0;
	if (failedTreeLevel)
		candTree_failed = fs->make<TTree>(theFileName+"_failed","Event Summary");
	myTree = new HZZ4lNtupleFactory(candTree, candTree_failed);
	const int nbins = 45;
	hCounter = fs->make<TH1F>("Counters", "Counters", nbins, 0., nbins);
	BookAllBranches();
	buildMELABranches();
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
MiniAnalyzer::endJob() 

{
	hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
	hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
	hCounter->SetBinContent(2 ,gen_ZZ4mu);
	hCounter->SetBinContent(3 ,gen_ZZ4e);
	hCounter->SetBinContent(4 ,gen_ZZ2mu2e);
	hCounter->SetBinContent(5 ,gen_ZZ2l2tau);
	hCounter->SetBinContent(6 ,gen_ZZ4mu_EtaAcceptance);
	hCounter->SetBinContent(7 ,gen_ZZ4mu_LeptonAcceptance);
	hCounter->SetBinContent(8 ,gen_ZZ2emu2tau);
	hCounter->SetBinContent(9 ,gen_ZZ4tau);
	hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
	hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
	hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
	hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
	hCounter->SetBinContent(19,gen_BUGGY);
	hCounter->SetBinContent(20,gen_Unknown);

	hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
	hCounter->SetBinContent(41,gen_sumGenMCWeight);
	hCounter->SetBinContent(42,gen_sumPUWeight);

	TH1 *h[1] ={ hCounter };
	for (int i = 0; i < 1; i++) {
		h[i]->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
		h[i]->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
		h[i]->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
		h[i]->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
		h[i]->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
		h[i]->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
		h[i]->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
		h[i]->GetXaxis()->SetBinLabel(8 ,"gen_ZZ2emu2tau");
		h[i]->GetXaxis()->SetBinLabel(9 ,"gen_ZZ4tau");
		h[i]->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
		h[i]->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
		h[i]->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
		h[i]->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
		h[i]->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
		h[i]->GetXaxis()->SetBinLabel(20,"gen_Unknown");

		h[i]->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
		h[i]->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
		h[i]->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");
	}

	delete myTree;

	return;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

void MiniAnalyzer::FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
		const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
	GenZ1Mass= pZ1.M();
	GenZ1Pt= pZ1.Pt();
	GenZ1Eta= pZ1.Eta();
	GenZ1Phi= pZ1.Phi();
	GenZ1Flav= Z1Id;

	GenZ2Mass= pZ2.M();
	GenZ2Pt= pZ2.Pt();
	GenZ2Eta= pZ2.Eta();
	GenZ2Phi= pZ2.Phi();
	GenZ2Flav= Z2Id;

	return;
}

void MiniAnalyzer::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
		const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2,
		const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
{
	GenLep1Pt=Lep1.Pt();
	GenLep1Eta=Lep1.Eta();
	GenLep1Phi=Lep1.Phi();
	GenLep1Id=Lep1Id;

	GenLep2Pt=Lep2.Pt();
	GenLep2Eta=Lep2.Eta();
	GenLep2Phi=Lep2.Phi();
	GenLep2Id=Lep2Id;

	GenLep3Pt=Lep3.Pt();
	GenLep3Eta=Lep3.Eta();
	GenLep3Phi=Lep3.Phi();
	GenLep3Id=Lep3Id;

	GenLep4Pt=Lep4.Pt();
	GenLep4Eta=Lep4.Eta();
	GenLep4Phi=Lep4.Phi();
	GenLep4Id=Lep4Id;

	//can comment this back in if Gen angles are needed for any reason...
	//TUtil::computeAngles(zzanalysis::tlv(Lep1), Lep1Id, zzanalysis::tlv(Lep2), Lep2Id, zzanalysis::tlv(Lep3), Lep3Id, zzanalysis::tlv(Lep4), Lep4Id, Gencosthetastar, GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, GenphistarZ1);

	return;
}

void MiniAnalyzer::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
	if (AssocLeps.size() >= 1) {
		GenAssocLep1Pt =AssocLeps.at(0)->p4().Pt();
		GenAssocLep1Eta=AssocLeps.at(0)->p4().Eta();
		GenAssocLep1Phi=AssocLeps.at(0)->p4().Phi();
		GenAssocLep1Id =AssocLeps.at(0)->pdgId();
	}
	if (AssocLeps.size() >= 2) {
		GenAssocLep2Pt =AssocLeps.at(1)->p4().Pt();
		GenAssocLep2Eta=AssocLeps.at(1)->p4().Eta();
		GenAssocLep2Phi=AssocLeps.at(1)->p4().Phi();
		GenAssocLep2Id =AssocLeps.at(1)->pdgId();
	}

	return;
}


void MiniAnalyzer::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
	GenHMass=pH.M();
	GenHPt=pH.Pt();
	GenHRapidity=pH.Rapidity();

	HqTMCweight=w;

	return;
}

void MiniAnalyzer::BookAllBranches(){
	//Event variables
	myTree->Book("RunNumber",RunNumber, failedTreeLevel >= minimalFailedTree);
	myTree->Book("EventNumber",EventNumber, failedTreeLevel >= minimalFailedTree);
	myTree->Book("LumiNumber",LumiNumber, failedTreeLevel >= minimalFailedTree);
	myTree->Book("NRecoMu",NRecoMu, failedTreeLevel >= fullFailedTree);
	myTree->Book("NRecoEle",NRecoEle, failedTreeLevel >= fullFailedTree);
	myTree->Book("Nvtx",Nvtx, failedTreeLevel >= fullFailedTree);
	myTree->Book("NObsInt",NObsInt, failedTreeLevel >= fullFailedTree);
	myTree->Book("NTrueInt",NTrueInt, failedTreeLevel >= fullFailedTree);

	myTree->Book("PFMET",PFMET, failedTreeLevel >= fullFailedTree);
	myTree->Book("nCleanedJets",nCleanedJets, failedTreeLevel >= fullFailedTree);
	myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30, failedTreeLevel >= fullFailedTree);
	myTree->Book("trigWord",trigWord, failedTreeLevel >= minimalFailedTree);
	myTree->Book("evtPassMETFilter",evtPassMETTrigger, failedTreeLevel >= minimalFailedTree);
	myTree->Book("ZZMass",ZZMass, false);
	myTree->Book("ZZsel",ZZsel, false);
	myTree->Book("ZZPt",ZZPt, false);
	myTree->Book("ZZEta",ZZEta, false);
	myTree->Book("ZZPhi",ZZPhi, false);
	myTree->Book("CRflag",CRflag, false);
	myTree->Book("Z1Mass",Z1Mass, false);
	myTree->Book("Z1Pt",Z1Pt, false);
	myTree->Book("Z1Flav",Z1Flav, false);


	//Z2 variables
	myTree->Book("Z2Mass",Z2Mass, false);
	myTree->Book("Z2Pt",Z2Pt, false);
	myTree->Book("Z2Flav",Z2Flav, false);
	myTree->Book("costhetastar",costhetastar, false);
	myTree->Book("helphi",helphi, false);
	myTree->Book("helcosthetaZ1",helcosthetaZ1, false);
	myTree->Book("helcosthetaZ2",helcosthetaZ2, false);
	myTree->Book("phistarZ1",phistarZ1, false);
	myTree->Book("phistarZ2",phistarZ2, false);
	myTree->Book("xi",xi, false);
	myTree->Book("xistar",xistar, false);


	myTree->Book("LepPt",LepPt, false);
	myTree->Book("LepEta",LepEta, false);
	myTree->Book("LepPhi",LepPhi, false);
	myTree->Book("LepLepId",LepLepId, false);


	//Jet variables
	myTree->Book("JetPt",JetPt, failedTreeLevel >= fullFailedTree);
	myTree->Book("JetEta",JetEta, failedTreeLevel >= fullFailedTree);
	myTree->Book("JetPhi",JetPhi, failedTreeLevel >= fullFailedTree);
	myTree->Book("JetMass",JetMass, failedTreeLevel >= fullFailedTree);

	myTree->Book("DiJetMass",DiJetMass, false);
	myTree->Book("DiJetDEta",DiJetDEta, false);
	myTree->Book("nExtraLep",nExtraLep, false);
	myTree->Book("nExtraZ",nExtraZ, false);
	myTree->Book("ExtraLepPt",ExtraLepPt, false);
	myTree->Book("ExtraLepEta",ExtraLepEta, false);
	myTree->Book("ExtraLepPhi",ExtraLepPhi, false);
	myTree->Book("ExtraLepLepId",ExtraLepLepId, false);

	myTree->Book("ZXFakeweight", ZXFakeweight, false);

	if (apply_K_NNLOQCD_ZZGG>0){
	  myTree->Book("KFactor_QCD_ggZZ_Nominal", KFactor_QCD_ggZZ_Nominal, failedTreeLevel >= minimalFailedTree);
	}
	if (apply_K_NLOEW_ZZQQB){
	  myTree->Book("KFactor_EW_qqZZ", KFactor_EW_qqZZ, failedTreeLevel >= minimalFailedTree);
	}
	if (apply_K_NNLOQCD_ZZQQB){
	  myTree->Book("KFactor_QCD_qqZZ_dPhi", KFactor_QCD_qqZZ_dPhi, failedTreeLevel >= minimalFailedTree);
	  myTree->Book("KFactor_QCD_qqZZ_M", KFactor_QCD_qqZZ_M, failedTreeLevel >= minimalFailedTree);
	  myTree->Book("KFactor_QCD_qqZZ_Pt", KFactor_QCD_qqZZ_Pt, failedTreeLevel >= minimalFailedTree);
	}

	myTree->Book("genFinalState", genFinalState, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genProcessId", genProcessId, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genHEPMCweight", genHEPMCweight, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genHEPMCweight_NNLO", genHEPMCweight_NNLO, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genHEPMCweight_POWHEGonly", genHEPMCweight_POWHEGonly, failedTreeLevel >= minimalFailedTree);
	myTree->Book("PUWeight", PUWeight, failedTreeLevel >= minimalFailedTree);
	myTree->Book("dataMCWeight", dataMCWeight, false);
	myTree->Book("trigEffWeight", trigEffWeight, false);
	myTree->Book("overallEventWeight", overallEventWeight, false);
	myTree->Book("HqTMCweight", HqTMCweight, failedTreeLevel >= minimalFailedTree);
	myTree->Book("xsec", xsection, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genxsec", genxsection, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genBR", genbranchingratio, failedTreeLevel >= minimalFailedTree);
	myTree->Book("genExtInfo", genExtInfo, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenHMass", GenHMass, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenHPt", GenHPt, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenHRapidity", GenHRapidity, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenZ1Mass", GenZ1Mass, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ1Pt", GenZ1Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ1Phi", GenZ1Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ1Flav", GenZ1Flav, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenZ2Mass", GenZ2Mass, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ2Pt", GenZ2Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ2Phi", GenZ2Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenZ2Flav", GenZ2Flav, failedTreeLevel >= minimalFailedTree);
	myTree->Book("GenLep1Pt", GenLep1Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep1Eta", GenLep1Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep1Phi", GenLep1Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep1Id", GenLep1Id, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep2Pt", GenLep2Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep2Eta", GenLep2Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep2Phi", GenLep2Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep2Id", GenLep2Id, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep3Pt", GenLep3Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep3Eta", GenLep3Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep3Phi", GenLep3Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep3Id", GenLep3Id, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep4Pt", GenLep4Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep4Eta", GenLep4Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep4Phi", GenLep4Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenLep4Id", GenLep4Id, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep1Pt", GenAssocLep1Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep1Eta", GenAssocLep1Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep1Phi", GenAssocLep1Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep1Id", GenAssocLep1Id, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep2Pt", GenAssocLep2Pt, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep2Eta", GenAssocLep2Eta, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep2Phi", GenAssocLep2Phi, failedTreeLevel >= fullFailedTree);
	myTree->Book("GenAssocLep2Id", GenAssocLep2Id, failedTreeLevel >= fullFailedTree);


	myTree->Book("LHEMotherPz", LHEMotherPz, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEMotherE", LHEMotherE, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEMotherId", LHEMotherId, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEDaughterPt", LHEDaughterPt, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEDaughterEta", LHEDaughterEta, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEDaughterPhi", LHEDaughterPhi, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEDaughterMass", LHEDaughterMass, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEDaughterId", LHEDaughterId, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEAssociatedParticlePt", LHEAssociatedParticlePt, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEAssociatedParticleEta", LHEAssociatedParticleEta, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEAssociatedParticlePhi", LHEAssociatedParticlePhi, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEAssociatedParticleMass", LHEAssociatedParticleMass, failedTreeLevel >= LHEFailedTree);
	myTree->Book("LHEAssociatedParticleId", LHEAssociatedParticleId, failedTreeLevel >= LHEFailedTree);

	myTree->Book("LHEPDFScale", LHEPDFScale, failedTreeLevel >= minimalFailedTree);

	// MELA branches are booked under buildMELA
}

void MiniAnalyzer::addweight(float &weight, float weighttoadd) {
	weight += weighttoadd;
}

void MiniAnalyzer::buildMELABranches(){
	/***********************/
	/***********************/
	/**   Reco branches   **/
	/***********************/
	/***********************/
	for (unsigned int it=0; it<recoMElist.size(); it++){
		MELAOptionParser* me_opt = new MELAOptionParser(recoMElist.at(it));
		if (recoMElist.at(it).find("Copy")!=string::npos) recome_copyopts.push_back(me_opt);
		else recome_originalopts.push_back(me_opt);
	}
	// Resolve original options
	for (unsigned int it=0; it<recome_originalopts.size(); it++){
		MELAOptionParser* me_opt = recome_originalopts.at(it);
		myTree->BookMELABranches(me_opt, false, 0);
	}
	// Resolve copy options
	for (unsigned int it=0; it<recome_copyopts.size(); it++){
		MELAOptionParser* me_opt = recome_copyopts.at(it);
		MELAOptionParser* original_opt=0;
		// Find the original options
		for (unsigned int ih=0; ih<recome_originalopts.size(); ih++){
			if (me_opt->testCopyAlias(recome_originalopts.at(ih)->getAlias())){
				original_opt = recome_originalopts.at(ih);
				break;
			}
		}
		if (original_opt==0) continue;
		else me_opt->pickOriginalOptions(original_opt);
		myTree->BookMELABranches(me_opt, false, 0);
	}

	/**********************/
	/**********************/
	/**   LHE branches   **/
	/**********************/
	/**********************/
	for (unsigned int it=0; it<lheMElist.size(); it++){
		MELAOptionParser* lheme_opt;
		// First find out if the option has a copy specification
		// These copy options will be evaulated in a separate loop
		if (lheMElist.at(it).find("Copy")!=string::npos){
			lheme_opt = new MELAOptionParser(lheMElist.at(it));
			lheme_copyopts.push_back(lheme_opt);
			continue;
		}

		// Create a hypothesis for each option
		MELAHypothesis* lheme_hypo = new MELAHypothesis(&mela, lheMElist.at(it));
		lheme_units.push_back(lheme_hypo);

		lheme_opt = lheme_hypo->getOption();
		if (lheme_opt->isAliased()) lheme_aliased_units.push_back(lheme_hypo);

		// Create a computation for each hypothesis
		MELAComputation* lheme_computer = new MELAComputation(lheme_hypo);
		lheme_computers.push_back(lheme_computer);

		// Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
		GMECHelperFunctions::addToMELACluster(lheme_computer, lheme_clusters);

		// Create the necessary branches for each computation
		myTree->BookMELABranches(lheme_opt, true, lheme_computer);
	}
	// Resolve copy options
	for (unsigned int it=0; it<lheme_copyopts.size(); it++){
		MELAOptionParser* lheme_opt = lheme_copyopts.at(it);
		MELAHypothesis* original_hypo=0;
		MELAOptionParser* original_opt=0;
		// Find the original options
		for (unsigned int ih=0; ih<lheme_aliased_units.size(); ih++){
			if (lheme_opt->testCopyAlias(lheme_aliased_units.at(ih)->getOption()->getAlias())){
				original_hypo = lheme_aliased_units.at(ih);
				original_opt = original_hypo->getOption();
				break;
			}
		}
		if (original_opt==0) continue;
		else lheme_opt->pickOriginalOptions(original_opt);
		// Create a new computation for the copy options
		MELAComputation* lheme_computer = new MELAComputation(original_hypo);
		lheme_computer->setOption(lheme_opt);
		lheme_computers.push_back(lheme_computer);

		// The rest is the same story...
		// Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
		GMECHelperFunctions::addToMELACluster(lheme_computer, lheme_clusters);

		// Create the necessary branches for each computation
		myTree->BookMELABranches(lheme_opt, true, lheme_computer);
	}
	// Loop over the computations to add any contingencies to aliased hypotheses
	for (unsigned int it=0; it<lheme_computers.size(); it++) lheme_computers.at(it)->addContingencies(lheme_aliased_units);

	if (DEBUG_MB){
		std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
		for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->Print();
		for (unsigned int icl=0; icl<lheme_clusters.size(); icl++) cout << "LHE ME cluster " << lheme_clusters.at(icl)->getName() << " is present in " << lheme_clusters.size() << " clusters with #Computations = " << lheme_clusters.at(icl)->getComputations()->size() << endl;
	}
}

void MiniAnalyzer::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG)
{
	//Initialize a new candidate into the tree
	//myTree->createNewCandidate(); // this doesn't do anything anymore

	//Reinitialize the per-candidate vectors (necessary because in CRs we can store more than 1 candidate per event)
	LepPt.clear();
	LepEta.clear();
	LepPhi.clear();
	LepLepId.clear();

	ExtraLepPt.clear();
	ExtraLepEta.clear();
	ExtraLepPhi.clear();
	ExtraLepLepId.clear();

	CRflag = CRFLAG;

	//Fill the info on the Higgs candidate
	ZZMass = cand.p4().mass();

	ZZPt  = cand.p4().pt();
	ZZEta = cand.p4().eta();
	ZZPhi = cand.p4().phi();


	DiJetMass  = cand.userFloat("DiJetMass");
	DiJetDEta  = cand.userFloat("DiJetDEta");

	//Fill the angular variables
	helcosthetaZ1 = cand.userFloat("costheta1");
	helcosthetaZ2 = cand.userFloat("costheta2");
	helphi       = cand.userFloat("phi");
	costhetastar = cand.userFloat("costhetastar");
	phistarZ1      = cand.userFloat("phistar1");
	//phistarZ2      = cand.userFloat("phistar2");
	xi            = cand.userFloat("xi");
	xistar        = cand.userFloat("xistar");

	// Get MELA probabilities
	pushRecoMELABranches(cand);


	//Z1 and Z2 variables
	const reco::Candidate* Z1;
	const reco::Candidate* Z2;
	vector<const reco::Candidate*> leptons;
	vector<const reco::Candidate*> fsrPhot;
	vector<short> fsrIndex;
	vector<string> labels;

	Z1   = cand.daughter("Z1");
	Z2   = cand.daughter("Z2");
	userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex);

	Z1Mass = Z1->mass();
	Z1Pt =   Z1->pt();
	Z1Flav =  getPdgId(Z1->daughter(0)) * getPdgId(Z1->daughter(1));

	Z2Mass = Z2->mass();
	Z2Pt =   Z2->pt();
	Z2Flav = getPdgId(Z2->daughter(0)) * getPdgId(Z2->daughter(1));


	Int_t sel = 0;

	// Precomputed selections
	bool candPass70Z2Loose = cand.userFloat("Z2Mass") &&
		cand.userFloat("MAllComb") &&
		cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
		ZZMass>70.;
	bool candPassFullSel70 = cand.userFloat("SR");
	bool candPassFullSel   = cand.userFloat("FullSel");
	bool candIsBest = cand.userFloat("isBestCand");
	bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut
	//std::cout<<"best "<<candIsBest<<"\t"<<candPass70Z2Loose<<"\t"<<candPassFullSel70<<"\t"<<candPassFullSel<<"\t"<<passMz_zz<<std::endl;

	if (candIsBest) {
		//    sel = 10; //FIXME see above
		if (candPass70Z2Loose) sel=70;
		if (candPassFullSel70){ // includes MZ2 > 12
			sel = 90;
			if (candPassFullSel){
				sel=100;
				if (passMz_zz) sel = 120;
			}
		}
	}



	if (!(evtPass)) {sel = -sel;} // avoid confusion when we write events which do not pass trigger/skim
	ZZsel = sel;


	for (unsigned int i=0; i<leptons.size(); ++i){

		//Fill the info on the lepton candidates
		LepPt .push_back( leptons[i]->pt() );
		LepEta.push_back( leptons[i]->eta() );
		LepPhi.push_back( leptons[i]->phi() );
		int id =  leptons[i]->pdgId();
		if(id == 22 && (i == 1 || i == 3)) id=-22; //FIXME this assumes a standard ordering of leptons.
		LepLepId.push_back( id );

	}
	//Fill the info on categorization
	nExtraLep = cand.userFloat("nExtraLep"); // Why is this still a float at this point?
	nExtraZ=cand.userFloat("nExtraZ");

	//Fill the info on the extra leptons
	for (int iExtraLep=1; iExtraLep<=(int)nExtraLep; iExtraLep++){
		TString extraString;extraString.Form("ExtraLep%d",iExtraLep);
		if (cand.hasUserCand(extraString.Data())){
			//for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
			reco::CandidatePtr candPtr=cand.userCand(extraString.Data());
			ExtraLepPt.push_back(candPtr->pt());
			ExtraLepEta.push_back(candPtr->eta());
			ExtraLepPhi.push_back(candPtr->phi());
			ExtraLepLepId.push_back(candPtr->pdgId());
		}

	}

	//Compute the data/MC weight
	dataMCWeight = 1.;
	trigEffWeight = 1.;

	dataMCWeight = 1.; 

	//Store an overall event weight (which is normalized by gen_sumWeights)
	overallEventWeight = genHEPMCweight * dataMCWeight * trigEffWeight;

}
void MiniAnalyzer::FillJet(const reco::GenJet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());

}

float MiniAnalyzer::EvalSpline(TSpline3* const& sp, float xval){
  double xmin = sp->GetXmin();
  double xmax = sp->GetXmax();
  double res=0;
  if (xval<xmin){
    res=sp->Eval(xmin);
    double deriv=sp->Derivative(xmin);
    res += deriv*(xval-xmin);
  }
  else if (xval>xmax){
    res=sp->Eval(xmax);
    double deriv=sp->Derivative(xmax);
    res += deriv*(xval-xmax);
  }
  else res=sp->Eval(xval);
  return res;
}

void MiniAnalyzer::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps){
  KFactor_QCD_ggZZ_Nominal=1;
  KFactor_QCD_qqZZ_dPhi=1;
  KFactor_QCD_qqZZ_M=1;
  KFactor_QCD_qqZZ_Pt=1;
  KFactor_EW_qqZZ=1;

  GenEventInfoProduct  genInfoP = *(genInfo.product());
  if (apply_K_NNLOQCD_ZZGG>0 && apply_K_NNLOQCD_ZZGG!=3){
    if (spkfactor_ggzz_nnlo[0]!=0) KFactor_QCD_ggZZ_Nominal = MiniAnalyzer::EvalSpline(spkfactor_ggzz_nnlo[0], GenHMass);
    if (apply_K_NNLOQCD_ZZGG==2){
      if (spkfactor_ggzz_nlo[0]!=0){
        float divisor = MiniAnalyzer::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
        KFactor_QCD_ggZZ_Nominal /= divisor;
      }
      else{
        KFactor_QCD_ggZZ_Nominal=0;
      }
    }
  }
  else if (apply_K_NNLOQCD_ZZGG==3){
    if (spkfactor_ggzz_nlo[0]!=0) KFactor_QCD_ggZZ_Nominal = MiniAnalyzer::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
  }

  if (genFinalState!=BUGGY){
    if (genZLeps.size()==4) {
      // Calculate NNLO/NLO QCD K factors for qqZZ
      if (apply_K_NNLOQCD_ZZQQB){
        bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
        KFactor_QCD_qqZZ_dPhi = kfactor_qqZZ_qcd_dPhi(fabs(GenZ1Phi-GenZ2Phi), (sameflavor) ? 1 : 2);
        KFactor_QCD_qqZZ_M    = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2) / kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 1);
        KFactor_QCD_qqZZ_Pt   = kfactor_qqZZ_qcd_Pt(GenHPt, (sameflavor) ? 1 : 2);
      }
      // Calculate NLO EWK K factors for qqZZ
      if (apply_K_NLOEW_ZZQQB){
        TLorentzVector GENZ1Vec, GENZ2Vec, GENZZVec;
        GENZ1Vec.SetPtEtaPhiM(GenZ1Pt, GenZ1Eta, GenZ1Phi, GenZ1Mass);
        GENZ2Vec.SetPtEtaPhiM(GenZ2Pt, GenZ2Eta, GenZ2Phi, GenZ2Mass);
        GENZZVec = GENZ1Vec + GENZ2Vec;

        KFactor_EW_qqZZ = EwkCorrections::getEwkCorrections(pruned, ewkTable, genInfoP, GENZ1Vec, GENZ2Vec);
      }
    }
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);

