// -*- C++ -*-
//
// Package:    PythiaAnalysis/MiniAnalyzer
// Class:      SmearJet
// 
/**\class SmearJet SmearJet.cc PythiaAnalysis/MiniAnalyzer/plugins/SmearJet.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Meng Xiao
//         Created:  Mon, 30 Jul 2018 16:39:23 GMT
//
//


// system include files
#include <memory>
#include "TLorentzVector.h"
#include "TRandom3.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"


#include <vector>
using namespace reco;
using namespace edm;
using namespace std;

//
// class declaration
//

class SmearJet : public edm::stream::EDProducer<> {
	public:
		explicit SmearJet(const edm::ParameterSet&);
		~SmearJet();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginStream(edm::StreamID) override;
		virtual void produce(edm::Event&, const edm::EventSetup&) override;
		virtual void endStream() override;
		edm::EDGetTokenT<edm::View<reco::Candidate> > muonToken_;
		edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;
		TRandom3* rgen_;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

		// ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SmearJet::SmearJet(const edm::ParameterSet& iConfig):
	muonToken_(consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"))),
	genJetToken_(consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genjet")))
{
	rgen_ = new TRandom3();
	produces<reco::GenJetCollection>();


}


SmearJet::~SmearJet()
{

	// do anything here that needs to be done at destruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
	void
SmearJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	Handle<edm::View<const reco::GenJet>> genjet;
	iEvent.getByToken(genJetToken_,genjet);
	Handle<edm::View<reco::Candidate>> lepton;
	iEvent.getByToken(muonToken_,lepton);
	auto outputMuons = std::make_unique<reco::GenJetCollection>();

	for (unsigned int i=0;i<genjet->size();i++){
	//std::cout<< "Comning to Jets"<<std::endl;
//	for( View<Candidate>::const_iterator jet= genjet->begin(); jet!= genjet->end(); ++ jet) {

		//reco::GenJet *mu = (reco::GenJet*)&*jet;

		reco::GenJet jet = genjet->at(i);
		reco::GenJet *mu = jet.clone(); 

		double pt = jet.pt();
		double newpt= rgen_->Gaus (pt, pt*0.15);
		//std::cout << "jet "<<pt<<"\t"<<newpt<<std::endl;
		bool overlap = false;
		if( newpt>25 && abs( jet.eta())<4.7){
			TLorentzVector p4;
			p4.SetPtEtaPhiM(newpt, jet.eta(), jet.phi(), jet.mass());
		  for( View<reco::Candidate>::const_iterator lep = lepton->begin(); lep != lepton->end(); ++ lep ){	
			  reco::Candidate::LorentzVector lepp4 =lep->p4();
			if(ROOT::Math::VectorUtil::DeltaR(lepp4.Vect(), jet.p4().Vect())< 0.3){
				overlap= true;
				break;
			}
		  }
		  if(!overlap){
			mu->setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu->mass()));	
//			std::cout<< "pushing jets"<<std::endl;
			outputMuons->push_back(*mu);
//			std::cout<< mu.pt()<<std::endl;
		  }
		}
	}
//	std::cout<< outputMuons->size()<<std::endl; 
		iEvent.put(std::move(outputMuons));


	}

	// ------------ method called once each stream before processing any runs, lumis or events  ------------
	void
		SmearJet::beginStream(edm::StreamID)
		{
		}

	// ------------ method called once each stream after processing all runs, lumis and events  ------------
	void
		SmearJet::endStream() {
		}

	// ------------ method called when starting to processes a run  ------------
	/*
	   void
	   SmearJet::beginRun(edm::Run const&, edm::EventSetup const&)
	   {
	   }
	   */

	// ------------ method called when ending the processing of a run  ------------
	/*
	   void
	   SmearJet::endRun(edm::Run const&, edm::EventSetup const&)
	   {
	   }
	   */

	// ------------ method called when starting to processes a luminosity block  ------------
	/*
	   void
	   SmearJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
	   {
	   }
	   */

	// ------------ method called when ending the processing of a luminosity block  ------------
	/*
	   void
	   SmearJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
	   {
	   }
	   */

	// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
	void
		SmearJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			//The following says we do not know what parameters are allowed so do no validation
			// Please change this to state exactly what you do use, even if it is no parameters
			edm::ParameterSetDescription desc;
			desc.setUnknown();
			descriptions.addDefault(desc);
		}

	//define this as a plug-in
	DEFINE_FWK_MODULE(SmearJet);
