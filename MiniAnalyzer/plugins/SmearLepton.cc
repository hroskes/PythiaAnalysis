// -*- C++ -*-
//
// Package:    PythiaAnalysis/MiniAnalyzer
// Class:      SmearLepton
// 
/**\class SmearLepton SmearLepton.cc PythiaAnalysis/MiniAnalyzer/plugins/SmearLepton.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Meng Xiao
//         Created:  Mon, 30 Jul 2018 16:39:47 GMT
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

#include "DataFormats/Candidate/interface/Candidate.h"
#include <vector>

//
// class declaration
//
using namespace reco;
using namespace edm;
using namespace std;


class SmearLepton : public edm::stream::EDProducer<> {
	public:
		explicit SmearLepton(const edm::ParameterSet&);
		~SmearLepton();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginStream(edm::StreamID) override;
		virtual void produce(edm::Event&, const edm::EventSetup&) override;
		virtual void endStream() override;
                edm::EDGetTokenT<edm::View<reco::Candidate> > muonToken_;
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
SmearLepton::SmearLepton(const edm::ParameterSet& iConfig):
	muonToken_(consumes<edm::View<reco::Candidate> >( iConfig.getParameter<edm::InputTag>("src")))

{
	rgen_ = new TRandom3();
	//
	  produces<reco::CandidateCollection>();


}


SmearLepton::~SmearLepton()
{

	// do anything here that needs to be done at destruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
	void
SmearLepton::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	Handle<edm::View<reco::Candidate>> pruned;
	iEvent.getByToken(muonToken_,pruned);
	auto outputMuons = std::make_unique<reco::CandidateCollection >();


	for( View<Candidate>::const_iterator par = pruned->begin(); par != pruned->end(); ++ par ) {
		int pdgid = par->pdgId();
		int status= par->status();
		if( (abs(pdgid)==11 || abs(pdgid)==13) && status==1){
//				reco::Candidate* mu = (reco::Candidate*)&*par; 
				reco::Candidate* mu = par->clone(); 
				TLorentzVector p4;
				double pt = par->pt();
				double newpt= rgen_->Gaus (pt, pt*0.02);
				if( abs(pdgid)==13 && newpt>5 && abs( par->eta())<2.4 ){
                     int hasmother=(par->mother()!=0);
                      int motherid =hasmother? par->mother()->pdgId():-1;
                      if(motherid == 23 || motherid==pdgid || motherid==25){
					p4.SetPtEtaPhiM(newpt, par->eta(), par->phi(), par->mass());
					mu->setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu->mass()));	
//					std::cout<< pt<<"\t"<<mu->pt()<<std::endl;
					outputMuons->push_back(mu);
}
				}
				else if (abs(pdgid)==11 && newpt>7 && abs( par->eta())<2.4){
                     int hasmother=(par->mother()!=0);
                      int motherid =hasmother? par->mother()->pdgId():-1;
                      if(motherid == 23 || motherid==pdgid || motherid==25){
					p4.SetPtEtaPhiM(newpt, par->eta(), par->phi(), par->mass());
					mu->setP4(reco::Particle::PolarLorentzVector(p4.Pt(), p4.Eta(), p4.Phi(), mu->mass()));	
//					std::cout<< pt<<"\t"<<mu->pt()<<std::endl;
					outputMuons->push_back(mu);
				}

			}
		}
}
	  iEvent.put(std::move(outputMuons));


}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
	void
SmearLepton::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SmearLepton::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   SmearLepton::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   SmearLepton::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   SmearLepton::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   SmearLepton::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SmearLepton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SmearLepton);
