#include <PythiaAnalysis/MiniAnalyzer/interface/DaughterDataHelpers.h>

#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>


using namespace std;
using namespace reco;

void userdatahelpers::embedDaughterData(pat::CompositeCandidate& cand) {

  for (unsigned i = 0; i<cand.numberOfDaughters(); ++i) {
    const reco::Candidate* d = cand.daughter(i);

    if(d->hasMasterClone()) d = d->masterClone().get();
    
    // We need the concrete object to access the method userFloat(). 
    // (A more general solution would be to creat a StringObjectFunction on the fly for each 
    // entry in userFloatNames(). That's maybe too time consuming (to be checked))
    if (const pat::CompositeCandidate* cc = dynamic_cast<const pat::CompositeCandidate*>(d)) {
      embedDaughterData(cand, i, cc);      
    } else if (const pat::Muon* mu = dynamic_cast<const pat::Muon*>(d)) {
      embedDaughterData(cand, i, mu);
    } else if (const pat::Electron* ele = dynamic_cast<const pat::Electron*>(d)) {
      embedDaughterData(cand, i, ele);
    } else if (const pat::Photon* ele = dynamic_cast<const pat::Photon*>(d)) {
      embedDaughterData(cand, i, ele);
    } else {
      //edm::LogError("") << "DaughterDataEmbedder: Unsupported daughter type";
    }
  }
}


float userdatahelpers::getUserFloat(const reco::Candidate* c, const char* name){
  if(c->hasMasterClone()) c = c->masterClone().get();
  if (const pat::Muon* mu = dynamic_cast<const pat::Muon*>(c)) {
    return mu->userFloat(name);
  } else if (const pat::Electron* ele = dynamic_cast<const pat::Electron*>(c)) {
    return ele->userFloat(name);
  } else if (const pat::Photon* ele = dynamic_cast<const pat::Photon*>(c)) {
    return ele->userFloat(name);
  } else if (const pat::CompositeCandidate* cc = dynamic_cast<const pat::CompositeCandidate*>(c)) {
    return cc->userFloat(name);
  }
  edm::LogError("") << "userdatahelpers::getUserFloat: Unsupported daughter type";

  return 0;
}

int userdatahelpers::hasUserFloat(const reco::Candidate* c, const char* name){
  if(c->hasMasterClone()) c = c->masterClone().get();
  if (const pat::Muon* mu = dynamic_cast<const pat::Muon*>(c)) {
    return mu->hasUserFloat(name);
  } else if (const pat::Electron* ele = dynamic_cast<const pat::Electron*>(c)) {
    return ele->hasUserFloat(name);
  } else if (const pat::Photon* ele = dynamic_cast<const pat::Photon*>(c)) {
    return ele->hasUserFloat(name);
  } else if (const pat::CompositeCandidate* cc = dynamic_cast<const pat::CompositeCandidate*>(c)) {
    return cc->hasUserFloat(name);
  }
  edm::LogError("") << "userdatahelpers::hasUserFloat: Unsupported daughter type";

  return -1;
}


void 
userdatahelpers::getSortedLeptons(const pat::CompositeCandidate& cand, vector<const Candidate*>& leptons, vector<string>& labels, vector<const Candidate*>& fsrPhotons, std::vector<short>& fsrIndex, bool is4l) {

  if (is4l) { // Regular 4 lepton SR/CR
    // Pointers to Z, sorted by mass
    vector<const Candidate*> Zs   = {cand.daughter("Z1"), cand.daughter("Z2")};

    //Pointer to leptons (Z11,Z12,Z21,Z22, to be sorted by charge)
    leptons = {Zs[0]->daughter(0), Zs[0]->daughter(1), Zs[1]->daughter(0), Zs[1]->daughter(1)};

    // Set prefixes for ZZCand userFloats
    string Z1Label = "d0.";
    string Z2Label = "d1.";
    if (cand.daughter("Z1")==cand.daughter(1)) swap(Z1Label,Z2Label);
    labels = {Z1Label+"d0.",Z1Label +"d1.", Z2Label+"d0.",Z2Label+"d1."};

    vector<unsigned> lOrder = {0,1,2,3};

    // Sort leptons by charge so that the order is Z1Lp, Z1Ln, Z2Lp, Z2Ln;
    // for TLEs, assume they are opposite-sign to the other lepton.
    // do nothing for the same-sign collections used for CRs
    bool need_swap = false;

    if(abs(leptons[0]->pdgId()) == 22 || abs(leptons[1]->pdgId()) == 22) {
        int non_TLE_index = -1;
        if(abs(leptons[0]->pdgId()) != 22) non_TLE_index = 0;
        if(abs(leptons[1]->pdgId()) != 22) non_TLE_index = 1;   
        if(non_TLE_index == -1) {
	  edm::LogError("") << "Found a Z candidate made of two TLE, this should never happen!";
	  abort();
	}
        if(leptons[non_TLE_index]->charge() < 0 && non_TLE_index == 0) need_swap = true; 
    } else {
      if (leptons[0]->charge() < 0 && leptons[0]->charge()*leptons[1]->charge()<0) {
        need_swap = true;
      }
    }
    if(need_swap) {
        swap(leptons[0],leptons[1]);
        swap(labels[0],labels[1]);
        swap(lOrder[0],lOrder[1]);
    }

    need_swap = false;
    if(abs(leptons[2]->pdgId()) == 22 || abs(leptons[3]->pdgId()) == 22) {
        int non_TLE_index = -1;
        if(abs(leptons[2]->pdgId()) != 22) non_TLE_index = 2;
        if(abs(leptons[3]->pdgId()) != 22) non_TLE_index = 3;   
        if(non_TLE_index == -1) {
	  edm::LogError("") << "Found a Z candidate made of two TLE, this should never happen!";
	  abort();
	}
        if(leptons[non_TLE_index]->charge() < 0 && non_TLE_index == 2) need_swap = true; 
    } else {
      if(leptons[2]->charge() < 0 && leptons[2]->charge()*leptons[3]->charge()<0) {        
        need_swap = true;
      }
    }
    if(need_swap) {
        swap(leptons[2],leptons[3]);
        swap(labels[2],labels[3]);
        swap(lOrder[2],lOrder[3]);
    }

    // Collect FSR
    for (unsigned iZ=0; iZ<2; ++iZ) {
      for (unsigned ifsr=2; ifsr<Zs[iZ]->numberOfDaughters(); ++ifsr) {
	const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Zs[iZ]->daughter(ifsr));
	int ilep = iZ*2+fsr->userFloat("leptIdx");
	fsrPhotons.push_back(fsr);
	fsrIndex.push_back(lOrder[ilep]);
      }
    }

  } else { // Z+l
    const Candidate* Z1 = cand.daughter(0); // the Z    
    leptons = {Z1->daughter(0), Z1->daughter(1), cand.daughter(1)};
    labels = {"d0.d0.","d0.d1.","d1."};
    vector<unsigned> lOrder = {0,1,2};

    if (leptons[0]->charge() < 0 && leptons[0]->charge()*leptons[1]->charge()<0) {
      swap(leptons[0],leptons[1]);
      swap(labels[0],labels[1]);
      swap(lOrder[0],lOrder[1]);
    }

    for (unsigned ifsr=2; ifsr<Z1->numberOfDaughters(); ++ifsr) { //FIXME this will not pick FSR from l
	const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(Z1->daughter(ifsr));
	int ilep = fsr->userFloat("leptIdx");
	fsrPhotons.push_back(fsr);
	fsrIndex.push_back(lOrder[ilep]);
    }
  }
}


void 
userdatahelpers::getSortedZLeptons(const pat::CompositeCandidate& cand, vector<const Candidate*>& leptons, vector<string>& labels, vector<const Candidate*>& fsrPhotons, std::vector<short>& fsrIndex) {

  // Pointer to leptons, to be sorted by charge, in order Lp, Ln
  leptons = {cand.daughter(0), cand.daughter(1)};

  labels = {"d0.","d1."};
  vector<unsigned> lOrder = {0,1};

  if (leptons[0]->charge() < 0 && leptons[0]->charge()*leptons[1]->charge()<0) {
    swap(leptons[0],leptons[1]);
    swap(labels[0],labels[1]);
    swap(lOrder[0],lOrder[1]);
  }
     
  // Collect FSR
  for (unsigned ifsr=2; ifsr<cand.numberOfDaughters(); ++ifsr) {
    const pat::PFParticle* fsr = static_cast<const pat::PFParticle*>(cand.daughter(ifsr));
    int ilep = fsr->userFloat("leptIdx");
    fsrPhotons.push_back(fsr);
    fsrIndex.push_back(lOrder[ilep]);
  }

}

