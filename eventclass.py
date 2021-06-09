import contextlib, itertools, methodtools, more_itertools, ROOT
from JHUGenMELA.MELA.lhefile import LHEEvent, LHEFileBase, LHEFile_Hwithdecay

class NumpyImport(object):
  def __getattr__(self, attr):
    global np
    import numpy as np
    return getattr(np, attr)
np = NumpyImport()

class LHEEvent_Hwithdecay_smear(LHEEvent):
  #smearing is not actually implemented yet
  @classmethod
  def extracteventparticles(cls, lines, isgen):
    assert not isgen
    daughters, mothers, associated = [], [], []
    ids = [None]
    mother1s = [None]
    mother2s = [None]
    for line in lines:
      id, status, mother1, mother2 = (int(_) for _ in line.split()[0:4])
      ids.append(id)
      if (1 <= abs(id) <= 6 or abs(id) == 21) and not isgen:
        line = line.replace(str(id), "0", 1)  #replace the first instance of the jet id with 0, which means unknown jet
      mother1s.append(mother1)
      mother2s.append(mother2)
      if status == -1:
        mothers.append(line)
      elif status == 1 and (1 <= abs(id) <= 6 or 11 <= abs(id) <= 16 or abs(id) in (21, 22)):
        while True:
          if mother1 != mother2 or mother1 is None:
            associated.append(line)
            break
          if ids[mother1] in (25, 39):
            daughters.append(line)
            break
          mother2 = mother2s[mother1]
          mother1 = mother1s[mother1]

    if not isgen: mothers = None
    return daughters, associated, mothers

class LHEFile_Hwithdecay_smear(LHEFileBase):
  lheeventclass = LHEEvent_Hwithdecay_smear
  @property
  def finalstateparticles(self):
    return itertools.chain(self.daughters, self.associated)

class Event(object):
  def __init__(self, i, gen, reco):
    self.__i = i
    self.__gen = gen
    self.__reco = reco

  @property
  def RunNumber(self): return 1
  @property
  def LumiNumber(self): return 1
  @property
  def EventNumber(self): return self.__i
  @property
  def NRecoMu(self): return sum(1 for id, p in self.__reco.finalstateparticles if abs(id) == 13)
  @property
  def NRecoEle(self): return sum(1 for id, p in self.__reco.finalstateparticles if abs(id) == 11)
  @property
  def Nvtx(self): return 0
  @property
  def NObsInt(self): return 0
  @property
  def NTrueInt(self): return 0
  @property
  def PFMET(self): return sum((p for id, p in self.__reco.finalstateparticles if abs(id) in {12, 14, 16}), ROOT.TLorentzVector()).Pt()
  @property
  def nCleanedJets(self): return sum(1 for id, p in self.__reco.finalstateparticles if abs(id) == 0)
  @property
  def nCleanedJetsPt30(self): return sum(1 for id, p in self.__reco.finalstateparticles if abs(id) == 0 and p.Pt() >= 30)
  @property
  def evtPassMETFilter(self): return 0
  @property
  def trigWord(self): return 0
  @methodtools.lru_cache()
  @property
  def ZZp4(self): return sum((p for id, p in self.__reco.daughters), ROOT.TLorentzVector())
  @methodtools.lru_cache()
  @property
  def ZZMass(self): return self.ZZp4.M()
  @property
  def ZZsel(self): return 0
  @property
  def ZZPt(self): return self.ZZp4.Pt()
  @property
  def ZZEta(self): return self.ZZp4.Eta()
  @property
  def ZZPhi(self): return self.ZZp4.Phi()
  @property
  def CRflag(self): return 0
  @methodtools.lru_cache()
  @property
  def sortedleptons(self):
    possibleZs = [sorted(pair, key=lambda x: x.first) for pair in itertools.combinations(self.__reco.daughters, 2) if abs(pair[0].first) in {11, 13} and sum(p.first for p in pair) == 0]
    Z1pair = min(possibleZs, key=lambda x: abs(sum((p for id, p in x), ROOT.TLorentzVector()).M()-125))
    l1p, l1m = Z1pair
    Z2pair, = {(l2p, l2m) for l2p, l2m in possibleZs if l2p is not l1p and l2m is not l1m}
    return Z1pair[0], Z1pair[1], Z2pair[0], Z2pair[1]
  @methodtools.lru_cache()
  @property
  def Z1(self): return self.sortedleptons[:2]
  @methodtools.lru_cache()
  @property
  def Z1p4(self): return sum((p for id, p in self.Z1), ROOT.TLorentzVector())
  @methodtools.lru_cache()
  @property
  def Z2(self): return self.sortedleptons[2:]
  @methodtools.lru_cache()
  @property
  def Z2p4(self): return sum((p for id, p in self.Z2), ROOT.TLorentzVector())
  @methodtools.lru_cache()
  @property
  def Z1Mass(self): return self.Z1p4.M()
  @property
  def Z1Pt(self): return self.Z1p4.Pt()
  @property
  def Z1Flav(self): return np.product([id for id, p in self.Z1])
  @methodtools.lru_cache()
  @property
  def Z2Mass(self): return self.Z2p4.M()
  @property
  def Z2Pt(self): return self.Z2p4.Pt()
  @property
  def Z2Flav(self): return np.product([id for id, p in self.Z2])
  @methodtools.lru_cache()
  @property
  def decayangles(self):
    angles = self.__reco.computeDecayAngles()
    np.testing.assert_almost_equal(angles.qH, self.ZZMass, decimal=5)
    np.testing.assert_almost_equal(angles.m1, self.Z1Mass, decimal=5)
    np.testing.assert_almost_equal(angles.m2, self.Z2Mass, decimal=5)
    return angles
  @property
  def costhetastar(self): return self.decayangles.costhetastar
  @property
  def helphi(self): return self.decayangles.Phi
  @property
  def helcosthetaZ1(self): return self.decayangles.costheta1
  @property
  def helcosthetaZ2(self): return self.decayangles.costheta2
  @property
  def phistarZ1(self): return self.decayangles.Phi1
  @property
  def phistarZ2(self): return 0
  @property
  def xi(self): return self.ZZp4.Phi()
  @property
  def xistar(self):
    boosted = ROOT.TLorentzVector(self.Z1p4)
    boosted.Boost(-self.ZZp4.BoostVector())
    return boosted.Phi()
  @property
  def LepPt(self):
    return [p.Pt() for id, p in self.sortedleptons]
  @property
  def LepEta(self):
    return [p.Eta() for id, p in self.sortedleptons]
  @property
  def LepPhi(self):
    return [p.Phi() for id, p in self.sortedleptons]
  @property
  def LepLepId(self):
    return [id for id, p in self.sortedleptons]

  @methodtools.lru_cache()
  @property
  def sortedextraleptons(self):
    return sorted(((id, p) for id, p in self.__reco.associated if abs(id) in (11, 13)), key=lambda x: x.second.Pt())
  @property
  def ExtraLepPt(self): return [p.Pt() for id, p in self.sortedextraleptons]
  @property
  def ExtraLepEta(self): return [p.Eta() for id, p in self.sortedextraleptons]
  @property
  def ExtraLepPhi(self): return [p.Phi() for id, p in self.sortedextraleptons]
  @property
  def ExtraLepLepId(self): return [id for id, p in self.sortedextraleptons]
  @property
  def nExtraLep(self):
    return len(self.sortedextraleptons)
  @property
  def nExtraZ(self):
    nep = nem = nmp = nmm = 0
    for id, p in self.sortedextraleptons:
      if id == +11: nem += 1
      if id == -11: nep += 1
      if id == +13: nmm += 1
      if id == -13: nmp += 1
    return min(nem, nep) + min(nmm, nmp)

  @methodtools.lru_cache()
  @property
  def sortedjets(self):
    return sorted(((id, p) for id, p in self.__reco.associated if id == 0), key=lambda x: x.second.Pt())
  @property
  def JetPt(self): return [p.Pt() for id, p in self.sortedjets]
  @property
  def JetEta(self): return [p.Eta() for id, p in self.sortedjets]
  @property
  def JetPhi(self): return [p.Phi() for id, p in self.sortedjets]
  @property
  def JetMass(self): return [p.M() for id, p in self.sortedjets]
  @property
  def DiJetMass(self):
    try:
      (id1, p1), (id2, p2) = self.sortedjets[:2]
    except ValueError:
      return -99
    return (p1+p2).M()
  @property
  def DiJetDEta(self):
    try:
      (id1, p1), (id2, p2) = self.sortedjets[:2]
    except ValueError:
      return -99
    return p1.Eta()-p2.Eta()

  @classmethod
  def branches(cls):
    return [
      "RunNumber",
      "EventNumber",
      "LumiNumber",
      "NRecoMu",
      "NRecoEle",
      "Nvtx",
      "NObsInt",
      "NTrueInt",
      "PFMET",
      "nCleanedJets",
      "nCleanedJetsPt30",
      "trigWord",
      "evtPassMETFilter",
      "ZZMass",
      "ZZsel",
      "ZZPt",
      "ZZEta",
      "ZZPhi",
      "CRflag",
      "Z1Mass",
      "Z1Pt",
      "Z1Flav",
      "Z2Mass",
      "Z2Pt",
      "Z2Flav",
      "costhetastar",
      "helphi",
      "helcosthetaZ1",
      "helcosthetaZ2",
      "phistarZ1",
      "phistarZ2",
      "xi",
      "xistar",
      "LepPt",
      "LepEta",
      "LepPhi",
      "LepLepId",
      "JetPt",
      "JetEta",
      "JetPhi",
      "JetMass",
      "DiJetMass",
      "DiJetDEta",
      "nExtraLep",
      "nExtraZ",
      "ExtraLepPt",
      "ExtraLepEta",
      "ExtraLepPhi",
      "ExtraLepLepId",
    ] or [
      "ZXFakeweight",
      "KFactor_QCD_ggZZ_Nominal",
      "genFinalState",
      "genProcessId",
      "genHEPMCweight",
      "genHEPMCweight_NNLO",
      "genHEPMCweight_POWHEGonly",
      "PUWeight",
      "dataMCWeight",
      "trigEffWeight",
      "overallEventWeight",
      "HqTMCweight",
      "xsec",
      "genxsec",
      "genBR",
      "genExtInfo",
      "GenHMass",
      "GenHPt",
      "GenHRapidity",
      "GenZ1Mass",
      "GenZ1Pt",
      "GenZ1Phi",
      "GenZ1Flav",
      "GenZ2Mass",
      "GenZ2Pt",
      "GenZ2Phi",
      "GenZ2Flav",
      "GenLep1Pt",
      "GenLep1Eta",
      "GenLep1Phi",
      "GenLep1Id",
      "GenLep2Pt",
      "GenLep2Eta",
      "GenLep2Phi",
      "GenLep2Id",
      "GenLep3Pt",
      "GenLep3Eta",
      "GenLep3Phi",
      "GenLep3Id",
      "GenLep4Pt",
      "GenLep4Eta",
      "GenLep4Phi",
      "GenLep4Id",
      "GenAssocLep1Pt",
      "GenAssocLep1Eta",
      "GenAssocLep1Phi",
      "GenAssocLep1Id",
      "GenAssocLep2Pt",
      "GenAssocLep2Eta",
      "GenAssocLep2Phi",
      "GenAssocLep2Id",
      "LHEMotherPz",
      "LHEMotherE",
      "LHEMotherId",
      "LHEDaughterPt",
      "LHEDaughterEta",
      "LHEDaughterPhi",
      "LHEDaughterMass",
      "LHEDaughterId",
      "LHEAssociatedParticlePt",
      "LHEAssociatedParticleEta",
      "LHEAssociatedParticlePhi",
      "LHEAssociatedParticleMass",
      "LHEAssociatedParticleId",
      "LHEPDFScale",
      "pConst_GG_SIG_ghg2_1_ghz1_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",
      "p_GG_SIG_ghg2_1_ghz2_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz4_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
      "p_GG_SIG_ghg2_1_ghza2_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghza4_1_JHUGen",
      "p_GG_SIG_ghg2_1_gha2_1_JHUGen",
      "p_GG_SIG_ghg2_1_gha4_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen",
      "p_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal",
      "p_JVBF_SIG_ghv1_1_JHUGen_JECNominal",
      "pConst_JQCD_SIG_ghg2_1_JHUGen_JECNominal",
      "p_JQCD_SIG_ghg2_1_JHUGen_JECNominal",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghza2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghza4_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_gha2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_gha4_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghza2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_ghza4_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_gha2_1_JHUGen_JECNominal",
      "p_JJVBF_SIG_ghv1_1_gha4_1_JHUGen_JECNominal",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal",
      "p_JJQCD_SIG_ghg4_1_JHUGen_JECNominal",
      "p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal",
      "pConst_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal",
      "p_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal",
      "pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz4_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal",
      "p_HadZH_SIG_ghza2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghza4_1_JHUGen_JECNominal",
      "p_HadZH_SIG_gha2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_gha4_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghza2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_ghza4_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_gha2_1_JHUGen_JECNominal",
      "p_HadZH_SIG_ghz1_1_gha4_1_JHUGen_JECNominal",
      "pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw1_1_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw2_1_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw4_1_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal",
      "p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal",
      "p_ttHUndecayed_SIG_kappa_1_JHUGen_JECNominal",
      "p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECNominal",
      "p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECNominal",
      "p_bbH_SIG_kappa_1_JHUGen_JECNominal",
      "pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp",
      "p_JVBF_SIG_ghv1_1_JHUGen_JECUp",
      "pConst_JQCD_SIG_ghg2_1_JHUGen_JECUp",
      "p_JQCD_SIG_ghg2_1_JHUGen_JECUp",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv4_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECUp",
      "p_JJVBF_SIG_ghza2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghza4_1_JHUGen_JECUp",
      "p_JJVBF_SIG_gha2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_gha4_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghza2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_ghza4_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_gha2_1_JHUGen_JECUp",
      "p_JJVBF_SIG_ghv1_1_gha4_1_JHUGen_JECUp",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECUp",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECUp",
      "p_JJQCD_SIG_ghg4_1_JHUGen_JECUp",
      "p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECUp",
      "pConst_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECUp",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECUp",
      "p_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECUp",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECUp",
      "pConst_HadZH_SIG_ghz1_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECUp",
      "p_HadZH_SIG_ghz2_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz4_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECUp",
      "p_HadZH_SIG_ghza2_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghza4_1_JHUGen_JECUp",
      "p_HadZH_SIG_gha2_1_JHUGen_JECUp",
      "p_HadZH_SIG_gha4_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghza2_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_ghza4_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_gha2_1_JHUGen_JECUp",
      "p_HadZH_SIG_ghz1_1_gha4_1_JHUGen_JECUp",
      "pConst_HadWH_SIG_ghw1_1_JHUGen_JECUp",
      "p_HadWH_SIG_ghw1_1_JHUGen_JECUp",
      "p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECUp",
      "p_HadWH_SIG_ghw2_1_JHUGen_JECUp",
      "p_HadWH_SIG_ghw4_1_JHUGen_JECUp",
      "p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECUp",
      "p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECUp",
      "p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECUp",
      "p_ttHUndecayed_SIG_kappa_1_JHUGen_JECUp",
      "p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECUp",
      "p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECUp",
      "p_bbH_SIG_kappa_1_JHUGen_JECUp",
      "pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn",
      "p_JVBF_SIG_ghv1_1_JHUGen_JECDn",
      "pConst_JQCD_SIG_ghg2_1_JHUGen_JECDn",
      "p_JQCD_SIG_ghg2_1_JHUGen_JECDn",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv4_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECDn",
      "p_JJVBF_SIG_ghza2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghza4_1_JHUGen_JECDn",
      "p_JJVBF_SIG_gha2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_gha4_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghza2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_ghza4_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_gha2_1_JHUGen_JECDn",
      "p_JJVBF_SIG_ghv1_1_gha4_1_JHUGen_JECDn",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECDn",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECDn",
      "p_JJQCD_SIG_ghg4_1_JHUGen_JECDn",
      "p_JJQCD_SIG_ghg2_1_ghg4_1_JHUGen_JECDn",
      "pConst_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECDn",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECDn",
      "p_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECDn",
      "p_JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECDn",
      "pConst_HadZH_SIG_ghz1_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECDn",
      "p_HadZH_SIG_ghz2_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz4_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECDn",
      "p_HadZH_SIG_ghza2_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghza4_1_JHUGen_JECDn",
      "p_HadZH_SIG_gha2_1_JHUGen_JECDn",
      "p_HadZH_SIG_gha4_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghza2_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_ghza4_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_gha2_1_JHUGen_JECDn",
      "p_HadZH_SIG_ghz1_1_gha4_1_JHUGen_JECDn",
      "pConst_HadWH_SIG_ghw1_1_JHUGen_JECDn",
      "p_HadWH_SIG_ghw1_1_JHUGen_JECDn",
      "p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECDn",
      "p_HadWH_SIG_ghw2_1_JHUGen_JECDn",
      "p_HadWH_SIG_ghw4_1_JHUGen_JECDn",
      "p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECDn",
      "p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECDn",
      "p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECDn",
      "p_ttHUndecayed_SIG_kappa_1_JHUGen_JECDn",
      "p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECDn",
      "p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECDn",
      "p_bbH_SIG_kappa_1_JHUGen_JECDn",
      "p_LepZH_SIG_ghz1_1_JHUGen",
      "p_LepZH_SIG_ghz1prime2_1E4_JHUGen",
      "p_LepZH_SIG_ghz2_1_JHUGen",
      "p_LepZH_SIG_ghz4_1_JHUGen",
      "p_LepZH_SIG_ghza1prime2_1E4_JHUGen",
      "p_LepZH_SIG_ghza2_1_JHUGen",
      "p_LepZH_SIG_ghza4_1_JHUGen",
      "p_LepZH_SIG_gha2_1_JHUGen",
      "p_LepZH_SIG_gha4_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghz2_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghz4_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghza2_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_ghza4_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_gha2_1_JHUGen",
      "p_LepZH_SIG_ghz1_1_gha4_1_JHUGen",
      "p_LepWH_SIG_ghw1_1_JHUGen",
      "p_LepWH_SIG_ghw1prime2_1E4_JHUGen",
      "p_LepWH_SIG_ghw2_1_JHUGen",
      "p_LepWH_SIG_ghw4_1_JHUGen",
      "p_LepWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen",
      "p_LepWH_SIG_ghw1_1_ghw2_1_JHUGen",
      "p_LepWH_SIG_ghw1_1_ghw4_1_JHUGen",
      "p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen",
      "p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen",
      "p_INDEPENDENT_SIG_gZPz1_1_JHUGen",
      "p_INDEPENDENT_SIG_gZPz2_1_JHUGen",
      "p_GG_SIG_gXg1_1_gXz1_1_JHUGen",
      "p_GG_SIG_gXg2_1_gXz2_1_JHUGen",
      "p_GG_SIG_gXg3_1_gXz3_1_JHUGen",
      "p_GG_SIG_gXg4_1_gXz4_1_JHUGen",
      "p_GG_SIG_gXg1_1_gXz5_1_JHUGen",
      "p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen",
      "p_GG_SIG_gXg1_1_gXz6_1_JHUGen",
      "p_GG_SIG_gXg1_1_gXz7_1_JHUGen",
      "p_GG_SIG_gXg5_1_gXz8_1_JHUGen",
      "p_GG_SIG_gXg5_1_gXz9_1_JHUGen",
      "p_GG_SIG_gXg5_1_gXz10_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz1_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz2_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz3_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz4_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz5_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz6_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz7_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz8_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz9_1_JHUGen",
      "p_QQB_SIG_XqqLR_1_gXz10_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz1_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz2_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz3_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz4_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz5_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz1_1_gXz5_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz6_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz7_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz8_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz9_1_JHUGen",
      "p_INDEPENDENT_SIG_gXz10_1_JHUGen",
      "pConst_GG_SIG_kappaTopBot_1_ghz1_1_MCFM",
      "p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM",
      "p_GG_BSI_kappaTopBot_1_ghz1_1_MCFM",
      "p_GG_BSI_kappaTopBot_1_ghz1_i_MCFM",
      "pConst_GG_BKG_MCFM",
      "p_GG_BKG_MCFM",
      "pConst_QQB_BKG_MCFM",
      "p_QQB_BKG_MCFM",
      "pConst_ZJJ_BKG_MCFM",
      "p_ZJJ_BKG_MCFM",
      "p_JJEW_SIG_ghv1_1_MCFM_JECNominal",
      "p_JJEW_BSI_ghv1_1_MCFM_JECNominal",
      "p_JJEW_BSI_ghv1_i_MCFM_JECNominal",
      "p_JJEW_BKG_MCFM_JECNominal",
      "pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal",
      "p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal",
      "p_JJVBF_S_BSI_ghv1_1_MCFM_JECNominal",
      "p_JJVBF_S_BSI_ghv1_i_MCFM_JECNominal",
      "p_JJVBF_SIG_ghv1_1_MCFM_JECNominal",
      "p_JJVBF_BSI_ghv1_1_MCFM_JECNominal",
      "p_JJVBF_BSI_ghv1_i_MCFM_JECNominal",
      "pConst_JJVBF_BKG_MCFM_JECNominal",
      "p_JJVBF_BKG_MCFM_JECNominal",
      "pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal",
      "p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal",
      "p_HadZH_S_BSI_ghz1_1_MCFM_JECNominal",
      "p_HadZH_S_BSI_ghz1_i_MCFM_JECNominal",
      "p_HadZH_SIG_ghz1_1_MCFM_JECNominal",
      "p_HadZH_BSI_ghz1_1_MCFM_JECNominal",
      "p_HadZH_BSI_ghz1_i_MCFM_JECNominal",
      "pConst_HadZH_BKG_MCFM_JECNominal",
      "p_HadZH_BKG_MCFM_JECNominal",
      "pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal",
      "p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal",
      "p_HadWH_S_BSI_ghw1_1_MCFM_JECNominal",
      "p_HadWH_S_BSI_ghw1_i_MCFM_JECNominal",
      "pConst_HadWH_BKG_MCFM_JECNominal",
      "p_HadWH_BKG_MCFM_JECNominal",
      "pConst_JJQCD_BKG_MCFM_JECNominal",
      "p_JJQCD_BKG_MCFM_JECNominal",
      "p_JJEW_SIG_ghv1_1_MCFM_JECUp",
      "p_JJEW_BSI_ghv1_1_MCFM_JECUp",
      "p_JJEW_BSI_ghv1_i_MCFM_JECUp",
      "p_JJEW_BKG_MCFM_JECUp",
      "pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECUp",
      "p_JJVBF_S_SIG_ghv1_1_MCFM_JECUp",
      "p_JJVBF_S_BSI_ghv1_1_MCFM_JECUp",
      "p_JJVBF_S_BSI_ghv1_i_MCFM_JECUp",
      "p_JJVBF_SIG_ghv1_1_MCFM_JECUp",
      "p_JJVBF_BSI_ghv1_1_MCFM_JECUp",
      "p_JJVBF_BSI_ghv1_i_MCFM_JECUp",
      "pConst_JJVBF_BKG_MCFM_JECUp",
      "p_JJVBF_BKG_MCFM_JECUp",
      "pConst_HadZH_S_SIG_ghz1_1_MCFM_JECUp",
      "p_HadZH_S_SIG_ghz1_1_MCFM_JECUp",
      "p_HadZH_S_BSI_ghz1_1_MCFM_JECUp",
      "p_HadZH_S_BSI_ghz1_i_MCFM_JECUp",
      "p_HadZH_SIG_ghz1_1_MCFM_JECUp",
      "p_HadZH_BSI_ghz1_1_MCFM_JECUp",
      "p_HadZH_BSI_ghz1_i_MCFM_JECUp",
      "pConst_HadZH_BKG_MCFM_JECUp",
      "p_HadZH_BKG_MCFM_JECUp",
      "pConst_HadWH_S_SIG_ghw1_1_MCFM_JECUp",
      "p_HadWH_S_SIG_ghw1_1_MCFM_JECUp",
      "p_HadWH_S_BSI_ghw1_1_MCFM_JECUp",
      "p_HadWH_S_BSI_ghw1_i_MCFM_JECUp",
      "pConst_HadWH_BKG_MCFM_JECUp",
      "p_HadWH_BKG_MCFM_JECUp",
      "pConst_JJQCD_BKG_MCFM_JECUp",
      "p_JJQCD_BKG_MCFM_JECUp",
      "p_JJEW_SIG_ghv1_1_MCFM_JECDn",
      "p_JJEW_BSI_ghv1_1_MCFM_JECDn",
      "p_JJEW_BSI_ghv1_i_MCFM_JECDn",
      "p_JJEW_BKG_MCFM_JECDn",
      "pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECDn",
      "p_JJVBF_S_SIG_ghv1_1_MCFM_JECDn",
      "p_JJVBF_S_BSI_ghv1_1_MCFM_JECDn",
      "p_JJVBF_S_BSI_ghv1_i_MCFM_JECDn",
      "p_JJVBF_SIG_ghv1_1_MCFM_JECDn",
      "p_JJVBF_BSI_ghv1_1_MCFM_JECDn",
      "p_JJVBF_BSI_ghv1_i_MCFM_JECDn",
      "pConst_JJVBF_BKG_MCFM_JECDn",
      "p_JJVBF_BKG_MCFM_JECDn",
      "pConst_HadZH_S_SIG_ghz1_1_MCFM_JECDn",
      "p_HadZH_S_SIG_ghz1_1_MCFM_JECDn",
      "p_HadZH_S_BSI_ghz1_1_MCFM_JECDn",
      "p_HadZH_S_BSI_ghz1_i_MCFM_JECDn",
      "p_HadZH_SIG_ghz1_1_MCFM_JECDn",
      "p_HadZH_BSI_ghz1_1_MCFM_JECDn",
      "p_HadZH_BSI_ghz1_i_MCFM_JECDn",
      "pConst_HadZH_BKG_MCFM_JECDn",
      "p_HadZH_BKG_MCFM_JECDn",
      "pConst_HadWH_S_SIG_ghw1_1_MCFM_JECDn",
      "p_HadWH_S_SIG_ghw1_1_MCFM_JECDn",
      "p_HadWH_S_BSI_ghw1_1_MCFM_JECDn",
      "p_HadWH_S_BSI_ghw1_i_MCFM_JECDn",
      "pConst_HadWH_BKG_MCFM_JECDn",
      "p_HadWH_BKG_MCFM_JECDn",
      "pConst_JJQCD_BKG_MCFM_JECDn",
      "p_JJQCD_BKG_MCFM_JECDn",
      "p_m4l_SIG",
      "p_m4l_BKG",
      "p_m4l_SIG_ScaleDown",
      "p_m4l_BKG_ScaleDown",
      "p_m4l_SIG_ResDown",
      "p_m4l_BKG_ResDown",
      "p_m4l_SIG_ScaleUp",
      "p_m4l_BKG_ScaleUp",
      "p_m4l_SIG_ResUp",
      "p_m4l_BKG_ResUp",
      "p_HadZH_mavjj_true_JECNominal",
      "p_HadWH_mavjj_true_JECNominal",
      "p_HadZH_mavjj_JECNominal",
      "p_HadWH_mavjj_JECNominal",
      "p_HadZH_mavjj_true_JECUp",
      "p_HadWH_mavjj_true_JECUp",
      "p_HadZH_mavjj_JECUp",
      "p_HadWH_mavjj_JECUp",
      "p_HadZH_mavjj_true_JECDn",
      "p_HadWH_mavjj_true_JECDn",
      "p_HadZH_mavjj_JECDn",
      "p_HadWH_mavjj_JECDn",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal_BestDJJ",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal_BestDJJ",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_BestDJJ",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal_BestDJJ",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECUp_BestDJJ",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECUp_BestDJJ",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECUp_BestDJJ",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECUp_BestDJJ",
      "pConst_JJVBF_SIG_ghv1_1_JHUGen_JECDn_BestDJJ",
      "p_JJVBF_SIG_ghv1_1_JHUGen_JECDn_BestDJJ",
      "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECDn_BestDJJ",
      "p_JJQCD_SIG_ghg2_1_JHUGen_JECDn_BestDJJ",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz2_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghz4_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza2_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_ghza4_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_gha2_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1_1_gha4_i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza1prime2_1E4i_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_ghz4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz2_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghz4_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza2_1_ghza4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza2_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza2_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza4_1_gha2_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_ghza4_1_gha4_1_JHUGen",
      "p_Gen_GG_SIG_ghg2_1_gha2_1_gha4_1_JHUGen",
    ]
