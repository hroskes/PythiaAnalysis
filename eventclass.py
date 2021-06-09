import contextlib, itertools, methodtools, more_itertools, ROOT
from JHUGenMELA.MELA.lhefile import LHEEvent, LHEFileBase, LHEFile_Hwithdecay

class NumpyImport(object):
  def __getattr__(self, attr):
    global np
    import numpy as np
    return getattr(np, attr)
np = NumpyImport()

class TVarImport(object):
  def __getattr__(self, attr):
    global TVar
    from JHUGenMELA.MELA.mela import TVar
    return getattr(TVar, attr)
TVar = TVarImport()

class ProbabilityLine(object):
  def __init__(self, name, production, process, matrixelement, alias=None, options={}, couplings={}, defaultme=None, cluster=None, forceincomingflavors=None, ispm4l=None, supermelasyst=None, ispmavjjtrue=None, ispmavjj=None):
    self.name = name
    self.alias = alias
    self.production = production
    self.process = process
    self.matrixelement = matrixelement
    self.options = options
    self.couplings = couplings
    self.defaultme = defaultme
    self.cluster = cluster
    self.forceincomingflavors = forceincomingflavors
    self.ispm4l = ispm4l
    self.supermelasyst = supermelasyst
    self.ispmavjjtrue = ispmavjjtrue
    self.ispmavjj = ispmavjj
  @classmethod
  def kwargsfromline(cls, line, alllines):
    kwargs = {k.lower(): v for k, v in (_.split(":") for _ in line.split(" "))}
    if "copy" in kwargs:
      copyline, = {line for line in alllines if "Name:"+kwargs["copy"]+" " in line}
      copykwargs = cls.kwargsfromline(copyline, alllines)
      copykwargs.update(kwargs)
      kwargs = copykwargs
      del kwargs["copy"]

    for thing in "production", "process", "matrixelement":
      if isinstance(kwargs[thing], str):
        kwargs[thing] = getattr(TVar, kwargs[thing])

    for thing in "options", "couplings":
      if thing in kwargs:
        kwargs[thing] = {k: v for k, v in (_.split("=") for _ in kwargs[thing].split(";"))}

    for name, val in kwargs.get("couplings", {}).items():
      real, imag = val.split(",")
      real = float(real)
      imag = float(imag)
      kwargs["couplings"][name] = real + imag*1j

    if kwargs.get("alias", None) == "<Name>":
      kwargs["alias"] = kwargs["name"]

    if "defaultme" in kwargs:
      kwargs["defaultme"] = float(kwargs["defaultme"])

    return kwargs

  @classmethod
  def fromline(cls, line, alllines):
    return cls(**cls.kwargsfromline(line, alllines))


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
  def __init__(self, i, gen, reco, xsec, genxsec, genBR):
    self.__i = i
    self.__gen = gen
    self.__reco = reco
    self.__xsec = xsec
    self.__genxsec = genxsec
    self.__genBR = genBR

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
    possibleZs = [sorted(pair, key=lambda x: x[0]) for pair in itertools.combinations(self.__reco.daughters, 2) if abs(pair[0][0]) in {11, 13} and sum(p[0] for p in pair) == 0]
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
    return sorted(((id, p) for id, p in self.__reco.associated if abs(id) in (11, 13)), key=lambda x: x[1].Pt())
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
    return sorted(((id, p) for id, p in self.__reco.associated if id == 0), key=lambda x: x[1].Pt())
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

  @property
  def ZXFakeweight(self): return 0
  @property
  def dataMCWeight(self): return 1
  @property
  def trigEffWeight(self): return 1
  @property
  def overallEventWeight(self): return self.genHEPMCweight
  @property
  def HqTMCweight(self): return 1
  @property
  def PUWeight(self): return 1
  @property
  def genHEPMCweight_NNLO(self): return self.genHEPMCweight
  @property
  def genHEPMCweight_POWHEGonly(self): return self.genHEPMCweight
  @property
  def trigEffWeight(self): return 1
  @property
  def genFinalState(self):
    return {
      169*169: 0,
      121*121: 1,
      169*121: 2,
    }[self.GenZ1Flav*self.GenZ2Flav]
  @property
  def genProcessId(self): return 0
  @property
  def genHEPMCweight(self): return self.__gen.weight
  @property
  def KFactor_QCD_ggZZ_Nominal(self): return 1.25752

  @methodtools.lru_cache()
  @property
  def Gensortedleptons(self):
    possibleZs = [sorted(pair, key=lambda x: x[0]) for pair in itertools.combinations(self.__gen.daughters, 2) if abs(pair[0][0]) in {11, 13} and sum(p[0] for p in pair) == 0]
    Z1pair = min(possibleZs, key=lambda x: abs(sum((p for id, p in x), ROOT.TLorentzVector()).M()-125))
    l1p, l1m = Z1pair
    Z2pair, = {(l2p, l2m) for l2p, l2m in possibleZs if l2p is not l1p and l2m is not l1m}
    return Z1pair[0], Z1pair[1], Z2pair[0], Z2pair[1]

  @methodtools.lru_cache()
  @property
  def GenHp4(self): return sum((p for id, p in self.__gen.daughters), ROOT.TLorentzVector())
  @property
  def GenHMass(self): return self.GenHp4.M()
  @property
  def GenHPt(self): return self.GenHp4.Pt()
  @property
  def GenHRapidity(self): return self.GenHp4.Rapidity()

  @methodtools.lru_cache()
  @property
  def GenZ1(self): return self.Gensortedleptons[:2]
  @methodtools.lru_cache()
  @property
  def GenZ1p4(self): return sum((p for id, p in self.GenZ1), ROOT.TLorentzVector())
  @methodtools.lru_cache()
  @property
  def GenZ2(self): return self.Gensortedleptons[2:]
  @methodtools.lru_cache()
  @property
  def GenZ2p4(self): return sum((p for id, p in self.GenZ2), ROOT.TLorentzVector())
  @methodtools.lru_cache()
  @property
  def GenZ1Mass(self): return self.GenZ1p4.M()
  @property
  def GenZ1Pt(self): return self.GenZ1p4.Pt()
  @property
  def GenZ1Phi(self): return self.GenZ1p4.Phi()
  @methodtools.lru_cache()
  @property
  def GenZ1Flav(self): return np.product([id for id, p in self.GenZ1])
  @methodtools.lru_cache()
  @property
  def GenZ2Mass(self): return self.GenZ2p4.M()
  @property
  def GenZ2Pt(self): return self.GenZ2p4.Pt()
  @property
  def GenZ2Phi(self): return self.GenZ2p4.Phi()
  @methodtools.lru_cache()
  @property
  def GenZ2Flav(self): return np.product([id for id, p in self.GenZ2])
  @methodtools.lru_cache()
  @property
  def Gendecayangles(self):
    angles = self.__gen.computeDecayAngles()
    np.testing.assert_almost_equal(angles.qH, self.GenHMass, decimal=5)
    np.testing.assert_almost_equal(angles.m1, self.GenZ1Mass, decimal=5)
    np.testing.assert_almost_equal(angles.m2, self.GenZ2Mass, decimal=5)
    return angles
  @property
  def Gencosthetastar(self): return self.Gendecayangles.costhetastar
  @property
  def Genhelphi(self): return self.Gendecayangles.Phi
  @property
  def GenhelcosthetaZ1(self): return self.Gendecayangles.costheta1
  @property
  def GenhelcosthetaZ2(self): return self.Gendecayangles.costheta2
  @property
  def GenphistarZ1(self): return self.Gendecayangles.Phi1
  @property
  def GenphistarZ2(self): return 0

  @methodtools.lru_cache()
  @property
  def GenLep1(self): return self.Gensortedleptons[0]
  @methodtools.lru_cache()
  @property
  def GenLep1p4(self): return self.GenLep1[1]
  @property
  def GenLep1Pt(self): return self.GenLep1p4.Pt()
  @property
  def GenLep1Eta(self): return self.GenLep1p4.Eta()
  @property
  def GenLep1Phi(self): return self.GenLep1p4.Phi()
  @property
  def GenLep1Id(self): return self.GenLep1[0]

  @methodtools.lru_cache()
  @property
  def GenLep2(self): return self.Gensortedleptons[1]
  @methodtools.lru_cache()
  @property
  def GenLep2p4(self): return self.GenLep2[1]
  @property
  def GenLep2Pt(self): return self.GenLep2p4.Pt()
  @property
  def GenLep2Eta(self): return self.GenLep2p4.Eta()
  @property
  def GenLep2Phi(self): return self.GenLep2p4.Phi()
  @property
  def GenLep2Id(self): return self.GenLep2[0]

  @methodtools.lru_cache()
  @property
  def GenLep3(self): return self.Gensortedleptons[2]
  @methodtools.lru_cache()
  @property
  def GenLep3p4(self): return self.GenLep3[1]
  @property
  def GenLep3Pt(self): return self.GenLep3p4.Pt()
  @property
  def GenLep3Eta(self): return self.GenLep3p4.Eta()
  @property
  def GenLep3Phi(self): return self.GenLep3p4.Phi()
  @property
  def GenLep3Id(self): return self.GenLep3[0]

  @methodtools.lru_cache()
  @property
  def GenLep4(self): return self.Gensortedleptons[3]
  @methodtools.lru_cache()
  @property
  def GenLep4p4(self): return self.GenLep4[1]
  @property
  def GenLep4Pt(self): return self.GenLep4p4.Pt()
  @property
  def GenLep4Eta(self): return self.GenLep4p4.Eta()
  @property
  def GenLep4Phi(self): return self.GenLep4p4.Phi()
  @property
  def GenLep4Id(self): return self.GenLep4[0]

  zerovector = ROOT.TLorentzVector(0, 0, 0, 0)

  @methodtools.lru_cache()
  @property
  def Gensortedextraleptons(self):
    return sorted(((id, p) for id, p in self.__reco.associated if abs(id) in (11, 13)), key=lambda x: x[1].Pt())
  @methodtools.lru_cache()
  @property
  def GenAssocLep1(self):
    try:
      return self.Gensortedextraleptons[0]
    except IndexError:
      return 0, self.zerovector
  @methodtools.lru_cache()
  @property
  def GenAssocLep1p4(self): return self.GenAssocLep1[1]
  @property
  def GenAssocLep1Pt(self): return self.GenAssocLep1p4.Pt()
  @property
  def GenAssocLep1Eta(self): return self.GenAssocLep1p4.Eta()
  @property
  def GenAssocLep1Phi(self): return self.GenAssocLep1p4.Phi()
  @property
  def GenAssocLep1Id(self): return self.GenAssocLep1[0]
  @methodtools.lru_cache()
  @property
  def GenAssocLep2(self):
    try:
      return self.Gensortedextraleptons[1]
    except IndexError:
      return 0, self.zerovector
  @methodtools.lru_cache()
  @property
  def GenAssocLep2p4(self): return self.GenAssocLep2[1]
  @property
  def GenAssocLep2Pt(self): return self.GenAssocLep2p4.Pt()
  @property
  def GenAssocLep2Eta(self): return self.GenAssocLep2p4.Eta()
  @property
  def GenAssocLep2Phi(self): return self.GenAssocLep2p4.Phi()
  @property
  def GenAssocLep2Id(self): return self.GenAssocLep2[0]

  @methodtools.lru_cache()
  @property
  def LHEDaughters(self): return self.__gen.daughters
  @property
  def LHEDaughterMass(self): return [p.M() for id, p in self.LHEDaughters]
  @property
  def LHEDaughterPt(self): return [p.Pt() for id, p in self.LHEDaughters]
  @property
  def LHEDaughterEta(self): return [p.Eta() for id, p in self.LHEDaughters]
  @property
  def LHEDaughterPhi(self): return [p.Phi() for id, p in self.LHEDaughters]
  @property
  def LHEDaughterId(self): return [id for id, p in self.LHEDaughters]

  @methodtools.lru_cache()
  @property
  def LHEAssociatedParticles(self): return self.__gen.associated
  @property
  def LHEAssociatedParticleMass(self): return [p.M() for id, p in self.LHEAssociatedParticles]
  @property
  def LHEAssociatedParticlePt(self): return [p.Pt() for id, p in self.LHEAssociatedParticles]
  @property
  def LHEAssociatedParticleEta(self): return [p.Eta() for id, p in self.LHEAssociatedParticles]
  @property
  def LHEAssociatedParticlePhi(self): return [p.Phi() for id, p in self.LHEAssociatedParticles]
  @property
  def LHEAssociatedParticleId(self): return [id for id, p in self.LHEAssociatedParticles]

  @methodtools.lru_cache()
  @property
  def LHEMothers(self): return self.__gen.mothers
  @property
  def LHEMotherPz(self): return [p.Pz() for id, p in self.LHEMothers]
  @property
  def LHEMotherE(self): return [p.E() for id, p in self.LHEMothers]
  @property
  def LHEMotherId(self): return [id for id, p in self.LHEMothers]

  @property
  def LHEPDFScale(self): return 1
  @property
  def genExtInfo(self): return -1

  @property
  def xsec(self): return self.__xsec
  @property
  def genxsec(self): return self.__genxsec
  @property
  def genBR(self): return self.__genBR

  #p_GG_SIG_gXg1_1_gXz1_1_JHUGen = p_GG_SIG_gXg2_1_gXz2_1_JHUGen = p_GG_SIG_gXg3_1_gXz3_1_JHUGen = p_GG_SIG_gXg4_1_gXz4_1_JHUGen = p_GG_SIG_gXg1_1_gXz5_1_JHUGen = p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen = p_GG_SIG_gXg1_1_gXz6_1_JHUGen = p_GG_SIG_gXg1_1_gXz7_1_JHUGen = p_GG_SIG_gXg5_1_gXz8_1_JHUGen = p_GG_SIG_gXg5_1_gXz9_1_JHUGen = p_GG_SIG_gXg5_1_gXz10_1_JHUGen = -999

  @classmethod
  def branches(cls, cjlstprocess):
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
    ] or [
      "pConst_GG_SIG_ghg2_1_ghz1_1_JHUGen",
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

  @methodtools.lru_cache()
  @classmethod
  def recoprobabilities(cls):
    from pyFragments.RecoProbabilities import theRecoProbabilities
    probs = [ProbabilityLine.fromline(line, theRecoProbabilities) for line in theRecoProbabilities]
    probs = [prob for prob in probs if prob.name.startswith("GG_SIG") and "MCFM" not in prob.name]
    return probs

  @property
  def gen(self): return self.__gen
  @property
  def reco(self): return self.__reco

for prob in Event.recoprobabilities():
  def f(self, process=prob.process, production=prob.production, me=prob.matrixelement, couplings=prob.couplings):
    reco = self.reco
    reco.setProcess(process, me, production)
    for k, v in couplings.items():
      setattr(reco, k, v)
    return reco.computeP()

  f.__name__ = "p_"+prob.name
  setattr(Event, f.__name__, property(f))
