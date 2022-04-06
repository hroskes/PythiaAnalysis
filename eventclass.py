import contextlib, csv, functools, itertools, methodtools, more_itertools, pathlib2, random, ROOT
from JHUGenMELA.MELA import mela
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

def fspath(path):
  if isinstance(path, str): return path
  return path.__fspath__()

thisfolder = pathlib2.Path(__file__).parent

class ProbabilityLine(object):
  def __init__(self, name, production, process, matrixelement, alias=None, couplings={}, defaultme=None, cluster=None, forceincomingflavors=None, ispm4l=None, supermelasyst=None, ispmavjjtrue=None, ispmavjj=None, addpconst=False, addpaux=False, subtractp=[], maxnumerator=None, maxdenominator=None, isgen=None, nobranch=None, dividep=None):
    self.name = name
    self.alias = alias
    self.production = production
    self.process = process
    self.matrixelement = matrixelement
    self.couplings = couplings
    self.defaultme = defaultme
    self.cluster = cluster
    self.forceincomingflavors = forceincomingflavors
    self.ispm4l = ispm4l
    self.supermelasyst = supermelasyst
    self.ispmavjjtrue = ispmavjjtrue
    self.ispmavjj = ispmavjj
    self.addpconst = addpconst
    self.addpaux = addpaux
    self.subtractp = subtractp
    self.maxnumerator = maxnumerator
    self.maxdenominator = maxdenominator
    self.isgen = isgen
    self.nobranch = nobranch
    self.dividep = dividep

  @classmethod
  def kwargsfromline(cls, line, alllines):
    kwargs = {k.lower(): v for k, v in (_.split(":") for _ in line.split(" "))}
    if "copy" in kwargs:
      copyline, = {line for line in alllines if "Name:"+kwargs["copy"]+" " in line}
      copykwargs = cls.kwargsfromline(copyline, alllines)
      copykwargs.update(kwargs)
      kwargs = copykwargs
      del kwargs["copy"]

    for thing in "production", "process", "matrixelement", "supermelasyst":
      if isinstance(kwargs.get(thing, None), str):
        kwargs[thing] = getattr(TVar, kwargs[thing])

    for thing in "options", "couplings":
      if thing in kwargs:
        kwargs[thing] = {k: v for k, v in (_.split("=") for _ in kwargs[thing].split(";"))}

    for name, val in kwargs.pop("options", {}).items():
      name = name.lower()
      assert name not in kwargs, name
      kwargs[name] = val

    for name, val in kwargs.get("couplings", {}).items():
      if name == "separateWWZZcouplings":
        del kwargs["couplings"][name]
        kwargs["couplings"]["differentiate_HWW_HZZ"] = bool(int(val))
      else:
        real, imag = val.split(",")
        real = float(real)
        imag = float(imag)
        kwargs["couplings"][name] = real + imag*1j

    if kwargs.get("alias", None) == "<Name>":
      kwargs["alias"] = kwargs["name"]

    if "defaultme" in kwargs:
      kwargs["defaultme"] = float(kwargs["defaultme"])

    for thing in "addpaux", "addpconst", "ispm4l", "isgen", "nobranch":
      if thing in kwargs:
        kwargs[thing] = bool(int(kwargs[thing]))

    for thing in "subtractp",:
      if thing in kwargs:
        kwargs[thing] = ["p_"+_ for _ in kwargs[thing].split(",")]

    return kwargs

  @classmethod
  def fromline(cls, line, alllines):
    return cls(**cls.kwargsfromline(line, alllines))

class LHEEvent_gen(LHEEvent):
  @classmethod
  def extracteventparticles(cls, lines, isgen):
    assert isgen
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

    if not daughters: daughters = dummydaughters() #silence TUtil::ConvertVectorFormat warning

    if not isgen: mothers = None
    return daughters, associated, mothers

ptsigmas = {
  11: 2.399/6,
  13: 2.169/6,
  0: 18./6,
}

minpts = {
  11: 7,
  13: 5,
  0: 30,
}

maxetas = {
  11: 2.5,
  13: 2.4,
  0: 4.7,
}

for dct in ptsigmas, minpts, maxetas:
  for _ in range(6): dct[_] = dct[0]
  for _ in list(dct.keys()): dct[-_] = dct[_]
del dct, _

def dummydaughters():
  global __dummydaughters
  try:
    return __dummydaughters
  except NameError:
    __dummydaughters = [
      mela.SimpleParticle_t(15, 0, 0, 0, 0),
      mela.SimpleParticle_t(-15, 0, 0, 0, 0),
    ]
  return dummydaughters()


class LHEEvent_reco(LHEEvent):
  #smearing is not actually implemented yet
  @classmethod
  def extracteventparticles(cls, lines, isgen):
    assert not isgen
    random.seed(hash(lines[0]))
    daughters, mothers, associated, leptons, taus, jets = [], [], [], {11: [], -11: [], 13: [], -13: []}, [], []
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
      elif status == 1 and (0 <= abs(id) <= 6 or 11 <= abs(id) <= 16 or abs(id) in (21, 22)):
        if abs(id) in (11, 13):
          leptons[id].append(line)
        elif abs(id) in (1, 2, 3, 4, 5, 21):
          jets.append(line)
        elif abs(id) in (12, 14, 16):
          pass
        elif abs(id) == 15:
          taus.append(line)
        else:
          assert False, id

    for lst in list(leptons.values()) + [jets]:
      lst[:] = [mela.SimpleParticle_t(_) for _ in lst]
      for i, (id, p) in reversed(list(enumerate(lst[:]))):
        ptsigma = ptsigmas[id]
        p.SetPtEtaPhiM(random.gauss(p.Pt(), ptsigma), p.Eta(), p.Phi(), p.M())
        if p.Pt() < minpts[id] or abs(p.Eta()) > maxetas[id]:
          del lst[i]

    associated += jets

    leptons[1] = leptons[11] + leptons[13]
    leptons[-1] = leptons[-11] + leptons[-13]

    if min(len(leptons[11]), len(leptons[-11])) + min(len(leptons[13]), len(leptons[-13])) < 2:
      daughters = []
    elif len(leptons[1]) == 2 == len(leptons[-1]):
      daughters = leptons[1] + leptons[-1]
    else:
      combinations = [[lep1, lep2, lep3, lep4] for (lep1, lep3), (lep2, lep4) in itertools.product(itertools.combinations(leptons[1], 2), itertools.combinations(leptons[-1], 2)) if lep1.first+lep2.first==0 and lep3.first+lep4.first==0]
      daughters = min(combinations, key=lambda x: abs(sum((lep[1] for lep in x), ROOT.TLorentzVector()).M()-125))
      associated += [_ for _ in leptons[1]+leptons[-1] if _ not in daughters]

    if len(daughters) == 4:
      possibleZs = [sorted(pair, key=lambda x: x[0]) for pair in itertools.combinations(daughters, 2) if abs(pair[0][0]) in {11, 13} and sum(p[0] for p in pair) == 0]
      Z1pair = min(possibleZs, key=lambda x: abs(sum((p for id, p in x), ROOT.TLorentzVector()).M()-91.1876))
      l1p, l1m = Z1pair
      Z2pair, = {(l2p, l2m) for l2p, l2m in possibleZs if l2p is not l1p and l2m is not l1m}
      daughters = [Z1pair[0], Z1pair[1], Z2pair[0], Z2pair[1]]

    if not daughters: daughters = dummydaughters() #silence TUtil::ConvertVectorFormat warning

    if not isgen: mothers = None
    return daughters, associated, mothers

class LHEFile_reco(LHEFileBase):
  lheeventclass = LHEEvent_reco
  @property
  def finalstateparticles(self):
    return itertools.chain(self.daughters, self.associated)

class LHEFile_gen(LHEFileBase):
  lheeventclass = LHEEvent_gen

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
  def PFMET(self):
    return -999 #below line doesn't work because we ignore taus
    return -sum((p for id, p in self.__reco.finalstateparticles), ROOT.TLorentzVector()).Pt()
  @methodtools.lru_cache()
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
    return self.__reco.daughters
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
    np.testing.assert_almost_equal(angles.qH, self.ZZMass, decimal=2)
    np.testing.assert_almost_equal(angles.m1, self.Z1Mass, decimal=2)
    np.testing.assert_almost_equal(angles.m2, self.Z2Mass, decimal=2)
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
  @methodtools.lru_cache()
  @property
  def VBFangles(self):
    if self.nCleanedJets < 2: return None
    angles = self.__reco.computeVBFAngles()
    return angles
  @property
  def costhetastarVBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.costhetastar
  @property
  def costheta1VBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.costheta1
  @property
  def costheta2VBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.costheta2
  @property
  def PhiVBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.Phi
  @property
  def Phi1VBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.Phi1
  @property
  def Q2V1VBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.Q2V1
  @property
  def Q2V2VBF(self):
    angles = self.VBFangles
    if angles is None: return -999
    return angles.Q2V2
  @methodtools.lru_cache()
  @property
  def VHangles(self):
    if self.nCleanedJets < 2: return None
    angles = self.__reco.computeVHAngles(TVar.Had_ZH)
    return angles
  @property
  def costhetastarVH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.costhetastar
  @property
  def costheta1VH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.costheta1
  @property
  def costheta2VH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.costheta2
  @property
  def PhiVH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.Phi
  @property
  def Phi1VH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.Phi1
  @property
  def mVstarVH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.mVstar
  @property
  def mVVH(self):
    angles = self.VHangles
    if angles is None: return -999
    return angles.mV
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
    }.get(self.GenZ1Flav*self.GenZ2Flav, -1)
  @property
  def genProcessId(self): return 0
  @property
  def genHEPMCweight(self): return self.__gen.weight
  @property
  def KFactor_QCD_ggZZ_Nominal(self): return 1.25752

  @methodtools.lru_cache()
  @property
  def Gensortedleptons(self):
    possibleZs = [sorted(pair, key=lambda x: x[0]) for pair in itertools.combinations(self.__gen.daughters, 2) if sum(p[0] for p in pair) == 0]
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

  @classmethod
  def branches(cls, cjlstprocess):
    result = [
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
      "costhetastarVBF",
      "costheta1VBF",
      "costheta2VBF",
      "PhiVBF",
      "Phi1VBF",
      "Q2V1VBF",
      "Q2V2VBF",
      "costhetastarVH",
      "costheta1VH",
      "costheta2VH",
      "PhiVH",
      "Phi1VH",
      "mVstarVH",
      "mVVH",
    ]
    if cjlstprocess.startswith("ggH"): result += [
      "KFactor_QCD_ggZZ_Nominal",
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
    else: result += [
      "p_Gen_Dec_SIG_ghz1_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz2_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghz4_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza2_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_ghza4_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_gha2_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1_1_gha4_i_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_ghz4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz2_1_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghz4_1_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_Dec_SIG_ghza1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghza1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza2_1_ghza4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza2_1_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghza2_1_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_ghza4_1_gha2_1_JHUGen",
      "p_Gen_Dec_SIG_ghza4_1_gha4_1_JHUGen",
      "p_Gen_Dec_SIG_gha2_1_gha4_1_JHUGen",
    ]
    if cjlstprocess.startswith("ZH"): result += [
      "p_Gen_ZH_SIG_ghz1_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz1prime2_1E4i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz2_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghz4_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza1prime2_1E4i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza2_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_ghza4_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_gha2_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1_1_gha4_i_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_ghz4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz2_1_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghz4_1_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_ZH_SIG_ghza1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghza1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza2_1_ghza4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza2_1_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghza2_1_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_ghza4_1_gha2_1_JHUGen",
      "p_Gen_ZH_SIG_ghza4_1_gha4_1_JHUGen",
      "p_Gen_ZH_SIG_gha2_1_gha4_1_JHUGen",
    ]
    if cjlstprocess.startswith("WH"): result += [
      "p_Gen_WH_SIG_ghw1_1_JHUGen",
      "p_Gen_WH_SIG_ghw1prime2_1E4_JHUGen",
      "p_Gen_WH_SIG_ghw2_1_JHUGen",
      "p_Gen_WH_SIG_ghw4_1_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw1prime2_1E4i_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw2_1_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw2_i_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw4_1_JHUGen",
      "p_Gen_WH_SIG_ghw1_1_ghw4_i_JHUGen",
      "p_Gen_WH_SIG_ghw1prime2_1E4_ghw2_1_JHUGen",
      "p_Gen_WH_SIG_ghw1prime2_1E4_ghw4_1_JHUGen",
      "p_Gen_WH_SIG_ghw2_1_ghw4_1_JHUGen",
    ]
    if cjlstprocess.startswith("VBF") and not cjlstprocess.endswith("0L1Zg_M125"): result += [
      "p_Gen_VBF_SIG_ghv1_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghv2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_JHUGen",
      "p_Gen_VBF_SIG_ghv4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghw1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghw2_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghv1_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghw2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghw2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_ghw4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw2_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghw4_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_gha2_1_gha4_1_JHUGen",
    ]
    if cjlstprocess.startswith("VBF") and cjlstprocess.endswith("0L1Zg_M125"): result += [
      "p_Gen_VBF_SIG_ghz1_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghz1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghz2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghz4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz2_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza1prime2_1E4_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghz4_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_ghza2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza1prime2_1E4_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_ghza4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza2_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_gha2_1_JHUGen",
      "p_Gen_VBF_SIG_ghza4_1_gha4_1_JHUGen",
      "p_Gen_VBF_SIG_gha2_1_gha4_1_JHUGen",
    ]
    if cjlstprocess.startswith("HJJ"): result += [
      "p_Gen_HJJ_SIG_ghg2_1_JHUGen",
      "p_Gen_HJJ_SIG_ghg2_1_ghg4_1_JHUGen",
      "p_Gen_HJJ_SIG_ghg4_1_JHUGen",
      "p_Gen_HJJ_SIG_ghg2_1_ghg4_i_JHUGen",
      "p_Gen_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal",
      "p_Gen_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal",
      "p_Gen_JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal",
    ]
    return result

  p_GG_BSI_kappaTopBot_1_ghz1_i_MCFM = -999
  p_ttHUndecayed_SIG_kappa_1_JHUGen_JECNominal = -999
  p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECNominal = -999
  p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECNominal = -999
  p_bbH_SIG_kappa_1_JHUGen_JECNominal = -999
  p_ttHUndecayed_SIG_kappa_1_JHUGen_JECUp = -999
  p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECUp = -999
  p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECUp = -999
  p_bbH_SIG_kappa_1_JHUGen_JECUp = -999
  p_ttHUndecayed_SIG_kappa_1_JHUGen_JECDn = -999
  p_ttHUndecayed_SIG_kappatilde_1_JHUGen_JECDn = -999
  p_ttHUndecayed_SIG_kappa_1_kappatilde_1_JHUGen_JECDn = -999
  p_bbH_SIG_kappa_1_JHUGen_JECDn = -999

  p_Gen_JJQCD_InitialQQ_SIG_ghg2_1_JHUGen_JECNominal = -999
  p_Gen_JJQCD_InitialQQ_SIG_ghg4_1_JHUGen_JECNominal = -999
  p_Gen_JJQCD_InitialQQ_SIG_ghg2_1_ghg4_1_JHUGen_JECNominal = -999

  @methodtools.lru_cache()
  @classmethod
  def recoprobabilities(cls):
    from pyFragments.RecoProbabilities import theRecoProbabilities
    probs = [ProbabilityLine.fromline(line, theRecoProbabilities) for line in theRecoProbabilities]
    probs = [prob for prob in probs]
    return probs

  @property
  def gen(self): return self.__gen
  @property
  def reco(self): return self.__reco

_2jetproductions = {TVar.Had_WH, TVar.Had_ZH, TVar.JJVBF, TVar.JJEW, TVar.JJQCD, TVar.Had_WH_S, TVar.Had_ZH_S, TVar.JJVBF_S, TVar.JJEW_S, TVar.JJQCD_S, TVar.Had_WH_TU, TVar.Had_ZH_TU, TVar.JJVBF_TU, TVar.JJEW_TU, TVar.JJQCD_TU}
_1jetproductions = {TVar.JQCD}
_lepVHproductions = {TVar.Lep_WH, TVar.Lep_ZH, TVar.Lep_WH_S, TVar.Lep_ZH_S, TVar.Lep_WH_TU, TVar.Lep_ZH_TU}
_prodproductions = _2jetproductions | _1jetproductions | _lepVHproductions | {TVar.GammaH, TVar.ttH, TVar.bbH}
for prob in Event.recoprobabilities():
  def f(self, process=prob.process, production=prob.production, me=prob.matrixelement, couplings=prob.couplings, addpconst=prob.addpconst, addpaux=prob.addpaux, ispm4l=prob.ispm4l, supermelasyst=prob.supermelasyst, isprod=prob.production in _prodproductions and prob.matrixelement != TVar.MCFM, isproddec=prob.production in _prodproductions and prob.matrixelement == TVar.MCFM, need2jets="J2JEC" in str(prob.cluster), need1jet="J1JEC" in str(prob.cluster), islepVH=prob.cluster in ("LepZH", "LepWH"), subtractp=prob.subtractp):
    reco = self.reco
    if islepVH: return -999, -999, -999
    if need2jets and self.nCleanedJets < 2: return -999, -999, -999
    if need1jet and self.nCleanedJets != 1: return -999, -999, -999
    reco.setProcess(process, me, production)
    for k, v in couplings.items():
      setattr(reco, k, v)
    if ispm4l:
      prob = reco.computePM4l(supermelasyst)
    elif isproddec:
      prob = reco.computeProdDecP()
    elif isprod:
      prob = reco.computeProdP()
    else:
      prob = reco.computeP()
    for name in subtractp:
      prob -= getattr(self, name)
    pconst = reco.getConstant() if addpconst else None
    paux = reco.getPAux() if addpaux else None
    return prob, pconst, paux

  f.__name__ = prob.name
  if not hasattr(Event, f.__name__): setattr(Event, f.__name__, methodtools.lru_cache()(property(f)))

  def p(self, allname=f.__name__): return getattr(self, allname)[0]
  p.__name__ = "p_"+prob.name
  if not hasattr(Event, p.__name__): setattr(Event, p.__name__, methodtools.lru_cache()(property(p)))

  if prob.addpconst:
    def pconst(self, allname=f.__name__): return getattr(self, allname)[1]
    pconst.__name__ = "pConst_"+prob.name
    if not hasattr(Event, pconst.__name__): setattr(Event, pconst.__name__, methodtools.lru_cache()(property(pconst)))

  if prob.addpaux:
    def paux(self, allname=f.__name__): return getattr(self, allname)[2]
    paux.__name__ = "pAux_"+prob.name
    if not hasattr(Event, paux.__name__): setattr(Event, paux.__name__, methodtools.lru_cache()(property(paux)))

@methodtools.lru_cache()
def CJLSTrow(cjlstprocess):
  with open(fspath(thisfolder/"samples_2018_MC.csv")) as f:
    reader = csv.DictReader(f)
    for row in reader:
      if row["identifier"].strip() == cjlstprocess:
        return row
    else:
      assert False, cjlstprocess

@methodtools.lru_cache()
def eventsubclass(cjlstprocess):
  class EventSubClass(Event):
    @methodtools.lru_cache()
    @classmethod
    def genprobabilities(cls):
      fragment = CJLSTrow(cjlstprocess)["::pyFragments"].split(";")[1].replace(".py", "")
      theLHEProbabilities = __import__("pyFragments."+fragment, fromlist="theLHEProbabilities").theLHEProbabilities
      probs = [ProbabilityLine.fromline(line, theLHEProbabilities) for line in theLHEProbabilities]
      probs = [prob for prob in probs]
      return probs

  for prob in EventSubClass.genprobabilities():
    def f(self, process=prob.process, production=prob.production, me=prob.matrixelement, couplings=prob.couplings, isprod=prob.production in _prodproductions and prob.matrixelement != TVar.MCFM, isproddec=prob.production in _prodproductions and prob.matrixelement == TVar.MCFM, dividep=prob.dividep):
      gen = self.gen
      if production == TVar.Had_ZH and abs(gen.associated[0][0]) >= 11:
        production = TVar.Lep_ZH
      if production == TVar.Had_WH and abs(gen.associated[0][0]) >= 11:
        production = TVar.Lep_WH
      gen.setProcess(process, me, production)
      for k, v in couplings.items():
        setattr(gen, k, v)
      if isproddec:
        prob = gen.computeProdDecP()
      elif isprod:
        prob = gen.computeProdP()
      else:
        prob = gen.computeP()
      if dividep is not None:
        try:
          prob /= getattr(self, "p_Gen_"+dividep)
        except ZeroDivisionError:
          prob = 0
      return prob

    f.__name__ = "p_Gen_"+prob.name
    if not hasattr(Event, f.__name__): setattr(Event, f.__name__, methodtools.lru_cache()(property(f)))

  EventSubClass.__name__ = "Event"+cjlstprocess
  return EventSubClass
