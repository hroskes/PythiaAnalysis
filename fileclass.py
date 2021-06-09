import abc, argparse, contextlib2, csv, methodtools, more_itertools, pathlib2

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--lhefilename", type=pathlib2.Path, required=True)
  p.add_argument("--cjlstprocess", required=True)
  p.add_argument("--outfilename", type=pathlib2.Path, required=True)
  p.add_argument("--overwrite", action="store_true")
  p.add_argument("--cjlstfolder", type=pathlib2.Path, default=pathlib2.Path("/work-zfs/lhc/GENtrees/210601_2018MC_photons"))
  args = p.parse_args()

thisfolder = pathlib2.Path(__file__).parent

from eventclass import Event, LHEFile_Hwithdecay, LHEFile_Hwithdecay_smear

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

class TFile(object):
  def __init__(self, filename, *args, **kwargs):
    self.__filename = filename
    self.__args = args
    self.__deleteifbad = kwargs.pop("deleteifbad", False)
    self.__entered = False
    assert not kwargs, kwargs
  def __enter__(self):
    import ROOT
    self.__write = False
    self.__bkpdirectory = ROOT.gDirectory.GetDirectory(ROOT.gDirectory.GetPath())
    self.__f = ROOT.TFile.Open(fspath(self.__filename), *self.__args)
    self.__entered = True
    if not self.__f:
      raise IOError(str(self.__filename)+" is a null pointer, see above for details.")
    if self.__f.IsZombie():
      self.__exit__()
      raise IOError(str(self.__filename)+" is a zombie, see above for details.")

    try:
      openoption = self.__args[0].upper()
    except IndexError:
      openoption = ""

    self.__write = {
      "": False,
      "READ": False,
      "NEW": True,
      "CREATE": True,
      "RECREATE": True,
      "UPDATE": True,
    }[openoption]

    return self.__f

  def __exit__(self, *errorstuff):
    try:
      if self.__write and (not any(errorstuff) or not self.__deleteifbad):
        self.__f.Write()
      self.__f.Close()
      self.__bkpdirectory.cd()
    except:
      if self.__write and self.__deleteifbad:
        self.__filename.unlink()
      raise

    if self.__write and self.__deleteifbad and any(errorstuff):
      self.__filename.unlink()

class Branch(object):
  __metaclass__ = abc.ABCMeta
  def __init__(self, name, typ):
    self.__name = name
    self.__type = typ

  @abc.abstractproperty
  def thingforsetbranchaddress(self): pass

  def attachtotree(self, t):
    t.SetBranchAddress(self.__name, self.thingforsetbranchaddress)
    target = type(getattr(t, self.__name))
    if target != self.__type:
      if target == int and self.__type == np.short:
        pass #ok
      elif target == float and self.__type == np.float32:
        pass #ok
      else:
        raise ValueError("Wrong type for {}: {}, should be {}".format(self.__name, self.__type, type(getattr(t, self.__name))))

  @abc.abstractmethod
  def setthingfortree(self, event): pass

  def setbranchvalue(self, event):
    newvalue = self.setthingfortree(event)
    self.lastsetbranchvalue = newvalue

  @property
  def name(self): return self.__name

class DummyBranch(Branch):
  """
  Branch that doesn't actually go into the tree but can be used to calculate things
  """
  def __init__(self, name):
    return super(DummyBranch, self).__init__(name, "DUMMY")
  @property
  def thingforsetbranchaddress(self): return None
  def attachtotree(self, t): pass

class NormalBranch(Branch):
  def __init__(self, name, typ):
    super(NormalBranch, self).__init__(name, typ)
    self.__array = np.array([0], typ)
  @property
  def thingforsetbranchaddress(self): return self.__array
  def setthingfortree(self, event):
    self.__array[0] = getattr(event, self.name)
    return self.__array[0]

class VectorBranch(Branch):
  def __init__(self, name, typ):
    import ROOT
    vectortyp = ROOT.vector(typ)
    super(VectorBranch, self).__init__(name, vectortyp)
    self.__vector = vectortyp(0)
  @property
  def thingforsetbranchaddress(self): return self.__vector
  def setthingfortree(self, event):
    self.__vector.clear()
    for value in getattr(event, self.name):
      self.__vector.push_back(value)
    return self.__vector

class CJLHEFile(contextlib2.ExitStack):
  def __init__(self, lhefilename, cjlstfolder, cjlstprocess, outfilename, overwrite=False):
    super(CJLHEFile, self).__init__()
    self.__lhefilename = lhefilename
    self.cjlstfolder = cjlstfolder
    self.cjlstprocess = cjlstprocess
    self.__outfilename = outfilename
    if overwrite:
      try:
        outfilename.unlink()
      except OSError:
        pass

    self.__branches = self.branches()

    branchnames = {branch.name for branch in self.__branches}
    targetbranchnames = set(Event.branches())
    assert branchnames == targetbranchnames, branchnames ^ targetbranchnames
    bad = {name for name in branchnames if not hasattr(Event, name)}
    assert not bad, bad

  def branches(self):
    float32 = np.float32
    return [
      NormalBranch("RunNumber", int),
      NormalBranch("LumiNumber", int),
      NormalBranch("EventNumber", long),
      NormalBranch("NRecoMu", int),
      NormalBranch("NRecoEle", int),
      NormalBranch("Nvtx", int),
      NormalBranch("NObsInt", int),
      NormalBranch("NTrueInt", float),
      NormalBranch("PFMET", float),
      NormalBranch("nCleanedJets", int),
      NormalBranch("nCleanedJetsPt30", int),
      NormalBranch("trigWord", int),
      NormalBranch("evtPassMETFilter", int),
      NormalBranch("CRflag", int),
      NormalBranch("ZZMass", float32),
      NormalBranch("ZZsel", int),
      NormalBranch("ZZPt", float32),
      NormalBranch("ZZEta", float32),
      NormalBranch("ZZPhi", float32),
      NormalBranch("Z1Flav", int),
      NormalBranch("Z1Mass", float32),
      NormalBranch("Z1Pt", float32),
      NormalBranch("Z2Flav", int),
      NormalBranch("Z2Mass", float32),
      NormalBranch("Z2Pt", float32),
      NormalBranch("costhetastar", float32),
      NormalBranch("helphi", float32),
      NormalBranch("helcosthetaZ1", float32),
      NormalBranch("helcosthetaZ2", float32),
      NormalBranch("phistarZ1", float32),
      NormalBranch("phistarZ2", float32),
      NormalBranch("xi", float32),
      NormalBranch("xistar", float32),
      VectorBranch("LepPt", float),
      VectorBranch("LepEta", float),
      VectorBranch("LepPhi", float),
      VectorBranch("LepLepId", "short"),
      VectorBranch("JetPt", float),
      VectorBranch("JetEta", float),
      VectorBranch("JetPhi", float),
      VectorBranch("JetMass", float),
      NormalBranch("DiJetMass", float32),
      NormalBranch("DiJetDEta", float32),
      NormalBranch("nExtraLep", int),
      NormalBranch("nExtraZ", int),
      VectorBranch("ExtraLepPt", float),
      VectorBranch("ExtraLepEta", float),
      VectorBranch("ExtraLepPhi", float),
      VectorBranch("ExtraLepLepId", "short"),
      NormalBranch("genFinalState", int),
      NormalBranch("trigEffWeight", float32),
      NormalBranch("dataMCWeight", float32),
      NormalBranch("KFactor_QCD_ggZZ_Nominal", float32),
      NormalBranch("genHEPMCweight_NNLO", float32),
      NormalBranch("HqTMCweight", float32),
      NormalBranch("ZXFakeweight", float32),
      NormalBranch("PUWeight", float32),
      NormalBranch("genHEPMCweight", float32),
      NormalBranch("overallEventWeight", float32),
      NormalBranch("genHEPMCweight_POWHEGonly", float32),
      NormalBranch("genProcessId", int),
      NormalBranch("xsec", float32),
      NormalBranch("genxsec", float32),
      NormalBranch("genBR", float32),
      NormalBranch("genExtInfo", int),
      NormalBranch("GenHMass", float32),
      NormalBranch("GenHPt", float32),
      NormalBranch("GenHRapidity", float32),
      NormalBranch("GenZ1Mass", float32),
      NormalBranch("GenZ1Pt", float32),
      NormalBranch("GenZ1Phi", float32),
      NormalBranch("GenZ1Flav", float32),
      NormalBranch("GenZ2Mass", float32),
      NormalBranch("GenZ2Pt", float32),
      NormalBranch("GenZ2Phi", float32),
      NormalBranch("GenZ2Flav", float32),
      NormalBranch("GenLep1Pt", float32),
      NormalBranch("GenLep1Eta", float32),
      NormalBranch("GenLep1Phi", float32),
      NormalBranch("GenLep1Id", int),
      NormalBranch("GenLep2Pt", float32),
      NormalBranch("GenLep2Eta", float32),
      NormalBranch("GenLep2Phi", float32),
      NormalBranch("GenLep2Id", int),
      NormalBranch("GenLep3Pt", float32),
      NormalBranch("GenLep3Eta", float32),
      NormalBranch("GenLep3Phi", float32),
      NormalBranch("GenLep3Id", int),
      NormalBranch("GenLep4Pt", float32),
      NormalBranch("GenLep4Eta", float32),
      NormalBranch("GenLep4Phi", float32),
      NormalBranch("GenLep4Id", int),
      NormalBranch("GenAssocLep1Pt", float32),
      NormalBranch("GenAssocLep1Eta", float32),
      NormalBranch("GenAssocLep1Phi", float32),
      NormalBranch("GenAssocLep1Id", int),
      NormalBranch("GenAssocLep2Pt", float32),
      NormalBranch("GenAssocLep2Eta", float32),
      NormalBranch("GenAssocLep2Phi", float32),
      NormalBranch("GenAssocLep2Id", int),
      VectorBranch("LHEMotherPz", float),
      VectorBranch("LHEMotherE", float),
      VectorBranch("LHEMotherId", "short"),
      VectorBranch("LHEDaughterPt", float),
      VectorBranch("LHEDaughterEta", float),
      VectorBranch("LHEDaughterPhi", float),
      VectorBranch("LHEDaughterMass", float),
      VectorBranch("LHEDaughterId", "short"),
      VectorBranch("LHEAssociatedParticlePt", float),
      VectorBranch("LHEAssociatedParticleEta", float),
      VectorBranch("LHEAssociatedParticlePhi", float),
      VectorBranch("LHEAssociatedParticleMass", float),
      VectorBranch("LHEAssociatedParticleId", "short"),
      NormalBranch("LHEPDFScale", float32),
    ]

  @property
  def cjlstfilename(self):
    return self.cjlstfolder/self.cjlstprocess/"ZZ4lAnalysis.root"

  def __enter__(self):
    super(CJLHEFile, self).__enter__()
    self.__outfile = self.enter_context(TFile(self.__outfilename, "CREATE", deleteifbad=True))
    ZZTree = self.__outfile.mkdir("ZZTree")
    with TFile(self.cjlstfilename) as CJLSTfile:
      ZZTree.cd()
      t = CJLSTfile.Get("ZZTree/candTree")
      self.__t = t.CloneTree(0, "fast")
    for branch in self.__branches:
      branch.attachtotree(self.__t)

    self.__gen = self.enter_context(LHEFile_Hwithdecay(fspath(self.__lhefilename), isgen=True))
    self.__reco = self.enter_context(LHEFile_Hwithdecay_smear(fspath(self.__lhefilename), isgen=False))

    self.__nentries = 0
    with open(fspath(self.__lhefilename)) as f:
      for line in f:
        if line.strip() == "<event>":
          self.__nentries += 1

    return self

  @methodtools.lru_cache()
  @property
  def xsecs(self):
    with open(fspath(thisfolder/"samples_2018_MC.csv")) as f:
      reader = csv.DictReader(f)
      for row in reader:
        if row["identifier"].strip() == self.cjlstprocess:
          break
      else:
        assert False, self.cjlstprocess

      xsec = float(row["crossSection=-1"])*float(row["BR=1"])
      variables = {k: v for k, v in (_.split("=") for _ in row["::variables"].split(";"))}
      genxsec = variables["GENXSEC"]
      genBR = variables["GENBR"]
      return xsec, genxsec, genBR
  @methodtools.lru_cache()
  @property
  def xsec(self): return self.xsecs[0]
  @methodtools.lru_cache()
  @property
  def genxsec(self): return self.xsecs[1]
  @methodtools.lru_cache()
  @property
  def genBR(self): return self.xsecs[2]
      

  def __iter__(self):
    for i, (gen, reco) in enumerate(more_itertools.more.zip_longest(self.__gen, self.__reco)):
      yield Event(i=i, gen=gen, reco=reco, xsec=self.xsec, genxsec=self.genxsec, genBR=self.genBR)

  def run(self):
    with self:
      for i, event in enumerate(self):
        for branch in self.__branches:
          branch.setbranchvalue(event)
        self.__t.Fill()

        if (i+1)%10000 == 0 or (i+1) == self.__nentries:
          print i+1, "/", self.__nentries
          for branch in self.__branches:
            assert getattr(self.__t, branch.name) == branch.lastsetbranchvalue, (branch.name, getattr(self.__t, branch.name), branch.lastsetbranchvalue)

def main(**kwargs):
  CJLHEFile(**kwargs).run()

if __name__ == "__main__":
  main(**args.__dict__)
