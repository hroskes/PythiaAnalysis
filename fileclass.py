import abc, argparse, contextlib2, more_itertools, pathlib2

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--lhefilename", type=pathlib2.Path, required=True)
  p.add_argument("--cjlstfilename", type=pathlib2.Path, required=True)
  p.add_argument("--outfilename", type=pathlib2.Path, required=True)
  p.add_argument("--overwrite", action="store_true")
  args = p.parse_args()

from eventclass import Event, LHEFile_Hwithdecay_smear

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

class JetVectorBranch(Branch):
  def __init__(self, name, typ):
    vectortyp = ROOT.vector(typ)
    super(JetVectorBranch, self).__init__(name, vectortyp)
    self.__vector = vectortyp(0)
  @property
  def thingforsetbranchaddress(self): return self.__vector
  def setthingfortree(self, event):
    self.__vector.clear()
    for value in getattr(event, self.name):
      self.__vector.push_back(value)
    return self.__vector

def branches():
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
  ]

class CJLHEFile(contextlib2.ExitStack):
  def __init__(self, lhefilename, cjlstfilename, outfilename, overwrite=False):
    super(CJLHEFile, self).__init__()
    self.__lhefilename = lhefilename
    self.__cjlstfilename = cjlstfilename
    self.__outfilename = outfilename
    if overwrite:
      try:
        outfilename.unlink()
      except OSError:
        pass

    self.__branches = branches()

    branchnames = {branch.name for branch in self.__branches}
    targetbranchnames = set(Event.branches())
    assert branchnames == targetbranchnames, branchnames ^ targetbranchnames

  def __enter__(self):
    super(CJLHEFile, self).__enter__()
    self.__outfile = self.enter_context(TFile(self.__outfilename, "CREATE", deleteifbad=True))
    with TFile(self.__cjlstfilename) as CJLSTfile:
      self.__outfile.cd()
      t = CJLSTfile.Get("ZZTree/candTree")
      self.__t = t.CloneTree(0, "fast")
    for branch in self.__branches:
      branch.attachtotree(self.__t)

    self.__gen = self.enter_context(LHEFile_Hwithdecay_smear(fspath(self.__lhefilename), isgen=True))
    self.__reco = self.enter_context(LHEFile_Hwithdecay_smear(fspath(self.__lhefilename), isgen=False))

    self.__nentries = 0
    with open(fspath(self.__lhefilename)) as f:
      for line in f:
        if line.strip() == "<event>":
          self.__nentries += 1

    return self

  def __iter__(self):
    for i, (gen, reco) in enumerate(more_itertools.more.zip_longest(self.__gen, self.__reco)):
      yield Event(i=i, gen=gen, reco=reco)

  def run(self):
    with self:
      for i, event in enumerate(self):
        for branch in self.__branches:
          branch.setbranchvalue(event)
        self.__t.Fill()

        if (i+1)%10000 == 0 or (i+1) == self.__nentries:
          print i+1, "/", self.__nentries

def main(**kwargs):
  CJLHEFile(**kwargs).run()

if __name__ == "__main__":
  main(**args.__dict__)
