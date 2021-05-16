#!/usr/bin/env python

import shutil, csv, pprint

def getdataset(identifier):
  for filename in "/afs/cern.ch/work/h/hroskes/public/CJLST/CMSSW_10_2_15/src/ZZAnalysis/AnalysisStep/test/prod/samples_2018_MC.csv", "/afs/cern.ch/work/h/hroskes/public/CJLST/CMSSW_10_2_15/src/ZZAnalysis/AnalysisStep/test/prod/samples_2018_MC_anomalous.csv":
    with open(filename) as f:
      for row in csv.DictReader(f):
        if row["identifier"] == identifier:
          return row["dataset"]
  assert False, identifier

with open("samples_2017_MC.csv") as f, open("samples_2017_MC.csv") as f2, open("samples_2017_MC_new.csv", "w") as newf:
  next(f) #header
  reader = csv.DictReader(f2)
  writer = csv.DictWriter(newf, fieldnames=reader.fieldnames, lineterminator="\n")
  writer.writeheader()
  for line in f:
    print line
    if not line.strip(): newf.write(line); continue
    row = next(reader)
    if not line.startswith("#"):
      row["dataset"] = getdataset(row["identifier"])
    writer.writerow(row)

shutil.move("samples_2017_MC_new.csv", "samples_2017_MC.csv")
