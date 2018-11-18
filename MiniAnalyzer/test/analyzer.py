import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall17MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/94X_mc2017_realistic_v11_ext1-v1/00000/DAAEAB58-9D0E-E811-B5A0-7CD30ACE141A.root'
)
		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('ZZ4lAnalysis.root')
		)

process.softLeptons= cms.EDProducer("SmearLepton",
		src = cms.InputTag("prunedGenParticles")
		)

process.cleanJets= cms.EDProducer("SmearJet",
		genjet = cms.InputTag("slimmedGenJets"),
		src = cms.InputTag("softLeptons")
		)

process.bareZCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
		checkCharge = cms.bool(True),
		cut = cms.string("mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())"),
		decay = cms.string('softLeptons@+ softLeptons@-')
		)

process.ZCand = cms.EDProducer("ZCandidateFiller",
		FSRMode = cms.string('RunII'),
		bestZAmong = cms.string(" mass > 40 && mass < 120"),
		src = cms.InputTag("bareZCand")
		)

process.bareZZCand = cms.EDProducer("CandViewShallowCloneCombiner",
#process.bareZZCand = cms.EDProducer("CandViewCombiner",
		checkCharge = cms.bool(True),
		cut = cms.string('deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02 &&deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(0).eta, daughter(1).daughter(0).phi)>0.02 &&deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).daughter(1).eta, daughter(1).daughter(1).phi)>0.02'),
		#cut = cms.string('1>0'),
		decay = cms.string('ZCand ZCand')
		)

process.ZZCand = cms.EDProducer("ZZCandidateFiller",
		ZRolesByMass = cms.bool(True),
		bestCandAmong = cms.PSet(
			isBestCand = cms.string("daughter(\'Z1\').mass>40 && daughter(\'Z1\').mass<120&&daughter(\'Z2\').mass>4  && daughter(\'Z2\').mass<120&&userFloat(\'mLL4\')>4&&userFloat(\'pt1\')>20 && userFloat(\'pt2\')>10&&mass>70&&userFloat(\'passSmartMLL\')&&daughter(\'Z2\').mass>12")
			),
		bestCandComparator = cms.string('byBestKD'),
		flags = cms.PSet(
			FullSel = cms.string("daughter(\'Z1\').mass>40 && daughter(\'Z1\').mass<120&&daughter(\'Z2\').mass>4  && daughter(\'Z2\').mass<120&&userFloat(\'mLL4\')>4&&userFloat(\'pt1\')>20 && userFloat(\'pt2\')>10&&mass>70&&userFloat(\'passSmartMLL\')&&daughter(\'Z2\').mass>12&&mass>100"),
			FullSel70 = cms.string("daughter(\'Z1\').mass>40 && daughter(\'Z1\').mass<120&&daughter(\'Z2\').mass>4  && daughter(\'Z2\').mass<120&&userFloat(\'mLL4\')>4&&userFloat(\'pt1\')>20 && userFloat(\'pt2\')>10&&mass>70&&userFloat(\'passSmartMLL\')&&daughter(\'Z2\').mass>12"),
			MAllComb = cms.string("userFloat(\'mLL4\')>4"),
			SR = cms.string("daughter(\'Z1\').mass>40 && daughter(\'Z1\').mass<120&&daughter(\'Z2\').mass>4  && daughter(\'Z2\').mass<120&&userFloat(\'mLL4\')>4&&userFloat(\'pt1\')>20 && userFloat(\'pt2\')>10&&mass>70&&userFloat(\'passSmartMLL\')&&daughter(\'Z2\').mass>12"),
			Z2Mass = cms.string("daughter(\'Z2\').mass>4  && daughter(\'Z2\').mass<120")
			),
		recoProbabilities = cms.vstring( 
			),
		src = cms.InputTag("bareZZCand"),
		)




process.ZZTree= cms.EDAnalyzer('MiniAnalyzer',
		genjet = cms.InputTag("cleanJets"),
		pruned = cms.InputTag("prunedGenParticles"),
		VVDecayMode = cms.int32(0),
		VVMode = cms.int32(1),
		setup = cms.int32(LEPTON_SETUP),
		failedTreeLevel = cms.int32(2),
		fileName = cms.untracked.string('candTree'),
		sampleName = cms.string(SAMPLENAME),
    		CandCollection = cms.untracked.string('ZZCand'),
		skipEmptyEvents = cms.bool(True),
    GenBR = cms.double(GENBR),
    GenXSEC = cms.double(GENXSEC),
    		xsec = cms.double(XSEC),
		lheProbabilities = cms.vstring(
			),
		recoProbabilities = cms.vstring( 
			),
		Apply_K_NNLOQCD_ZZGG = cms.int32(globals().get("APPLY_K_NNLOQCD_ZZGG", 0)),
		Apply_K_NNLOQCD_ZZQQB = cms.bool(globals().get("APPLY_K_NNLOQCD_ZZQQB", False)),
		Apply_K_NLOEW_ZZQQB = cms.bool(globals().get("APPLY_K_NLOEW_ZZQQB", False)),

		)


process.Candidates = cms.Path(process.softLeptons+process.cleanJets+process.bareZCand+process.ZCand+process.bareZZCand+process.ZZCand)

process.trees = cms.EndPath(process.ZZTree)
		#process.p = cms.Path(process.ZZTree)






