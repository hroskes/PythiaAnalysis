import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			#        '/store/mc/RunIISummer15GS/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8_V2/GEN-SIM/MCRUN2_71_V1-v1/110000/02995DC1-B191-E711-AB1E-001EC9AF64B6.root '
			#'/store/mc/RunIIFall17MiniAOD/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/60000/421F9A3D-A1D7-E711-8629-0CC47A13CCFC.root'
			#'/store/mc/RunIISummer16MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJ_NNLOPS_JHUgenV702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/44EEBBEB-E0DE-E611-A786-0CC47A07F9FE.root'
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
		bestZAmong = cms.string(" mass > 12 && mass < 120"),
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
		recoProbabilities = cms.vstring( ('Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen Options:AddPConst=1',
			'Name:GG_SIG_ghg2_1_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:ZZGG MatrixElement:JHUGen',
			'Name:GG_SIG_ghg2_1_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:ZZGG MatrixElement:JHUGen'
			)),
		src = cms.InputTag("bareZZCand"),
		)




process.ZZTree= cms.EDAnalyzer('MiniAnalyzer',
		genjet = cms.InputTag("cleanJets"),
		pruned = cms.InputTag("prunedGenParticles"),
		VVDecayMode = cms.int32(0),
		VVMode = cms.int32(1),
		failedTreeLevel = cms.int32(2),
		fileName = cms.untracked.string('candTree'),
		sampleName = cms.string('HJJ0PM_M125'),
    		CandCollection = cms.untracked.string('ZZCand'),
		skipEmptyEvents = cms.bool(True),
    GenBR = cms.double(0.0002745),
    GenXSEC = cms.double(29.99),
    		xsec = cms.double(0.01333521),
		lheProbabilities = cms.vstring(
			'Name:SampleDecayHypothesisJHUGen Alias:<Name> Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0 isGen:1 NoBranch:1',
			'Name:Dec_SIG_ghz1_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz1=1,0 Options:DivideP=SampleDecayHypothesisJHUGen isGen:1',
			'Name:Dec_SIG_ghz4_1_JHUGen Process:SelfDefine_spin0 Production:ZZINDEPENDENT MatrixElement:JHUGen Couplings:ghz4=1,0 Options:DivideP=SampleDecayHypothesisJHUGen isGen:1',
			'Name:SampleProductionHypothesisJHUGen Alias:<Name> Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg2=1,0;ghz1=1,0 Cluster:BestLOAssociatedVBF isGen:1 NoBranch:1',
			'Name:HJJ_SIG_ghg2_1_JHUGen Alias:SampleProductionHypothesisJHUGen Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg2=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:BestLOAssociatedVBF isGen:1',
			'Name:HJJ_SIG_ghg4_1_JHUGen Process:SelfDefine_spin0 Production:JJQCD MatrixElement:JHUGen Couplings:ghg4=1,0 Options:DivideP=SampleProductionHypothesisJHUGen Cluster:BestLOAssociatedVBF isGen:1',
			),
		recoProbabilities = cms.vstring( ('Name:GG_SIG_ghg2_1_ghz1_1_JHUGen Alias:<Name> Process:HSMHiggs Production:ZZGG MatrixElement:JHUGen Options:AddPConst=1',
			'Name:GG_SIG_ghg2_1_ghz2_1_JHUGen Alias:<Name> Process:H0hplus Production:ZZGG MatrixElement:JHUGen',
			'Name:GG_SIG_ghg2_1_ghz4_1_JHUGen Alias:<Name> Process:H0minus Production:ZZGG MatrixElement:JHUGen'
			))
		)


process.Candidates = cms.Path(process.softLeptons+process.cleanJets+process.bareZCand+process.ZCand+process.bareZZCand+process.ZZCand)

process.trees = cms.EndPath(process.ZZTree)
		#process.p = cms.Path(process.ZZTree)






