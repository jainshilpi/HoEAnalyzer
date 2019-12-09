import FWCore.ParameterSet.Config as cms

process = cms.Process("HoE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load("Geometry.HcalEventSetup.CaloTowerTopology_cfi")
process.load("Configuration.Geometry.GeometryExtended2017_cff")
process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 500

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')     # MC

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
                   ####samples i made (Single Ele gun )
                  '/store/user/swmukher/SingleElectronPt2to100/SingleEle_MINIAOD_v3/191125_090555/0000/miniaod_99.root',
                  '/store/user/swmukher/SingleElectronPt2to100/SingleEle_MINIAOD_v3/191125_090555/0000/miniaod_98.root'
                )
                            )

process.demo = cms.EDAnalyzer('HoEAnalyzer',
       electrons = cms.InputTag('slimmedElectrons'),
       pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
       hbheInput = cms.InputTag("reducedEgamma" ,  "reducedHBHEHits" ,  "PAT"),
                              )
process.TFileService = cms.Service("TFileService", fileName = cms.string('test_MINIAOD_2018.root'))

process.p = cms.Path(process.demo)
