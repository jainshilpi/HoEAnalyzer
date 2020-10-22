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
#process.load("Configuration.Geometry.GeometryExtended2017_cff")
#process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")
process.load("Configuration.Geometry.GeometryExtended2021Reco_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '') # 2018 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2023_realistic_v3', '') # 2023 MC
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '') # 2021 MC

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                ####samples i made (Single Ele gun )
                                #'/store/user/swmukher/SingleElectronPt2to100/SingleEle_MINIAOD_v3/191125_090555/0000/miniaod_99.root',
                                #'/store/user/swmukher/SingleElectronPt2to100/SingleEle_MINIAOD_v3/191125_090555/0000/miniaod_98.root'
                                # Run3 file to check hcal depths
                                #'/store/mc/Run3Summer19MiniAOD/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v2/260000/FE3A7D24-F46E-8744-B92B-F6115FD395A0.root'
                                # example 2021 file
                                #'/store/mc/Run3Summer19MiniAOD/QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/2C20FC5F-B670-4C46-BCC9-02EDF6EFA5F6.root'
                                #'/store/mc/Run3Winter20DRMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8/MINIAODSIM/DRFlatPU30to80_110X_mcRun3_2021_realistic_v6-v2/50000/902C30AC-3F0B-F646-B8C1-B5819EE94912.root'
                                #'/store/mc/Run3Winter20DRMiniAOD/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/230000/96ABF29F-33FF-754F-B524-A22A9E47C212.root'

                            )
)

#we now setup VID for HEEP V7.0, we will do this with the handy function below
#which creates the standard process.egmGsfElectronIDSequence for VID which we
#add to our path
#from HEEP.VID.tools import setupVIDForHEEPV70
#setupVIDForHEEPV70(process,useMiniAOD=useMiniAOD)

#process.demo = cms.EDAnalyzer('EleHoEAnalyzer',
process.demo = cms.EDAnalyzer('FlatEleHoEAnalyzer',
                              genEventSrc = cms.InputTag('generator'),
                              electrons = cms.InputTag('slimmedElectrons'),
                              pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                              #hbheInput = cms.InputTag("reducedEgamma" ,  "reducedHBHEHits" ,  "PAT"),
                              hbheInput = cms.InputTag("reducedEgamma" ,  "reducedHBHEHits" ),
                              ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                              eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                              esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                              genParticleSrc       = cms.InputTag("prunedGenParticles"),
                              #genParticleSrc       = cms.InputTag("genParticles"),
                              Run2_2018_ = cms.bool(False),
                              rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),

                              #trkIsoMap=cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
                              #vid=cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                              #vidBitmap=cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70Bitmap")
                              output_file = cms.string("flat_electron.root"),

                              
                          )
#process.TFileService = cms.Service("TFileService", fileName = cms.string('electron.root'))

process.p = cms.Path(process.demo)
