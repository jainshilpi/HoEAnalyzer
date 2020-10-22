if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

from multiprocessing import Process

from CRABClient.UserUtilities import config
config = config()

### CRAB Common
#config.General.workArea = '/afs/cern.ch/work/a/afiqaize/egm/hcal_depth_200122/CMSSW_10_6_8/src/HoE/HoEAnalyzer/test/crab/'
#config.General.workArea = '/afs/cern.ch/work/a/afiqaize/egm/hcal_depth_200122/CMSSW_11_0_1/src/HoE/HoEAnalyzer/test/crab/'
config.General.workArea = 'crab_HoE'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.JobType.psetName = 'flat_ele_cfg.py'
config.JobType.psetName = 'ele_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 2500
config.Data.totalUnits = -1
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.publication = False
#config.Data.outLFNDirBase = '/store/group/phys_egamma/shilpi/hcal_depth_200701/'
#config.Data.outLFNDirBase = '/store/group/phys_egamma/shilpi/hcal_depth_200715/'
config.Data.outLFNDirBase = '/store/group/phys_egamma/shilpi/hcal_depth_200818/'

config.Site.storageSite = 'T2_CH_CERN'

config.Data.ignoreLocality = False
#config.Site.whitelist = ['T2_CH_CERN','T2_IT_Pisa','T2_UK_London_IC','T2_HU_Budapest','T2_DE_DESY']

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

### -------------------------------------- ###
### CRAB Request Specific
### -------------------------------------- ###

job_req = {
    ## 106X
    #"dy":               "/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"tt":               "/TTToSemiLeptonic_TuneCP5_14TeV-powheg-pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"yj":               "/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_15_30":     "/QCD_Pt_15to30_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_30_50":     "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_50_80":     "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_80_120":    "/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_120_170":   "/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_170_300":   "/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_300_470":   "/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_470_600":   "/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_600_800":   "/QCD_Pt_600to800_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_800_1000":  "/QCD_Pt_800to1000_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_1000_1400": "/QCD_Pt_1000to1400_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_1400_1800": "/QCD_Pt_1400to1800_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_1800_2400": "/QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_2400_3200": "/QCD_Pt_2400to3200_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",
    #"qcd_pt_3200":      "/QCD_Pt_3200toInf_TuneCP5_14TeV_pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM",

    ## 110X
#    "dy":               "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Run3Winter20DRMiniAOD-DRFlatPU30to80_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "tt":               "/TTbar_14TeV_TuneCP5_Pythia8/Run3Winter20DRMiniAOD-FlatPU0to80ALCARECO_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "qcd3000":          "/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8/Run3Winter20DRMiniAOD-DRFlatPU30to80_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "qcd7000":          "/QCD_Pt_15to7000_TuneCP5_Flat_14TeV_pythia8/Run3Winter20DRMiniAOD-FlatPU0to80ALCARECO_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "yy":               "/DoublePhoton_FlatPt-5To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "zp":               "/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/MINIAODSIM",
    "ttMtt1000ToInf":               "/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6_ext1-v1/MINIAODSIM",

    ## 110X private
    #"pu0":              "/DoubleEle-FlatPt2To250-PU0_13TeV-110X_mcRun3_2021_realistic_v6/afiqaize-MINIAODSIM-28ea11753f803e89e7fd3e375196acc9/USER"
}

for req, dat in job_req.items():
    config.General.requestName = req
    config.Data.inputDataset = dat
    #config.Data.secondaryInputDataset = ''
    #config.JobType.pyCfgParams = ['doTrigger=True', 'isMC=True']
    config.JobType.outputFiles = ['flat_electron.root']
    #config.Data.lumiMask = ''
    #config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_DESY']

    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
