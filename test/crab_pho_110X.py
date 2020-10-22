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
config.JobType.psetName = 'pho_cfg.py'
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
config.Data.allowNonValidInputDataset = True
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


#    "qcd3000":          "/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8/Run3Winter20DRMiniAOD-DRFlatPU30to80_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "zp":               "/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/MINIAODSIM",
#    "doublePho5To300": "/DoublePhoton_FlatPt-5To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v2/MINIAODSIM",
#    "doublePho300To1000": "/DoublePhoton_FlatPt-300To1000/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v1/MINIAODSIM",
#    "doublePho1000To1500": "/DoublePhoton_FlatPt-1000To1500/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v1/MINIAODSIM"
    
    ## 110X private
    #"pu0":              "/DoubleEle-FlatPt2To250-PU0_13TeV-110X_mcRun3_2021_realistic_v6/afiqaize-MINIAODSIM-28ea11753f803e89e7fd3e375196acc9/USER"
    "doublePho0to200NOPU" : "/PhotonFlatPt0To200/Run3Winter20DRMiniAOD-NoPU_110X_mcRun3_2021_realistic_v6-v3/MINIAODSIM",
    "doublePho200to500NOPU" : "/PhotonFlatPt200To500/Run3Winter20DRMiniAOD-NoPU_110X_mcRun3_2021_realistic_v6-v3/MINIAODSIM"
    
}

for req, dat in job_req.items():
    config.General.requestName = req
    config.Data.inputDataset = dat
    #config.Data.secondaryInputDataset = ''
    #config.JobType.pyCfgParams = ['doTrigger=True', 'isMC=True']
    config.JobType.outputFiles = ['flat_photon.root']
    #config.Data.lumiMask = ''
    #config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_DESY']

    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
