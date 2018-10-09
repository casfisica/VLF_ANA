from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TASK'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_Data_cfg.py'
#config.JobType.maxMemoryMB = 2500

config.Data.inputDataset = 'DATASAMPLE'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 5000
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/jruizalv/VLF_ANA/OUTPUTDIR/'
config.Data.outLFNDirBase ='/store/user/csalazar/TEST_VLF_ANA/OUTPUTDIR/'
config.Data.publication = False
config.Data.outputDatasetTag = 'TASK'

config.Site.storageSite = 'T2_CH_CERN'
