from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DeepNtuples_ttbar'
# config.General.requestName = 'DeepNtuples_qcd'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'DeepNtuplizer.py'
config.JobType.inputFiles = ['QGL_cmssw8020_v2.db']

config.section_("Data")
# config.Data.inputDataset = '/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8/PhaseIIMTDTDRAutumn18MiniAOD-PU200_103X_upgrade2023_realistic_v2-v1/MINIAODSIM'
config.Data.inputDataset = '/TTbar_14TeV_TuneCP5_Pythia8/PhaseIIMTDTDRAutumn18MiniAOD-PU200_103X_upgrade2023_realistic_v2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'CRAB3_tutorial_May2015_MC_analysis'

config.section_("Site")
config.Site.storageSite = 'T1_DE_KIT_Disk'

config.section_("User")
config.User.voGroup = 'dcms'

