from CRABClient.UserUtilities import config
config = config()

# sample = "TTToHadronic_noPU"
sample = "TTToHadronic_PU200"

config.General.requestName       = sample
config.General.workArea          = "/afs/cern.ch/user/e/etsai/workspace/SecondaryVertexing_CMSSW_14_1_0_pre0/src/DeepNTuples/DeepNtuplizer/production"
config.General.transferOutputs   = True
config.General.instance          = "prod"

config.JobType.pluginName        = "Analysis"
config.JobType.psetName          = "DeepNtuplizer.py"
config.JobType.inputFiles        = ["QGL_cmssw8020_v2.db"]
config.JobType.maxMemoryMB       = 3000
config.JobType.maxJobRuntimeMin  = 60

config.Data.splitting            = "FileBased"
config.Data.unitsPerJob          = 1
config.Data.outLFNDirBase        = "/store/user/etsai"
config.Data.publication          = False
config.Data.outputDatasetTag     = sample
config.Data.userInputFiles       = open(sample + ".list").readlines()
config.Data.outputPrimaryDataset = "DeepNTuples"

config.Site.storageSite          = "T2_CH_CERN"
config.Site.whitelist            = ["T2_CH_CERN"]
