import FWCore.ParameterSet.Config as cms

process = cms.Process('ParticleFlowDQMOffline')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v4'

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# my analyzer
process.load('DQMOffline.ParticleFlow.runBasic_cfi')

# back to original script
with open('fileList_2.log') as f:
    lines = f.readlines()

#Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(lines))

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))


process.p = cms.Path(process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain*process.PFAnalyzer)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

## Schedule definition
process.schedule = cms.Schedule(
    process.p,
    process.DQMoutput_step
    )







