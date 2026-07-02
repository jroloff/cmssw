# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 --filein file:/store/mc/RunIISummer16DR80Premix/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RAW-RECO/HighMET-PUMoriond17_skim_test_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/70000/E4E4E436-B57C-E711-8312-0CC47AA992B2.root --fileout file:output_custom_PF.root --conditions auto:phase1_2024_realistic --step RAW2DIGI,L1Reco,RECO --datatier RECO --eventcontent RECO --era Run3 --python_filename myCustomReco_cfg.py --no_exec -n 100 --mc
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://xrootd-cms.infn.it//store/mc/Run3Winter26Digi/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/GEN-SIM-RAW/150X_mcRun3_2026_realistic_v4-v2/100000/0066f01c-a58b-4417-bf6f-df4d0058940e.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:output_custom_PF.root'),
    outputCommands = process.RECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# -----------------------------------------------------------------------
# Load the standard PF cluster prerequisites.
# -----------------------------------------------------------------------
#process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cff")
#process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterHCAL_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowEHCluster_cff")
#process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")

#process.pfEHRecoTask = cms.Task(
#    process.pfEHClusterProducer,
#)

# -----------------------------------------------------------------------
# Output: write the new collections alongside the existing RECO content.
# -----------------------------------------------------------------------
process.RECOoutput.outputCommands.extend(cms.untracked.vstring(
        "keep recoPFCandidates_*_*_*",
        "keep recoPFEHClusters_*_*_*",
        "keep recoPFEHClusters_*_*_*",
        "keep recoPFBlocks_*_*_*",
    ),)

# Schedule
# -----------------------------------------------------------------------

#process.reconstruction_step.associate(process.pfEHRecoTask)
#process.pfEHRecoPath = cms.Path(cms.Task(process.pfEHRecoTask))
#process.schedule.append(process.pfEHRecoPath)


process.load("RecoParticleFlow.PFClusterProducer.particleFlowEHCluster_cff")

# What label does the task actually hold?
for m in process.particleFlowEHTask:
    print(m)

# Is the module in the task the same object as what's on the process?
print(process.pfEHClusterProducer)
print(id(process.pfEHClusterProducer))
#process.schedule.associate(process.particleFlowEHTask)
#process.schedule.append(process.particleFlowEHTask)

process.reconstruction_step.associate(process.particleFlowEHTask)
process.pfEHClusterProducerPath = cms.Path(cms.Task(process.pfEHClusterProducer))
process.schedule.append(process.pfEHClusterProducerPath)
#process.schedule.append(process.particleFlowEHSequence)


process.p = cms.Path(
    process.particleFlowEHCluster
    )


#print(process.schedule)
#process.options.wantSummary = cms.untracked.bool(True)

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



