import FWCore.ParameterSet.Config as cms
jetAnalysis = cms.PSet(
    observables     = cms.vstring('pt,50.,0.,300.', 'hcalE,50,0,300', 'eOverP,50,0,5'),
    chargedPFObservables     = cms.vstring('pt,50.,0.,300.'),
    neutralPFObservables     = cms.vstring('pt,50.,0.,300.'),
    cutList     = cms.vstring('eta,-2,2', 'nTrkInBlock,0.5,1.5'),
    jetCutList     = cms.vstring('pt,150,10000'),
)
