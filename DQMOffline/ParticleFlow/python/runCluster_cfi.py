import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

PFClusterAnalyzer = DQMEDAnalyzer("PFClusterAnalyzer",
    pfJetCollection        = cms.InputTag("ak4PFJetsCHS"),
    pfCandidates             = cms.InputTag("particleFlow"),
    PVCollection             = cms.InputTag("offlinePrimaryVertices"),
    recoPFClusters           = cms.InputTag("particleFlowClusterECAL"),

    TriggerResultsLabel        = cms.InputTag("TriggerResults::HLT"),
    TriggerNames = cms.vstring("HLT_PFJet450"),
    eventSelection = cms.string("nocut"),



    pfAnalysis = cms.PSet(
      # Bins of NPV for plots
      NPVBins = cms.vdouble(0,100),

      # A list of observables for which plots should be made.
      # The format should be a list of semicolon-separated values.
      # The first is the observable name, corresponding to a key in m_funcMap 
      # in PFAnalysis. The second is the TLatex string that serves as the x-axis
      # title (which must not include a semicolon).
      # The last values are the bins. If three values are given, then the values,
      # in order, are the number of bins, the lowest, and the highest values.
      # If any other number is given, this is just a list of bins for the histogram.
      observables     = cms.vstring('pt;p_{T,PFC};50.;0.;350.', 
                                    'eta;#eta;50;-5;5',
                                    'phi;#phi;50;-3.14;3.14',
                                    'energy;E;50;0;300',
                                    'time;time;60;0;30',
                                    'timeError;timeError;100;0;1',
                                    'depth;depth;100;0;10',
                                    'layer;layer;30;-15;15',
                                   ),

      

      # This is a list of multidimensional cuts that are applied for the plots.
      # In the case of multiple bins, every combination of bin is tested.
      # The format should be a list of semicolon-separated values.
      # The first is the observable name, corresponding to a key in m_funcMap 
      # in PFAnalysis. The last values are the bins, following the same
      # conventions as the observables.
      #
      # Since we may want to test multiple sets of cuts simultaneously, 
      # these are separated by '[' 
      # For example, for 
      # cutList     = cms.vstring('[pt;1;0;10000]'),
      # there is one histogram made for PFCs with 0 < pT < 10000.
      # Similarly, for cutList     = cms.vstring('[pt;1;0;10000]', '[pt;1;05;10000][eta;1;-5;5]'),
      # there is one histogram made for PFCs with 0 < pT < 10000,
      # and one histogram made for PFCs with 5 < pT < 10000 and -5 < eta < 5.
     
      cutList     = cms.vstring(
                                '[pt;1;0;10000]',
                                #'[pt;0;1;2;4;6;10;20;40;60;100][abseta;0;1.5;2.0;2.5;2.8;2.85;2.9;2.95;3]',
                                #'[pt;0;2;5;10;20;50;100;1000]',
                                #'[pt;1;0;10000][abseta;0;1;2;2.5;2.6;2.7;2.8;2.9;3;3.5;4.0;4.5]',
                                #'[pt;1;0;10000][abseta;0;1;1.5;2;2.5;3;3.5;4.0]',
                               ),

    )


)
