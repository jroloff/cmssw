#ifndef PFClusterAnalyzer_H
#define PFClusterAnalyzer_H

/** \class JetMETAnalyzer
 *
 *  DQM PF candidate analysis monitoring
 *
 *  \author J. Roloff - Brown University
 *
 */

#include <memory>
#include <fstream>
#include <utility>
#include <string>
#include <cmath>
#include <map>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
class PFClusterAnalyzer;

class PFClusterAnalyzer : public DQMEDAnalyzer {
public:
  /// Constructor
  PFClusterAnalyzer(const edm::ParameterSet&);

  /// Destructor
  ~PFClusterAnalyzer() override;

  /// Inizialize parameters for histo binning
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  /// Initialize run-based parameters
  void dqmBeginRun(const edm::Run&, const edm::EventSetup&) override;

private:
  struct binInfo;
  // A map between an observable name and a function that obtains that observable from a  PFCluster.
  // This allows us to construct more complicated observables easily, and have it more configurable
  // in the config file.
  std::map<std::string, std::function<std::vector<double>(const reco::PFCluster)>> m_funcMap;
  std::map<std::string, std::function<double(const reco::PFClusterCollection)>>
      m_eventFuncMap;

  std::map<std::string,
           std::function<double(const std::vector<reco::PFClusterRef> pfCands, const reco::PFJet)>>
      m_jetWideFuncMap;

  std::map<std::string, std::function<double(const reco::PFCluster, const reco::PFJet)>> m_pfInJetFuncMap;
  std::map<std::string, std::function<double(const reco::PFJet)>> m_jetFuncMap;

  std::map<std::string, std::function<bool(const edm::Handle<std::vector<reco::PFJet> >& pfJets)>> m_eventSelectionMap;

  binInfo getBinInfo(std::string);

  // Book MonitorElements
  void bookMESetSelection(std::string, DQMStore::IBooker&);


  std::vector<int> getPFBin(const reco::PFCluster pfCand, int i);
  //int getJetBin(const reco::PFJet jetCand, int i);


  int getBinNumber(double binVal, std::vector<double> bins);
  std::vector<int> getBinNumbers(std::vector<std::vector<double> > binVal, std::vector<std::vector<double>> bins);
  std::vector<double> getBinList(std::string binString);

  std::vector<std::string> getAllSuffixes(std::vector<std::string> observables,
                                          std::vector<std::vector<double>> binnings);
  std::string stringWithDecimals(int bin, std::vector<double> bins);

  std::string getSuffix(std::vector<int> binList,
                        std::vector<std::string> observables,
                        std::vector<std::vector<double>> binnings);


  // Various functions designed to get information from a PF canddidate
  static std::vector<double> getPt(const reco::PFCluster pfCand) { return {pfCand.pt()}; }
  static std::vector<double> getEnergy(const reco::PFCluster pfCand) { return {pfCand.energy()}; }
  static std::vector<double> getEta(const reco::PFCluster pfCand) { return {pfCand.eta()}; }
  static std::vector<double> getAbsEta(const reco::PFCluster pfCand) { return {std::abs(pfCand.eta())}; }
  static std::vector<double> getPhi(const reco::PFCluster pfCand) { return {pfCand.phi()}; }
  static std::vector<double> getTime(const reco::PFCluster pfCand) { return {pfCand.time()}; }
  static std::vector<double> getTimeError(const reco::PFCluster pfCand) { return {pfCand.timeError()}; }
  static std::vector<double> getDepth(const reco::PFCluster pfCand) { return {pfCand.depth()}; }
  static std::vector<double> getLayer(const reco::PFCluster pfCand) { return {double(pfCand.layer())}; }


  static std::vector<double> getRecHitTimes(const reco::PFCluster pfCand){
    std::vector<double> times;
    const std::vector<reco::PFRecHitFraction> recHits= pfCand.recHitFractions();
    for(auto &rechit : recHits){
      times.push_back(rechit.recHitRef()->time()); 
 

    }
    return times;

  }



  bool passesTriggerSelection(const edm::Handle<edm::TriggerResults>& triggerResults, const edm::TriggerNames& triggerNames, const std::vector<std::string> triggerOptions){

    // Hack to make it pass the lowest unprescaled HLT?
    Int_t JetHiPass = 0;

    const unsigned int nTrig(triggerNames.size());
    for (unsigned int i = 0; i < nTrig; ++i) {
      for (unsigned int j = 0; j < triggerOptions.size(); ++j) {
        if(triggerOptions[j] == "") {
          JetHiPass=1;
          break;
        }
        if (triggerNames.triggerName(i).find(triggerOptions[j]) != std::string::npos && triggerResults->accept(i)) {
          JetHiPass = 1;
          break;
        }
      }
      if(JetHiPass) break;
    }

    if (!JetHiPass)
      return false;
    return true;
  }


  static bool passesNoCutSelection(const edm::Handle<std::vector<reco::PFJet> >& pfJets) {
    return true;
  }

  static bool passesDijetSelection(const edm::Handle<std::vector<reco::PFJet> >& pfJets) {
    if (pfJets->size() < 2)
      return false;
    if (pfJets->at(0).pt() < 450)
      return false;
    if (pfJets->at(0).pt() / pfJets->at(1).pt() > 2)
      return false;
    
    return true;
  } 


  edm::EDGetTokenT<reco::PFClusterCollection> thePfCandidateCollection_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetsToken_;
  edm::InputTag srcWeights;

  edm::EDGetTokenT<edm::ValueMap<float>> weightsToken_;
  edm::ValueMap<float> const* weights_;

  edm::EDGetTokenT<GenEventInfoProduct> tok_ew_;

  edm::InputTag theTriggerResultsLabel_;
  edm::InputTag vertexTag_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  std::string m_selection;


  std::vector<std::vector<std::string>> m_allSuffixes;

  // The directory where the output is stored
  std::string m_directory;

  // All of the histograms, stored as a map between the histogram name and the histogram
  std::map<std::string, MonitorElement*> map_of_MEs;

  //std::map<reco::PFCluster::ParticleType, std::string> m_particleTypeName;



  //check later if we need only one set of parameters
  edm::ParameterSet parameters_;

  typedef std::vector<std::string> vstring;
  typedef std::vector<double> vDouble;

  vstring m_triggerOptions;
  // Information on which observables to make histograms for.
  // In the config file, this should come as a comma-separated list of
  // the observable name, the number of bins for the histogram, and
  // the lowest and highest values for the histogram.
  // The observable name should have an entry in m_funcMap to define how
  // it can be retrieved from a PFCluster.
  vstring m_observables;
  vstring m_observableNames;


  // Information on what cuts should be applied to PFClusters that are
  // being monitored. In the config file, this should come as a comma-separated list of
  // the observable name, and the lowest and highest values for the histogram.
  // The observable name should have an entry in m_funcMap to define how
  // it can be retrieved from a PFCluster.
  vstring m_cutList;
  std::vector<std::vector<std::string>> m_fullCutList;
  std::vector<std::vector<std::vector<double>>> m_binList;

  vDouble m_npvBins;

  std::vector<std::string> m_pfNames;
};

struct PFClusterAnalyzer::binInfo {
  std::string observable;
  std::string axisName;
  int nBins;
  double binMin;
  double binMax;
};

DEFINE_FWK_MODULE(PFClusterAnalyzer);
#endif
