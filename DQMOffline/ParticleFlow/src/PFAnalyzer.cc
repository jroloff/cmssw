/** \class PFAnalyzer
 *
 *  DQM ParticleFlow analysis monitoring
 *
 *  \author J. Roloff - Brown University
 *
 */

#include "DQMOffline/ParticleFlow/interface/PFAnalyzer.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"

#include <string>
#include <cmath>

using namespace edm;
using namespace reco;

// ***********************************************************
PFAnalyzer::PFAnalyzer(const edm::ParameterSet& pSet){
  m_directory = "ParticleFlow";
  parameters_ = pSet.getParameter<edm::ParameterSet>("jetAnalysis2");

  thePfCandidateCollection_ = consumes<reco::PFCandidateCollection>(pSet.getParameter<edm::InputTag>("pfCandidates"));
  pfJetsToken_ = consumes<reco::PFJetCollection>(pSet.getParameter<edm::InputTag>("pfJetCollection"));

  m_observables = parameters_.getParameter<vstring>("observables");

  // TODO we might want to define some observables that are only for PFCs in jets,
  // and that might take a jet as an input to the function as well

  m_cutList = parameters_.getParameter<vstring>("cutList");
  m_jetCutList = parameters_.getParameter<vstring>("jetCutList");

  // Link observable strings to the static functions defined in the header file
  // Many of these are quite trivial, but this enables a simple way to include a 
  // variety of observables on-the-fly.
  // These should not be deleted, even if you do not currently plan to make
  // plots with them.
  m_funcMap["pt"] = getPt;
  m_funcMap["eta"] = getEta;
  m_funcMap["phi"] = getPhi;
  m_funcMap["hcalE"] = getHCalEnergy;
  m_funcMap["eOverP"] = getEoverP;
  m_funcMap["nTrkInBlock"] = getNTracksInBlock;

  // Link jet observables to static functions in the header file.
  // This is very similar to m_funcMap, but for jets instead.
  m_jetFuncMap["pt"] = getJetPt;

  // Convert the cutList strings into real cuts that can be applied
  // The format should be three comma separated values
  // with the first number being the name of the observable
  // (corresponding to a key in m_funcMap),
  // the second being the minimum value, and the last being the max.
  for(unsigned int i=0; i<m_cutList.size(); i++){
    size_t pos = m_cutList[i].find(",");
    std::string observableName = m_cutList[i].substr(0, pos);
    m_cutList[i].erase(0, pos + 1);

    pos = m_cutList[i].find(",");
    m_cutMins.push_back( atof(m_cutList[i].substr(0, pos).c_str()));
    m_cutList[i].erase(0, pos + 1);

    m_cutMaxes.push_back(atof(m_cutList[i].c_str()));
    m_cutList[i] = observableName;
  }

  // Convert the jetCutList strings into real cuts that can be applied
  // The format should be three comma separated values,
  // with the first number being the name of the observable
  // (corresponding to a key in m_jetFuncMap),
  // the second being the minimum value, and the last being the max value.
  for(unsigned int i=0; i<m_jetCutList.size(); i++){
    size_t pos = m_jetCutList[i].find(",");
    std::string observableName = m_jetCutList[i].substr(0, pos);
    m_jetCutList[i].erase(0, pos + 1);

    pos = m_jetCutList[i].find(",");
    m_jetCutMins.push_back( atof(m_jetCutList[i].substr(0, pos).c_str()));
    m_jetCutList[i].erase(0, pos + 1);

    m_jetCutMaxes.push_back(atof(m_jetCutList[i].c_str()));
    m_jetCutList[i] = observableName;
  }
}

// ***********************************************************
PFAnalyzer::~PFAnalyzer() {
  LogTrace("PFAnalyzer") << "[PFAnalyzer] Saving the histos";
}

// ***********************************************************
void PFAnalyzer::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {
  ibooker.setCurrentFolder(m_directory);

  // TODO: Make it possible to use an arbitrary list of bins instead of evenly space bins
  //
  // Books a histogram for each histogram in the config file.
  // The format for the observables should be four comma separated values,
  // with the first being the observable name (corresponding to one of 
  // the keys in m_funcMap), the second being the number of bins,
  // and the last two being the min and max value for the histogram respectively.
  for(unsigned int i=0; i<m_observables.size(); i++){
    size_t pos = m_observables[i].find(",");
    std::string observableName = m_observables[i].substr(0, pos);
    m_observables[i].erase(0, pos + 1);
    
    pos = m_observables[i].find(",");
    int nBins = atoi(m_observables[i].substr(0, pos).c_str());
    m_observables[i].erase(0, pos + 1);
    
    pos = m_observables[i].find(",");
    float binMin = atof(m_observables[i].substr(0, pos).c_str());
    m_observables[i].erase(0, pos + 1);
    
    float binMax = atof(m_observables[i].c_str());
    m_observables[i] = observableName;


    // For each observable, we make a couple histograms based on a few generic categorizations.
    // In all cases, the PFCs that go into these histograms must pass the PFC selection from m_cutList.

    // This histogram is all PFCs that pass the basic selection
    std::string histName = Form("allPFC_%s", observableName.c_str());
    MonitorElement* mHist = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHist));

    // This histogram is all neutral PFCs that pass the basic selection
    histName = Form("neutralHadPFC_%s", observableName.c_str());
    MonitorElement* mHistNeutral = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistNeutral));

    // This histogram is all neutral PFCs that pass the basic selection
    histName = Form("chargedHadPFC_%s", observableName.c_str());
    MonitorElement* mHistCharged = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistCharged));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("electronPFC_%s", observableName.c_str());
    MonitorElement* mHistElectron = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistElectron));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("muonPFC_%s", observableName.c_str());
    MonitorElement* mHistMuon = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistMuon));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("gammaPFC_%s", observableName.c_str());
    MonitorElement* mHistGamma = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistGamma));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("hadHFPFC_%s", observableName.c_str());
    MonitorElement* mHistHadHF = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistHadHF));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("emHFPFC_%s", observableName.c_str());
    MonitorElement* mHistEMHF = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistEMHF));



    // These histograms are for PFCs passing the basic selection, and which are matched to jets
    // that pass the jet selection
    histName = Form("allPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistInJet));

    // This histogram is all neutral PFCs that pass the basic selection
    histName = Form("neutralHadPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistNeutralInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistNeutralInJet));

    // This histogram is all neutral PFCs that pass the basic selection
    histName = Form("chargedHadPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistChargedInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistChargedInJet));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("electronPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistElectronInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistElectronInJet));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("muonPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistMuonInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistMuonInJet));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("gammaPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistGammaInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistGammaInJet));

    // This histogram is electron PFCs that pass the basic selection
    histName = Form("hadHFPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistHadHFInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistHadHFInJet));
    
    // This histogram is electron PFCs that pass the basic selection
    histName = Form("emHFPFC_jetMatched_%s", observableName.c_str());
    MonitorElement* mHistEMHFInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHistEMHFInJet));
  }
}

void PFAnalyzer::bookMESetSelection(std::string DirName, DQMStore::IBooker& ibooker) {
  ibooker.setCurrentFolder(DirName);
}

// ***********************************************************
void PFAnalyzer::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
}

bool PFAnalyzer::passesPFCuts(const reco::PFCandidate pfCand){
  for(unsigned int i=0; i<m_cutList.size(); i++){
    if(m_funcMap[m_cutList[i]](pfCand) < m_cutMins[i] || m_funcMap[m_cutList[i]](pfCand) >= m_cutMaxes[i])  return false;
  }

  return true; 
}

bool PFAnalyzer::passesJetCuts(const reco::PFJet jetCand){
  for(unsigned int i=0; i<m_jetCutList.size(); i++){
    if(m_jetFuncMap[m_jetCutList[i]](jetCand) < m_jetCutMins[i] || m_jetFuncMap[m_jetCutList[i]](jetCand) >= m_jetCutMaxes[i])  return false;
  }

  return true;
}

// ***********************************************************
void PFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Retrieve the PFCs
  edm::Handle<reco::PFCandidateCollection> pfCollection;
  iEvent.getByToken(thePfCandidateCollection_, pfCollection);
  if (!pfCollection.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF candidate \n";
    return;
  }

  edm::Handle<PFJetCollection> pfJets;
  iEvent.getByToken(pfJetsToken_, pfJets);
  if (!pfJets.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF jets \n";
    return;
  }

  // TODO Perform some event selection
  // Probably we want to define a few different options for how the selection will work
  // if(!passSelection(iEvent)) return;


  for (reco::PFCandidateCollection::const_iterator recoPF = pfCollection->begin(); recoPF != pfCollection->end(); ++recoPF) {
    if(!passesPFCuts(*recoPF)) continue;
    bool isCharged = false;

    // Eventually, we might want the hist name to include the cuts that we are applying, 
    // so I am keepking it as a separate string for now, even though it is redundant.
    // Make plots of all observables
    for(unsigned int i=0; i<m_observables.size(); i++){
      std::string histName = Form("%s", m_observables[i].c_str());
      map_of_MEs[m_directory + "/allPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));

      switch(recoPF->particleId()){
        case PFCandidate::ParticleType::h :
          map_of_MEs[m_directory + "/chargedHadPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::h0 :
          map_of_MEs[m_directory + "/neutralHadPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::e :
          map_of_MEs[m_directory + "/electronPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::mu :
          map_of_MEs[m_directory + "/muonPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::gamma :
          map_of_MEs[m_directory + "/gammaPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::h_HF :
          map_of_MEs[m_directory + "/hadTowerPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        case PFCandidate::ParticleType::egamma_HF :
          map_of_MEs[m_directory + "/emTowerPFC_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
          break;
        default:
           break;
      }
    }
  }


  // Make plots of all observables, this time for PF candidates within jets
  for (reco::PFJetCollection::const_iterator cjet = pfJets->begin(); cjet != pfJets->end(); ++cjet) {
    if(!passesJetCuts(*cjet)) continue;
    std::vector<reco::PFCandidatePtr> pfConstits = cjet->getPFConstituents();

    for(auto recoPF:pfConstits){
      if(!passesPFCuts(*recoPF)) continue;

      for(unsigned int i=0; i<m_observables.size(); i++){
        std::string histName = Form("%s", m_observables[i].c_str());
        map_of_MEs[m_directory + "/allPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));

        switch(recoPF->particleId()){
          case PFCandidate::ParticleType::h :
            map_of_MEs[m_directory + "/chargedHadPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case PFCandidate::ParticleType::h0 :
            map_of_MEs[m_directory + "/neutralHadPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
             break;
          case PFCandidate::ParticleType::e :
            map_of_MEs[m_directory + "/electronPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case PFCandidate::ParticleType::mu :
            map_of_MEs[m_directory + "/muonPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case PFCandidate::ParticleType::gamma :
            map_of_MEs[m_directory + "/gammaPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case PFCandidate::ParticleType::h_HF :
            map_of_MEs[m_directory + "/hadTowerPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          case PFCandidate::ParticleType::egamma_HF :
            map_of_MEs[m_directory + "/emTowerPFC_jetMatched_" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
            break;
          default:
             break;
        }
      }
    }
  }

}
