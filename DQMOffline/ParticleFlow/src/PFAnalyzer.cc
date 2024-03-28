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
  m_chargedPFObservables = parameters_.getParameter<vstring>("chargedPFObservables");
  m_neutralPFObservables = parameters_.getParameter<vstring>("neutralPFObservables");


  m_cutList = parameters_.getParameter<vstring>("cutList");
  m_jetCutList = parameters_.getParameter<vstring>("jetCutList");

  // Link observable strings to the static functions defined in the header file
  // Many of these are quite trivial, but this enables a simple way to include a 
  // variety of observables on-the-fly
  m_funcMap["pt"] = getPt;
  m_funcMap["eta"] = getEta;
  m_funcMap["phi"] = getPhi;
  m_funcMap["hcalE"] = getHCalEnergy;
  m_funcMap["eOverP"] = getEoverP;
  m_funcMap["nTrkInBlock"] = getNTracksInBlock;

  m_jetFuncMap["pt"] = getJetPt;

  // Convert the cutList strings into real cuts that can be applied
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

    std::string histName = Form("%s", observableName.c_str());
    MonitorElement* mPt = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mPt));

    histName = Form("jetMatched_%s", observableName.c_str());
    MonitorElement* mPtInJet = ibooker.book1D(histName, ";PFC p_{T} [GeV];", nBins, binMin, binMax);
    map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mPtInJet));
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

/*
  //Vertex information
  Handle<VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);
  if (!vertexHandle.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid vertex collection \n";
    return;
  }
*/

  /*
  int numPV = 0;
  if (vertexHandle.isValid()) {
    VertexCollection vertexCollection = *(vertexHandle.product());
    numPV = vertexCollection.size();
  }

  */

  edm::Handle<PFJetCollection> pfJets;
  iEvent.getByToken(pfJetsToken_, pfJets);
  if (!pfJets.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF jets \n";
    return;
  }

  // Perform some event selection?

  // Apply cuts to jets for the matching studies
  // It's probably better to slim this down so we don't have to check this for every jet with every PFC
  

  for (reco::PFCandidateCollection::const_iterator recoPF = pfCollection->begin(); recoPF != pfCollection->end(); ++recoPF) {
    if(!passesPFCuts(*recoPF)) continue;
    // Make plots of all observables
    for(unsigned int i=0; i<m_observables.size(); i++){
      std::string histName = Form("%s", m_observables[i].c_str());
      map_of_MEs[m_directory + "/" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
    }
    for (reco::PFJetCollection::const_iterator cjet = pfJets->begin(); cjet != pfJets->end(); ++cjet) {
      if(!passesJetCuts(*cjet)) continue;
      if (deltaR(cjet->eta(), cjet->phi(), recoPF->eta(), recoPF->phi()) > 0.2) continue;

      for(unsigned int i=0; i<m_observables.size(); i++){
        std::string histName = Form("jetMatched_%s", m_observables[i].c_str());
        map_of_MEs[m_directory + "/" + histName]->Fill(m_funcMap[m_observables[i]](*recoPF));
      }
    }

/*
    bool isCharged = false;
    bool isJetMatched = false;
    // Make plots for charged PFCs
    if(isCharged){
    }
    // Make plots for neutral PFCs
    else{
    }

    // Make plots for PFCs matched to jets (with some selection)
    if(isJetMatched){
    }
*/
  }



}
