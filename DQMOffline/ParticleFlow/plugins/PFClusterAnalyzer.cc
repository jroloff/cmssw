/** \class PFClusterAnalyzer
 *
 *  DQM ParticleFlow analysis monitoring
 *
 *  \author J. Roloff - Brown University
 *
 */

#include "DQMOffline/ParticleFlow/plugins/PFClusterAnalyzer.h"
#include <iostream>

// ***********************************************************
PFClusterAnalyzer::PFClusterAnalyzer(const edm::ParameterSet& pSet) {
  m_directory = "ParticleFlow";
  parameters_ = pSet.getParameter<edm::ParameterSet>("pfAnalysis");

  thePfCandidateCollection_ = consumes<reco::PFClusterCollection>(pSet.getParameter<edm::InputTag>("recoPFClusters"));
  pfJetsToken_ = consumes<reco::PFJetCollection>(pSet.getParameter<edm::InputTag>("pfJetCollection"));

  theTriggerResultsLabel_ = pSet.getParameter<edm::InputTag>("TriggerResultsLabel");
  m_selection = pSet.getParameter<std::string>("eventSelection");

  triggerResultsToken_ = consumes<edm::TriggerResults>(edm::InputTag(theTriggerResultsLabel_));
  m_triggerOptions = pSet.getParameter<vstring>("TriggerNames");

  vertexTag_ = pSet.getParameter<edm::InputTag>("PVCollection");
  vertexToken_ = consumes<std::vector<reco::Vertex>>(edm::InputTag(vertexTag_));

  tok_ew_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  m_observables = parameters_.getParameter<vstring>("observables");
  m_npvBins = parameters_.getParameter<vDouble>("NPVBins");

  // List of cuts applied to PFCs that we want to plot
  m_cutList = parameters_.getParameter<vstring>("cutList");

  m_eventSelectionMap["dijet"] = &passesDijetSelection;
  m_eventSelectionMap["nocut"] = &passesNoCutSelection;

  // Link observable strings to the static functions defined in the header file
  // Many of these are quite trivial, but this enables a simple way to include a
  // variety of observables on-the-fly.
  m_funcMap["pt"] = &getPt;
  m_funcMap["energy"] = getEnergy;
  m_funcMap["eta"] = getEta;
  m_funcMap["abseta"] = getAbsEta;
  m_funcMap["phi"] = getPhi;
  m_funcMap["time"] = getTime;
  m_funcMap["timeError"] = getTimeError;
  m_funcMap["depth"] = getDepth;
  m_funcMap["layer"] = getLayer;


  // Convert the cutList strings into real cuts that can be applied
  // The format should be three comma separated values
  // with the first number being the name of the observable
  // (corresponding to a key in m_funcMap),
  // the second being the minimum value, and the last being the max.
  for (unsigned int i = 0; i < m_cutList.size(); i++) {
    m_fullCutList.push_back(std::vector<std::string>());
    while (m_cutList[i].find("]") != std::string::npos) {
      size_t pos = m_cutList[i].find("]");
      m_fullCutList[i].push_back(m_cutList[i].substr(1, pos));
      m_cutList[i].erase(0, pos + 1);
    }
  }

  for (unsigned int i = 0; i < m_fullCutList.size(); i++) {
    m_binList.push_back(std::vector<std::vector<double>>());
    for (unsigned int j = 0; j < m_fullCutList[i].size(); j++) {
      size_t pos = m_fullCutList[i][j].find(";");
      std::string observableName = m_fullCutList[i][j].substr(0, pos);
      m_fullCutList[i][j].erase(0, pos + 1);

      m_binList[i].push_back(getBinList(m_fullCutList[i][j]));
      m_fullCutList[i][j] = observableName;
    }
  }

}

// ***********************************************************
PFClusterAnalyzer::~PFClusterAnalyzer() { LogTrace("PFClusterAnalyzer") << "[PFClusterAnalyzer] Saving the histos"; }

// ***********************************************************
void PFClusterAnalyzer::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {
  ibooker.setCurrentFolder(m_directory);

  for (unsigned int i = 0; i < m_fullCutList.size(); i++) {
    m_allSuffixes.push_back(getAllSuffixes(m_fullCutList[i], m_binList[i]));
  }


  for (unsigned int npv = 0; npv < m_npvBins.size() - 1; npv++) {
    std::string npvString = Form("npv_%.0f_%.0f", m_npvBins[npv], m_npvBins[npv + 1]);
    // TODO: Make it possible to use an arbitrary list of bins instead of evenly space bins?
    // It is not clear if this is straightforward to do with these classes and CMSSW.
    // If it is, it should be an easy change to the code.
    //
    //
    // Books a histogram for each histogram in the config file.
    // The format for the observables should be four comma separated values,
    // with the first being the observable name (corresponding to one of
    // the keys in m_funcMap), the second being the number of bins,
    // and the last two being the min and max value for the histogram respectively.
    for (unsigned int i = 0; i < m_observables.size(); i++) {
      std::string cObservable = m_observables[i];
      PFClusterAnalyzer::binInfo obsInfo = getBinInfo(cObservable);

      if (npv == 0)
        m_observableNames.push_back(obsInfo.observable);

      for (unsigned int j = 0; j < m_allSuffixes.size(); j++) {
        for (unsigned int n = 0; n < m_allSuffixes[j].size(); n++) {
          // Loop over all of the different types of PF candidates
            // For each observable, we make a couple histograms based on a few generic categorizations.
            // In all cases, the PFCs that go into these histograms must pass the PFC selection from m_cutList.
            std::string histName = Form("allPFC_%s%s_%s",
                                        obsInfo.observable.c_str(),
                                        m_allSuffixes[j][n].c_str(),
                                        npvString.c_str());
            MonitorElement* mHist = ibooker.book1D(
                histName, Form(";%s;", obsInfo.axisName.c_str()), obsInfo.nBins, obsInfo.binMin, obsInfo.binMax);
            map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHist));
          }

      }
    }

  }

  std::string histName = Form("NPV");
  MonitorElement* mHist = ibooker.book1D(histName, Form(";%s;", "N_PV"), 100, 0, 100);
  map_of_MEs.insert(std::pair<std::string, MonitorElement*>(m_directory + "/" + histName, mHist));
}

PFClusterAnalyzer::binInfo PFClusterAnalyzer::getBinInfo(std::string observableString) {
  PFClusterAnalyzer::binInfo binningDetails;

  size_t pos = observableString.find(";");
  binningDetails.observable = observableString.substr(0, pos);
  observableString.erase(0, pos + 1);

  std::vector<double> binList = getBinList(observableString);
  pos = observableString.find(";");
  binningDetails.axisName = observableString.substr(0, pos);
  observableString.erase(0, pos + 1);

  pos = observableString.find(";");
  binningDetails.nBins = atoi(observableString.substr(0, pos).c_str());
  observableString.erase(0, pos + 1);

  pos = observableString.find(";");
  binningDetails.binMin = atof(observableString.substr(0, pos).c_str());
  observableString.erase(0, pos + 1);

  binningDetails.binMax = atof(observableString.c_str());

  return binningDetails;
}

void PFClusterAnalyzer::bookMESetSelection(std::string DirName, DQMStore::IBooker& ibooker) {
  ibooker.setCurrentFolder(DirName);
}

// ***********************************************************
void PFClusterAnalyzer::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}


// How many significant digits do we need to save for the values to be distinct?
std::string PFClusterAnalyzer::stringWithDecimals(int bin, std::vector<double> bins) {
  double diff = bins[bin + 1] - bins[bin];
  double sigFigs = log10(diff);

  // We only want to save as many significant digits as we need to.
  // Currently, we might lose some information, so we should think about
  // if we want to specify more digits
  if (sigFigs >= 1) {
    return Form("%.0f_%.0f", bins[bin], bins[bin + 1]);
  }

  int nDecimals = int(-1 * sigFigs) + 1;
  // We do not want to use decimals since these can mess up histogram retrieval in some cases.
  // Instead, we use a 'p' to indicate the decimal.
  double newDigit = std::abs((bins[bin] - int(bins[bin])) * pow(10, nDecimals));
  double newDigit2 = std::abs((bins[bin + 1] - int(bins[bin + 1])) * pow(10, nDecimals));
  std::string signStringLow = "";
  std::string signStringHigh = "";
  if (bins[bin] < 0)
    signStringLow = "m";
  if (bins[bin + 1] < 0)
    signStringHigh = "m";

  int higherDigitsLow = (bins[bin]>0)?floor(bins[bin]):ceil(bins[bin]);
  int higherDigitsHigh = (bins[bin+1]>0)?floor(bins[bin+1]):ceil(bins[bin+1]);

  return Form("%s%dp%.0f_%s%dp%.0f",
              signStringLow.c_str(),
              std::abs(higherDigitsLow),
              newDigit,
              signStringHigh.c_str(),
              std::abs(higherDigitsHigh),
              newDigit2);
}

std::vector<double> PFClusterAnalyzer::getBinList(std::string binString) {
  std::vector<double> binList;

  while (binString.find(";") != std::string::npos) {
    size_t pos = binString.find(";");
    binList.push_back(atof(binString.substr(0, pos).c_str()));
    binString.erase(0, pos + 1);
  }
  binList.push_back(atof(binString.c_str()));

  if (binList.size() == 3) {
    int nBins = int(binList[0]);
    double minVal = binList[1];
    double maxVal = binList[2];
    binList.clear();

    for (int i = 0; i <= nBins; i++) {
      binList.push_back(minVal + i * (maxVal - minVal) / nBins);
    }
  }

  return binList;
}


std::vector<std::string> PFClusterAnalyzer::getAllSuffixes(std::vector<std::string> observables,
                                                    std::vector<std::vector<double>> binnings) {
  int nTotalBins = 1;
  std::vector<int> nBins;
  for (unsigned int i = 0; i < binnings.size(); i++) {
    nTotalBins = (binnings[i].size() - 1) * nTotalBins;
    nBins.push_back(binnings[i].size() - 1);
  }

  std::vector<std::vector<int>> binList;

  for (int i = 0; i < nTotalBins; i++) {
    binList.push_back(std::vector<int>());
  }

  int factor = nTotalBins;
  int otherFactor = 1;
  for (unsigned int i = 0; i < binnings.size(); i++) {
    factor = factor / nBins[i];

    for (int k = 0; k < factor; k++) {
      for (int j = 0; j < nBins[i]; j++) {
        for (int m = 0; m < otherFactor; m++) {
          int binNumber = k*nBins[i] + j * otherFactor + m;
          binList[binNumber].push_back(j);
        }
      }
    }
    otherFactor = otherFactor * nBins[i];
  }

  std::vector<std::string> allSuffixes;
  for (int i = 0; i < nTotalBins; i++) {
    allSuffixes.push_back(getSuffix(binList[i], observables, binnings));
  }

  return allSuffixes;
}

// Get a unique string corresponding to the selection cuts
std::string PFClusterAnalyzer::getSuffix(std::vector<int> binList,
                                  std::vector<std::string> observables,
                                  std::vector<std::vector<double>> binnings) {
  std::string suffix = "";
  for (unsigned int i = 0; i < binList.size(); i++) {
    if (binList[i] < 0)
      return "";
    std::string digitString = stringWithDecimals(binList[i], binnings[i]);

    suffix = Form("%s_%s_%s", suffix.c_str(), observables[i].c_str(), digitString.c_str());
  }

  return suffix;
}

int PFClusterAnalyzer::getBinNumber(double binVal, std::vector<double> bins) {
  if (binVal < bins[0])
    return -1;
  for (unsigned int i = 0; i < bins.size(); i++) {
    if (binVal < bins[i])
      return i - 1;
  }

  return -1;
}

std::vector<int> PFClusterAnalyzer::getBinNumbers(std::vector<std::vector<double> > binVal, std::vector<std::vector<double>> bins) {
  std::vector<int> newbins;
  
  for(unsigned int k=0; k<binVal.size(); k++){
  std::vector<int> cbins;
  std::vector<int> nBins;
  for (unsigned int i = 0; i < binVal[k].size(); i++) {
    int cbin = getBinNumber(binVal[k][i], bins[i]);
    if (cbin < 0)
      return newbins;
    nBins.push_back(bins[i].size() - 1);
    cbins.push_back(cbin);
  }

  int bin = 0;
  int factor = 1;
  for (unsigned int i = 0; i < binVal[k].size(); i++) {
    bin += cbins[i] * factor;
    factor = factor * nBins[i];
  }
  newbins.push_back(bin);
  }

  return newbins;
}

std::vector<int> PFClusterAnalyzer::getPFBin(const reco::PFCluster pfBlock, int i) {
  std::vector<std::vector<double> > binVals;
  for (unsigned int j = 0; j < m_fullCutList[i].size(); j++) {
    binVals.push_back(m_funcMap[m_fullCutList[i][j]](pfBlock));
  }

  return getBinNumbers(binVals, m_binList[i]);
}


// ***********************************************************
void PFClusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const edm::Handle<GenEventInfoProduct> genEventInfo = iEvent.getHandle(tok_ew_);
  double eventWeight = 1;
  if (genEventInfo.isValid()) {
    eventWeight = genEventInfo->weight();
  }


  //Vertex information
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);

  if (!vertexHandle.isValid()) {
    LogDebug("") << "PFClusterAnalyzer: Could not find vertex collection" << std::endl;
  }
  int numPV = 0;

  if (vertexHandle.isValid()) {
    reco::VertexCollection vertex = *(vertexHandle.product());
    for (reco::VertexCollection::const_iterator v = vertex.begin(); v != vertex.end(); ++v) {
      if (v->isFake())
        continue;
      if (v->ndof() < 4)
        continue;
      if (fabs(v->z()) > 24.0)
        continue;
      ++numPV;
    }
  }

  int npvBin = getBinNumber(numPV, m_npvBins);
  if (npvBin < 0)
    return;
  std::string npvString = Form("npv_%.0f_%.0f", m_npvBins[npvBin], m_npvBins[npvBin + 1]);

  // **** Get the TriggerResults container
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  if(!triggerResults.isValid()){
      edm::LogError("PFClusterAnalyzer") << "invalid trigger result \n";
      return;
  }
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);


  edm::Handle<reco::PFClusterCollection> pfCollection;
  iEvent.getByToken(thePfCandidateCollection_, pfCollection);
  if (!pfCollection.isValid()) {
    edm::LogError("PFClusterAnalyzer") << "invalid collection: PF clusters \n";
    return;
  }

  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByToken(pfJetsToken_, pfJets);
  if (!pfJets.isValid()) {
    edm::LogError("PFAnalyzer") << "invalid collection: PF jets \n";
    return;
  }

  if(!passesTriggerSelection(triggerResults, triggerNames, m_triggerOptions)){
    return;
  }


  for (reco::PFClusterCollection::const_iterator recoPF = pfCollection->begin(); recoPF != pfCollection->end();
       ++recoPF) {
    for (unsigned int j = 0; j < m_fullCutList.size(); j++) {
      std::vector<int> binNumbers = getPFBin(*recoPF, j);
      for(unsigned int k=0; k<binNumbers.size(); k++){
      int binNumber = binNumbers[k];
      if (binNumber < 0)
        continue;
      if (binNumber >= int(m_allSuffixes[j].size())) {
        continue;
      }
      std::string binString = m_allSuffixes[j][binNumber];

      // Eventually, we might want the hist name to include the cuts that we are applying,
      // so I am keepking it as a separate string for now, even though it is redundant.
      // Make plots of all observables
      for (unsigned int i = 0; i < m_observables.size(); i++) {
        std::string histName = Form("%s%s_%s", m_observableNames[i].c_str(), binString.c_str(), npvString.c_str());
        std::vector<double> valsX = m_funcMap[m_observableNames[i]](*recoPF);
        for(unsigned int m=0; m<valsX.size(); m++){
          map_of_MEs[m_directory + "/allPFC_" + histName]->Fill(valsX[m], eventWeight);
        }
      }
    }
    }
  }

}
