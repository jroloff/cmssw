#include "DQMOffline/ParticleFlow/interface/DQMExample_Step3.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

//
// -------------------------------------- Constructor
// --------------------------------------------
//
DQMExample_Step3::DQMExample_Step3(const edm::ParameterSet &ps) {
  edm::LogInfo("DQMExample_Step3") << "Constructor  DQMExample_Step3::DQMExample_Step3 " << std::endl;

  // Get parameters from configuration file
  numMonitorName_ = ps.getParameter<std::string>("numMonitorName");
  denMonitorName_ = ps.getParameter<std::string>("denMonitorName");
}

//
// -- Destructor
//
DQMExample_Step3::~DQMExample_Step3() {
  edm::LogInfo("DQMExample_Step3") << "Destructor DQMExample_Step3::~DQMExample_Step3 " << std::endl;
}

//
// -------------------------------------- beginJob
// --------------------------------------------
//
void DQMExample_Step3::beginJob() { edm::LogInfo("DQMExample_Step3") << "DQMExample_Step3::beginJob " << std::endl; }
//
// -------------------------------------- get and book in the endJob
// --------------------------------------------
//
void DQMExample_Step3::dqmEndJob(DQMStore::IBooker &ibooker_, DQMStore::IGetter &igetter_) {
  // create and cd into new folder
  ibooker_.setCurrentFolder("What_I_do_in_the_client/Ratio");

  // get available histograms
  MonitorElement *numerator = igetter_.get(numMonitorName_);
  MonitorElement *denominator = igetter_.get(denMonitorName_);

  if (!numerator || !denominator) {
    edm::LogError("DQMExample_Step3") << "MEs not found!" << std::endl;
    return;
  }

  // book new histogram
  h_ptRatio = ibooker_.book1D("ptRatio", "pt ratio pf matched objects", 50, 0., 100.);
  h_ptRatio->setAxisTitle("pt [GeV]");

  for (int iBin = 1; iBin < numerator->getNbinsX(); ++iBin) {
    if (denominator->getBinContent(iBin) == 0)
      h_ptRatio->setBinContent(iBin, 0.);
    else
      h_ptRatio->setBinContent(iBin, numerator->getBinContent(iBin) / denominator->getBinContent(iBin));
  }
}

//
// -------------------------------------- get in the endLumi if needed
// --------------------------------------------
//
void DQMExample_Step3::dqmEndLuminosityBlock(DQMStore::IBooker &ibooker_,
                                             DQMStore::IGetter &igetter_,
                                             edm::LuminosityBlock const &iLumi,
                                             edm::EventSetup const &iSetup) {
  edm::LogInfo("DQMExample_Step3") << "DQMExample_Step3::endLumi " << std::endl;
}
