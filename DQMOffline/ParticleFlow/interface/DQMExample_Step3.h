#ifndef DQMExample_Step3_H
#define DQMExample_Step3_H

// Framework
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

// DQM
#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

class DQMExample_Step3 : public DQMEDHarvester {
public:
  DQMExample_Step3(const edm::ParameterSet &ps);
  ~DQMExample_Step3() override;

protected:
  void beginJob() override;
  void dqmEndLuminosityBlock(DQMStore::IBooker &,
                             DQMStore::IGetter &,
                             edm::LuminosityBlock const &,
                             edm::EventSetup const &) override;  // performed in the endLumi
  void dqmEndJob(DQMStore::IBooker &,
                 DQMStore::IGetter &) override;  // performed in the endJob

private:
  // private variables

  // variables from config file
  std::string numMonitorName_;
  std::string denMonitorName_;

  // Histograms
  MonitorElement *h_ptRatio;
};

#endif
