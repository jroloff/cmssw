// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

// user include files
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDigi/interface/HcalCalibrationEventTypes.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtfeWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/Provenance/interface/Provenance.h"

#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

namespace HcalMinbias {}

// constructors and destructor
class AnalyzerMinbias : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit AnalyzerMinbias(const edm::ParameterSet&);
  ~AnalyzerMinbias() override = default;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

private:
  void analyzeHcal(const HcalRespCorrs* myRecalib,
                   const HBHERecHitCollection& HithbheNS,
                   const HBHERecHitCollection& HithbheMB,
                   const HFRecHitCollection& HithfNS,
                   const HFRecHitCollection& HithfMB,
                   int algoBit,
                   bool fill);

  struct myInfo {
    double theMB0, theMB1, theMB2, theMB3, theMB4;
    double theNS0, theNS1, theNS2, theNS3, theNS4;
    double theDif0, theDif1, theDif2, runcheck;
    myInfo() {
      theMB0 = theMB1 = theMB2 = theMB3 = theMB4 = 0;
      theNS0 = theNS1 = theNS2 = theNS3 = theNS4 = 0;
      theDif0 = theDif1 = theDif2 = runcheck = 0;
    }
  };

  // ----------member data ---------------------------
  const std::string fOutputFileName;
  const bool runNZS_, theRecalib_, ignoreL1_;
  std::string hcalfile_;
  std::ofstream* myout_hcal;
  TFile* hOutputFile_;
  TTree* myTree_;
  TH1D *h_Noise[4], *h_Signal[4];
  double rnnum_;

  // Root tree members
  double rnnumber;
  int mydet, mysubd, depth, iphi, ieta, cells, trigbit;
  float phi, eta;
  float mom0_MB, mom1_MB, mom2_MB, mom3_MB, mom4_MB, occup;
  float mom0_Noise, mom1_Noise, mom2_Noise, mom3_Noise, mom4_Noise;
  float mom0_Diff, mom1_Diff, mom2_Diff, mom3_Diff, mom4_Diff;

  std::map<std::pair<int, HcalDetId>, myInfo> myMap_;
  const edm::EDGetTokenT<HBHERecHitCollection> tok_hbherecoMB_, tok_hbherecoNoise_;
  const edm::EDGetTokenT<HORecHitCollection> tok_horecoMB_, tok_horecoNoise_;
  const edm::EDGetTokenT<HFRecHitCollection> tok_hfrecoMB_, tok_hfrecoNoise_;
  const edm::EDGetTokenT<HBHERecHitCollection> tok_hbheNormal_;
  const edm::EDGetTokenT<L1GlobalTriggerObjectMapRecord> tok_hltL1GtMap_;
  const edm::ESGetToken<HcalRespCorrs, HcalRespCorrsRcd> tok_respCorr_;
};

AnalyzerMinbias::AnalyzerMinbias(const edm::ParameterSet& iConfig)
    : fOutputFileName(iConfig.getUntrackedParameter<std::string>("HistOutFile")),
      runNZS_(iConfig.getUntrackedParameter<bool>("RunNZS", true)),
      theRecalib_(iConfig.getParameter<bool>("Recalib")),
      ignoreL1_(iConfig.getUntrackedParameter<bool>("IgnoreL1", true)),
      tok_hbherecoMB_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInputMB"))),
      tok_hbherecoNoise_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInputNoise"))),
      tok_horecoMB_(consumes<HORecHitCollection>(iConfig.getParameter<edm::InputTag>("hoInputMB"))),
      tok_horecoNoise_(consumes<HORecHitCollection>(iConfig.getParameter<edm::InputTag>("hoInputNoise"))),
      tok_hfrecoMB_(consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("hfInputMB"))),
      tok_hfrecoNoise_(consumes<HFRecHitCollection>(iConfig.getParameter<edm::InputTag>("hfInputNoise"))),
      tok_hbheNormal_(consumes<HBHERecHitCollection>(edm::InputTag("hbhereco"))),
      tok_hltL1GtMap_(consumes<L1GlobalTriggerObjectMapRecord>(edm::InputTag("hltL1GtObjectMap"))),
      tok_respCorr_(esConsumes<HcalRespCorrs, HcalRespCorrsRcd>()) {
  usesResource(TFileService::kSharedResource);
  // get name of output file with histogramms
  // get token names of modules, producing object collections
}

void AnalyzerMinbias::beginJob() {
  std::string det[4] = {"HB", "HE", "HO", "HF"};
  char name[80], title[80];
  for (int subd = 0; subd < 4; ++subd) {
    sprintf(name, "Noise_%s", det[subd].c_str());
    sprintf(title, "Energy Distribution for Noise in %s", det[subd].c_str());
    h_Noise[subd] = new TH1D(name, title, 100, -10., 10.);
    sprintf(name, "Signal_%s", det[subd].c_str());
    sprintf(title, "Energy Distribution for Signal in %s", det[subd].c_str());
    h_Signal[subd] = new TH1D(name, title, 100, -10., 10.);
  }

  edm::Service<TFileService> fs;
  hOutputFile_ = new TFile(fOutputFileName.c_str(), "RECREATE");
  myTree_ = fs->make<TTree>("RecJet", "RecJet Tree");
  myTree_->Branch("mydet", &mydet, "mydet/I");
  myTree_->Branch("mysubd", &mysubd, "mysubd/I");
  myTree_->Branch("cells", &cells, "cells");
  myTree_->Branch("depth", &depth, "depth/I");
  myTree_->Branch("ieta", &ieta, "ieta/I");
  myTree_->Branch("iphi", &iphi, "iphi/I");
  myTree_->Branch("eta", &eta, "eta/F");
  myTree_->Branch("phi", &phi, "phi/F");
  myTree_->Branch("mom0_MB", &mom0_MB, "mom0_MB/F");
  myTree_->Branch("mom1_MB", &mom1_MB, "mom1_MB/F");
  myTree_->Branch("mom2_MB", &mom2_MB, "mom2_MB/F");
  myTree_->Branch("mom3_MB", &mom3_MB, "mom3_MB/F");
  myTree_->Branch("mom4_MB", &mom4_MB, "mom4_MB/F");
  myTree_->Branch("mom0_Noise", &mom0_Noise, "mom0_Noise/F");
  myTree_->Branch("mom1_Noise", &mom1_Noise, "mom1_Noise/F");
  myTree_->Branch("mom2_Noise", &mom2_Noise, "mom2_Noise/F");
  myTree_->Branch("mom3_Noise", &mom3_Noise, "mom3_Noise/F");
  myTree_->Branch("mom4_Noise", &mom4_Noise, "mom4_Noise/F");
  myTree_->Branch("mom0_Diff", &mom0_Diff, "mom0_Diff/F");
  myTree_->Branch("mom1_Diff", &mom1_Diff, "mom1_Diff/F");
  myTree_->Branch("mom2_Diff", &mom2_Diff, "mom2_Diff/F");
  myTree_->Branch("occup", &occup, "occup/F");
  myTree_->Branch("trigbit", &trigbit, "trigbit/I");
  myTree_->Branch("rnnumber", &rnnumber, "rnnumber/D");

  myMap_.clear();
}

//  EndJob
//
void AnalyzerMinbias::endJob() {
  int ii = 0;
  for (std::map<std::pair<int, HcalDetId>, myInfo>::const_iterator itr = myMap_.begin(); itr != myMap_.end(); ++itr) {
#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("AnalyzerMinimumBias") << "Fired trigger bit number " << itr->first.first;
#endif
    myInfo info = itr->second;
    if (info.theMB0 > 0) {
      mom0_MB = info.theMB0;
      mom1_MB = info.theMB1;
      mom2_MB = info.theMB2;
      mom3_MB = info.theMB3;
      mom4_MB = info.theMB4;
      mom0_Noise = info.theNS0;
      mom1_Noise = info.theNS1;
      mom2_Noise = info.theNS2;
      mom3_Noise = info.theNS3;
      mom4_Noise = info.theNS4;
      mom0_Diff = info.theDif0;
      mom1_Diff = info.theDif1;
      mom2_Diff = info.theDif2;
      rnnumber = info.runcheck;
      trigbit = itr->first.first;
      mysubd = itr->first.second.subdet();
      depth = itr->first.second.depth();
      ieta = itr->first.second.ieta();
      iphi = itr->first.second.iphi();
#ifdef EDM_ML_DEBUG
      edm::LogVerbatim("AnalyzerMinimumBias") << " Result=  " << trigbit << " " << mysubd << " " << ieta << " " << iphi
                                     << " mom0  " << mom0_MB << " mom1 " << mom1_MB << " mom2 " << mom2_MB << " mom3 "
                                     << mom3_MB << " mom4 " << mom4_MB << " mom0_Noise " << mom0_Noise << " mom1_Noise "
                                     << mom1_Noise << " mom2_Noise " << mom2_Noise << " mom3_Noise " << mom3_Noise
                                     << " mom4_Noise " << mom4_Noise << " mom0_Diff " << mom0_Diff << " mom1_Diff "
                                     << mom1_Diff << " mom2_Diff " << mom2_Diff;
#endif
      myTree_->Fill();
      ii++;
    }
  }
  cells = ii;
#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("AnalyzerMinimumBias") << "cells " << cells;
#endif
  hOutputFile_->Write();
  hOutputFile_->cd();
  myTree_->Write();
  for (int i = 0; i < 4; i++) {
    h_Noise[i]->Write();
    h_Signal[i]->Write();
  }
  hOutputFile_->Close();
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void AnalyzerMinbias::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  rnnum_ = (float)iEvent.run();
  const HcalRespCorrs* myRecalib = nullptr;
  if (theRecalib_) {
    myRecalib = &iSetup.getData(tok_respCorr_);
  }  // theRecalib

  const edm::Handle<HBHERecHitCollection> hbheNormal = iEvent.getHandle(tok_hbheNormal_);
  if (!hbheNormal.isValid()) {
    edm::LogVerbatim("AnalyzerMinimumBias") << " hbheNormal failed";
  } else {
    edm::LogVerbatim("AnalyzerMinimumBias") << " The size of the normal collection " << hbheNormal->size();
  }

  const edm::Handle<HBHERecHitCollection> hbheNS = iEvent.getHandle(tok_hbherecoNoise_);
  if (!hbheNS.isValid()) {
    edm::LogWarning("AnalyzerMinimumBias") << "HcalCalibAlgos: Error! can't get hbheNoise product!";
    return;
  }
  const HBHERecHitCollection HithbheNS = *(hbheNS.product());
  edm::LogVerbatim("AnalyzerMinimumBias") << "HBHE NS size of collection " << HithbheNS.size();
  if (runNZS_ && HithbheNS.size() != 5184) {
    edm::LogWarning("AnalyzerMinimumBias") << "HBHE NS problem " << rnnum_ << " size " << HithbheNS.size();
    return;
  }

  const edm::Handle<HBHERecHitCollection> hbheMB = iEvent.getHandle(tok_hbherecoMB_);
  if (!hbheMB.isValid()) {
    edm::LogWarning("AnalyzerMinimumBias") << "HcalCalibAlgos: Error! can't get hbhe product!";
    return;
  }
  const HBHERecHitCollection HithbheMB = *(hbheMB.product());
  edm::LogVerbatim("AnalyzerMinimumBias") << "HBHE MB size of collection " << HithbheMB.size();
  if (runNZS_ && HithbheMB.size() != 5184) {
    edm::LogWarning("AnalyzerMinimumBias") << "HBHE problem " << rnnum_ << " size " << HithbheMB.size();
    return;
  }

  const edm::Handle<HFRecHitCollection> hfNS = iEvent.getHandle(tok_hfrecoNoise_);
  if (!hfNS.isValid()) {
    edm::LogWarning("AnalyzerMinimumBias") << "HcalCalibAlgos: Error! can't get hfNoise product!";
    return;
  }
  const HFRecHitCollection HithfNS = *(hfNS.product());
  edm::LogVerbatim("AnalyzerMinimumBias") << "HF NS size of collection " << HithfNS.size();
  if (runNZS_ && HithfNS.size() != 1728) {
    edm::LogWarning("AnalyzerMinimumBias") << "HF NS problem " << rnnum_ << " size " << HithfNS.size();
    return;
  }

  const edm::Handle<HFRecHitCollection> hfMB = iEvent.getHandle(tok_hfrecoMB_);
  if (!hfMB.isValid()) {
    edm::LogWarning("AnalyzerMinimumBias") << "HcalCalibAlgos: Error! can't get hf product!";
    return;
  }
  const HFRecHitCollection HithfMB = *(hfMB.product());
  edm::LogVerbatim("AnalyzerMinimumBias") << "HF MB size of collection " << HithfMB.size();
  if (runNZS_ && HithfMB.size() != 1728) {
    edm::LogWarning("AnalyzerMinimumBias") << "HF problem " << rnnum_ << " size " << HithfMB.size();
    return;
  }

  if (ignoreL1_) {
    analyzeHcal(myRecalib, HithbheNS, HithbheMB, HithfNS, HithfMB, 1, true);
  } else {
    const edm::Handle<L1GlobalTriggerObjectMapRecord> gtObjectMapRecord = iEvent.getHandle(tok_hltL1GtMap_);
    if (gtObjectMapRecord.isValid()) {
      const std::vector<L1GlobalTriggerObjectMap>& objMapVec = gtObjectMapRecord->gtObjectMap();
      int ii(0);
      bool ok(false), fill(true);
      for (std::vector<L1GlobalTriggerObjectMap>::const_iterator itMap = objMapVec.begin(); itMap != objMapVec.end();
           ++itMap, ++ii) {
        bool resultGt = (*itMap).algoGtlResult();
        if (resultGt == 1) {
          ok = true;
          int algoBit = (*itMap).algoBitNumber();
          analyzeHcal(myRecalib, HithbheNS, HithbheMB, HithfNS, HithfMB, algoBit, fill);
          fill = false;
          std::string algoNameStr = (*itMap).algoName();
#ifdef EDM_ML_DEBUG
          edm::LogVerbatim("AnalyzerMinimumBias") << "Trigger[" << ii << "] " << algoNameStr << " bit " << algoBit << " entered";
#endif
        }
      }
      if (!ok)
        edm::LogVerbatim("AnalyzerMinimumBias") << "No passed L1 Triggers";
    }
  }
}

void AnalyzerMinbias::analyzeHcal(const HcalRespCorrs* myRecalib,
                                  const HBHERecHitCollection& HithbheNS,
                                  const HBHERecHitCollection& HithbheMB,
                                  const HFRecHitCollection& HithfNS,
                                  const HFRecHitCollection& HithfMB,
                                  int algoBit,
                                  bool fill) {
  // Noise part for HB HE
  std::map<std::pair<int, HcalDetId>, myInfo> tmpMap;
  tmpMap.clear();

  for (HBHERecHitCollection::const_iterator hbheItr = HithbheNS.begin(); hbheItr != HithbheNS.end(); hbheItr++) {
    // Recalibration of energy
    float icalconst = 1.;
    DetId mydetid = hbheItr->id().rawId();
    if (theRecalib_)
      icalconst = myRecalib->getValues(mydetid)->getValue();

    HBHERecHit aHit(hbheItr->id(), hbheItr->energy() * icalconst, hbheItr->time());
    double energyhit = aHit.energy();

    DetId id = (*hbheItr).detid();
    HcalDetId hid = HcalDetId(id);
    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    if (itr1 == myMap_.end()) {
      myInfo info;
      myMap_[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr1->second.theNS0++;
    itr1->second.theNS1 += energyhit;
    itr1->second.theNS2 += (energyhit * energyhit);
    itr1->second.theNS3 += (energyhit * energyhit * energyhit);
    itr1->second.theNS4 += (energyhit * energyhit * energyhit * energyhit);
    itr1->second.runcheck = rnnum_;
    if (fill)
      h_Noise[hid.subdet() - 1]->Fill(energyhit);

    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));
    if (itr2 == tmpMap.end()) {
      myInfo info;
      tmpMap[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr2->second.theNS0++;
    itr2->second.theNS1 += energyhit;
    itr2->second.theNS2 += (energyhit * energyhit);
    itr2->second.theNS3 += (energyhit * energyhit * energyhit);
    itr2->second.theNS4 += (energyhit * energyhit * energyhit * energyhit);
    itr2->second.runcheck = rnnum_;

  }  // HBHE_NS

  // Signal part for HB HE

  for (HBHERecHitCollection::const_iterator hbheItr = HithbheMB.begin(); hbheItr != HithbheMB.end(); hbheItr++) {
    // Recalibration of energy
    float icalconst = 1.;
    DetId mydetid = hbheItr->id().rawId();
    if (theRecalib_)
      icalconst = myRecalib->getValues(mydetid)->getValue();

    HBHERecHit aHit(hbheItr->id(), hbheItr->energy() * icalconst, hbheItr->time());
    double energyhit = aHit.energy();

    DetId id = (*hbheItr).detid();
    HcalDetId hid = HcalDetId(id);

    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));

    if (itr1 == myMap_.end()) {
      myInfo info;
      myMap_[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr1->second.theMB0++;
    itr1->second.theDif0 = 0;
    itr1->second.theMB1 += energyhit;
    itr1->second.theMB2 += (energyhit * energyhit);
    itr1->second.theMB3 += (energyhit * energyhit * energyhit);
    itr1->second.theMB4 += (energyhit * energyhit * energyhit * energyhit);
    itr1->second.runcheck = rnnum_;
    float mydiff = 0.0;
    if (itr2 != tmpMap.end()) {
      mydiff = energyhit - (itr2->second.theNS1);
      itr1->second.theDif0++;
      itr1->second.theDif1 += mydiff;
      itr1->second.theDif2 += (mydiff * mydiff);
      if (fill)
        h_Signal[hid.subdet() - 1]->Fill(mydiff);
    }
  }  // HBHE_MB

  // HF

  for (HFRecHitCollection::const_iterator hbheItr = HithfNS.begin(); hbheItr != HithfNS.end(); hbheItr++) {
    // Recalibration of energy
    float icalconst = 1.;
    DetId mydetid = hbheItr->id().rawId();
    if (theRecalib_)
      icalconst = myRecalib->getValues(mydetid)->getValue();

    HFRecHit aHit(hbheItr->id(), hbheItr->energy() * icalconst, hbheItr->time());
    double energyhit = aHit.energy();
    // Remove PMT hits
    if (fabs(energyhit) > 40.)
      continue;
    DetId id = (*hbheItr).detid();
    HcalDetId hid = HcalDetId(id);

    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));

    if (itr1 == myMap_.end()) {
      myInfo info;
      myMap_[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr1->second.theNS0++;
    itr1->second.theNS1 += energyhit;
    itr1->second.theNS2 += (energyhit * energyhit);
    itr1->second.theNS3 += (energyhit * energyhit * energyhit);
    itr1->second.theNS4 += (energyhit * energyhit * energyhit * energyhit);
    itr1->second.runcheck = rnnum_;
    if (fill)
      h_Noise[hid.subdet() - 1]->Fill(energyhit);

    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));
    if (itr2 == tmpMap.end()) {
      myInfo info;
      tmpMap[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr2->second.theNS0++;
    itr2->second.theNS1 += energyhit;
    itr2->second.theNS2 += (energyhit * energyhit);
    itr2->second.theNS3 += (energyhit * energyhit * energyhit);
    itr2->second.theNS4 += (energyhit * energyhit * energyhit * energyhit);
    itr2->second.runcheck = rnnum_;

  }  // HF_NS

  // Signal part for HF

  for (HFRecHitCollection::const_iterator hbheItr = HithfMB.begin(); hbheItr != HithfMB.end(); hbheItr++) {
    // Recalibration of energy
    float icalconst = 1.;
    DetId mydetid = hbheItr->id().rawId();
    if (theRecalib_)
      icalconst = myRecalib->getValues(mydetid)->getValue();
    HFRecHit aHit(hbheItr->id(), hbheItr->energy() * icalconst, hbheItr->time());

    double energyhit = aHit.energy();
    // Remove PMT hits
    if (fabs(energyhit) > 40.)
      continue;

    DetId id = (*hbheItr).detid();
    HcalDetId hid = HcalDetId(id);

    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    std::map<std::pair<int, HcalDetId>, myInfo>::iterator itr2 = tmpMap.find(std::pair<int, HcalDetId>(algoBit, hid));

    if (itr1 == myMap_.end()) {
      myInfo info;
      myMap_[std::pair<int, HcalDetId>(algoBit, hid)] = info;
      itr1 = myMap_.find(std::pair<int, HcalDetId>(algoBit, hid));
    }
    itr1->second.theMB0++;
    itr1->second.theDif0 = 0;
    itr1->second.theMB1 += energyhit;
    itr1->second.theMB2 += (energyhit * energyhit);
    itr1->second.theMB3 += (energyhit * energyhit * energyhit);
    itr1->second.theMB4 += (energyhit * energyhit * energyhit * energyhit);
    itr1->second.runcheck = rnnum_;
    float mydiff = 0.0;
    if (itr2 != tmpMap.end()) {
      mydiff = energyhit - (itr2->second.theNS1);
      itr1->second.theDif0++;
      itr1->second.theDif1 += mydiff;
      itr1->second.theDif2 += (mydiff * mydiff);
      if (fill)
        h_Signal[hid.subdet() - 1]->Fill(mydiff);
    }
  }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AnalyzerMinbias);
