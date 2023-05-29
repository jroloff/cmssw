#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimTracker/TrackTriggerAssociation/interface/StubAssociation.h"
#include "L1Trigger/TrackTrigger/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/LayerEncoding.h"
#include "L1Trigger/TrackerTFP/interface/KalmanFilterFormats.h"

#include <TProfile.h>
#include <TH1F.h>

#include <vector>
#include <deque>
#include <set>
#include <cmath>
#include <numeric>
#include <sstream>

using namespace std;
using namespace edm;
using namespace tt;
using namespace trackerTFP;

namespace trklet {

  /*! \class  trklet::AnalyzerKFout
   *  \brief  Class to analyze hardware like structured Track Collection generated by Kalman Filter output module
   *  \author Thomas Schuh
   *  \date   2021, Aug
   */
  class AnalyzerKFout : public one::EDAnalyzer<one::WatchRuns, one::SharedResources> {
  public:
    AnalyzerKFout(const ParameterSet& iConfig);
    void beginJob() override {}
    void beginRun(const Run& iEvent, const EventSetup& iSetup) override;
    void analyze(const Event& iEvent, const EventSetup& iSetup) override;
    void endRun(const Run& iEvent, const EventSetup& iSetup) override {}
    void endJob() override;

  private:
    //
    void associate(const set<TTTrackRef>& ttTracks, const StubAssociation* ass, set<TPPtr>& tps, int& sum) const;

    // ED input token of accepted Stubs
    EDGetTokenT<StreamsTrack> edGetTokenAccepted_;
    // ED input token of lost Tracks
    EDGetTokenT<StreamsTrack> edGetTokenLost_;
    // ED input token of TTStubRef to TPPtr association for tracking efficiency
    EDGetTokenT<StubAssociation> edGetTokenSelection_;
    // ED input token of TTStubRef to recontructable TPPtr association
    EDGetTokenT<StubAssociation> edGetTokenReconstructable_;
    // Setup token
    ESGetToken<Setup, SetupRcd> esGetTokenSetup_;
    // DataFormats token
    ESGetToken<DataFormats, DataFormatsRcd> esGetTokenDataFormats_;
    // stores, calculates and provides run-time constants
    const Setup* setup_ = nullptr;
    // enables analyze of TPs
    bool useMCTruth_;
    //
    int nEvents_ = 0;

    // Histograms

    TProfile* prof_;
    TProfile* profChannel_;
    TH1F* hisChannel_;

    // printout
    stringstream log_;
  };

  AnalyzerKFout::AnalyzerKFout(const ParameterSet& iConfig) : useMCTruth_(iConfig.getParameter<bool>("UseMCTruth")) {
    usesResource("TFileService");
    // book in- and output ED products
    const string& label = iConfig.getParameter<string>("LabelKFout");
    const string& branchAccepted = iConfig.getParameter<string>("BranchAcceptedTracks");
    const string& branchLost = iConfig.getParameter<string>("BranchLostTracks");
    edGetTokenAccepted_ = consumes<StreamsTrack>(InputTag(label, branchAccepted));
    edGetTokenLost_ = consumes<StreamsTrack>(InputTag(label, branchLost));
    if (useMCTruth_) {
      const auto& inputTagSelecttion = iConfig.getParameter<InputTag>("InputTagSelection");
      const auto& inputTagReconstructable = iConfig.getParameter<InputTag>("InputTagReconstructable");
      edGetTokenSelection_ = consumes<StubAssociation>(inputTagSelecttion);
      edGetTokenReconstructable_ = consumes<StubAssociation>(inputTagReconstructable);
    }
    // book ES products
    esGetTokenSetup_ = esConsumes<Setup, SetupRcd, Transition::BeginRun>();
    // log config
    log_.setf(ios::fixed, ios::floatfield);
    log_.precision(4);
  }

  void AnalyzerKFout::beginRun(const Run& iEvent, const EventSetup& iSetup) {
    // helper class to store configurations
    setup_ = &iSetup.getData(esGetTokenSetup_);
    // book histograms
    Service<TFileService> fs;
    TFileDirectory dir;
    dir = fs->mkdir("KFout");
    prof_ = dir.make<TProfile>("Counts", ";", 9, 0.5, 9.5);
    prof_->GetXaxis()->SetBinLabel(1, "Stubs");
    prof_->GetXaxis()->SetBinLabel(2, "Tracks");
    prof_->GetXaxis()->SetBinLabel(3, "Lost Tracks");
    prof_->GetXaxis()->SetBinLabel(4, "Matched Tracks");
    prof_->GetXaxis()->SetBinLabel(5, "All Tracks");
    prof_->GetXaxis()->SetBinLabel(6, "Found TPs");
    prof_->GetXaxis()->SetBinLabel(7, "Found selected TPs");
    prof_->GetXaxis()->SetBinLabel(8, "Lost TPs");
    prof_->GetXaxis()->SetBinLabel(9, "All TPs");
    // channel occupancy
    constexpr int maxOcc = 180;
    const int numChannels = setup_->tfpNumChannel();
    hisChannel_ = dir.make<TH1F>("His Channel Occupancy", ";", maxOcc, -.5, maxOcc - .5);
    profChannel_ = dir.make<TProfile>("Prof Channel Occupancy", ";", numChannels, -.5, numChannels - .5);
  }

  void AnalyzerKFout::analyze(const Event& iEvent, const EventSetup& iSetup) {
    // read in kf products
    Handle<StreamsTrack> handleAccepted;
    iEvent.getByToken<StreamsTrack>(edGetTokenAccepted_, handleAccepted);
    Handle<StreamsTrack> handleLost;
    iEvent.getByToken<StreamsTrack>(edGetTokenLost_, handleLost);
    // read in MCTruth
    const StubAssociation* selection = nullptr;
    const StubAssociation* reconstructable = nullptr;
    if (useMCTruth_) {
      Handle<StubAssociation> handleSelection;
      iEvent.getByToken<StubAssociation>(edGetTokenSelection_, handleSelection);
      selection = handleSelection.product();
      prof_->Fill(9, selection->numTPs());
      Handle<StubAssociation> handleReconstructable;
      iEvent.getByToken<StubAssociation>(edGetTokenReconstructable_, handleReconstructable);
      reconstructable = handleReconstructable.product();
    }
    // analyze kfout products and associate found tracks with reconstrucable TrackingParticles
    set<TPPtr> tpPtrs;
    set<TPPtr> tpPtrsSelection;
    set<TPPtr> tpPtrsLost;
    int allMatched(0);
    int allTracks(0);
    for (int region = 0; region < setup_->numRegions(); region++) {
      int nStubsRegion(0);
      int nTracksRegion(0);
      int nLostRegion(0);
      for (int channel = 0; channel < setup_->tfpNumChannel(); channel++) {
        const int index = region * setup_->tfpNumChannel() + channel;
        const StreamTrack& accepted = handleAccepted->at(index);
        const StreamTrack& lost = handleLost->at(index);
        hisChannel_->Fill(accepted.size());
        profChannel_->Fill(channel, accepted.size());
        set<TTTrackRef> tracks;
        for (const FrameTrack& frame : accepted)
          if (frame.first.isNonnull())
            tracks.insert(frame.first);
        nTracksRegion += tracks.size();
        nStubsRegion += accumulate(tracks.begin(), tracks.end(), 0, [](int& sum, const TTTrackRef& ttTrackRef) {
          return sum += (int)ttTrackRef->getStubRefs().size();
        });
        set<TTTrackRef> tracksLost;
        for (const FrameTrack& frame : lost)
          if (frame.first.isNonnull())
            tracksLost.insert(frame.first);
        nLostRegion += tracksLost.size();
        allTracks += tracks.size();
        int tmp(0);
        associate(tracks, selection, tpPtrsSelection, tmp);
        associate(tracksLost, selection, tpPtrsLost, tmp);
        associate(tracks, reconstructable, tpPtrs, allMatched);
      }
      prof_->Fill(1, nStubsRegion);
      prof_->Fill(2, nTracksRegion);
      prof_->Fill(3, nLostRegion);
    }
    deque<TPPtr> tpPtrsRealLost;
    set_difference(tpPtrsLost.begin(), tpPtrsLost.end(), tpPtrs.begin(), tpPtrs.end(), back_inserter(tpPtrsRealLost));
    prof_->Fill(4, allMatched);
    prof_->Fill(5, allTracks);
    prof_->Fill(6, tpPtrs.size());
    prof_->Fill(7, tpPtrsSelection.size());
    prof_->Fill(8, tpPtrsRealLost.size());
    nEvents_++;
  }

  void AnalyzerKFout::endJob() {
    if (nEvents_ == 0)
      return;
    // printout SF summary
    const double totalTPs = prof_->GetBinContent(9);
    const double numStubs = prof_->GetBinContent(1);
    const double numTracks = prof_->GetBinContent(2);
    const double numTracksLost = prof_->GetBinContent(3);
    const double totalTracks = prof_->GetBinContent(5);
    const double numTracksMatched = prof_->GetBinContent(4);
    const double numTPsAll = prof_->GetBinContent(6);
    const double numTPsEff = prof_->GetBinContent(7);
    const double numTPsLost = prof_->GetBinContent(8);
    const double errStubs = prof_->GetBinError(1);
    const double errTracks = prof_->GetBinError(2);
    const double errTracksLost = prof_->GetBinError(3);
    const double fracFake = (totalTracks - numTracksMatched) / totalTracks;
    const double fracDup = (numTracksMatched - numTPsAll) / totalTracks;
    const double eff = numTPsEff / totalTPs;
    const double errEff = sqrt(eff * (1. - eff) / totalTPs / nEvents_);
    const double effLoss = numTPsLost / totalTPs;
    const double errEffLoss = sqrt(effLoss * (1. - effLoss) / totalTPs / nEvents_);
    const vector<double> nums = {numStubs, numTracks, numTracksLost};
    const vector<double> errs = {errStubs, errTracks, errTracksLost};
    const int wNums = ceil(log10(*max_element(nums.begin(), nums.end()))) + 5;
    const int wErrs = ceil(log10(*max_element(errs.begin(), errs.end()))) + 5;
    log_ << "                        KFout SUMMARY                        " << endl;
    //log_ << "number of stubs       per TFP = " << setw(wNums) << numStubs << " +- " << setw(wErrs) << errStubs << endl;
    log_ << "number of tracks      per TFP = " << setw(wNums) << numTracks << " +- " << setw(wErrs) << errTracks
         << endl;
    log_ << "number of lost tracks per TFP = " << setw(wNums) << numTracksLost << " +- " << setw(wErrs) << errTracksLost
         << endl;
    log_ << "          tracking efficiency = " << setw(wNums) << eff << " +- " << setw(wErrs) << errEff << endl;
    log_ << "     lost tracking efficiency = " << setw(wNums) << effLoss << " +- " << setw(wErrs) << errEffLoss << endl;
    log_ << "                    fake rate = " << setw(wNums) << fracFake << endl;
    log_ << "               duplicate rate = " << setw(wNums) << fracDup << endl;
    log_ << "=============================================================";
    LogPrint("L1Trigger/TrackerTFP") << log_.str();
  }

  //
  void AnalyzerKFout::associate(const set<TTTrackRef>& ttTracks,
                                const StubAssociation* ass,
                                set<TPPtr>& tps,
                                int& sum) const {
    for (const TTTrackRef& ttTrack : ttTracks) {
      const vector<TTStubRef>& ttStubRefs = ttTrack->getStubRefs();
      const vector<TPPtr>& tpPtrs = ass->associate(ttStubRefs);
      if (tpPtrs.empty())
        continue;
      sum++;
      copy(tpPtrs.begin(), tpPtrs.end(), inserter(tps, tps.begin()));
    }
  }

}  // namespace trklet

DEFINE_FWK_MODULE(trklet::AnalyzerKFout);