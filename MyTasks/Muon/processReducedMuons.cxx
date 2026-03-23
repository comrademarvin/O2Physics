#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Configurable.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;

struct processReducedMuons {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  AnalysisCompositeCut* muonCuts;

  void init(InitContext&) {
    histos.add("hPt", "Muon pT", HistType::kTH1F, {{100, 0, 100}});
    histos.add("hEta", "Muon eta", HistType::kTH1F, {{100, -5.0, -1.5}});
    histos.add("hPosZ", "Event z-vertex position", HistType::kTH1F, {{100, -20.0, 20.0}});

    // add muons cuts of interest
    muonCuts = new AnalysisCompositeCut(true);
    muonCuts->AddCut(o2::aod::dqcuts::GetCompositeCut("muonQualityCuts"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("MCHMID"));
    //muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonHighPt2"));

    // Initialize VarManager
    VarManager::SetDefaultVarNames();
  }

  void processMuonStandalone(MyMuonTracks const& muons) {
    for (auto& muon : muons) {
      // Reset the variable values
      VarManager::ResetValues();

      // Fill the variables for the current muon
      VarManager::fgValues[VarManager::kPt] = muon.pt();
      VarManager::fgValues[VarManager::kEta] = muon.eta();
      VarManager::fgValues[VarManager::kMuonRAtAbsorberEnd] = muon.rAtAbsorberEnd();
      VarManager::fgValues[VarManager::kMuonPDca] = muon.pDca();
      VarManager::fgValues[VarManager::kMuonChi2] = muon.chi2();
      VarManager::fgValues[VarManager::kMuonChi2MatchMCHMID] = muon.chi2MatchMCHMID();
      VarManager::fgValues[VarManager::kMuonTrackType] = muon.trackType();

      // Apply the cut
      if (muonCuts->IsSelected(VarManager::fgValues)) {
        histos.fill(HIST("hPt"), muon.pt());
        histos.fill(HIST("hEta"), muon.eta());
      }
    }
  };

  void processMuonAssoc(aod::ReducedMuonsAssoc const& assocs, aod::ReducedEvents const& events, MyMuonTracks const& muons) {
    for (auto& assoc : assocs) {
      //auto muon = muons[assoc.reducedmuonId()];
      auto event = assoc.template reducedevent_as<aod::ReducedEvents>();

      histos.fill(HIST("hPosZ"), event.posZ());
    }
  };

  PROCESS_SWITCH(processReducedMuons, processMuonStandalone, "Run muon selection on stabdalone muon", false);
  PROCESS_SWITCH(processReducedMuons, processMuonAssoc, "Run event+muon on associated muons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<processReducedMuons>(cfgc)};
}