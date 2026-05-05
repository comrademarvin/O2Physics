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

#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define datamodel by joining relevant tables (with MC labels)
using MyMuonTracksLabelled = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;

// define datamodel bit mask for the variable manager to know which variables to fill for the muon objects
constexpr static uint32_t MuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct processReducedMuonsMC {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  AnalysisCompositeCut* muonCuts;

  void init(InitContext&) {
    // define histograms of interest
    histos.add("hPt", "Muon pT", HistType::kTH1F, {{100, 0, 100}});
    histos.add("hEta", "Muon eta", HistType::kTH1F, {{100, -5.0, -1.5}});
    histos.add("DCAxy", "Muon DCAxy", HistType::kTH1F, {{100, 0, 10}});

    // Initialize variable manager for applying cuts (DQ framework thing)
    VarManager::SetDefaultVarNames();

    // add muon cuts of interest (re-applying some to be safe)
    muonCuts = new AnalysisCompositeCut(true);
    muonCuts->AddCut(o2::aod::dqcuts::GetCompositeCut("muonQualityCuts"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("matchedMchMid"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonHighPt5"));

    // provides the list of required variables so that the variable manager knows what to fill
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
  }

  void process(MyMuonTracksLabelled const& muons, aod::ReducedMCTracks const& mcTracks)
  {
    for (auto& muon : muons) {
      // Reset the variable values and fill current muon variables
      VarManager::ResetValues();
      VarManager::FillTrack<MuonFillMap>(muon);

      // check with has_reducedMCTrack, because the label can be missing
      if (!muon.has_reducedMCTrack()) {
        continue;
      }
      auto mcMuon = mcTracks.rawIteratorAt(muon.reducedMCTrackId()); // access the associated MC track

      // mcReducedFlags is the MC-signal bit mask produced by tableMakerMC_withAssoc, when passing the cfgMCsignals in the config
      const uint16_t mcFlags = muon.mcReducedFlags();
      // for now, assume only configured exactly ONE MC signal in tableMakerMC_withAssoc config (more can be added later)
      // in that case, bit 0 encodes the decision for that signal
      if ((mcFlags & static_cast<uint16_t>(0x1u)) == 0u) {
        continue;
      }

      // re-apply muon cuts to be on the safe side
      if (!muonCuts->IsSelected(VarManager::fgValues)) {
        continue;
      }

      // fill histograms
      histos.fill(HIST("hPt"), muon.pt());
      histos.fill(HIST("hEta"), muon.eta());

      // compute and fill DCAxy (offset)
      float DCAxy = std::sqrt(std::pow(muon.fwdDcaX(), 2) + std::pow(muon.fwdDcaY(), 2));
      histos.fill(HIST("DCAxy"), DCAxy);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<processReducedMuonsMC>(cfgc)};
}