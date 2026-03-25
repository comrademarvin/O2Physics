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
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;

//constexpr static uint32_t EventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t MuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

struct processReducedMuons {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  AnalysisCompositeCut* muonCuts;
  AnalysisCompositeCut* eventCuts;

  void init(InitContext&) {
    // standalone muon
    histos.add("hPt", "Muon p_{T} distribution;p_{T} (GeV/c);Counts", HistType::kTH1F, {{100, 0, 40}});
    histos.add("hPt_before_cuts", "Muon p_{T} distribution;p_{T} (GeV/c);Counts", HistType::kTH1F, {{100, 0, 80}});
    histos.add("hEta", "Muon eta distribution;#eta;Counts", HistType::kTH1F, {{100, -5.0, -1.5}});
    histos.add("hEta_before_cuts", "Muon eta distribution;#eta;Counts", HistType::kTH1F, {{100, -5.0, -1.5}});
    histos.add("hPhi", "Muon Phi Angle Distribution;#phi;Counts", HistType::kTH1F, {{100, -4.0, 4.0}});
    histos.add("hChi", "Muon track #chi^{2} distribution;#chi^{2};Counts", HistType::kTH1F, {{100, 0, 20}});
    histos.add("hChiMatchMCHMID", "Muon track #chi^{2} of MCH+MID match distribution;#chi^{2} (MCH+MID);Counts", HistType::kTH1F, {{100, 0, 20}});
    histos.add("hDCAxy", "Muon DCA_{x}+DCA_{y} distribution;DCA_{x} (cm);DCA_{y};Counts", HistType::kTH2F, {{100, -10.0, 10.0},{100, -10.0, 10.0}});

    // standalone events
    histos.add("hPosZ", "Event Z-vertex distribution;PosZ (cm);Counts", HistType::kTH1F, {{100, -15.0, 15.0}});
    histos.add("hCentFT0C", "Run 3 centrality from FT0C multiplicity", HistType::kTH1F, {{100, 0.0, 100.0}});

    // associated muon+events
    histos.add("hPosZ_assoc", "Associated Event z-vertex position", HistType::kTH1F, {{100, -20.0, 20.0}});
    histos.add("hMuonsPerEvent", "Number of standalone muons per event", HistType::kTH1D, {{10, 0, 10}});
    histos.add("hMuonAmbiguity", "Number of standalone muons per recosntructed ID", HistType::kTH1D, {{5, 0, 5}});

    // Initialize VarManager
    VarManager::SetDefaultVarNames();

    // add muons cuts of interest
    muonCuts = new AnalysisCompositeCut(true);
    muonCuts->AddCut(o2::aod::dqcuts::GetCompositeCut("muonQualityCuts"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("MCHMID"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonLowPt"));
    //muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonHighPt2"));

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void processEvents(MyEvents const& events) {
    for (auto& event : events) {
      histos.fill(HIST("hPosZ"), event.posZ());
      histos.fill(HIST("hCentFT0C"), event.centFT0C());
    }
  }

  void processMuonStandalone(MyMuonTracks const& muons) {
    for (auto& muon : muons) {
      // Reset the variable values and fill current muon variables
      VarManager::ResetValues();
      VarManager::FillTrack<MuonFillMap>(muon);

      // Fill histograms before cuts
      histos.fill(HIST("hPt_before_cuts"), muon.pt());
      histos.fill(HIST("hEta_before_cuts"), muon.eta());

      // Apply the cut
      if (muonCuts->IsSelected(VarManager::fgValues)) {
        histos.fill(HIST("hPt"), muon.pt());
        histos.fill(HIST("hEta"), muon.eta());
        histos.fill(HIST("hPhi"), muon.phi());
        histos.fill(HIST("hChi"), muon.chi2());
        histos.fill(HIST("hChiMatchMCHMID"), muon.chi2MatchMCHMID());
        histos.fill(HIST("hDCAxy"), muon.fwdDcaX(), muon.fwdDcaY());
      }
    }
  };

  void processMuonAssoc(aod::ReducedMuonsAssoc const& assocs, MyEvents const& events, MyMuonTracks const& muons) {
    std::unordered_map<int, int> eventMuonCount; // Map to count muons per event
    std::unordered_map<int, int> muonAmbiguityCount; // Map to count ambiguous muons

    for (auto& assoc : assocs) {
      auto muon = assoc.template reducedmuon_as<MyMuonTracks>();
      auto event = assoc.template reducedevent_as<MyEvents>();

      // Reset the variable values and fill current muon variables
      VarManager::ResetValues();
      VarManager::FillTrack<MuonFillMap>(muon);

      // Apply the muon cuts
      if (muonCuts->IsSelected(VarManager::fgValues)) {
        if (muon.isAmbiguous() > 0) {
          std::cout << "Ambiguous track, skipping..." << std::endl;
          continue; // Skip ambiguous tracks
        }

        histos.fill(HIST("hPosZ_assoc"), event.posZ());

        // Increment the muon count for the current event
        // If the event does not exist in the map, it will be added with a default value of 0
        eventMuonCount[event.globalIndex()]++;
        muonAmbiguityCount[muon.globalIndex()]++;
      }
    }

    // Fill the histogram for the number of muons per event
    for (const auto& [eventID, muonCount] : eventMuonCount) {
      histos.fill(HIST("hMuonsPerEvent"), muonCount);
    }
    for (const auto& [muonID, muonCount] : muonAmbiguityCount) {
      histos.fill(HIST("hMuonAmbiguity"), muonCount);
    }
  };

  PROCESS_SWITCH(processReducedMuons, processEvents, "Run analysis on standalone events", true);
  PROCESS_SWITCH(processReducedMuons, processMuonStandalone, "Run analysis on stabdalone muon", true);
  PROCESS_SWITCH(processReducedMuons, processMuonAssoc, "Run event+muon analysis on associated muons", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<processReducedMuons>(cfgc)};
}