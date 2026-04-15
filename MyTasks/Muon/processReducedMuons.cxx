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
  AnalysisCut* standaloneMuonCut;

  void init(InitContext&) {
    // standalone muon
    histos.add("hPt", "Muon p_{T} distribution after cuts;p_{T} (GeV/c);Counts", HistType::kTH1F, {{100, 0, 20}}, true);
    histos.add("hPt_before_cuts", "Muon p_{T} distribution before cuts;p_{T} (GeV/c);Counts", HistType::kTH1F, {{100, 0, 20}}, true);
    histos.add("hEta", "Muon eta distribution after cuts;#eta;Counts", HistType::kTH1F, {{100, -5.0, -1.5}}, true);
    histos.add("hEta_before_cuts", "Muon eta distribution before cuts;#eta;Counts", HistType::kTH1F, {{100, -5.0, -1.5}}, true);
    histos.add("hPhi", "Muon phi angle distribution;#phi;Counts", HistType::kTH1F, {{100, -4.0, 4.0}}, true);
    histos.add("hChi", "Muon track #chi^{2} distribution;#chi^{2};Counts", HistType::kTH1F, {{100, 0, 20}}, true);
    histos.add("hChiMatchMCHMID", "Muon track #chi^{2} of MCH+MID match distribution;#chi^{2} (MCH+MID);Counts", HistType::kTH1F, {{100, 0, 20}}, true);
    histos.add("hDCAxy", "Muon DCA_{X}+DCA_{Y} distribution;DCA_{X} (cm);DCA_{Y} (cm);Counts", HistType::kTH2F, {{100, -50.0, 50.0},{100, -50.0, 50.0}});

    // standalone events
    histos.add("hPosZ", "Event Z-vertex distribution;Vtx Z (cm);Counts", HistType::kTH1F, {{100, -20.0, 20.0}}, true);
    histos.add("hCentFT0C", "Run 3 centrality from FT0C multiplicity", HistType::kTH1F, {{100, 0.0, 100.0}}, true);
    histos.add("hMultTPC", "Track multiplicity from TPC;Multiplicity;Counts", HistType::kTH1D, {{100, 0, 15000}}, true);

    // associated muon+events
    histos.add("hPosZ_assoc", "Associated muon events Z-vertex position;Vtx Z (cm);Counts", HistType::kTH1F, {{100, -20.0, 20.0}}, true);
    histos.add("hMuonsPerEvent", "Number of standalone muons per event;Muons per Event;Counts", HistType::kTH1D, {{10, 0, 10}}, true);
    histos.add("hMuonAmbiguity", "Number of events per standalone muon;Events per Muon;Counts", HistType::kTH1D, {{4, 0, 4}});

    // Initialize VarManager
    VarManager::SetDefaultVarNames();

    // add muons cuts of interest
    muonCuts = new AnalysisCompositeCut(true);
    muonCuts->AddCut(o2::aod::dqcuts::GetCompositeCut("muonQualityCuts"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("MCHMID"));
    muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonLowPt"));
    //muonCuts->AddCut(o2::aod::dqcuts::GetAnalysisCut("muonHighPt2"));

    standaloneMuonCut = o2::aod::dqcuts::GetAnalysisCut("MCHMID");

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void processEvents(MyEvents const& events) {
    for (auto& event : events) {
      histos.fill(HIST("hPosZ"), event.posZ());
      histos.fill(HIST("hCentFT0C"), event.centFT0C());
      histos.fill(HIST("hMultTPC"), event.multTPC());
    }
  }

  void processMuonStandalone(MyMuonTracks const& muons) {
    for (auto& muon : muons) {
      // Reset the variable values and fill current muon variables
      VarManager::ResetValues();
      VarManager::FillTrack<MuonFillMap>(muon);

      // Fill standalone muon histograms before cuts
      if (standaloneMuonCut->IsSelected(VarManager::fgValues)) {
        histos.fill(HIST("hPt_before_cuts"), muon.pt());
        histos.fill(HIST("hEta_before_cuts"), muon.eta());
      }

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