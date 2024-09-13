#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct wMuonFwdEfficiency {
    // Histogram registry: an object to hold your histograms
    HistogramRegistry histos{"histos", {},
    OutputObjHandlingPolicy::AnalysisObject};

    void init(InitContext const&)
    {
        // define axes you want to use
        const AxisSpec axisCounter{1, 0, +1, ""};
        const AxisSpec axisEta{10, -4.0, -2.5, "#eta"};
        const AxisSpec axisPt{15, 0.0, +80.0, "p_{T} (GeV/c)"};

        // create histograms
        histos.add("eventCounterReco", "eventCounterReco", kTH1F, {axisCounter});
        histos.add("eventCounterSim", "eventCounterSim", kTH1F, {axisCounter});
        histos.add("yPtHist", "yPtHist", kTH2F, {axisPt, axisEta});
        histos.add("yPtTruthHist", "yPtTruthHist", kTH2F, {axisPt, axisEta});
        histos.add("PtTruthHist", "PtTruthHist", kTH1F, {axisPt});
    }

    using muonTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

    void processReco(aod::Collision const& collision, muonTracks const& tracks, aod::McParticles const&)
    {
        histos.fill(HIST("eventCounterReco"), 0.5);
        for (auto& track : tracks) {
            if(track.has_mcParticle()){
                auto mcParticle = track.mcParticle();
                auto statusCode = mcParticle.getGenStatusCode();
                if (abs(mcParticle.pdgCode())==13 && statusCode >= 21 && statusCode <= 34) {
                    histos.fill(HIST("yPtHist"), mcParticle.pt(), mcParticle.eta());
                    LOGF(info, "Properties of reco W muon: pt %f; eta %f; collision %d", mcParticle.pt(), mcParticle.eta(), mcParticle.mcCollisionId());
                }
            }
        }
    }
    PROCESS_SWITCH(wMuonFwdEfficiency, processReco, "W->muon reconstructed information", true);

    //Filter<Tracks> etaFilter = track::eta < -2.5 && track::eta > -4;

    void processSim(aod::McParticles const& mcParticles)
    {
        histos.fill(HIST("eventCounterSim"), 0.5);
        for (const auto& mcParticle : mcParticles) {
            // investigate the hardest process
            auto statusCode = mcParticle.getGenStatusCode();
            if (statusCode >= 11 && statusCode <= 34) {
                LOGF(info, "Hardest process particle has pdg code %d, from collision %d", mcParticle.pdgCode(), mcParticle.mcCollisionId());
            }

            if (abs(mcParticle.pdgCode())==13 && statusCode >= 21 && statusCode <= 34) { // muon from the hardest process
                if (mcParticle.eta() > -4.0 && mcParticle.eta() < -2.5) { // forward region
                    histos.fill(HIST("yPtTruthHist"), mcParticle.pt(), mcParticle.eta());
                }
                histos.fill(HIST("PtTruthHist"), mcParticle.pt());

                // auto const& muMother = mcParticle.mothers_first_as<aod::McParticles>();
                // LOGF(info, "First mother: %d has pdg code %d", muMother.globalIndex(), muMother.pdgCode());
            }
        }
    }
    PROCESS_SWITCH(wMuonFwdEfficiency, processSim, "W->muon simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<wMuonFwdEfficiency>(cfgc)};
}