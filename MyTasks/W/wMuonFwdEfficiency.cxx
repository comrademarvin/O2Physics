#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <iostream>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct wMuonFwdEfficiency {
    // Histogram registry: an object to hold your histograms
    HistogramRegistry histos{"histos", {},
    OutputObjHandlingPolicy::AnalysisObject};

    // store id of selected tracks to account for ambiguous tracks
    std::vector<int64_t> selectedTracksID;

    void init(InitContext const&)
    {
        // define axes you want to use
        const AxisSpec axisCounter{1, 0, +1, ""};
        const AxisSpec axisEta{10, -4.0, -2.5, "#eta"};
        const AxisSpec axisPt{20, 0.0, +80.0, "p_{T} (GeV/c)"};
        const AxisSpec axisDeltaPt{24, 0.0, +6.0, "|p_{T}^{true} - p_{T}^{reco}|(GeV/c)"};

        // create histograms
        histos.add("eventCounterReco", "eventCounterReco", kTH1F, {axisCounter});
        histos.add("eventCounterSim", "eventCounterSim", kTH1F, {axisCounter});
        histos.add("yPtRecoHist", "yPtRecoHist", kTH2F, {axisPt, axisEta});
        histos.add("yPtTruthHist", "yPtTruthHist", kTH2F, {axisPt, axisEta});
        histos.add("PtRecoHist", "PtRecoHist", kTH1F, {axisPt});
        histos.add("PtTruthHist", "PtTruthHist", kTH1F, {axisPt});
        histos.add("PtResolution", "PtResolution", kTH1F, {axisDeltaPt});
    }

    using muonTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

    void processReco(aod::Collision const& collision, muonTracks const& tracks, aod::McParticles const&)
    {
        histos.fill(HIST("eventCounterReco"), 0.5);
        for (auto& track : tracks) {
            if(track.has_mcParticle()){
                auto mcParticle = track.mcParticle();
                auto statusCode = abs(mcParticle.getGenStatusCode());
                if (abs(mcParticle.pdgCode())==13) {
                    auto muMother = mcParticle.mothers_first_as<aod::McParticles>();
                    auto muMotherPDG = abs(muMother.pdgCode());
                    auto muID = mcParticle.globalIndex();
                    bool hasWmother = false;
                    // std::cout << "==== Mu from W decay chain: 13(" << mcParticle.eta() << "," << mcParticle.pt() << ")"
                    //             << " <- " << muMotherPDG << "(" << muMother.eta() << "," << muMother.pt() << ")";

                    if (muMotherPDG != 24) {
                        // check if W mother in decay chain
                        auto mcPart(muMother);
                        auto mcPartPDG = abs(mcPart.pdgCode());

                        while (mcPart.has_mothers()) {
                            mcPart = *(mcPart.mothers_first_as<aod::McParticles>());
                            mcPartPDG = abs(mcPart.pdgCode());
                            //std::cout << " <- " << mcPartPDG << "(" << mcPart.eta() << "," << mcPart.pt() << ")";

                            if (mcPartPDG == 24) {
                                hasWmother = true;
                                break;
                            }
                        }
                    } else {
                        hasWmother = true;
                    }
                    // if (hasWmother) std::cout << " (PICKED)";
                    // std::cout << std::endl;

                    if (hasWmother) {
                        // check if duplicate (ambiguous track)
                        int occuranceCount = count(selectedTracksID.begin(), selectedTracksID.end(), muID);

                        if (occuranceCount < 1) {
                            histos.fill(HIST("PtRecoHist"), mcParticle.pt());
                            histos.fill(HIST("yPtRecoHist"), mcParticle.pt(), mcParticle.eta());
                            histos.fill(HIST("PtResolution"), abs(mcParticle.pt()-track.pt()));
                            selectedTracksID.emplace(selectedTracksID.end(), muID); // add track ID to selected tracks
                        }
                    }
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
            //auto statusCode = abs(mcParticle.getGenStatusCode());
            // if (statusCode <= 34) {
            //     LOGF(info, "Hardest process particle has pdg code %d, from collision %d", mcParticle.pdgCode(), mcParticle.mcCollisionId());
            // }

            if (abs(mcParticle.pdgCode())==13) {
                auto muMother = mcParticle.mothers_first_as<aod::McParticles>();
                auto muMotherPDG = abs(muMother.pdgCode());

                if (muMotherPDG == 24) {
                    //if (mcParticle.eta() > -4.0 && mcParticle.eta() < -2.5) { // muon from W in forward region
                    histos.fill(HIST("yPtTruthHist"), mcParticle.pt(), mcParticle.eta());
                    histos.fill(HIST("PtTruthHist"), mcParticle.pt());
                    //}
                }
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