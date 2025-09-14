#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct wMuonFwdEfficiency {
    // Histogram registry: an object to hold your histograms
    HistogramRegistry histos{"histos", {},
    OutputObjHandlingPolicy::AnalysisObject};

    // definition of reconstructed info to be saved for processing
    struct RecoTrackInfo {
        int64_t trackID;
        uint8_t trackType;
        float chi2;
        float chi2MatchMCHMID;
        float chi2MatchMCHMFT;
        float matchScoreMCHMFT;
        int32_t matchMFTTrackId;
        float pDca;
        float pt;
        float eta;
        float ptResolution;
    };

    // Map to group reconstructed tracks by their associated simulated track ID
    std::unordered_map<int64_t, std::vector<RecoTrackInfo>> trackGroups;

    // output textfile for muon tracks
    std::ofstream muonTracksOut;

    void init(InitContext const&)
    {
        // open output file and write header
        muonTracksOut.open("muonTracks.csv");
        if (!muonTracksOut.is_open()) {
            LOGF(fatal, "Failed to open muonTracks.csv for writing");
        }
        muonTracksOut << "trackID,trackType,phi,tgl,signed1Pt,nClusters,pDCA,rAtAbsorberEnd,sign,chi2,chi2MatchMCHMID,chi2MatchMCHMFT,trackTime,eta,pt,p" << std::endl;

        // define axes you want to use
        const AxisSpec axisCounter{1, 0, +1, ""};
        const AxisSpec axisEta{20, -4.0, -2.5, "#eta"};
        const AxisSpec axisPt{20, 0.0, +80.0, "p_{T} (GeV/c)"};
        const AxisSpec axisDeltaPt{40, 0.0, +20.0, "|p_{T}^{true} - p_{T}^{reco}|(GeV/c)"};
        const AxisSpec axisChi2{50, 0.0, +10.0, "#chi^{2}"};
        const AxisSpec axisChi2Global{50, 0.0, +100.0, "#chi^{2}"};
        const AxisSpec axisDCA{100, 0.0, +1000.0, "pDCA"};
        const AxisSpec axisTrackType{5, 0, 5, "Track Type"};

        // create histograms
        histos.add("eventCounterReco", "eventCounterReco", kTH1F, {axisCounter});
        histos.add("eventCounterSim", "eventCounterSim", kTH1F, {axisCounter});
        histos.add("yPtRecoHist", "yPtRecoHist", kTH2F, {axisPt, axisEta});
        histos.add("yPtTruthHist", "yPtTruthHist", kTH2F, {axisPt, axisEta});
        histos.add("PtRecoHist", "PtRecoHist", kTH1F, {axisPt});
        histos.add("PtTruthHist", "PtTruthHist", kTH1F, {axisPt});
        histos.add("PtResolution", "PtResolution", kTH1F, {axisDeltaPt});
        histos.add("chi2", "chi2", kTH1F, {axisChi2});
        histos.add("chi2_global_track", "chi2_global_track", kTH1F, {axisChi2Global});
        histos.add("chi2MatchMCHMID", "chi2MatchMCHMID", kTH1F, {axisChi2});
        histos.add("chi2MatchMCHMFT", "chi2MatchMCHMFT", kTH1F, {axisChi2});
        histos.add("pDCA", "pDCA", kTH1F, {axisDCA});
        histos.add("trackType", "trackType", kTH1D, {axisTrackType});
    }

    using muonTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

    //void processReco(aod::Collision const& collision, muonTracks const& tracks, aod::McParticles const&) // run with collisions
    void processReco(muonTracks const& tracks, aod::McParticles const&) // run without collisions
    {
        histos.fill(HIST("eventCounterReco"), 0.5);

        for (auto& track : tracks) {
            if(track.has_mcParticle()) { // check if track has associated MC particle and has a valid chi2 match
                auto mcParticle = track.mcParticle();
                int64_t recoTrackID = track.globalIndex();    // Reconstructed track ID

                if (abs(mcParticle.pdgCode())==13) { // check if associated MC track is a muon
                    int64_t mcTrackID = mcParticle.globalIndex(); // Simulated track ID

                    // check if muon is from W decay
                    auto muMother = mcParticle.mothers_first_as<aod::McParticles>();
                    auto muMotherPDG = abs(muMother.pdgCode());
                    
                    bool hasWmother = false;
                    //std::cout << "==== Mu decay chain: 13 <- " << muMotherPDG;

                    if (muMotherPDG != 24) {
                        // check if W mother in decay chain
                        auto mcPart(muMother);
                        auto mcPartPDG = abs(mcPart.pdgCode());

                        while (mcPart.has_mothers()) {
                            mcPart = *(mcPart.mothers_first_as<aod::McParticles>());
                            mcPartPDG = abs(mcPart.pdgCode());
                            //std::cout << " <- " << mcPartPDG;

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
                        auto muTrackType = static_cast<int64_t>(track.trackType());
                        auto muChi2 = track.chi2();
                        auto muChi2MatchMCHMID = track.chi2MatchMCHMID();
                        auto muChi2MatchMCHMFT = track.chi2MatchMCHMFT();
                        auto muMatchScoreMCHMFT = track.matchScoreMCHMFT();
                        auto muMatchMFTTrackId = track.matchMFTTrackId();
                        auto muDca = track.pDca();
                        auto muPt = track.pt();
                        auto muEta = track.eta();
                        auto ptResolution = abs(mcParticle.pt() - track.pt());
        
                        // Store the reconstructed track information in the map for additional processing, associated with the simulated track ID
                        RecoTrackInfo trackInfo = {recoTrackID, muTrackType, muChi2, muChi2MatchMCHMID, muChi2MatchMCHMFT, muMatchScoreMCHMFT, muMatchMFTTrackId,
                            muDca, muPt, muEta, ptResolution};
                        trackGroups[mcTrackID].push_back(trackInfo);

                        // save all muon tracks to the output file
                        muonTracksOut << recoTrackID << "," << muTrackType << "," << track.phi() << "," << track.tgl() << "," << track.signed1Pt() << ","
                                    << static_cast<int64_t>(track.nClusters()) << "," << muDca << "," << track.rAtAbsorberEnd() << ","
                                    << static_cast<int64_t>(track.sign()) << "," << muChi2 << "," << muChi2MatchMCHMID << "," << muChi2MatchMCHMFT << ","
                                    << track.trackTime() << "," << muEta << "," << muPt << "," << track.p() << std::endl;

                        histos.fill(HIST("trackType"), muTrackType);
                        
                        // basic cuts before plotting
                        if (muTrackType == 3 && muChi2MatchMCHMID > 0) { // only look at standalone tracks for now
                            // write to histograms
                            histos.fill(HIST("PtRecoHist"), muPt);
                            histos.fill(HIST("yPtRecoHist"), muPt, muEta);
                            histos.fill(HIST("PtResolution"), ptResolution);
                            histos.fill(HIST("chi2"), muChi2);
                            histos.fill(HIST("chi2MatchMCHMID"), muChi2MatchMCHMID);
                            histos.fill(HIST("pDCA"), muDca);
                        }

                        // global muon tracks
                        if (muTrackType == 0 && muChi2MatchMCHMFT > 0) {
                            histos.fill(HIST("chi2_global_track"), muChi2);
                            histos.fill(HIST("chi2MatchMCHMFT"), muMatchScoreMCHMFT);
                        }
                    }
                }
            }
        }

        // additional track processing to check the multiple reconstructed tracks associated with the same simulated track
        for (auto& [mcTrackID, recoTracks] : trackGroups) {
            // prcess tracks to choose track with lowest chi2 (this is now pointless?)
            // int trackIndex = 0;
            // int chosenIndex = trackIndex;
            // float minChi2 = 100000.0;
            // for (const auto& trackInfo : recoTracks) {
            //     if (trackInfo.chi2 < minChi2) {
            //         minChi2 = trackInfo.chi2;
            //         chosenIndex = trackIndex; // remember the index of the track with the lowest chi2
            //     }
            //     trackIndex++;
            // }

            // choose when to print the information
            bool printRecoTracks = false;
            if (printRecoTracks) {
                // printing
                if (recoTracks.size() > 1) {
                    // Output a heading for the simulated track
                    LOGF(info, "Simulated track ID %ld has multiple associated reconstructed tracks:", mcTrackID);

                    // List the associated reconstructed tracks
                    for (const auto& trackInfo : recoTracks) {
                        LOGF(info, "RecoTrack ID = %ld, Track Type = %d, chi2 = %.2f, chi2MatchMCHMID = %.2f, chi2MatchMCHMFT = %.2f, matchScoreMCHMFT = %.2f, MFT Track ID = %d",
                        trackInfo.trackID, trackInfo.trackType, trackInfo.chi2, trackInfo.chi2MatchMCHMID, trackInfo.chi2MatchMCHMFT, trackInfo.matchScoreMCHMFT, trackInfo.matchMFTTrackId);
                    }
                }
                //LOGF(info, "Chosen track has index %d with track ID %ld", chosenIndex, recoTracks[chosenIndex].trackID);
            }
        }
    }
    PROCESS_SWITCH(wMuonFwdEfficiency, processReco, "W->muon reconstructed information", true);

    //Filter<Tracks> etaFilter = track::eta < -2.5 && track::eta > -4;

    void processSim(aod::McCollision const& collision, aod::McParticles const& mcParticles)
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