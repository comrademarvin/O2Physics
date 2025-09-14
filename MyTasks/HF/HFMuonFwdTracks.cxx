#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct HFMuonFwdTracks {
    // Histogram registry: an object to hold your histograms
    HistogramRegistry histos{"histos", {},
    OutputObjHandlingPolicy::AnalysisObject};

    // output textfile for muon tracks
    std::ofstream muonTracksOut;

    void init(InitContext const&)
    {
        // open output file and write header
        muonTracksOut.open("muonTracks_HF.csv");
        if (!muonTracksOut.is_open()) {
            LOGF(fatal, "Failed to open muonTracks.csv for writing");
        }
        muonTracksOut << "trackID,trackType,phi,tgl,signed1Pt,nClusters,pDCA,rAtAbsorberEnd,sign,chi2,chi2MatchMCHMID,chi2MatchMCHMFT,trackTime,eta,pt,p" << std::endl;

        // define axes you want to use
        const AxisSpec axisCounter{1, 0, +1, ""};
        const AxisSpec axisEta{10, -4.0, -2.5, "#eta"};
        const AxisSpec axisPt{10, 0.0, +20.0, "p_{T} (GeV/c)"};
        const AxisSpec axisTrackType{5, 0, 5, "Track Type"};

        // create histograms
        histos.add("eventCounterReco", "eventCounterReco", kTH1F, {axisCounter});
        histos.add("muPtHistReco", "muPtHistReco", kTH1F, {axisPt});
        histos.add("trackType", "trackType", kTH1D, {axisTrackType});
    }

    Int_t GetFlavour(Int_t pdgCode)
    {
        //
        // return the flavour of a particle
        // input: pdg code of the particle
        // output: Int_t
        //         3 in case of strange (open and hidden)
        //         4 in case of charm (")
        //         5 in case of beauty (")
        //
        Int_t pdg = TMath::Abs(pdgCode);
        // Resonance
        if (pdg > 100000)
            pdg %= 100000;
        if (pdg > 10000)
            pdg %= 10000;
        // meson ?
        if (pdg > 10)
            pdg /= 100;
        // baryon ?
        if (pdg > 10)
            pdg /= 10;
        return pdg;
    }

    using muonTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

    void process(muonTracks const& tracks, aod::McParticles const&)
    {
        histos.fill(HIST("eventCounterReco"), 0.5);
        for (auto& track : tracks) {
            if(track.has_mcParticle()){
                auto mcParticle = track.mcParticle();
                int64_t recoTrackID = track.globalIndex();    // Reconstructed track ID
                //auto statusCode = mcParticle.getGenStatusCode();
                if (abs(mcParticle.pdgCode())==13) {
                    int64_t mcTrackID = mcParticle.globalIndex(); // Simulated track ID

                    // check if muon is from HF decay
                    auto muMother = mcParticle.mothers_first_as<aod::McParticles>();
                    auto muMotherPDG = abs(muMother.pdgCode());
                    auto muMotherFlavour = GetFlavour(muMotherPDG);

                    bool hasHFmother = false;
                    std::cout << "==== Mu decay chain: 13 <- " << muMotherPDG;

                    if (muMotherFlavour != 4 && muMotherFlavour != 5) { // direct mother is not HF hadron
                        // check if HF hadron mother in decay chain
                        auto mcPart(muMother);
                        auto mcPartPDG = abs(mcPart.pdgCode());
                        auto mcPartFlavour = GetFlavour(mcPartPDG);

                        while (mcPart.has_mothers()) {
                            mcPart = *(mcPart.mothers_first_as<aod::McParticles>());
                            mcPartPDG = abs(mcPart.pdgCode());
                            mcPartFlavour = GetFlavour(mcPartPDG);
                            std::cout << " <- " << mcPartPDG;

                            if (mcPartFlavour == 4 || mcPartFlavour == 5) {
                                hasHFmother = true;
                                break;
                            }
                        }
                    } else {
                        hasHFmother = true;
                    }
                    if (hasHFmother) std::cout << " (PICKED)";
                    std::cout << std::endl;

                    if (hasHFmother) {
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

                        histos.fill(HIST("trackType"), muTrackType);

                        if (muTrackType == 3) histos.fill(HIST("muPtHistReco"), muPt); // look at pT for standalone tracks

                        // save all muon tracks to the output file
                        muonTracksOut << recoTrackID << "," << muTrackType << "," << track.phi() << "," << track.tgl() << "," << track.signed1Pt() << ","
                                    << static_cast<int64_t>(track.nClusters()) << "," << muDca << "," << track.rAtAbsorberEnd() << ","
                                    << static_cast<int64_t>(track.sign()) << "," << muChi2 << "," << muChi2MatchMCHMID << "," << muChi2MatchMCHMFT << ","
                                    << track.trackTime() << "," << muEta << "," << muPt << "," << track.p() << std::endl;
                    };

                    // check whether is prompt HF hadron
                    // auto mcPart(muMother);
                    // auto prevMcPart(muMother);
                    // auto mcPartPDG = abs(mcPart.pdgCode());
                    // int isPrompt = 1;
                    // std::cout << "==== Forward muon decay chain: mu";
                    // while (mcPart.has_mothers() && (abs(mcPart.getGenStatusCode()) > 80 || mcPart.getGenStatusCode() == 0)) { // print out mother chain
                    //     std::cout << " <- " << mcPartPDG;
                    //     prevMcPart = *(mcPart);
                    //     mcPart = *(mcPart.mothers_first_as<aod::McParticles>());
                    //     mcPartPDG = abs(mcPart.pdgCode());
                    // }
                    // if (div(abs(prevMcPart.pdgCode()), 100).quot != div(muMotherPDG, 100).quot) isPrompt = 0;
                    // std::cout << "; isPrompt = " << isPrompt << std::endl;
                }
            }
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFMuonFwdTracks>(cfgc)};
}