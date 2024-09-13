#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace o2;
using namespace o2::framework;

// output textfile for muon tracks
std::ofstream muonTracksOut("muontracks.txt");

struct HFMuonFwdTracks {
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
        const AxisSpec axisPt{10, 0.0, +20.0, "p_{T} (GeV/c)"};

        // create histograms
        histos.add("eventCounterReco", "eventCounterReco", kTH1F, {axisCounter});
        histos.add("muPtHistReco", "muPtHistReco", kTH1F, {axisPt});
        histos.add("muPtHistRecoD", "muPtHistRecoD", kTH1F, {axisPt});

        // set textfile header
        muonTracksOut << "ID,eta,pt,p,phi,motherPDG,nClusters,pDca,chi2,chi2MatchMCHMID,chi2MatchMCHMFT,isPrompt\n";
    }

    using muonTracks = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

    void process(aod::Collision const& collision, muonTracks const& tracks, aod::McParticles const&)
    {
        histos.fill(HIST("eventCounterReco"), 0.5);
        for (auto& track : tracks) {
            if(track.has_mcParticle()){
                auto mcParticle = track.mcParticle();
                //auto statusCode = mcParticle.getGenStatusCode();
                if (abs(mcParticle.pdgCode())==13) {
                    auto muID = mcParticle.globalIndex();
                    auto muEta = mcParticle.eta();
                    if (muEta >= -4.0 && muEta <= -2.5) histos.fill(HIST("muPtHistReco"), mcParticle.pt()); // forward eta region

                    // look for HF mother
                    auto const& muMother = mcParticle.mothers_first_as<aod::McParticles>();
                    auto muMotherPDG = abs(muMother.pdgCode());
                    auto motherStatusCode = muMother.getGenStatusCode();
                    if (abs(muMotherPDG) >= 411 && abs(muMotherPDG) <= 435) { // D meson parent - can be prompt/no-prompt
                        histos.fill(HIST("muPtHistRecoD"), mcParticle.pt());
                    }

                    // check if duplicate (ambiguous track)
                    int occuranceCount = count(selectedTracksID.begin(), selectedTracksID.end(), muID);

                    // write to textfile
                    if (occuranceCount < 1) {
                        auto muRecoEta = track.eta();
                        auto muRecoPt = track.pt();
                        auto muRecoP = track.p();
                        auto muRecoPhi = track.phi();
                        auto muCluster = track.nClusters();
                        auto muDCA = track.pDca();
                        auto muChi2 = track.chi2();
                        auto muChi2MCHMID = track.chi2MatchMCHMID();
                        auto muChi2MCHMFT = track.chi2MatchMCHMFT();
                        int isPrompt = 0;
                        if (abs(motherStatusCode) >= 81 && abs(motherStatusCode) <= 89) isPrompt = 1;
                        muonTracksOut << muID << "," << muRecoEta << "," << muRecoPt << "," << muRecoP << "," << muRecoPhi << "," << muMotherPDG << "," << std::to_string(muCluster) << "," << muDCA << "," << muChi2 << "," << muChi2MCHMID << "," << muChi2MCHMFT << "," << isPrompt << std::endl;
                        selectedTracksID.emplace(selectedTracksID.end(), muID);
                    }
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