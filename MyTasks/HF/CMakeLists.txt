o2physics_add_dpl_workflow(hf-muon-fwd-tracks
                    SOURCES HFMuonFwdTracks.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)