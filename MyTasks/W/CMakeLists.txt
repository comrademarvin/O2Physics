o2physics_add_dpl_workflow(w-muon-fwd-efficiency
                    SOURCES wMuonFwdEfficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)