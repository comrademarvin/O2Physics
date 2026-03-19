o2physics_add_dpl_workflow(muon-reduced-reader
                    SOURCES processReducedMuons.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(muon-reduced-reader-mc
                    SOURCES processReducedMuonsMC.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)