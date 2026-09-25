# Input data provenance

The four adjacent text files are public DUNE TDR simulation inputs, downloaded from the ancillary material associated with DUNE Collaboration, “Experiment Simulation Configurations Approximating DUNE TDR,” arXiv:2103.04797 (2021), https://doi.org/10.48550/arXiv.2103.04797. The ancillary flux is the G4LBNF v3r5p4 optimized-engineered-beam FHC far-detector flux. The DUNE flux documentation describes units and format at https://glaucus.crc.nd.edu/DUNEFluxes/documentation.html.

| Local file | Original ancillary basename | Use |
|---|---|---|
| `dune_tdr_fhc_fd_flux.txt` | `histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEFD_globes_flux.txt` | Unoscillated FHC far-detector flux histogram, neutrinos GeV\(^{-1}\) m\(^{-2}\) POT\(^{-1}\) |
| `dune_tdr_cc_cross_sections.dat` | `xsec_cc.dat` | GENIE charged-current cross-section table used by the GLoBES configuration |
| `dune_tdr_numu_fhc_selection_efficiency.txt` | `post_dis_FHC_numu_sig.txt` | Published FHC \(\nu_\mu\) selection efficiency vector |
| `dune_tdr_globes_config.glb` | `DUNE_GLoBES.glb` | Baseline, density, binning, exposure, and response configuration |

The values are carried into this repository as plain text to make the analysis reproducible. This project does not claim ownership of the DUNE inputs; users should consult the arXiv ancillary record and DUNE documentation for source and reuse terms. The flux file's 1297 km reference plane and the GLoBES oscillation baseline of 1284.9 km differ; the simulation records that discrepancy explicitly.
