# chiralpl

Simulates vibronic spectra of aggregates of organic chromophores using the Holstein Hamiltonian in the multiparticle basis spectra.
Designed for chiral aggregates.
See EXAMPLE.inp for an example input file.
Compilation requires intel fortran compiler ifx and the MKL library (LAPACK).

## Features
- [x] Arbitrary aggregate size
- [x] Point dipole approximation dipole-dipole coupling of N nearest neighbours
- [x] Simulation and plotting of absorption and PL spectra with up to 4 vibrational quanta in the ground state.
- [x] Properties calculated over configurational average of "diagonal disorder".
- [x] Parallelisation of configurational average using OpenMP
- [x] Simulation of CD spectra
- [ ] Simulation of CPL spectra (Nearly there!)
- [ ] Effects of charge-transfer interactions.
- [ ] Simulation of very large systems using sparse matrices.

### References
- Hestand, N. J.; Spano, F. C. Expanded Theory of H- and J-Molecular Aggregates: The Effects of Vibronic Coupling and Intermolecular Charge Transfer. Chem. Rev. 2018, 118 (15), 7069–7163. https://doi.org/10.1021/acs.chemrev.7b00581.
- Spano, F. C.; Meskers, S. C. J.; Hennebicq, E.; Beljonne, D. Using Circularly Polarized Luminescence to Probe Exciton Coherence in Disordered Helical Aggregates. The Journal of Chemical Physics 2008, 129 (2), 024704. https://doi.org/10.1063/1.2943647.
- Spano, F. C. Absorption and Emission in Oligo-Phenylene Vinylene Nanoaggregates: The Role of Disorder and Structural Defects. The Journal of Chemical Physics 2002, 116 (13), 5877–5891. https://doi.org/10.1063/1.1446034.
- Philpott, M. R. Theory of the Coupling of Electronic and Vibrational Excitations in Molecular Crystals and Helical Polymers. The Journal of Chemical Physics 1971, 55 (5), 2039–2054. https://doi.org/10.1063/1.1676371.

