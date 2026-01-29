# HBSS_Cont

**HBSS_Cont** is a Matlab toolbox for nonlinear frequency response analysis
based on Harmonic Balance in extended state space combined with pseudo-arc-length continuation.

The toolbox supports continuous and discrete-time state-space models, created analytically or identified experimentally with methods such as NSI (Nonlinear Subspace Identification) or NFR-ID. 

Stability analysis of periodic
solutions is also performed through Floquet multipliers, with automatic detection of
Fold, Period-Doubling, and Neimark–Sacker bifurcations.

---

## Main features

- Harmonic Balance formulation in extended state space
- Pseudo-arc-length continuation of nonlinear frequency response curves
- Support for theoretical and experimentally identified models
- Stability analysis via monodromy matrix and Floquet multipliers
- Automatic detection of bifurcation events

---

## Documentation

A detailed description of the theory, algorithms, and examples is provided in  `docs/HBSS_Cont_User_Guide.pdf`. 

---

## Getting started

Tested on Matlab R2025b.

Run any script in the `examples/` folder.
Each example automatically adds the `src/` directory to the MATLAB path.

---

## License

This software is intended for research purposes and is released under the GNU General Public License v3 (GPLv3). See the `LICENSE` file for details.


Disclaimer: this software is provided "as is", without warranty of any kind. The authors shall not be held liable for any inaccuracies in the results or damages resulting from the use of the code.

---

## Citation

If you use **HBSS_Cont** in academic work, please cite the following paper, which provides methodological description of the method:

**D. Anastasio, S. Marchesiello,** *Nonlinear frequency response curves estimation and stability analysis of randomly excited systems in the subspace framework.* **Nonlinear Dynamics**,  2023. DOI: [10.1007/s11071-023-08280-6](https://doi.org/10.1007/s11071-023-08280-6)

---

## References

The methodology implemented in this toolbox has been adopted and discussed in:

[1] D. Anastasio, S. Marchesiello, *Nonlinear frequency response curves estimation
and stability analysis of randomly excited systems
in the subspace framework*, Nonlinear Dynamics, 2023. DOI: [10.1007/s11071-023-08280-6](https://doi.org/10.1007/s11071-023-08280-6)

[2] D. Anastasio, S. Marchesiello, G. Kerschen, *Estimation of the periodic solutions of geometrically nonlinear structures by broadband excitation*, Proceedings of ISMA2024. [Link](http://past.isma-isaac.be/downloads/isma2024/proceedings/Contribution_417_proceeding_3.pdf) 

[3] D. Anastasio, G. Raze, G. Kerschen, *Frequency-domain system identification of nonlinear structures using experimental continuation data*, Journal of Vibration and Control, 2025. DOI: [10.1177/10775463251405337](https://doi.org/10.1177/10775463251405337) 

[4] D. Anastasio, G. Raze, S. Marchesiello, G. Kerschen, *Identification of primary and secondary resonances: experimental continuation and broadband data-based modeling*, Nonlinear Dynamics, 2026.DOI: [10.1007/s11071-025-11940-4](https://doi.org/10.1007/s11071-025-11940-4) 
