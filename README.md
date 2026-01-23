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

A detailed description of the theory, algorithms, and examples is provided in:

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

If you use **HBSS_Cont** in academic work, please cite:

D. Anastasio, S. Marchesiello, Nonlinear frequency response curves estimation
and stability analysis of randomly excited systems
in the subspace framework, Nonlinear Dynamics, 2023. doi: https://doi.org/10.1007/s11071-023-08280-6

---

## References

[1] D. Anastasio, S. Marchesiello, Nonlinear frequency response curves estimation
and stability analysis of randomly excited systems
in the subspace framework, Nonlinear Dynamics, 2023. doi: https://doi.org/10.1007/s11071-023-08280-6
