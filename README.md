# Review and Critical Analysis of The BROCODE Model: A Novel Mathematical Model for the Briggs-Rauscher Reaction

This repository contains the MATLAB code and computational models developed for the **Network Modelling** course project. The work consists of a critical review and a stochastic extension of the **BROCODE** model, a reduced mathematical framework designed to capture the dynamics of the oscillatory **Briggs-Rauscher (BR)** reaction.

## Authors
* **Federico Pontiroli** - University of Trento
* **Antonella Giorgio** - University of Trento
* **Denise Triacca** - University of Trento

---

## Project Description
The Briggs-Rauscher reaction is an oscillating chemical system that shifts periodically between a radical process (blue phase) and a non-radical process (amber/yellow phase). The BROCODE model aims to simplify the original complex chemical network (which involves up to 15 coupled differential equations) by reducing it to just 3 ordinary differential equations (ODEs) tracking key species: iodous acid, iodide, and molecular iodine.

In this study, we critically analyzed this reduction framework through three main steps:
1. **Deterministic Consistency:** Evaluation of the Quasi-Steady State Approximation (QSSA) applied to low-concentration intermediates, revealing a non-negligible sensitivity of the system dynamics to the selection of kinetic parameters.
2. **Stochastic Framework:** Extension of the model to a stochastic environment by implementing several Stochastic Simulation Algorithms (SSAs) to evaluate the effect of the reaction volume on the convergence between deterministic and stochastic descriptions.
3. **Antioxidant Assay:** Simulating the action of an antioxidant species (ferulic acid) scavenging hydroperoxyl and iodine dioxide radicals to observe and measure the inhibition of chemical oscillations.

---

## Implemented Code and Algorithms
The codebase compares different variants of stochastic simulation algorithms in terms of computational performance and numerical diagnostics:
* **DM** (Direct Method)
* **FRM** (First Reaction Method)
* **NRM** (Next Reaction Method)
* **RSSA / RSSA_dep** (Rejection Stochastic Simulation Algorithm optimized using a dependency graph)

### Directory Structure
* `src/models/`: Contains the definitions of kinetic models and reactant/product stoichiometric matrices (`vMinusBRr`, `vPlusBRr`).
* `src/utils/`: Utility functions such as `dependencygraph.m` (to map structural dependencies between reactions for optimized solvers) and `generateStochasticrates.m` (to convert deterministic rate constants into stochastic propensities).
* `scripts/benchmarks.m`: Script used to benchmark the computational performance of the algorithms (utilizes `parfor` loops for parallel simulation execution).


---

##References

* **[1]** Heath Dimsey, Larry Forbes, and Andrew Bassom. *"The BROCODE model: A novel mathematical model for the Briggs-Rauscher reaction"*. Journal of Mathematical Chemistry, 63 (May 2025), pp. 1806–1825. DOI: [10.21203/rs.3.rs-6765229/v1](https://doi.org/10.21203/rs.3.rs-6765229/v1).
* **[2]** L. Marchetti,Corrado Priami, Vo Hong
Thanh, et al. *Simulation algorithms for compu-
tational systems biology* , [Springer], [2017].
* **[3]** Rinaldo Cervellati et al. *"The Briggs-Rauscher reaction as a test to measure the activity of antioxidants"*. Helvetica Chimica Acta, 84(12), 2001, pp. 3533–3547.
* **[4]** Refer to the pdf for full references!
