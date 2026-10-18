# Biomathematics-374
Mathematical Modelling of Alzheimer’s Disease: A Nonlinear Dynamical System Model of Protein Accumulation and Neuroinflammation

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=nala-steward/Biomathematics-374&file=Code.m)

**Threshold-Dependent Neuroinflammation in a Nonlinear Model of Alzheimer's Disease**

A four-variable ordinary differential equation (ODE) model, analysed in MATLAB, that asks whether positive feedback between amyloid-β, tau, neuroinflammation and neuronal death can produce self-sustaining chronic inflammation, and where the tipping thresholds lie.

**Overview**

The classic amyloid cascade hypothesis treats Alzheimer's disease as a one-way process. This project instead models a feedback loop in which amyloid-β drives inflammation, inflammation drives neuronal death, and dying neurons release further amyloid-β. Saturating Hill functions are used for the nonlinear responses.

**State variables**
  - A:	Amyloid-β
  - T:	Tau
  - I:	Neuroinflammation
  - D: Dead neurons

**Key results**
  - **Three equilibria at baseline** (k_A1 = 1.3): a stable healthy-like state (A = 0.85, I = 0.10, D = 0.49), an unstable saddle that acts as a threshold, and a stable chronic state (A = 5.98, I = 0.99, D = 8.07).
  - **Two saddle-node bifurcations** at k_A1 ≈ 1.155 and 1.485. Between them the healthy and chronic states coexist, producing hysteresis.
  - Bistability requires sigmoidal responses. With n = m = 1 the system has a single equilibrium.
  - **Bistable region**: a diagonal band in the (k_A1, K_A) plane covering about 21% of the explored grid (331 of 1,600 combinations).
  - **Tipping**: starting from the healthy state, a transition to the chronic state occurred only when initial inflammation exceeded I(0) ≈ 2.04. None of the amyloid-β pulses tested triggered a transition.
  - **Sensitivity**: healthy-state amyloid-β depends mainly on amyloid clearance and production, followed by inflammation resolution. Tau parameters had little influence.

**What the code does**
  - **Equilibria and stability**: finds all equilibria at baseline parameters, computes the Jacobian and its eigenvalues, and classifies each point.
  - **Time course**s: simulates trajectories from zero initial conditions and from different starting values of A.
  - **Bifurcation analysis**: varies k_A1 to locate the saddle-node bifurcations and the hysteresis region.
  - **Parameter-plane map**: counts equilibria across a 40 × 40 grid of k_A1 and K_A.
  - **Critical perturbation tests**: finds the initial inflammation needed to tip the system into the chronic state.
  - **Local sensitivity analysis**: computes normalised sensitivities of the healthy-state amyloid level to ±1% changes in each parameter.

**Limitations**
  - Parameters are illustrative rather than fitted to data.
  - Direct amyloid-β and tau toxicity are switched off, so tau's minor role is partly a consequence of the model's assumptions.
  - Amyloid-tau coupling, a finite healthy-neuron population, ageing and random fluctuations are not included.
  - The model is well-mixed and deterministic, with no spatial spread.
  - The sensitivity analysis covers only the healthy equilibrium, not the bifurcation thresholds.
  - Grid-based root finding can miss closely spaced equilibria.

**Possible extensions**
  - Link tau dynamics to amyloid-β accumulation
  - Add a healthy-neuron population that declines as neurons die
  - Separate pro- and anti-inflammatory microglial states
  - Vary k_A1 and k_A2 independently
  - Test combined perturbations and perturbations of D
  - Run sensitivity analysis at the tipping threshold
  - Fit parameters to biomarker or protein-turnover data
  - Add spatial structure and brain connectivity

**References**
Fornari, S., Schäfer, A., Jucker, M., Goriely, A. and Kuhl, E. (2019) 'Prion-like spreading of Alzheimer's disease within the brain's connectome', Journal of the Royal Society Interface.
Hao, W. and Friedman, A. (2016) 'Mathematical model on Alzheimer's disease', BMC Systems Biology, 10, 108.

Add the other works cited in your report as needed.

Author

**Nala Steward, Jadean Govender, Stellenbosch University**
