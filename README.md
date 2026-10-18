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

Authors

**Nala Steward, Jadean Govender, Stellenbosch University**

**References**
Alzheimer’s Association, 2024. 2024 Alzheimer’s disease facts and figures. Alzheimer’s & Dementia 20, 3708–3821. https://doi.org/10.1002/alz.13809
Alzheimer’s Association, n.d. What is Alzheimer’s? [WWW Document]. Alzheimer’s Association. URL https://www.alz.org/alzheimers-dementia/what-is-alzheimers (accessed 10.12.26).
Alzheimer’s Disease International, n.d. ADI - Dementia statistics. URL https://www.alzint.org/about/dementia-facts-figures/dementia-statistics/ (accessed 10.11.26).
Betts, J.G., Young, K.A., Wise, J.A., Johnson, E., Poe, B., Kruse, D.H., Korol, O., Johnson, J.E., Womble, M., DeSaix, P., Betts, J.G., Young, K.A., Wise, J.A., Johnson, E., Poe, B., Kruse, D.H., Korol, O., Johnson, J.E., Womble, M., DeSaix, P., 2022. 12.1 Basic Structure and Function of the Nervous System - Anatomy and Physiology 2e | OpenStax [WWW Document]. URL https://openstax.org/books/anatomy-and-physiology-2e/pages/12-1-basic-structure-and-function-of-the-nervous-system (accessed 10.16.26).
Cleveland Clinic, 2023. Nervous System: What Does It Do? [WWW Document]. Cleveland Clinic. URL https://my.clevelandclinic.org/health/body/21202-nervous-system (accessed 10.13.26).
d‘Errico, P., Meyer-Luehmann, M., 2020. Mechanisms of Pathogenic Tau and Aβ Protein Spreading in Alzheimer’s Disease. Front Aging Neurosci 12, 265. https://doi.org/10.3389/fnagi.2020.00265
Doig, A.J., 2018. Positive Feedback Loops in Alzheimer’s Disease: The Alzheimer’s Feedback Hypothesis. Journal of Alzheimer’s Disease 66, 25–36. https://doi.org/10.3233/JAD-180583
Fornari, S., Schäfer, A., Jucker, M., Goriely, A., Kuhl, E., 2019. Prion-like spreading of Alzheimer’s disease within the brain’s connectome. J R Soc Interface 16, 20190356. https://doi.org/10.1098/rsif.2019.0356
Hao, W., Friedman, A., 2016. Mathematical model on Alzheimer’s disease. BMC Syst Biol 10, 108. https://doi.org/10.1186/s12918-016-0348-2
Hardy, J., 2006. Alzheimer’s disease: The amyloid cascade hypothesis: An update and reappraisal.
Hardy, J.A., Higgins, G.A., 1992. Alzheimer’s disease: the amyloid cascade hypothesis. Science 256, 184–186.
Harmon, A., 2014. Neurons (nerve cells) | Health and Medicine | Research Starters | EBSCOhost [WWW Document]. EBSCO. URL https://www.ebsco.com/ (accessed 10.16.26).
Harrar, S., 2024. Why Is It So Hard To Treat Alzheimer’s? [WWW Document]. AARP. URL https://www.aarp.org/health/conditions-treatments/why-is-there-no-cure-for-alzheimers/ (accessed 10.12.26).
Heneka, M.T., Carson, M.J., Khoury, J.E., Landreth, G.E., Brosseron, F., Feinstein, D.L., Jacobs, A.H., Wyss-Coray, T., Vitorica, J., Ransohoff, R.M., Herrup, K., Frautschy, S.A., Finsen, B., Brown, G.C., Verkhratsky, A., Yamanaka, K., Koistinaho, J., Latz, E., Halle, A., Petzold, G.C., Town, T., Morgan, D., Shinohara, M.L., Perry, V.H., Holmes, C., Bazan, N.G., Brooks, D.J., Hunot, S., Joseph, B., Deigendesch, N., Garaschuk, O., Boddeke, E., Dinarello, C.A., Breitner, J.C., Cole, G.M., Golenbock, D.T., Kummer, M.P., 2015. Neuroinflammation in Alzheimer’s disease. The Lancet Neurology 14, 388–405. https://doi.org/10.1016/S1474-4422(15)70016-5
Hippius, H., Neundörfer, G., 2003. The discovery of Alzheimer’s disease. Dialogues Clin Neurosci 5, 101–108. https://doi.org/10.31887/DCNS.2003.5.1/hhippius
Kaufman, S.K., Diamond, M.I., 2013. Prion-Like Propagation of Protein Aggregation and Related Therapeutic Strategies. Neurotherapeutics 10, 371–382. https://doi.org/10.1007/s13311-013-0196-3
Korczyn, A.D., 2012. Why have we failed to cure Alzheimer’s disease? - PubMed [WWW Document]. URL https://pubmed.ncbi.nlm.nih.gov/22258512/ (accessed 10.12.26).
Kuznetsov, Y.A., 2006. Saddle-node bifurcation. Scholarpedia 1, 1859. https://doi.org/10.4249/scholarpedia.1859
Liu, E., Zhang, Y., Wang, J.-Z., 2024. Updates in Alzheimer’s disease: from basic research to diagnosis and therapies. Transl Neurodegener 13, 45. https://doi.org/10.1186/s40035-024-00432-x
Long, C., Fritts, A., Broadway, J., Brawman-Mintzer, O., Mintzer, J., 2025. Neuroinflammation: A Driving Force in the Onset and Progression of Alzheimer’s Disease. Journal of Clinical Medicine 14, 331. https://doi.org/10.3390/jcm14020331
National Institute on Aging, 2024. What Happens to the Brain in Alzheimer’s Disease? [WWW Document]. National Institute on Aging. URL https://www.nia.nih.gov/health/alzheimers-causes-and-risk-factors/what-happens-brain-alzheimers-disease (accessed 10.12.26).
National Institute on Aging, n.d. What Is Alzheimer’s Disease? | National Institute on Aging [WWW Document]. URL https://www.nia.nih.gov/health/alzheimers-and-dementia/what-alzheimers-disease (accessed 10.11.26).
Nisbet, R.M., Polanco, J.-C., Ittner, L.M., Götz, J., 2015. Tau aggregation and its interplay with amyloid-β. Acta Neuropathol 129, 207–220. https://doi.org/10.1007/s00401-014-1371-2
Rasmussen, J., Langerman, H., 2019. Alzheimer’s Disease – Why We Need Early Diagnosis. Degener Neurol Neuromuscul Dis 9, 123–130. https://doi.org/10.2147/DNND.S228939
ScienceDirect, n.d. Bistability - an overview | ScienceDirect Topics [WWW Document]. URL https://www.sciencedirect.com/topics/mathematics/bistability (accessed 10.18.26a).
ScienceDirect, n.d. Hysteresis Phenomenon - an overview | ScienceDirect Topics [WWW Document]. URL https://www.sciencedirect.com/topics/engineering/hysteresis-phenomenon (accessed 10.18.26b).
Vaquer-Alicea, J., Diamond, M.I., 2019. Propagation of Protein Aggregation in Neurodegenerative Diseases. Annual Review of Biochemistry 88, 785–810. https://doi.org/10.1146/annurev-biochem-061516-045049
Walker, L.C., n.d. Aβ plaques. Free Neuropathol 1, 31. https://doi.org/10.17879/freeneuropathology-2020-3025
Webers, A., Heneka, M.T., Gleeson, P.A., 2020. The role of innate immune responses and neuroinflammation in amyloid accumulation and progression of Alzheimer’s disease. Immunology and Cell Biology 98, 28–41. https://doi.org/10.1111/imcb.12301
World Health Organization, n.d. 11 million lives lost each year: urgent action needed on neurological care [WWW Document]. URL https://www.who.int/news/item/14-10-2025-11-million-lives-lost-each-year-urgent-action-needed-on-neurological-care (accessed 10.11.26).
