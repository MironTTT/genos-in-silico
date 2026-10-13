# genos-in-silico

GENOS-IN-SILICO
A Local, Disease-Agnostic Platform for De Novo Drug Discovery with Integrated Multi-Scale Physiological Simulation
Document type: Project Overview
Version: 0.1.0 (Inception)
Status: Active
Date: 2026-10-13

Abstract
Genos-in-Silico is a fully local, offline desktop platform for the computational discovery of novel therapeutic molecules. It couples a generative deep-learning module for molecular design with a deterministic, physics- and physiology-grounded evaluation pipeline that transports each candidate molecule through a sequence of increasingly selective filters: physicochemical arbitration, multi-compartment pharmacokinetics, tissue-permeability estimation, and finally interaction with a mechanistic model of the target pathology at the cellular level.

The platform is disease-agnostic by design. Its first instantiation targets Alzheimer's disease through the amyloid-β pathway, but the generator, the physicochemical arbiter, the pharmacokinetic core, and the user interface are all written independently of any specific pathology. Retargeting the platform to a new disease requires replacing exactly four components — the target-tissue model, the pathology dynamics, the therapeutic mechanism, and the disease-specific constants — while everything else remains untouched.

The system runs entirely on commodity hardware: a single consumer GPU with 8 GB of VRAM for neural inference, a multi-core CPU for numerical simulation, and 16 GB of RAM. No cloud services, no external APIs, no large language models at runtime. Every computation is reproducible, auditable, and owned by the researcher.

# 1. Vision
Modern drug discovery is fragmented. Generative chemistry tools produce molecules but stop at validity. Pharmacokinetic simulators model drug disposition but cannot design drugs. Neural simulators reproduce tissue electrophysiology but never see a drug. Clinical trials integrate everything — but only after a decade of preclinical work and billions of dollars.

Genos-in-Silico closes this loop locally.

It is built on a simple conviction: a single researcher, on a single workstation, should be able to move from an empty screen to a ranked list of drug candidates with predicted pharmacokinetics, predicted target-tissue exposure, predicted therapeutic effect, and predicted toxicity — in minutes, not years. Not to replace wet-lab validation, but to make the in silico stage rigorous, transparent, and genuinely informative, so that the expensive experiments that follow are chosen with care.

The platform is not a black box. It is not a cloud service. It is not a wrapper around a proprietary model. It is a scientific instrument — a piece of software whose every equation, every constant, and every threshold is visible, editable, and reproducible.

# 2. Scientific Philosophy
Three principles shape every design decision.

First: physics before statistics. Where a process can be described by a mechanistic model — mass conservation, first-order kinetics, diffusion on a graph — the platform uses one. Machine learning is reserved for the parts of the problem where mechanism is genuinely unknown or intractable: the generation of novel molecular structures. Everything downstream of generation is deterministic and auditable.

Second: multi-scale integration over single-scale depth. The value of the platform does not come from modeling any single layer with maximum fidelity. It comes from connecting layers that are usually modeled in isolation — from the SMILES string of a molecule, through its concentration in plasma, through its permeability across the blood–brain barrier, to its effect on a single neuron's survival. This integration is the platform's scientific contribution.

Third: reproducibility as a first-class feature. Every simulation produces a complete log: the molecule, its descriptors, the patient parameters, the numerical solver settings, the time series of every state variable. Any result can be re-run bit-for-bit. Any parameter can be traced to its source — a literature reference, a user slider, or a model prediction.

# 3. What the Platform Does
## 3.1. Generates Novel Molecules
At the top of the pipeline sits a variational autoencoder (VAE) operating on SMILES strings. The model is deliberately compact — on the order of two to five million parameters — so that it fits comfortably within the memory budget of a consumer GPU. It is trained on the ZINC-250k dataset, a curated collection of drug-like molecules chosen for its balance of chemical diversity, cleanliness, and computational tractability.

The VAE learns a continuous latent representation of chemical space. Once trained, it can:

Sample novel molecules by drawing from the latent prior.

Interpolate between two known molecules, producing structurally intermediate candidates.

Mutate a molecule by perturbing its latent vector, exploring local chemical neighborhoods.

Cross over two parent molecules by mixing their latent representations, mimicking a form of generative recombination.

The output is a stream of candidate SMILES strings — chemically valid by construction, novel by sampling, and biased toward the drug-like region of chemical space by training.

## 3.2. Filters with Physicochemical Rigor
Not every generated molecule deserves to enter a physiological simulation. The physicochemical arbiter is a deterministic, rule-based filter that rejects molecules violating well-established drug-likeness criteria before they consume computational resources.

The arbiter evaluates:

Molecular weight — above 500 Da, oral bioavailability becomes unreliable.

Lipophilicity — outside a target LogP window, membrane permeation suffers or nonspecific binding dominates.

Aqueous solubility — insufficient solubility in plasma compromises distribution.

Toxicophores — SMARTS-pattern matches against a curated library of known structural alerts.

This stage is fast, conservative, and transparent. It does not use machine learning. It does not guess. It applies hard thresholds derived from decades of medicinal chemistry literature, and it reports its reasoning: every rejected molecule carries a machine-readable explanation. When the user wants to know why a compound was excluded, the answer is immediate and unambiguous.

## 3.3. Simulates Whole-Body Pharmacokinetics
Molecules that pass the arbiter enter a five-compartment physiologically based pharmacokinetic (PBPK) model. The compartments are:

Gut — site of absorption.

Blood — central distribution compartment.

Liver — primary site of metabolic clearance.

Kidney — primary site of renal elimination.

Brain — the target organ (in the Alzheimer's instantiation).

Each compartment is governed by an ordinary differential equation conserving mass. Drug flows between compartments according to organ-specific blood flow rates. Tissue-plasma partition coefficients — how much drug concentrates in each organ relative to plasma — are estimated from molecular descriptors via QSAR relationships. Clearance from liver and kidney follows first-order kinetics with organ-specific rate constants.

The model is parameterized with reference physiological values drawn from ICRP Publication 89. But the parameters are not fixed. The user controls them in real time through the interface:

"Old liver" reduces hepatic clearance and perfusion, simulating age-related metabolic decline.

"Diseased kidney" reduces renal elimination, simulating chronic kidney disease.

Dose determines the initial drug load in the gut compartment.

This is the platform's pharmacological realism: the same molecule produces different plasma profiles, different tissue exposures, and different toxicities depending on the patient being simulated. A candidate that works in a healthy thirty-year-old may fail catastrophically in a comorbid seventy-year-old, and the platform shows this explicitly.

## 3.4. Estimates Blood–Brain-Barrier Permeability
For central nervous system targets, the most common reason a promising molecule fails is that it never reaches the brain. The platform includes an explicit blood–brain-barrier permeability estimator based on molecular polar surface area and lipophilicity — the two descriptors most strongly correlated with passive BBB crossing in the literature.

The estimator returns a permeability coefficient between zero and one, which modulates the effective blood flow to the brain compartment in the PBPK model. A molecule with poor predicted permeability simply never accumulates in the brain, no matter how well it performs in the test tube. The user can additionally adjust a BBB integrity slider, modeling pathological barrier compromise — as occurs in neuroinflammation, aging, and several neurodegenerative conditions.

## 3.5. Models the Disease at the Cellular Level
The brain is represented as a graph of exactly one thousand spiking neurons, each modeled as a leaky integrate-and-fire unit, interconnected according to a NetworkX graph. The exact number — one thousand — is chosen to fit within the computational budget of the target hardware while remaining large enough to exhibit nontrivial network dynamics: propagating activity, cascading failures, and emergent rescue patterns.

In the Alzheimer's instantiation, the pathology is amyloid-β:

Initial state: one hundred neurons carry an amyloid-β load, representing an early symptomatic stage.

Propagation: Aβ diffuses along synaptic edges, spreading to neighboring neurons.

Synaptic degradation: edges exposed to Aβ above a threshold lose conductance and are pruned from the graph.

Apoptosis: neurons whose Aβ load exceeds a lethal threshold are marked dead and removed from the active network.

Therapeutic intervention: drug molecules that have crossed the blood–brain barrier reduce local Aβ concentration at rates determined by their concentration and binding properties.

Off-target toxicity: at high local concentrations, the drug damages healthy synapses and neurons — a trade-off that the platform makes visible, not hidden.

The result of each simulation is not a single score but a trajectory: the time series of neuronal survival, synaptic density, and drug concentration in each compartment, visualized in real time.

# 4. Disease Scalability
Genos-in-Silico is not an Alzheimer's platform. It is a platform whose first instantiation is Alzheimer's. The distinction matters architecturally.

To retarget the platform to a new disease — Parkinson's, Huntington's, a cancer, a cardiac arrhythmia, a metabolic disorder — exactly four components are replaced:

1. The target-tissue model. The spiking neural network of src/brain/ is replaced by whatever tissue model is appropriate: a cardiac electrophysiology model, a tumor-growth model, a hepatic metabolism graph. The interface contract — advance the model by one time step, given a drug concentration — remains identical.

2. The pathology dynamics. The rules describing how the disease evolves in the absence of treatment — amyloid-β accumulation in Alzheimer's, α-synuclein aggregation in Parkinson's, oncogenic signaling in cancer — are replaced. Again, the interface is fixed.

3. The therapeutic mechanism. How the drug acts on the pathology, and how it damages healthy tissue at high concentrations, is disease-specific and is replaced.

4. The disease-specific constants. Thresholds, rate constants, and initial conditions live in a single configuration module and are swapped wholesale.

Everything else — the generator, the physicochemical arbiter, the PBPK core, the GUI, the simulation orchestrator, the logging infrastructure — remains unchanged. A researcher who wants to model a new disease does not need to understand the VAE, the arbiter, or the PBPK integrator. They need to write one new tissue model, one new pathology module, one new therapeutic mechanism, and one new constants file.

This is the operational meaning of "disease-scalable."

# 5. Technical Architecture

## 5.1. Hardware Target
GPU: NVIDIA RTX 3050 (4 GB VRAM) — used for VAE inference only, in batched mode.

CPU: multi-core — used for all ODE integration, graph simulation, and descriptor computation.

RAM: 16 GB minimum.

Storage: ~5 GB for datasets, model weights, and logs.

Network: none required at runtime.

## 5.2. Software Stack
Layer	Library	Reason
GUI	PySide6	Native performance, dark scientific theming, mature
Plotting	pyqtgraph, matplotlib	Real-time 2D plots without rendering overhead
Cheminformatics	RDKit	Industry-standard molecular descriptors and SMILES handling
Deep learning	PyTorch	VAE training and inference, GPU acceleration
Numerical integration	SciPy (solve_ivp, LSODA)	Adaptive stiff ODE solver with rigorous error control
Graph modeling	NetworkX	Brain network construction and traversal
Validation	pydantic, dataclasses	Type-safe configuration and state objects
## 5.3. Design Principles
Immutability of molecular entities. A Molecule object, once constructed from a SMILES string, is frozen. It cannot be mutated. The arbiter either accepts or rejects — never modifies. This guarantees that a molecule's descriptor set is internally consistent throughout the pipeline.

Single source of physical truth. All physiological constants — organ volumes, blood flows, clearance rates — live in one configuration module. No magic numbers appear anywhere else in the codebase. A change to the physiological model is a change to one file.

Strict CPU/GPU separation. VAE inference runs on the GPU in batches. All ODE integration and graph simulation run on the CPU. No tensors cross into SciPy, no NumPy arrays cross into PyTorch. The boundary is explicit and enforced.

Time-scale separation. The brain model advances at millisecond resolution. The PBPK model advances at minute resolution. Coupling occurs through interpolation: the drug concentration at the brain boundary is sampled from the PBPK solution and fed into the neural simulation.

GUI as observer, not participant. The simulation core is unaware of Qt. It publishes events — ticks, neuron deaths, log lines — to an internal event bus. The GUI subscribes. This decoupling enables future headless batch screening: running thousands of simulations without a GUI, for high-throughput virtual screening campaigns.

## 5.4. Project Structure
The codebase is organized into a small number of clearly bounded modules:

src/core/ — domain entities: molecules, patients, simulation state.

src/generator/ — VAE model, tokenizer, sampling, training.

src/filter/ — physicochemical arbiter and descriptor wrappers.

src/pbpk/ — compartment definitions, ODE equations, solver, partition coefficients, BBB estimator.

src/brain/ — neural network, neuron model, amyloid dynamics, toxicity.

src/orchestrator/ — pipeline integration, worker threads, event bus.

src/gui/ — PySide6 windows, panels, widgets, theme.

config/ — physiological constants, patient profiles, settings.

data/ — datasets, model weights, SMARTS toxicophore patterns.

tests/ — validation suite including mass-balance checks.

docs/ — charter, journal, specifications.

Each directory has one job. Each file has one job. The boundary between them is sharp.

# 6. User Experience
The user interacts with the platform through a dark, scientific-themed desktop interface built on PySide6. The design is intentionally austere: no 3D molecular viewers, no animated proteins, no decorative graphics. Every pixel serves the simulation.

The main window is organized around a control panel, a log panel, and a plot panel.

The control panel exposes the simulation's parameters as live sliders and inputs:

Drug dose.

Disease stage.

Liver health.

Kidney health.

Blood–brain-barrier integrity.

Moving any slider immediately updates the underlying simulation. The user does not need to restart; the next integration step reflects the new parameters. This enables interactive exploration: what happens if the patient's kidneys decline mid-treatment? What if the disease progresses faster than expected? The platform answers in real time.

The log panel streams a minute-by-minute textual record of the simulation: drug concentrations, neuron counts, synaptic density, notable events. Every line is timestamped and reproducible.

The plot panel displays live 2D charts:

Concentration–time curves for each compartment.

Neuronal survival over time.

Synaptic density over time.

Therapeutic index — a composite view of efficacy versus toxicity.

The user is the observer, not the driver. The platform provides information; the human decides what to do with it. This is a deliberate stance: the platform is a scientific instrument, not an autonomous agent.

# 7. What the Platform Is Not
To be explicit about scope:

It is not a replacement for wet-lab experiments. No simulation, however mechanistically grounded, can substitute for empirical validation. The platform's purpose is to make the in silico stage more informative, so that the experiments chosen for the wet lab are more likely to succeed.

It is not a clinical decision tool. The patient models are abstractions, calibrated to reference physiology, not to individuals. No output of this platform should be used to guide patient care.

It is not a general-purpose molecular dynamics engine. It does not simulate atomic motion, protein folding, or binding free energies. It operates at the level of molecular descriptors and tissue concentrations.

It is not a cloud service. It runs on the user's machine, uses the user's GPU, and stores data on the user's disk. Nothing is transmitted anywhere.

# 8. Current Scope (v0.1)
The first version delivers:

A trained VAE generator operating on SMILES.

A deterministic physicochemical arbiter.

A five-compartment PBPK simulator with parametric organ health.

A blood–brain-barrier permeability estimator.

A 1,000-neuron spiking network model of amyloid-β-driven neurodegeneration.

A desktop GUI with real-time logs and 2D plots.

A test suite covering mass balance, descriptor consistency, and network dynamics.

Explicitly out of scope for v0.1:

Tau-protein pathology.

Multi-drug regimens.

Patient-specific imaging integration.

Cloud or distributed execution.

3D molecular visualization.

Protein structure prediction or molecular docking.

These are deferred, not rejected. The architecture accommodates them.

# 9. Success Criteria
The platform is considered successful at version 0.1 if:

The full pipeline — generation, arbitration, PBPK, BBB, neural simulation — executes end-to-end for a single candidate molecule in under five minutes on the target hardware.

The PBPK module conserves mass to a tolerance of 10⁻⁶ over simulated intervals when clearance is disabled — a necessary validation that the ODE system is correctly formulated.

The neural network exhibits monotonic neurodegeneration in the absence of intervention and dose-dependent rescue in its presence — a sanity check that the disease model behaves as biology dictates.

The GUI remains responsive at thirty frames per second or higher during simulation, even with all panels active.

A single researcher, unfamiliar with the codebase, can retarget the platform to a new disease by replacing the four designated modules, without touching the rest of the system.

The fifth criterion is the most important. It defines what "disease-scalable" means in practice.

# 10. Long-Term Trajectory
Version 0.1 is a foundation. The roadmap beyond it is deliberately open, but a few directions are natural:

Additional diseases. Parkinson's disease, Huntington's disease, and several cancers are already architecturally supported; only the biology needs writing.

Additional molecular modalities. Beyond small molecules: peptides, oligonucleotides, and targeted protein degraders.

Multi-drug regimens. Extending the PBPK model from single-compound to multi-compound disposition, with interaction terms.

Population simulation. Replacing single-patient simulation with Monte Carlo sampling over patient parameter distributions, yielding population-level pharmacokinetic and pharmacodynamic profiles.

Integration with external data. Importing patient-specific physiological parameters from clinical datasets, within the constraints of local-only execution.

Reproducibility infrastructure. Formal provenance tracking for every simulation, enabling published results to be re-executed exactly.

None of these are commitments. They are directions the architecture is designed to allow.

# 11. References
Delaney, J. S. (2004). ESOL: Estimating aqueous solubility directly from molecular structure. Journal of Chemical Information and Computer Sciences, 44(3), 1000–1005.

ICRP (2002). Basic anatomical and physiological data for use in radiological protection: reference values. ICRP Publication 89.

Clark, D. E. (1999). Rapid calculation of polar molecular surface area and its application to the prediction of transport phenomena. Journal of Pharmaceutical Sciences, 88(8), 807–814.

Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. Advanced Drug Delivery Reviews, 46(1–3), 3–26.

Abbott, N. J., Rönnbäck, L., & Hansson, E. (2006). Astrocyte–endothelial interactions at the blood–brain barrier. Nature Reviews Neuroscience, 7(1), 41–53.

End of project overview — v0.1.0.
