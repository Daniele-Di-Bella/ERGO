# ERGO: Modelling the Fate of Mutant Glycoproteins in the Endoplasmic Reticulum

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18483823.svg)](https://doi.org/10.5281/zenodo.18483823)

This repository contains a Python implementation that operationalizes the conceptual model described [in our work](https://www.biorxiv.org/content/10.64898/2025.11.30.691435v1), enabling researchers to simulate, visualize, and test hypotheses about glycoprotein dynamics in the endoplasmatic reticulum (ER).

## Introduction

The ER serves as the primary checkpoint where molecular chaperones (like calnexin/CNX), quality control lectins (like OS9), and glycan-modifying enzymes (mannosidases, glucosidases, UGGT) collectively determine whether a glycoprotein proceeds to the Golgi, remains for additional folding attempts, or gets targeted for ER-associated degradation (ERAD).

The fate of each glycoprotein is intimately tied to the structure of its N-glycans. Progressive mannose trimming by ER mannosidases (ERManI, EDEM1/2/3, EREM) generates distinct glycoforms — ranging from M9 (fully decorated) down to M5 (extensively trimmed) — that serve as molecular "timers" for quality control, because species lacking mannoses have a lower affinity for UGGT. Glucosylation by UGGT and deglucosylation by glucosidase II (GII) govern cycles of calnexin binding and release, providing opportunities for folding. Meanwhile, specific demannosylated glycoforms (like M8C, M7AC, M6) serve as recognition signals for ERAD lectins such as OS9, which recruit the degradation machinery.

This project presents a modeling framework that integrates:
- N-glycan trimming pathways across three mannose branches (A, B, C)
- Lectin binding dynamics (CNX for quality control cycles; OS9 for ERAD targeting)
- Enzymatic transformations (UGGT, GII, mannosidases, EREM)
- Competing fates (secretion to Golgi vs. ERAD-mediated degradation)

By providing an explicit, ODE-based representation of these processes, the framework allows researchers to explore how different components and parameters shape outcomes. Whether you're interested in testing specific mechanistic hypotheses, comparing model predictions with experimental data, or teaching ER glycobiology through interactive simulations, this code offers a flexible foundation for exploration.

## Model Overview

The model tracks 13 free **monoglycosylated** species in the ER lumen, along with their complexes with quality control lectins:

### N-Glycan Species
- Glucosylated forms: G1M9, G1M8B, G1M8C, G1M7BC
- Deglucosylated forms: M9, M8A, M8B, M8C, M7AB, M7AC, M7BC, M6, M5
  - Early glycoforms (M9, M8)
  - Extensively trimmed glycoforms (M8C, M7AC, M7BC, M6, M5)

### Key Processes Modeled

1. UGGT-catalyzed glucosylation: Adds glucose to M9, M8B, M8C, M7BC → enables CNX rebinding
2. GII-catalyzed deglucosylation: Removes glucose from G1-species
3. Mannose trimming across three branches:
   - Branch B: ERManI, EDEM2 (e.g., M9 → M8B → M7BC)
   - Branches A & C: EDEM1, EDEM3 (e.g., M9 → M8A → M7AB; M8C → M7AC → M6 → M5)
   - EREM: Direct Glc-Man cleavage on glucosylated species 
4. Lectin binding/unbinding:
   - CNX binds G1-species (quality control cycles)
   - OS9 binds extensively trimmed species (ERAD targeting)
5. Competing sinks:
   - Secretion: free glycoforms can exit to Golgi
   - ERAD: OS9-bound species are degraded

### Outputs

Running the simulation produces:
- Time-course plots showing the evolution of free, secreted, and degraded glycoprotein pools;
- Species-level dynamics for all 13 N-glycan structures;
- Comparative scenarios (e.g., UGGT inhibition, mannosidase inhibition, combined perturbations).

All simulations start from a single glycoform (M9, 1 µM) and track its processing over time (up to ~28 hours by default).

## Model Equations
**G1M9**, the mono-glucosylated glycoform is — in the _ALG6_ KO background — produced only by UGGT-mediated re-glucosylation of M9. Its free concentration is reduced both by ER GluII deglucosylation and by reversible binding to CNX. Because CNX is assumed to display identical association and dissociation kinetics toward all mono-glucosylated glycospecies \cite{Schrag-2001}, a single dissociation constant governs the equilibrium of each CNX–glycan complex:

$$
\begin{align}
K_{d\_\mathrm{CNX}}
&= \frac{[CNX] [G1M9]}{[CNX\_{G1M9}]} \\
&= \frac{[CNX] [G1M8B]}{[CNX\_{G1M8B}]} \\
&= \frac{[CNX] [G1M8C]}{[CNX\_{G1M8C}]} \\
&= \frac{[CNX] [G1M7BC]}{[CNX\_{G1M7BC}]}\\
&= \frac{k_{\mathrm{off\_CNX}}}{k_{\mathrm{on\_CNX}}}.
\end{align}
$$

The total CNX concentration therefore reads:

$$
\begin{align}
[CNX]_{\mathrm{tot}}
&= [CNX] \\
&+ [CNX\_{G1M9}] \\ 
&+ [CNX\_{G1M8B}] \\
&+ [CNX\_{G1M8C}] \\
&+ [CNX\_{G1M7BC}].
\end{align}
$$

ER ManI, the EDEMs and ER\ EM can digest G1M9. UGGT-mediated re-glucosylation slows ERAD, and CNX binding reduces the G1M9 free pool.  Collecting all contributions, the time evolution of G1M9 is:

$$
\begin{align}
\frac{d[G1M9]}{dt}
&= k_{UGGT_{G1M9}}[UGGT] [M9] \\
&+ k_{\mathrm{off\_CNX}}[CNX\_{G1M9}] \\
& - \Big(
k_{GluII_{M9}}[GluII] \\
&+ k_{ERManI_{G1M8B}}[ERManI] \\
&+ k_{EDEM1_{G1M8C}}[EDEM1] \\
&+ k_{EDEM3_{G1M8C}}[EDEM3] \\
&+ k_{EDEM2_{G1M8B}}[EDEM2] \\
&+ k_{ER\ EM_{M8A}}[ER\ EM]
\Big) [G1M9]  \\
& - k_{\mathrm{on\_CNX}}[CNX] [G1M9]. 
\end{align}
$$

**G1M8B**  is generated by ERManI or EDEM2 acting on G1M9, or by UGGT acting on M8B, and it is digested by GluII. The ER endo-mannosidase digests G1M8B to give M7AB. Binding/dissociating to/from CNX also subtracts/adds free G1M8B from the solution:

$$
\begin{align}
\frac{d[G1M8B]}{dt}&= \Big(k_{ERManI_{G1M8B}}[ERManI]\\
&+k_{EDEM2_{G1M8B}}[EDEM2]\Big)[G1M9]\\
&+k_{UGGT_{G1M8B}}[UGGT][M8B]\\
&-\Big(k_{GluII_{M8B}}[GluII]\\
&+k_{CNX_{on}}[CNX]\\
&+k_{ER\ EM_{M7AB}}[ER\ EM]\Big)[G1M8B]\\
&+k_{\mathrm{off\_CNX}}[CNX\_{G1M8B}]
\end{align}
$$

**G1M8C** is generated by EDEM1 or EDEM3 acting on G1M9, or by UGGT acting on M8C, it can be digested by EDEM2 and ERManI and by the ER EM. G1M8C can bind to and dissociate from CNX and OS9. For the latter equilibrium, we assume only one kinetic constant for all glycospecies for the binding to OS9. We write the dissociation equilibrium constant as $K_{d\_\mathrm{OS9}}$ with $k_{\mathrm{off\_{OS9}}}$ and $k_{\mathrm{on\_{OS9}}}$  rates:

$$
\begin{align}
K_{d\_\mathrm{OS9}}&=\frac{[OS9][M8C]}{[OS9\_{M8C}]}\\
&=\frac{[OS9][G1M8B]}{[OS9\_G1M8B]}\\
&=\frac{[OS9][G1M8C]}{[OS9\_G1M8C]}\\
&=\frac{[OS9][M7AC]}{[OS9\_{M7AC}]}\\
&=\frac{[OS9][M7BC]}{[OS9\_{M7BC}]}\\
&=\frac{[OS9][G1M7BC]}{[OS9\_G1M7BC]}\\
&=\frac{[OS9][M6]}{[OS9\_{M6}]}\\
&=\frac{[OS9][M5]}{[OS9\_{M5}]}\\
&=\frac{k_{\mathrm{off\_{OS9}}}}{k_{\mathrm{on\_{OS9}}}}
\end{align}
$$

with: 

$$
\begin{align}
    [OS9]_{tot}&=[OS9] \\
    &+[OS9\_{M8C}]\\
    &+[OS9\_G1M8B]\\
    &+[OS9\_G1M8C]\\
    &+[OS9\_{ }]\\
    &+[OS9\_{M7BC}]\\
    &+[OS9\_G1M7BC]\\
    &+[OS9\_{M6}]\\
    &+[OS9\_{M5}].
\end{align}
$$

In other words we assume that OS9 has the same on/off kinetics with all its ligand glycospecies and therefore the dissociation constant is one and the same for each OS9 complex with any of its ligand glycospecies. 

The ER endo-mannosidase (ER EM) digests G1M8C to give  .  Binding to CNX also subtracts free G1M8C from the solution. Overall, the rate of change of G1M8C concentration is:

$$
\begin{align}
\frac{d[G1M8C]}{dt}&=\Big(k_{EDEM1_{G1M8C}}[EDEM1]\\
&+k_{EDEM3_{G1M8C}}[EDEM3]\Big)[M9]\\
&+k_{UGGT_{G1M8C}}[UGGT][M8C]\\
&-\Big(k_{EDEM2_{G1M7BC}}[EDEM2]\\
&+k_{ERManI_{G1M7BC}}[ERManI]\\
&+k_{\mathrm{on\_{OS9}}}[OS9]\\
&+k_{GluII_{M8C}}[GluII]\\
&+k_{CNX_{on}}[CNX]\\
&+k_{ER\ EM_{ }}[ER\ EM]\Big)[G1M8C]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_G1M8C]\\
&+k_{\mathrm{off\_CNX}}[CNX\_{G1M8C}]
\end{align}
$$
    
**M9**, the high-mannose glycan, is the default _N_-glycan attached to the nascent mutant glycoprotein in the _ALG6_ KO background. It is the substrate of UGGT-mediated re-glucosylation, and of ERManI, EDEM2 and EDEM1/3 demannosylation. It can be generated back when GluII removes the terminal glucose from the mono-glucosylated G1M9 species previously produced by UGGT:

$$
\begin{align}
\frac{d[M9]}{dt}&=k_{GluII_{M9}}[GluII][G1M9]\\
&-\Big(k_{UGGT_{G1M9}}[UGGT]\\
&+k_{EDEM2_{M8B}}[EDEM2]\\
&+(k_{EDEM1_{M8A}}+k_{EDEM1_{M8C}})[EDEM1]\\
&+(k_{EDEM3_{M8A}}+k_{EDEM3_{M8C}})[EDEM3]\\
&+k_{ERManI_{M8B}}[ERManI]\Big)[M9]
\end{align}
$$

**M8A**, in an _ALG6_ KO background is generated  by:
* ER EM, from monoglucosylated species G1M9, G1M8B, G1M8C;
* ER exo-mannosidases acting on M9: not ER ManI or EDEM2 (which are specific for branch B), but EDEM1 and EDEM3.

The rate of change in the concentration of M8A can be calculated as:

$$
\begin{align}
\frac{d[M8A]}{dt}
%&k_{ER\ EM}Man\ IA*(G3M9+G2M9+G1M9)\\
&=k_{ER\ EM_{M8A}}[ER\ EM][G1M9]\\
&+\Big(k_{EDEM1_{M8A}}[EDEM1]\\
&+k_{EDEM3_{M8A}}[EDEM3]\Big)[M9]\\
&-\Big(k_{EDEM1_{ }^{C}}[EDEM1]\\
&+k_{EDEM3_{ }^{C}}[EDEM3]\\
&+k_{EDEM2_{M7AB}}[EDEM2]\\
&+k_{ERManI_{M7AB}}[ERManI]\Big)[M8A]
\end{align}
$$
 
**M8B** is generated by either ERManI or EDEM2 acting on M9, it is digested by EDEM1, EDEM3, and it is reglucosylated by UGGT:

$$
\begin{align}
\frac{d[M8B]}{dt}&=\Big(k_{ERManI_{M8B}}[ERManI]\\
&+k_{EDEM2_{M8B}}[EDEM2]\Big)[M9]\\
&-\Big((k_{EDEM1_{M7BC}}+k_{EDEM1_{M7AB}})[EDEM1]\\
& +(k_{EDEM3_{M7BC}}+k_{EDEM3_{M7AB}})[EDEM3]\\
&+k_{UGGT_{G1M8B}}[UGGT]\Big)[M8B]
\end{align}
$$

**M8C** is generated by either EDEM1 or EDEM3 acting on M9, or ER GluII acting on G1M8C and it is digested by EDEM2 or ER ERManI, reglucosylated by UGGT and it enters an equilibrium with OS9.  The rate of change of free M8C concentration is:

$$
\begin{align}
\frac{d[M8C]}{dt}&=\Big(k_{EDEM1_{M8C}}[EDEM1]\\
&+k_{EDEM3_{M8C}}[EDEM3]\Big)[M9]\\
&+k_{GluII_{M8C}}[GluII][G1M8C]\\
&-\Big(k_{EDEM2_{M7BC}}[EDEM2]\\
&+k_{ERManI_{M7BC}}[ERManI]\\
&+k_{EDEM1_{ }^{A}}[EDEM1]\\
&+k_{EDEM3_{ }^{A}}[EDEM3]\\
&+k_{\mathrm{on\_{OS9}}}[OS9]\\
&-k_{UGGT_{G1M8B}}[UGGT]\Big)[M8C]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_{M8C}]
\end{align}
$$
    
**M7AB** is generated by Man IB or EDEM2 acting on M8A, or by the ER EM endomannosidase complex acting on G1M8B, and it is digested by EDEM1 or EDEM3:

$$
\begin{align}
\frac{d[M7AB]}{dt}&=\Big(k_{EDEM2_{M7AB}}[EDEM2]\\
&+k_{ERManI_{M7AB}}[ERManI]\Big)[M8A]\\
&+k_{ER\ EM_{M7AB}}[ER\ EM][G1M8B]\\
&+\Big(k_{EDEM1_{M7AB}}[EDEM1]\\
&+k_{EDEM3_{M7AB}}[EDEM3]\Big)[M8B]\\
&-\Big(k_{EDEM1_{M6}^{C}}[EDEM1]\\
&+k_{EDEM3_{M6}^{C}}[EDEM3]\Big)[M7AB]
\end{align}
$$

**M7BC** is generated by Man IB or EDEM2 acting on M8C or by EDEM1 or EDEM3 on M8B, it is reglucosylated by UGGT and it is involved in the equilibrium binding to OS9:

$$
\begin{align}
\frac{d[M7BC]}{dt}&=\Big(k_{EDEM2_{M7BC}}[EDEM2]\\
&+k_{ERManI_{M7BC}}[ERManI]\Big)[M8C]\\
&+\Big(k_{EDEM1_{M7BC}}[EDEM1]\\
&+k_{EDEM3_{M7BC}}[EDEM3]\Big)[M8B]\\
&-\Big(k_{UGGT_{G1M7BC}}[UGGT]\\
&+k_{\mathrm{on\_{OS9}}}[OS9]\\
&+k_{EDEM1_{M6}^{A}}[EDEM1]\\
&+k_{EDEM3_{M6}^{A}}[EDEM3]\Big)[M7BC]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_{M7BC}]
\end{align}
$$

**G1M7BC** is generated by UGGT acting on M7BC, Man\ IB or EDEM2 acting on G1M8C or by EDEM1 or EDEM3 on G1M8B, it is deglucosylated by ER GluII and it is involved in the equilibria binding to CNX and OS9:

$$
\begin{align}
\frac{d[G1M7BC]}{dt}&=k_{UGGT_{G1M7BC}}[UGGT]*[M7BC]\\
&+\Big(k_{EDEM2_{M7BC}}[EDEM2]\\
&+k_{ERManI_{M7BC}}[ERManI]\Big)[G1M8C]\\
&+\Big(k_{EDEM1_{M7BC}}[EDEM1]\\
&+k_{EDEM3_{M7BC}}[EDEM3]\Big)[G1M8B]\\
&-\Big(k_{CNX_{on}}[CNX]\\
&+k_{GluII_{M7BC}}[GluII]\\
&+k_{\mathrm{on\_{OS9}}}[OS9]\Big)[G1M7BC]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_G1M7BC]\\
&+k_{CNX_{off}}[CNX\_{G1M7BC}]
\end{align}
$$

**M7AC** is generated by EDEM1 or EDEM3 acting on M8A or M8C, by the ER\ EM endomannosidase complex acting on G1M8C, it is digested by EDEM2 or ERManI and it is involved in the equilibrium binding to OS9:

$$
\begin{align}
\frac{d[M7AC]}{dt}&=\Big(k_{EDEM1_{M7AC}^{C}}[EDEM1]\\
&+k_{EDEM3_{M7AC}^{C}}[EDEM3]\Big)[M8A]\\
&+\Big(k_{EDEM1_{M7AC}^{A}}[EDEM1]\\
&+k_{EDEM3_{M7AC}^{A}}[EDEM3]\Big)[M8C]\\
&+k_{ER\ EM_{M7AC}}[ER\ EM][G1M8C]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_{M7AC}]\\
&-\Big(k_{EDEM2_{M6}}[EDEM2]\\
&+k_{ERManI_{M6}}[ERManI]\\
&+k_{\mathrm{on\_{OS9}}}[OS9]\Big)[M7AC]
\end{align}
$$

**M6** is generated by EDEM1,3 acting on M7BC or M7AB, by the ER EM endomannosidase complex acting on G1M7BC, and by ERManI or EDEM2 acting on M7AC; it is involved in the equilibrium binding to OS9:

$$
\begin{align}
\frac{d[M6]}{dt}&=k_{ER\ EM_{M6}}[ER\ EM][G1M7BC]\\
&+\Big(k_{EDEM2_{M6}}[EDEM2]\\
&+k_{ERManI_{M6}}[ERManI]\Big)[M7AC]\\
&+\Big(k_{EDEM1_{M6}^{C}}[EDEM1]\\
&+k_{EDEM3_{M6}^{C}}[EDEM3]\Big)[M7AB]\\
&+\Big(k_{EDEM1_{M6}^{A}}[EDEM1]\\
&+k_{EDEM3_{M6}^{A}}[EDEM3]\Big)[M7BC]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_{M6}]\\
&-k_{\mathrm{on\_{OS9}}}[OS9][M6]
\end{align}
$$

**M5** is generated by EDEM1,3 acting on M6, and it is involved in the equilibrium binding to OS9:

$$
\begin{align}
\frac{d[M5]}{dt}&=\Big(k_{EDEM1_{M5}}[EDEM1]\\
& +k_{EDEM3_{M5}}[EDEM3]\Big)[M6]\\
&+k_{\mathrm{off\_{OS9}}}[OS9\_{M5}]\\
&-k_{\mathrm{on\_{OS9}}}[OS9][M5]
\end{align}
$$

The sink ERAD term is proportional to the concentration of the glycospecies:OS9 complex, and the $K_{\mathrm{ERAD}}$ sink rate is one and the same for all seven OS9 complexes. ER GluII does not remove Glc from glucosylated OS9-bound substrates: as we said, CNX- and OS9-binding are mutually exclusive (a ternary complex of any glycospecies with both lectins is not envisaged).

1. 
$$
\begin{align}
\frac{d[OS9\_{M8C}]}{dt}&= k_{\mathrm{on\_{OS9}}}[OS9][M8C]\\
&-k_{ERAD}[OS9\_{M8C}]
\end{align}
$$

2.
$$
\begin{align}
\frac{d[OS9\_{M7AC}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][M7AC]\\
&-k_{ERAD}[OS9\_{M7AC}]\\
\end{align}
$$

3.
$$
\begin{align}
\frac{d[OS9\_{M7BC}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][M7BC]\\
&-k_{ERAD}[OS9\_{M7BC}]\\
\end{align}
$$

4. 
$$
\begin{align}
\frac{d[OS9\_{M6}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][M6]\\
&-k_{ERAD}[OS9\_{M6}]
\end{align}
$$

5.
$$
\begin{align}
\frac{d[OS9\_{M5}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][M5]\\
&-k_{ERAD}[OS9\_{M5}]
\end{align}
$$

6.
$$
\begin{align}
\frac{d[OS9\_{G1M8C}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][G1M8C]\\
&-k_{ERAD}[OS9\_{G1M8C}]
\end{align}
$$

8.
$$
\begin{align}
\frac{d[OS9\_{G1M7BC}]}{dt}&=k_{\mathrm{on\_{OS9}}}[OS9][G1M7BC]\\
&-k_{ERAD}[OS9\_{G1M7BC}]
\end{align}  
$$












## Model Assumptions and Limitations

As with any model, this framework makes specific simplifying assumptions:

1. **Well-mixed compartment**: The ER is treated as a single homogeneous space (no spatial gradients);
2. **Mass action kinetics**: All reactions follow simple bimolecular or pseudo-first-order kinetics;
3. **No explicit folding states**: Glycoproteins are distinguished only by their N-glycan structure, not by protein conformation; 
4. **Fixed enzyme pools**: Enzyme concentrations are constant;
5. **Single N-glycan per protein**: The model doesn't account for proteins with multiple glycosylation sites;
6. **Deterministic dynamics**: Uses ODEs rather than stochastic simulation.

## Code Structure

The implementation consists of five main sections:

### 1. Initialization
- Defines the 13 free glycan species and their lectin-binding properties
- Sets initial condition (M9 = 1.0 µM; all others = 0)
- Specifies enzyme concentrations, binding kinetics, and catalytic rate constants

### 2. ODE System (`ode_rhs` function)
- Computes the rate of change for all species: `dy/dt = f(y, params)`
- Handles:
  - Lectin binding/unbinding
  - Enzymatic transformations 
  - Secretion and ERAD sinks 
- Uses algebraic constraints for free lectin pools (CNX_free, OS9_free)

### 3. Solver (`run_case` function)
- Numerically integrates the ODEs using `scipy.integrate.solve_ivp`
- 
### 4. Plotting Functions
- `plot_aggregated`: Shows three aggregate curves (Free in ER, Secreted, Degraded)
- `plot_species`: Displays all 13 individual N-glycan trajectories
- Both use log-scale time axes and include t≈0 markers from initial conditions

### 5. Scenario Comparison
Four simulation conditions are run side-by-side:
- **A (Baseline)**: All enzymes active
- **B (UGGT inhibited)**: No glucosylation → CNX cycles disrupted
- **C (Mannosidases inhibited)**: No trimming → glycans remain in early forms
- **D (UGGT + Mannosidases inhibited)**: Both pathways blocked

Results are displayed in a 4×2 grid (left: aggregate; right: species-level).

## Getting Started

### Prerequisites
- Python 3.8+
- Required libraries: `numpy`, `scipy`, `matplotlib`

### Running in Google Colab (recommended)

You can run this code directly in Google Colab (no installation required):

1. Upload the `.py` file to Colab or copy-paste the code into a notebook
2. Install dependencies (usually `numpy` and `matplotlib` pre-installed): `!pip install scipy`
3. Run the cells sequentially
4. Modify parameters and re-run to explore different scenarios

### Running locally
#### Installation
```bash
# Clone the repository
git clone https://github.com/Daniele-Di-Bella/ER_glycoforms_fates_modelling.git
cd ER_glycoforms_fates_modelling

# Install dependencies
pip install numpy scipy matplotlib jupyterlab
```

#### Running the Simulation

Simply launch Jupyter Lab 
```bash
jupyter lab
```
and execute the script cell-by-cell in a Jupyter Notebook.

Output: A multi-panel figure displaying:
- Left column: Aggregate pools (Free, Secreted, Degraded) for each condition
- Right column: Individual N-glycan species trajectories

## Adapting the Framework

This model is designed as an open and extensible research tool. Here are some ways to customize it:

### 1. Change Enzyme Levels
Simulate overexpression or knockdown experiments:
```python
# Example: Double UGGT concentration
pmod = {"UGGT": 0.4}  # baseline was 0.2 µM
t, free_sum, SEC, degraded, sol_y = run_case(pmod)
```

### 2. Modify Kinetic Parameters
Explore how rate constants affect system behavior:
```python
# Example: Slow down mannosidase activity by 50%
pmod = {k: 0.5 * params[k] for k in params if k.startswith("k_ERManI_")}
```

### 3. Test Specific Inhibitors
Model pharmacological perturbations:
```python
# Inhibit only EDEM1/2/3 (leave ERManI and EREM active)
pmod = {k: 0.0 for k in params if "EDEM" in k and k.startswith("k_")}
```

### 4. Adjust Initial Conditions
Start from different glycoforms or mixtures:
```python
# Start from M8C instead of M9
y0_custom = np.zeros(len(full_ER))
y0_custom[name_to_index["M8C"]] = 1.0
```

### 5. Add New Reactions
Extend the reaction network in the `ode_rhs` function:
- Add new glycan species or intermediate states
- Include additional lectins or chaperones
- Model spatial compartmentalization (e.g., separate ER subdomains)

### 6. Compare with Experimental Data
Generate multiple scenarios, and see which scenario fits better your own measurements. 

## Reproducing Figures from the Paper

If you're interested in reproducing the specific results presented in our publication:

> Daniele Di Bella, Andrea Lia and Pietro Roversi (2026). Mathematical modelling of glycoprotein fate in the ER [version 1; peer review: 1 approved with reservations]. Wellcome Open Res, 11:185. https://doi.org/10.12688/wellcomeopenres.25558.1

The code as provided generates the Figures 4, 5 and 6 in the paper. To reproduce them run the script with default parameters. 

## Citation

If you use or adapt this framework in your research, please cite our paper:

**Formatted citation:**
```
Daniele Di Bella, Andrea Lia and Pietro Roversi (2026). Mathematical modelling of glycoprotein fate in the ER [version 1; peer review: 1 approved with reservations]. Wellcome Open Res, 11:185. https://doi.org/10.12688/wellcomeopenres.25558.1
```

**BibTeX:**
```bibtex
@article{DiBella2026glycoprotein,
  title={Mathematical modelling of glycoprotein fate in the ER},
  author={Di Bella, Daniele, Lia, Andrea and Roversi, Pietro},
  journal={Wellcome Open Res},
  volume={11},
  number={185},
  year={2026},
  link={https://doi.org/10.12688/wellcomeopenres.25558.1}
}
```

## License

This project is released under the **MIT License**, meaning you are free to:
- Use the code for any purpose (research, teaching, commercial applications)
- Modify and extend it for your own needs
- Distribute your modifications
- Incorporate it into larger projects

The only requirement is attribution (see LICENSE file for details).

We believe open science accelerates discovery—feel free to build upon this work!

## Questions, Issues, or Contributions?

We welcome feedback and collaboration! If you:
- Encounter bugs or numerical issues
- Have suggestions for model extensions
- Want to share how you've adapted the framework
- Need help interpreting results

Please:
- Open an issue on this GitHub repository
- Contact the corresponding authors at daniele.dibella@ibba.cnr.com, pietro.roversi@cnr.it

**Happy modelling! 🧬**
