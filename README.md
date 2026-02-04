<a href="https://doi.org/10.5281/zenodo.18483822"><img src="https://zenodo.org/badge/1081804198.svg" alt="DOI"></a>

# ERGO: Modelling the Fate of Mutant Glycoproteins in the Endoplasmic Reticulum

This repository contains a Python implementation that operationalizes the conceptual model described in our paper, enabling researchers to simulate, visualize, and test hypotheses about glycoprotein dynamics in the endoplasmatic reticulum (ER).

## Introduction

The ER serves as the primary checkpoint where molecular chaperones (like calnexin/CNX), quality control lectins (like OS9), and glycan-modifying enzymes (mannosidases, glucosidases, UGGT) collectively determine whether a glycoprotein proceeds to the Golgi, remains for additional folding attempts, or gets targeted for ER-associated degradation (ERAD).

The fate of each glycoprotein is intimately tied to the structure of its N-glycans. Progressive mannose trimming by ER mannosidases (ERManI, EDEM1/2/3, EREM) generates distinct glycoforms—ranging from M9 (fully decorated) down to M5 (extensively trimmed)—that serve as molecular "timers" for quality control, because species lacking mannoses have a lower affinity for UGGT. Glucosylation by UGGT and deglucosylation by glucosidase II (GII) govern cycles of calnexin binding and release, providing opportunities for folding. Meanwhile, specific demannosylated glycoforms (like M8C, M7AC, M6) serve as recognition signals for ERAD lectins such as OS9, which recruit the degradation machinery.

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

> Andrea Lia, Daniele Di Bella and Pietro Roversi (2025) *A Conceptual Framework to Model the Fate of Glycoproteins in the Endoplasmic Reticulum* [DOI or Preprint link]

The code as provided generates the Figure [X] in the paper. To reproduce run the script with default parameters. 

## Citation

If you use or adapt this framework in your research, please cite our paper:

**Formatted citation:**
```
Andrea Lia, Daniele Di Bella and Pietro Roversi (2025). Mathematical modelling of glycoprotein fate in the ER. [Journal Name], [Volume(Issue)], [pages]. DOI: [insert DOI]
```

**BibTeX:**
```bibtex
@article{Lia2024glycoprotein,
  title={Mathematical modelling of glycoprotein fate in the ER},
  author={Lia, Andrea and Di Bella, Daniele and Roversi, Pietro},
  journal={[Journal Name]},
  volume={[Volume]},
  number={[Issue]},
  pages={[pages]},
  year={[Year]},
  doi={[DOI]}
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
