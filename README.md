# In Silico Pharmacological Modelling: **_Daniellia oliveri_** and **_Flabellaria paniculata_** Bioactives vs. mPGES-1

## Project Overview
This project performs in silico molecular docking to investigate the interaction bioactives obtained from _Daniellia oliveri_ and _Flabellaria paniculata_ and mPGES-1 (microsomal Prostaglandin E Synthase-1). This research builds upon the phytochemical identification work of Prof. Sofidiya's group and Oladosu's group, aimed at understanding the competitive inhibitory potential of these bioactives against the mPGES-1 enzyme to block inflammatory $PGE_2$ biosynthesis.

## Objective
To physically model the bioactive 3D conformations bioactives obtained from **_Daniellia oliveri_** and **_Flabellaria paniculata_** and evaluate their spatial, electrostatic, and structural compatibility within the active site binding pocket of human **mPGES-1** (using structural biology data from the Protein Data Bank- PDB: **4AL0**).

## Biological Context
- Target: mPGES-1 (PDB: 4AL0). An inducible enzyme crucial in the inflammatory pathway (synthesizing PGE2). The script dynamically retains its essential GSH cofactor during preparation and detects the interfacial dimeric binding cleft.
- Compound Library:
- - Daniellia oliveri (Phenolics & Flavonoids: Caffeic acid, Rutin, Quercetin, etc.)
- - Flabellaria paniculata (Terpenoids & Sterols: Friedelin, Manool, alpha-Humulene, etc.)
- Controls: Standard NSAIDs (Ibuprofen, Diclofenac, Aspirin, Celecoxib).

## Methodology/ Pipeline Architecture
The workflow is driven by ```main.py```, which seamlessly links four automated scripts:
1. Ligand Preparation (```ligand_prep.py```): Uses RDKit to convert SMILES strings to 3D dimensions. Generates 50 conformers using ETKDGv3, minimizes energy via the MMFF94 force field, and performs thermodynamic and Graph Automorphism RMSD pruning to find the optimal 3D shape.
2. Receptor Preparation (```receptor_prep.py```): Uses BioPython to parse the raw PDB hierarchy. It dynamically isolates the functional chains, calculates the mathematical centroid of the native ligand for the docking grid, and purges crystallization artifacts while keeping structural cofactors.
3. Simulation Execution (```docking.py```): Utilizes Meeko and OpenBabel for rigid receptor and flexible ligand PDBQT conversions (with physiological pH 7.4 protonation). Executes the AutoDock Vina C++ binary asynchronously via Python's subprocess module.
4. Complex Generation (```generate_complex.py```): Extracts the top binding pose (Thermodynamic Minimum) from Vina's multi-model output, truncates proprietary syntax, and cleanly merges the ligand and receptor coordinates into a single .pdb file for immediate visualization in PyMOL/Chimera.

## Tech Stack
- Python (3.10+)
- RDKit, Biopython, Meeko
- AutoDock Vina
- PyMOL (Molecular Visualization)

## Key Results
The docking simulation successfully identified the binding mode of the bioactives within the mPGES-1 interfacial cleft.
| Compound             | Affinity_kcal_mol |
|----------------------|-------------------|
| Celecoxib            | -8.377            |
| Apigenin-7-glucoside | -7.977            |
| Quercetin            | -7.456            |
| Ibuprofen            | -7.417            |
| Kaempferol           | -7.145            |
| Catechin             | -6.883            |
| Friedelin            | -6.743            |
| Caffeic acid         | -6.725            |
| Rutin                | -6.664            |
| Manool               | -6.57             |
| Ferulic acid         | -6.502            |
| p-Coumaric acid      | -6.487            |
| beta-Elemene         | -6.336            |
| Diclofenac           | -6.26             |
| Beta-Sitosterol      | -6.215            |
| Aspirin              | -5.99             |
| alpha-Sabinene       | -5.436            |
| alpha-Humulene       | -4.449            |

## Visual
![rutin_docking_final.png](results/rutin_docking_final.png)

*Figure: 3D Visualization showing Apigenin-7glucoside (green) docked in the mPGES-1 pocket, interacting with key residues (Arg70, Arg126).*

## Project Structure
<details>
<summary>Click to expand the project tree structure</summary>

```bash
📁 mpges1-3d-danielliaxflabellaria/
    ├── 📁 data/
    │   └── 📁 raw/
    │       └── 📄 4AL0.pdb
    ├── 📁 notes/
    │   ├── 📄 Batch docking results summary_ terminal.txt
    │   ├── 📄 Faculty of Pharmacy Virtual Confere.md
    │   ├── 📄 prof sofidiyah's anti-inflammatory research.txt
    │   └── 📄 rutin3d Technical Report.docx
    ├── 📁 scripts/
    │   ├── 📄 docking.py
    │   ├── 📄 generate_complex.py
    │   ├── 📄 ligand_prep.py
    │   └── 📄 receptor_prep.py
    ├── 📄 .gitignore
    ├── 📄 LICENSE
    ├── 📄 main.py
    ├── 📄 README.md
    ├── 📄 requirements.txt
    ├── 📄 run_pipeline.ipynb
    └── 📄 vina.exe
```

</details>

## How to Run
1. Ensure all dependencies are installed
```commandline
pip install rdkit-pypi biopython meeko vina gemmi openbabel
```
2. Clone the repository
3. Ensure your raw target protein (4AL0.pdb) is located in ./data/raw/.
2. Run the full pipeline from your terminal
```commandline
python main.py
```
3. Once complete, open PyMOL and run the styling script
```text
@results/style_complex.pml
```
_Note: If running on your local system, you will also need the standalone AutoDock Vina executable (vina.exe for Windows or the vina binary for Linux/macOS) placed in your system path or project root._

## Output Structure
The pipeline automatically creates a structured, reproducible workspace:
- ```data/processed/ligand_prep/``` : Optimized .sdf conformers.
- ```data/processed/receptor_prep/``` : Cleaned, dehydrated .pdb receptor.
- ```data/processed/docking/``` : Prepared .pdbqt files and raw Vina logs/poses.
- ```results/``` :
- - ```master_docking_scores.csv``` (A ranked summary of all binding affinities).
- - ```*_mpges1_complex.pdb``` (Publication-ready merged complexes for 3D visualization).

## References
- Protein Data Bank (PDB ID: 4AL0)
- Forli Lab, Scripps Research (Meeko API)
- Prof. Sofidiya et al.

_Created by Daniella Ene-Obong, B.Pharm Candidate, University of Lagos_
