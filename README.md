# conformer-energetics

Tools to calculate conformer energetics of solvated drug-like molecules from MD trajectories.

## Contents

- `Conf/`: Example dataset and workflow for propylene glycol.
- `prob_map.py`: Computes OCCO dihedral angles and H–H distances, classifies gauche/trans conformers, counts transitions, and constructs 2D probability maps.
- `err_calc.py`: Checks ergodicity by analysing the number of gauche-trans transitions per ns (expects low standard deviation across trajectories). Also estimates uncertainties on basin populations by randomising frames, performing block averaging, and computing the standard error of the mean.

