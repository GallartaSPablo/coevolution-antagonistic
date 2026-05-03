This repository contains the source code to simulate and analyze the coevolutionary dynamics of antagonistic networks (exploiter-victim). The project is divided into two phases: the numerical simulation implemented in `C`, and the structural analysis and visualization implemented in `Python` (Jupyter Notebook).

## Project Structure

It is recommended to structure the repository as follows to facilitate reproducibility:

* `src/`: source codes for the C simulations.
  * `Coevolution_antagonistic_networks.c`: simulation for temporal tracking of traits and fitness.
* `notebooks/`: analysis notebooks.
  * `plot_figures_def.ipynb`: processing of the generated CSVs, topological analysis (NODF, intervality gaps), and figure generation.
* `data/` (not included in the repo): directory where the C programs export the data in `.csv` format (e.g., `/Files/Temporal/`, `VE_networks/`).

## Requirements and Dependencies

**For the numerical simulations (C):**
* GCC Compiler.
* Standard math library (`-lm`).

**For analysis and visualization (Python 3.11.9):**
* `numpy`, `pandas`, `csv`, `os`
* `networkx` (network handling and metrics)
* `matplotlib`, `seaborn` (plots, animations, and heatmaps)

##  How to run the code

### 1. Compiling the C models
To achieve maximum performance, the codes are compiled with the `-O3` optimization flag. Open your terminal and run:

\`\`\`bash
# Compile the model
gcc -o bin/Coevolution_antagonistic_networks src/Coevolution_antagonistic_WoL_networks.c -O3 -lm -Wall -Wextra
\`\`\`

### 2. Running the simulations
The compiled binaries require input parameters (Network ID, parameters $\xi_V$ and $\xi_E$ that control fitness, and control flags).

Example of running the temporal model:
\`\`\`bash
./bin/Coevolution_antagonistic_networks_temporal <net_id> <xi_d_V> <xi_d_E> <temporal_flag>
\`\`\`
*Note: the numerical results will be automatically saved in CSV files that will be read by the Notebook.*

### 3. Topological Analysis and Figures
Once the C simulations have finished and exported the CSVs:
1. Open Jupyter Notebook.
2. Run `plot_figures_def.ipynb`.
3. The notebook will automatically read the adjacency matrices (bipartite), and generate the figures for the paper.

## 📄 License
This project is distributed under the [GPLv3] license. See the `LICENSE` file for more details.
