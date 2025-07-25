# Building Lyapunov Functionals Algorithmically

Consider the hyperbolic partial differential equation (PDE)

<p align="center">
  <img src="PDEmain.png" alt="" style="width: 240px; display: block; margin: auto;">
</p>

This program attempts to build a Lyapunov functional for the above equation, showing that the energy of the
solution decays as time goes to infinity. The algorithm is described in detail in the paper
"Large-time asymptotics for hyperbolic systems with non-symmetric relaxation: An algorithmic approach"
by Timothée Crin-Barat, Lorenzo Liverani (myself), Ling-Yun Shou and Enrique Zuazua.

The general idea is to first decompose the matrix B into a skew-symmetric (Ba) and symmetric part (Bs):
<p align="center">
  <img src="Bdecomp.png" alt="" style="width: 340px; display: block; margin: auto;">
</p>
Then, one obtains the equation
<p align="center">
  <img src="PDEmod.png" alt="" style="width: 240px; display: block; margin: auto;">
</p>

At this point, one build the functional by iteratively summing together suitable scalar products of the form
⟨X^{k-1} U, X^{k-1} A U⟩ or ⟨X^{k-1} U, X^{k-1} Ba U⟩,
where X^{k-1} is a product of matrices of the form: Bs*(some product combination of A and Ba).

The whole procedure can be described as a tree search. For the full description, please refer to
the original paper.

If you use this code in your research, please cite:

```bibtex
@article{HypNonSym,
title = {Large-time asymptotics for hyperbolic systems with non-symmetric relaxation: An algorithmic approach},
journal = {Journal de Mathématiques Pures et Appliquées},
volume = {202},
pages = {103757},
year = {2025},
issn = {0021-7824},
doi = {https://doi.org/10.1016/j.matpur.2025.103757},
url = {https://www.sciencedirect.com/science/article/pii/S0021782425001011},
author = {Timothée Crin-Barat and Lorenzo Liverani and Ling-Yun Shou and Enrique Zuazua},
}
```

## Features

- Check the matrix condition at every node, to decide which direction to take (1: right; 0: mixed; -1: left)
- Visualize the obtained tree
- Output the LaTeX code to visualize clearly the functional (just copy-paste it into a LaTeX compiler)
- Several presets are available, including all of the examples in the original paper (Section 9).

## Missing Features

- The code currently generates only the functional for the decay in high frequencies (Section 5 of the paper). The functional for low frequencies will be added soon
- Each term of the functional has to be multiplied for a suitably small ε>0. This is not yet implemented
- A final predicted decay rate for the system will be added soon to the output, as well as a numerical check to see if the algorithm is optimal
- Some small additional features such as: checking the Inhomogeneous Kalman Rank condition, developing the whole algorithm for the mixed case, the Sugimoto system among the presets...

## Usage

```bash
python main.py
```
At this point, you can choose between:
1: Input your own A and B
2: Investigate random A and B, choosing first their size and the rank of Bs (this mode is highly unstable for the time being, due to missing features. I suggest NOT to use it)
3: Investigate a suitable preset


## Requirements

- Python 3.9

## File Structure

```
project/
├── main.py
├── tools/
│   ├── tree.py
│   ├── latex.py
│   ├── create_system.py
│   ├── functional.py
│   └── matrix.py
├── output.txt
└── README.md
```

## Output

The program generates output in `output.txt` containing:
- Matrix information
- Tree structure
- Lyapunov functional
- LaTeX output

For an example, please check the uploaded `output.txt` file (containing the output for the Timoshenko system)
