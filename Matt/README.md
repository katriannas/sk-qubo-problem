# Matt: O(N²) greedy algorithms for QUBO ground states

Two greedy algorithms. In both, every spin is fixed once along its local field, and afterwards only flips that lower the energy are made. Each folder has a step-by-step README, a plot, and a standalone C file.

| Folder | How the next spin is chosen | Error vs reference, N = 31 / 511 / 4095 |
|---|---|---|
| [`greedy/`](greedy) | field² ÷ remaining coupling strength; no guesses about unfixed spins | +0.04% / −0.95% / −0.34% |
| [`greedy_lookahead/`](greedy_lookahead) | same, plus predicted signs of unfixed spins | −0.02% / −1.28% / −0.65% |
| original greedy (for comparison) | largest \|field\| first | +0.56% / −0.88% / −0.34% |

Negative means the average energy is lower than the reference value.
