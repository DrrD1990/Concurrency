# Effect of concurrency on stochastic epidemic dydnamics in SIS and SIR models on networks

We provide here code to accompany the following paper:

R. Liu, M. Ogura, E. Fonseca dos Reis, N. Masuda. 
Effects of concurrency on epidemic spreading in Markovian temporal networks. 
Preprint on arXiv: [arXiv:2201.00754](https://arxiv.org/abs/2201.00754) 

## Jupyter notebook to create figures

```
figure_X.ipynb
```
contains the code to produce fig X in the manuscript, where $X = 2, 3, 4, 5, 6, 7$, and $8$.
The [notebooks/data/ folder](https://github.com/DrrD1990/concurrency/tree/main/notebooks/data) contains the numerical results based on which Figures 5-8 are produced.

## Python code for running stochastic SIS and SIR dynamics

The following Python codes are also included:

- SIS_Model_1prime.py -- to simulate the SIS model on networks with partnership model 1'.
- SIS_Model_2.py -- to simulate the SIS model on networks with partnership model 2.
- SIS_Model_1dprime.py -- to simulate the SIS model on networks with partnership model 1''.
- SIS_Model_3.py -- to simulate the SIS model on networks with partnership model 3.
- SIR_Model_1prime.py -- to simulate the SIR model on networks with partnership model 1'.
- SIR_Model_2.py -- to simulate the SIR model on networks with partnership model 2.
- SIR_Model_1dprime.py -- to simulate the SIR model on networks with partnership model 1''.
- SIR_Model_3.py -- to simulate the SIR model on networks with partnership model 3.
