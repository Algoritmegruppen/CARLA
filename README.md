# Community-Aware Recursive Layout Algorithm (CARLA)

This repository implements **CARLA**, a graph layout algorithm introduced in
our submission to the Graph Drawing Conference. CARLA is a force-directed
layout method based on Fruchterman–Reingold, enhanced with two key
contributions: (1) **recursive pre-positioning** based on community structure,
and (2) a simple **edge-weighting scheme** that strengthens intra-community
connections relative to inter-community ones, controlled by a tunable parameter
$\alpha \in (0,1]$. The layout process builds a multilevel visualization by
recursively laying out the quotient graph of communities, resulting in improved
readability and clearer visual separation of densely connected groups.

![](https://raw.githubusercontent.com/Algoritmegruppen/CARLA/refs/heads/main/dolphins.png)

The code supports layout generation, metric-based evaluation (edge crossings,
normalized edge length), and visual inspection of standard benchmark graphs
(e.g., karate, dolphin, football). It includes comparisons against
Fruchterman–Reingold (FR), ForceAtlas2 (FA2), and GRA (Huang et
al. 2020). Results show that CARLA improves layout clarity and structural
fidelity without sacrificing runtime efficiency. Dependencies: Python 3.13,
NetworkX 3.4.2. See the paper for implementation details and performance
analysis.

## Installing and running

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/Algoritmegruppen/CARLA
cd CARLA
```

Or, if using only NetworkX:

```bash
pip install networkx
```

Run the algorithm on an example graph:

```bash
python carla.py football.gexf
```
