# Sources

Checked 2026-11-25.

## Primary

- G. Mitchison and R. Durbin. Bounds on the learning capacity of some
  multi-layer networks. Biological Cybernetics 60 (1989) 345–356. Rigorous
  upper bounds for committee machines.
- E. Barkai, D. Hansel and H. Sompolinsky. Broken symmetries in multilayered
  perceptrons. Physical Review A 45 (1992) 4146. Replica analysis of the tree
  committee machine.
- A. Engel, H. M. Köhler, F. Tschepke, H. Vollmayr and A. Zippelius. Storage
  capacity and learning algorithms for two-layer neural networks. Physical
  Review A 45 (1992) 7590.
- M. Stojnic. Capacity of the treelike sign perceptrons neural networks with
  one hidden layer: RDT based upper bounds. arXiv:2312.08244 (2023). Proves the
  replica symmetric value is an upper bound.
- M. Stojnic. Lifted RDT based capacity analysis of the 1-hidden layer treelike
  sign perceptrons neural networks. arXiv:2312.08257 (2023). Best proved upper
  bound we found (3.43 for K = 3).

## Tools used in the proof

- T. M. Cover. Geometrical and statistical properties of systems of linear
  inequalities with applications in pattern recognition. IEEE Trans.
  Electronic Computers EC-14 (1965) 326–334.
- J. G. Wendel. A problem in geometric probability. Mathematica Scandinavica
  10 (1962) 109–111.
- Y. Gordon. On Milman's inequality and random subspaces which escape through a
  mesh in ℝⁿ. Geometric Aspects of Functional Analysis, LNM 1317 (1988).
- C. Thrampoulidis, S. Oymak and B. Hassibi. Regularized linear regression: a
  precise analysis of the estimation error. COLT 2015 (convex Gaussian min-max
  theorem).
- K. R. Davidson and S. J. Szarek. Local operator theory, random matrices and
  Banach spaces. Handbook of the Geometry of Banach Spaces, vol. 1 (2001).
- D. P. Bertsekas. Convex Optimization Theory. Athena Scientific (2009)
  (exact penalty for convex programs).

## Later work that retests, refines or disputes it

- D. Malzahn and A. Engel. Correlations between hidden units in multilayer
  neural networks and replica symmetry breaking. arXiv:cond-mat/9906201 (1999).
  One-step replica symmetry breaking for K = 3 committee machines.
- Stojnic's later fully lifted random duality papers (arXiv:2402.05696,
  arXiv:2402.05719) treat other activations or wide networks; none closes
  K = 3.

## Datasets

None. The result is mathematical.

## Searches run

- Firecrawl research paper search and web search, 2026-11-25, with these
  queries: "storage capacity tree committee machine sign perceptron one hidden
  layer bounds"; "treelike committee machine capacity rigorous lower bound";
  "lower bound storage capacity committee machine nonoverlapping receptive
  fields K=3 rigorous"; "memorization capacity lower bound treelike committee
  machine sign activations achievability algorithm"; "fully lifted random
  duality theory committee machine sign perceptron capacity"; "Exact capacity
  of the wide hidden layer treelike neural networks fully lifted RDT";
  "tree committee machine K=3 storage capacity one-step replica symmetry
  breaking value Barkai Hansel Sompolinsky 1992"; "Mitchison Durbin 1989
  committee machine capacity". Exa was unavailable (connection failure).
- Checked the introductions of arXiv:2312.08244 and arXiv:2312.08257 for
  previous lower bounds; they list upper bounds only.
