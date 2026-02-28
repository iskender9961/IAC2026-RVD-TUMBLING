# Reachability Example: Double Integrator

Pedagogical toy examples demonstrating **nominal**, **stochastic**, and **robust** reachability analysis on a 2D double integrator system.

## System

```
x_{k+1} = A x_k + B u_k (+w_k)

A = [1 1; 0 1],  B = [0.5; 1],  |u| ≤ 1
State: x = [position; velocity]
```

## Scripts

| Script | Description |
|--------|-------------|
| `nominal_backward.m` | Backward reachable set (safe-start region) to target, **no disturbance** |
| `nominal_forward.m` | Forward reachable set from target, **no disturbance** |
| `stochastic_backward.m` | Backward reachable set with **Gaussian noise** chance constraints (95%) |
| `stochastic_forward.m` | Forward reachable set with **Gaussian noise** tightening |
| `robust_backward.m` | Backward reachable set with **bounded worst-case** disturbance |
| `robust_forward.m` | Forward reachable set with **bounded worst-case** disturbance |
| `plot_overlay_comparison.m` | Overlay of all three backward sets showing nesting |
| `common_utils.m` | Shared parameters and colors |

## Forward vs Backward

- **Forward**: Starting from a known initial set, what states can be reached in N steps?
  - R_{k+1} = (A * R_k ⊕ B * U) ∩ X_state

- **Backward**: Which initial states can reach the target in N steps?
  - S_k = Pre(S_{k+1}) ∩ X_state,  where Pre(S) = {x : ∃ u ∈ U s.t. Ax + Bu ∈ S}

## Nominal vs Stochastic vs Robust

| Method | Disturbance | Guarantee | Size |
|--------|------------|-----------|------|
| **Nominal** | None | Deterministic | Largest |
| **Stochastic** | w ~ N(0,W) | P(safe) ≥ 95% | Middle |
| **Robust** | \|w\| ≤ w_max | For ALL w | Smallest |

The sets satisfy the **inclusion hierarchy**:

```
X_robust ⊆ X_stochastic ⊆ X_nominal
```

## What is Approximate vs Certified

- **Nominal**: Exact (certified) if dynamics are perfect and disturbance-free
- **Stochastic**: Certified inner approximation (Bonferroni + Gaussian quantile)
- **Robust**: Certified inner approximation (conservative row-by-row support function)

The stochastic and robust sets use conservative elimination of the "∃ u" quantifier,
so they are inner approximations of the true reachable sets. This means:
- Any state declared safe IS safe (under the respective assumptions)
- Some safe states may be missed (conservatism)

## Connection to Rendezvous/Docking (CWH)

These toy scripts use a 2D double integrator. The main project uses the full 3D
**Clohessy-Wiltshire-Hill (CWH)** relative dynamics with:
- 6-state model (3D position + velocity)
- Rotating LOS cone constraints (body-frame corridor)
- Directional per-constraint erosion model
- Synchronization range bound r < 2*a_max/ω²

The mathematical framework is identical:
- H-representation polytope propagation
- Support function tightening for ∃u and ∀w
- Bonferroni chance constraints for stochastic case

The CWH reachability modules in `reachability_nominal/`, `reachability_stochastic/`,
and `reachability_robust/` implement these concepts for the actual spacecraft problem.

## Running

From MATLAB (in the project root):
```matlab
cd Example
nominal_backward
stochastic_backward
robust_backward
plot_overlay_comparison
```

Figures are saved to `Example/figures/`.
