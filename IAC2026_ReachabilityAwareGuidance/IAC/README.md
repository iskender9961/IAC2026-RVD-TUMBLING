# IAC 2026 Demo: Tumbling-Target Approach and Docking

This folder contains a self-contained, reproducible Python demo scenario for:

- nonlinear truth propagation (ECI, two-body + J2),
- HCW-based receding-horizon guidance in LVLH,
- a time-varying line-of-sight docking corridor rotated by target tumbling,
- low-thrust spiral approach from up to 200 m initial separation,
- synchronized rotating stand-off hold at 0.5 m with anti-overshoot behavior,
- empirical safe-region mapping over initial conditions and thrust-to-mass ratios.

All code and generated artifacts stay under `IAC/`.

## Run

From repository root:

```bash
python IAC/run_all.py
```

Optional full sweep:

```bash
python IAC/run_all.py --mode full
```

Optional sweep controller mode:

```bash
python IAC/run_all.py --controller-mode qp

The default quick mode uses:

- controller sample time `Ts = 0.1 s`,
- LOS corridor parameters `x0 = 1.5 m`, `c_x = 1.0`,
- tumbling rate `1 deg/s`,
- thrust-to-mass sweep `{0.01, 0.05, 0.10, 0.20} m/s^2`.
```

## Generated Outputs

After running, the pipeline writes:

- `IAC/data/sweep_results.npz`
  - safe maps, dock maps, margins, docking times, grid and parameter axes
- `IAC/data/summary.json`
  - key stats and paths to generated artifacts
- `IAC/data/trajectory_case_*.npz`
  - illustrative trajectory time histories
- `IAC/figures/*.png`
  - safe-region maps
  - trajectory overlays with rotating corridor snapshots
  - feasibility margin over time
  - docking performance vs thrust-to-mass ratio
- `IAC/animations/approach_3d.mp4` or `.gif` (fallback)
- `IAC/animations/approach_topdown.mp4` or `.gif` (fallback)
- if video writers are unavailable: `IAC/animations/frames/*`

## Notes on Guarantees

- Current `compute_certified_region(...)` is an **empirical inner approximation** (binary erosion of sampled-safe grid cells).
- It is conservative in practice but **not** a formal proof-based certified set.
- The theorem-based certified-set API is kept as a planned extension and documented in `IAC/note.tex`.

## Folder Structure

```text
IAC/
  README.md
  run_all.py
  config.py
  dynamics.py
  constraints.py
  controller_mpc.py
  reachability.py
  plots.py
  animate.py
  note.tex
  data/
  figures/
  animations/
```
