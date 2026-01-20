# Physic_project_double_pendule — README
https://github.com/Algath/Physic_project_double_pendulum

This repository implements a double pendulum simulation, parameter estimation from tracker data, and basic analysis tools in Julia.

## Current project structure

```
/
├── data/                  # Optional input datasets (CSV, etc.)
├── assets/                # Tracking exports, videos and images
│   ├── comparison_data.csv
│   └── First_Video_2s.trk
├── imports.jl             # Centralized package imports
├── setup.jl               # One-time package installation script
├── pendule.jl             # Main simulation + optimization script
├── statistiques.jl        # Post-processing / statistical analysis
├── theorie.md             # Theory, derivations and equations
└── README.md              # This file
```

## Files overview

- `setup.jl`: Install required Julia packages (run once).
- `imports.jl`: Central place to `using` / `import` packages needed by scripts.
- `pendule.jl`: Loads tracker `.trk` data, extracts angles, defines the equations of motion, solves the ODEs and runs parameter optimization.
- `statistiques.jl`: Scripts for statistical analysis and comparisons between tracked and simulated trajectories.
- `theorie.md`: Theoretical background, derivations and comparisons (now contains the equations used in `pendule.jl`).
- `assets/`: Tracking exports (Tracker `.trk`), CSV comparison data and images used for documentation and plotting.
- `data/`: Place for additional input datasets if needed.

## Quick start

1. Install Julia (recommended 1.6+ or newer) and open a REPL in the project folder.
2. Run the setup script once to install dependencies:

```bash
julia --project=.
include("setup.jl")
```

3. In the REPL, load common imports and run the simulation:

```julia
include("imports.jl")
include("pendule.jl")
```

4. Run statistical analysis after simulation (optional):

```julia
include("statistiques.jl")
```

## Notes

- The main numeric model is implemented in `pendule.jl` using `DifferentialEquations.jl` and an optimization step with `Optim.jl` to estimate masses and damping.
- Tracker input is expected in `assets/First_Video_2s.trk`. Adjust `pendule.jl` if you have different file names or folder layout.

---

See `theorie.md` for detailed theory and equations used in the implementation.
