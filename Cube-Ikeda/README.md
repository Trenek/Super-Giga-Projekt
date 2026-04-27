# for now description not up to date, but 'Build & run' is right


# Cubic Ikeda DDE — Bifurcation & Poincaré Analysis

**Cubic Ikeda delay differential equation (DDE)**:

```
x'(t) = a · x(t−τ) · (1 − x(t−τ)²)
```

---

## Project structure

```
.
├── main.cpp                  # entry point — attractor visualisation, bifurcation diagram & Poincaré map
├── computations/
│   ├── solveOde.cpp/.h       # ODE solver helpers, solution curve extraction
│   └── toOdeFuncs.cpp/.h     # pseudospectral approximation matrix
├── plot_traj.py              # Python plotting of output/curve.csv
├── venv/                     # Python virtual environment
├── output/                   # generated files (created automatically)
└── Makefile
```

---
## Setup

For plotting the trajectory you need python libraries provided in requirements. \
For installation:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

If `capd-config` is not on your `PATH`, set `CAPDBINDIR` in the Makefile to the directory containing it:

```makefile
CAPDBINDIR = /path/to/capd/build/bin/
```

## Build & run

### Using Make 

```bash
make          # compile only
make run      # compile + run C++ binary
make plot     # compile + run + generate Python plots
make clean    # remove binary and object files
```
