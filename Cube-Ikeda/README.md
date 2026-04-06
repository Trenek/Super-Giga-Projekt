## Overview

The idea of this part of the project is to apprximate a solution to the cubic Ikeda equation.\
A simple pipeline of the code:

1. **C++ program** (`main.cpp`) generates data and writes it to a CSV file.
2. **Python script** (`plot_traj.py`) reads the CSV and generates a plot (`PNG`).

All generated files are stored in the `output/` folder.


## Setup

For plotting the trajectory you need python libraries provided in requirements. \
For installation:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```


## Run

To approximate and visualise the trajectory, it is enough to run
```bash
bash run.sh
```

## Output
Current version of the visualisation is provided in the `trajectory.png` file.
