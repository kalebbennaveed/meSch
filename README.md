# meSch: Multi-Agent Energy-Aware Scheduling for Task Persistence

[![IROS 2025](https://img.shields.io/badge/IROS-2025-blue)](https://www.iros.org/)

**meSch** is a scheduling algorithm for a team of autonomous robots that operate on long-term persistent tasks. The proposed framework accounts for the limited battery capacity of the robots and ensures that the robots return to charge their batteries **one at a time** at the single charging station.

[Demo Video](https://www.youtube.com/watch?v=8Hn8znYLPk0)

## Summary

The protocol is:
- **Applicable** to general nonlinear robot models under certain assumptions
- **Does not require** robots to be deployed at different times
- **Handles** robots with different discharge rates
- **Supports** mobile charging stations with state uncertainty

Feasibility of the algorithm for ensuring persistent charging is established under certain assumptions. Efficacy is validated through **simulation** and **hardware experiments**.

## Requirements

- [Julia](https://julialang.org/) 1.6 or later

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/meSch.git
   cd meSch
   ```

2. Start Julia and activate the project:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Usage

### Running Simulations

The code includes Jupyter notebooks for different experimental scenarios:

| Notebook | Description |
|----------|-------------|
| `notebooks/meSch-2Q.ipynb` | Two quadrotors with fixed charging station |
| `notebooks/meSch-3Q.ipynb` | Three quadrotors with fixed charging station |
| `notebooks/meSch-mobile.ipynb` | Two quadrotors with mobile charging station (rover) |
| `notebooks/benchmarking.ipynb` | Benchmarking experiments |

Run the notebooks from the project root. The first cell activates the project and installs dependencies.

### Programmatic Usage

```julia
# Load dependencies and experiment module
include("src/utils.jl")
include("src/ExpDynamicsLibrary.jl")
include("src/RefTrajectoryLibrary.jl")
include("src/ExpControllerLibrary.jl")
include("src/ExpmeSchFlylab_2.jl")  # or ExpmeSchFlylab_3, ExpmeSchDI2_mobile

# Run 2-quad simulation (see notebook for initial conditions)
sol, params = meSchExpDIFL.simulate(quad_ic, t0, t_max, discharge_rate)
```

## Project Structure

```
meSch/
├── src/
│   ├── meSch.jl              # Main module
│   ├── ExpmeSchFlylab_2.jl   # 2-quad meSch experiments
│   ├── ExpmeSchFlylab_3.jl   # 3-quad meSch experiments
│   ├── ExpmeSchDI2_mobile.jl # Mobile charging station experiments
│   ├── ExpDynamicsLibrary.jl # Robot dynamics (quadrotor, double-integrator)
│   ├── ExpControllerLibrary.jl # MPC and geometric controller
│   ├── RefTrajectoryLibrary.jl # Reference trajectories (Lissajous, figure-8)
│   ├── utils.jl              # Utilities (quaternions, etc.)
│   ├── Visualization.jl      # MeshCat visualization
│   └── utils/                # Mesh assets for visualization
├── notebooks/
│   ├── meSch-2Q.ipynb
│   ├── meSch-3Q.ipynb
│   ├── meSch-mobile.ipynb
│   └── benchmarking.ipynb
├── Project.toml
└── README.md
```

## Robot Model

- **High-level:** Double integrator model with Convex MPC for planning
- **Low-level:** Quadrotor dynamics with geometric controller tracking
- **Battery:** Linear discharge model with state of charge (SoC)

## License

See [LICENSE](LICENSE) for details.

---

## Publication & Citation

If you use this code in your research, please cite:

> **meSch: Multi-Agent Energy-Aware Scheduling for Task Persistence**  
> Kaleb Ben Naveed, An Dang, Rahul Kumar, Dimitra Panagou  
> *2025 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)*  
> Pages 21458–21465 | DOI: [10.1109/IROS60139.2025.11247605](https://doi.org/10.1109/IROS60139.2025.11247605)

**BibTeX:**
```bibtex
@INPROCEEDINGS{11247605,
  author={Ben Naveed, Kaleb and Dang, An and Kumar, Rahul and Panagou, Dimitra},
  booktitle={2025 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)},
  title={meSch: Multi-Agent Energy-Aware Scheduling for Task Persistence},
  year={2025},
  pages={21458-21465},
  doi={10.1109/IROS60139.2025.11247605}
}
```
