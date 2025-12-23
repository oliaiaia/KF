# Kalman Filter (KF + EKF) — C++ Implementation

## Description

This repository contains **C++ implementations of classical state estimation filters**:

* **Kalman Filter (KF)** — for linear dynamic systems
* **Extended Kalman Filter (EKF)** — for nonlinear systems with linear approximation

This project is useful for **state estimation tasks** in applications such as:
- Robotics (position, velocity, orientation estimation)
- Sensor data processing
- Navigation and SLAM
- Dynamic system control

## Key Features

* **Pure C++ implementation** with minimal dependencies
* **Comprehensive examples** for both filters:
  * KF for linear models
  * EKF for nonlinear models
* **Step-by-step usage examples** included in the project
* **CMake-based build system** for easy compilation
* **Data visualization** through Matplot++ integration

## Repository Structure

```
KF/
├── data/                    # Test datasets
├── include/                 # Header files
├── matplotplusplus/        # Local visualization library
├── src/                    # Filter implementation source code
├── CMakeLists.txt          # CMake configuration
├── example0.md             # Basic KF usage example
├── example1.md             # EKF usage example
└── example2.md             # Detailed step-by-step breakdown
```

## Building the Project

The project uses **CMake** for building. Follow these steps:

```bash
git clone https://github.com/oliaiaia/KF.git
cd KF
mkdir build && cd build
cmake ..
make
```

Executable files for examples and filter modules will be created in the build directory.

## Usage Examples

The repository includes **detailed markdown examples** (`example0.md`, `example1.md`, `example2.md`) that cover:

* How to **initialize the filters**
* How to **feed data** into the filters
* How to **retrieve state estimates**
* The **key differences** between KF and EKF

These examples serve as both **learning tutorials** and **practical templates** for real-world applications.

## Core Concepts (Reference)

### Kalman Filter (KF)

KF is applicable when the system is linear and noises are Gaussian:

**System equations:**
- State vector: `x_k`
- State transition: `x_{k+1} = A * x_k + B * u_k + w_k`
- Measurement: `z_k = H * x_k + v_k`

**Where:**
- `A` — state transition matrix (n×n)
- `B` — control input matrix (n×m)
- `H` — observation matrix (p×n)
- `w_k` ~ N(0, Q) — process noise (covariance Q)
- `v_k` ~ N(0, R) — measurement noise (covariance R)

**Prediction step:**
```
x̂_{k|k-1} = A * x̂_{k-1|k-1} + B * u_k
P_{k|k-1} = A * P_{k-1|k-1} * A^T + Q
```

**Update step:**
```
K_k = P_{k|k-1} * H^T * (H * P_{k|k-1} * H^T + R)^{-1}
x̂_{k|k} = x̂_{k|k-1} + K_k * (z_k - H * x̂_{k|k-1})
P_{k|k} = (I - K_k * H) * P_{k|k-1}
```

### Extended Kalman Filter (EKF)

EKF extends KF to handle **nonlinear systems**:

**Nonlinear system equations:**
- State transition: `x_{k+1} = f(x_k, u_k) + w_k`
- Measurement: `z_k = h(x_k) + v_k`

**Linearization using Jacobians:**
- State transition Jacobian: `F_k = ∂f/∂x |_{x=x̂_{k-1|k-1}}`
- Observation Jacobian: `H_k = ∂h/∂x |_{x=x̂_{k|k-1}}`

**Prediction step:**
```
x̂_{k|k-1} = f(x̂_{k-1|k-1}, u_k)
P_{k|k-1} = F_k * P_{k-1|k-1} * F_k^T + Q
```

**Update step:**
```
K_k = P_{k|k-1} * H_k^T * (H_k * P_{k|k-1} * H_k^T + R)^{-1}
x̂_{k|k} = x̂_{k|k-1} + K_k * (z_k - h(x̂_{k|k-1}))
P_{k|k} = (I - K_k * H_k) * P_{k|k-1}
```

## Requirements

* C++17 or higher
* CMake 3.10+
* Eigen library for matrix operations
* Matplot++ (included locally for visualization)

## License

[Specify your license here, e.g., MIT, Apache 2.0]