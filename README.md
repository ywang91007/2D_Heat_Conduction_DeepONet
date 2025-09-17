# 2D Heat Conduction with DeepONet

This repository contains an implementation of a **Deep Operator Network (DeepONet)** using **PyTorch** to solve the two-dimensional steady-state heat conduction problem under varying boundary conditions and geometric parameters. The project demonstrates how machine learning can be applied to accelerate the solution of partial differential equations (PDEs).

---

## Project Overview
- Implemented a **DeepONet** in PyTorch to predict the temperature distribution in 2D steady-state heat conduction problems.  
- Generated datasets by combining **analytical solutions** and initial conditions using MATLAB scripts.  
- Trained and validated the network with the generated data, comparing predictions against analytical solutions.  
- Visualized results to evaluate **accuracy, generalization**, and **computational efficiency**.  

---

## Features
- **Custom PyTorch dataset class** for handling training data.  
- Implementation of **DeepONet architecture** with branch and trunk networks.  
- **Training loop** with loss monitoring.  
- **Visualization** of predicted vs analytical temperature fields.  

---
## Results
- Predicted temperature fields closely match analytical solutions.
- DeepONet demonstrates strong generalization across varying boundary conditions and geometric configurations.
- Significant reduction in computation time compared to traditional solvers.
