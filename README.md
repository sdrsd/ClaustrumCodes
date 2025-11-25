ClaustrumCodes

Computational models of simulating neuronal responses and neural network dynamics within the claustral circuit.
This repository contains the MATLAB source code used to simulate neuronal responses and network dynamics in the claustrum. 
The models are constrained by empirically characterized inhibitory synaptic connectivity, as described in the associated research.


Repository Structure

The codebase is organized primarily by the figures they generate in the corresponding publication.

Main Simulation Scripts

These scripts are the entry points for the simulations. Running them will generate the data/plots for the specific panels in the paper.

Fig1a_Model_global_CC.m: Simulation of the global model dynamics for the CC configuration.

Fig1b_Model_global_CS.m: Simulation of the global model dynamics for the CS configuration.

Fig2a_Model_local_CC.m: Simulation of local circuit dynamics (CC configuration).

Fig2b_Model_local_CS.m: Simulation of local circuit dynamics (CS configuration).

Configuration & Utilities

Entire population settings.xlsx: A configuration file containing the physiological parameters and settings for the entire neuronal population used in the simulations.

Utility/: A folder containing helper functions and subroutines required by the main model scripts.

Getting Started

Prerequisites

MATLAB: The code is written entirely in MATLAB. Ensure you have a recent version installed.

Toolboxes: Check if standard toolboxes (like the Statistics and Machine Learning Toolbox) are required depending on your specific MATLAB installation.

Installation

Clone this repository to your local machine:
git clone [https://github.com/sdrsd/ClaustrumCodes.git](https://github.com/sdrsd/ClaustrumCodes.git)

Add the repository folder (and the Utility subfolder) to your MATLAB path.

Usage

To reproduce the results for a specific figure:

Open MATLAB.

Navigate to the ClaustrumCodes directory.

Run the desired script in the Command Window:

% Example: Run the simulation for Figure 1a

Fig1a_Model_global_CC

Data & Parameters

The simulation parameters are largely defined in Entire population settings.xlsx.
