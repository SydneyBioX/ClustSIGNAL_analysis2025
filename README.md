# ClustSIGNAL_analysis2025

This repository contains code for data analysis and figures reported in the manuscript "**ClustSIGNAL identifies cell types and subtypes using an1 adaptive smoothing approach for scalable spatial clustering**".

## Data analysis

-   **Simulation analysis** - see Simulation_scripts folder.

    -   Simulated_dataAnalysis: simulated data generation, parameter testing, benchmarking adaptive smoothing, stress testing - sparsity and segmentation errors, assessing clustering stability.

    -   simulation_functions: functions supporting simulation analysis.

-   Benchmarking analysis - see Benchmarking_scripts folder.

    -   Benchmarking_dataAnalysis: benchmarking runs 4 methods (ClustSIGNAL, BANKSY, BASS, SpatialPCA) on 4 real-world datasets (SeqFISH mouse embryo, Xenium breast cancer, CosMx lung cancer, MERFISH mouse hypothalamus).

    -   benchmarking_functions: wrapper functions for methods.

-   Main figures - see Figures folder.

    -   Fig2\_ simulations: simulation figure panels.

    -   Fig3_benchmarking: benchmarking figure panels.

    -   Fig4_dataAnalysis: data analysis figure panels, and related analyses.
