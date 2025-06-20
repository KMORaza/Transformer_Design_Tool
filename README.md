# Transformer Design Tool

The Transformer Design Tool provides a suite of modules for analyzing and designing transformers, focusing on electrical, magnetic, thermal, mechanical, and acoustic performance.

## Functioning

- **Main Application**: Launches the main window, providing access to all analysis modules.
- **Multi-Winding Transformer Design**: Designs transformers with multiple secondary windings, calculating turns, wire gauges, and flux densities.
- **Core and Copper Loss Estimation**: Computes hysteresis, eddy current, and copper losses based on transformer parameters.
- **Core Saturation Analysis**: Evaluates peak flux density and saturation status, visualizing the B-H curve.
- **Core Sizing**: Determines core area using empirical formulas.
- **Interleaving and Layering**: Calculates winding interleaving patterns and layer counts.
- **Harmonic and Frequency Response Analysis (HFRA)**: Simulates mechanical and acoustic responses under harmonic excitation.
- **Leakage and Stray Field Estimation**: Estimates leakage inductance and stray magnetic fields.
- **Magnetic Flux Density Visualization**: Displays flux density waveforms over time.
- **Winding Configuration Optimization**: Optimizes winding parameters to minimize leakage inductance.
- **Ratio and Wire Gauge Calculation**: Computes turns ratios and wire gauges based on currents and voltage regulation.
- **Temperature-Corrected Resistance Modeling**: Calculates winding resistances and copper losses at different temperatures.
- **Transient Thermal Simulation (FEM)**: Performs finite element method (FEM) simulations for transient temperature distribution.
- **Vibration and Acoustic Noise Modeling**: Analyzes vibrations and noise using FEM, focusing on magnetostriction and Lorentz forces.

The tool uses a centralized `Transformer` class to perform core calculations, integrating results across modules for consistency.

## Simulation Logic

Simulations are driven by user inputs (e.g., voltages, frequency, power, material properties) and rely on physics-based models implemented in the `Transformer` class. 

1. **Input Validation**: Ensures all inputs (voltages, frequency, power, efficiency, temperatures) are positive and within realistic ranges.
2. **Parameter Selection**: Uses transformer type, core material, and winding configuration to select appropriate coefficients (e.g., core loss, flux density limits, coupling coefficients).
3. **Core Sizing**: Calculates core area based on power and frequency using empirical formulas.
4. **Turns and Wire Gauge Calculation**: Determines primary and secondary turns and wire gauges based on voltage ratios and current densities.
5. **Loss Calculations**: Computes core losses (hysteresis and eddy current) using the Steinmetz equation and copper losses using temperature-corrected resistances.
6. **Magnetic Analysis**: Calculates flux densities, mutual inductances, leakage inductances, and stray fields.
7. **Thermal Simulation**: Uses FEM to model transient heat transfer, incorporating conduction, convection, and radiation.
8. **Vibration and Noise Simulation**: Employs FEM to analyze mechanical displacements and acoustic noise from magnetic forces.
9. **Result Visualization**: Displays results via text outputs and plots (e.g., flux density waveforms, temperature distributions, vibration maps).

Each module interfaces with the `Transformer` class, passing inputs and receiving a dictionary of results for display and plotting.

## Design Logic

The design logic focuses on creating a transformer that meets user-specified electrical and physical requirements while optimizing performance. 

- **Transformer Type and Core Selection**: Supports various transformer types (Power, Distribution, Isolation, Current, Potential, High-Frequency) and core configurations (Core-Type, Shell-Type). Parameters like core loss coefficient, maximum flux density, and coupling coefficient are adjusted accordingly.
- **Core Area Calculation**: Uses an empirical formula based on power, frequency, and transformer type to size the core, ensuring sufficient magnetic flux capacity.
- **Turns Ratio and Winding Design**: Calculates turns ratios with a regulation factor (1.05 for Power/Distribution, 1.02 otherwise) to account for voltage drops. Wire gauges are selected based on a current density of 2 A/mm².
- **Winding Configuration**: Supports Concentric or Interleaved configurations, affecting leakage inductance and stray fields. Interleaving patterns and layer counts are determined to optimize performance.
- **Loss Minimization**: Balances core and copper losses to achieve the specified efficiency, using material-specific Steinmetz parameters and temperature-corrected resistances.
- **Thermal Management**: Designs for safe operation by simulating temperature rise, ensuring hotspot temperatures remain within limits.
- **Mechanical Stability**: Minimizes vibrations and noise through material selection (e.g., Young’s modulus) and damping considerations.

The design integrates multiple constraints (e.g., saturation limits, fill factor, bobbin utilization) to produce a feasible transformer configuration.

## Code Structure

The codebase is modular, with distinct files for each functional component, ensuring maintainability and reusability. 

- **transformer.py**: The core `Transformer` class, containing the `calculate` method that performs all major computations (core sizing, turns, losses, inductances, flux densities, thermal and mechanical properties). It integrates helper classes for specific calculations.
- **CoreSizing.py**: Implements the `CoreSizing` class, calculating core area using an empirical formula based on power, frequency, and transformer type.
- **RatioAndWireGaugeCalc.py**: Contains the `RatioAndWireGaugeCalc` class, computing turns ratios and selecting wire gauges from an AWG table based on current density.
- **InterleavingAndLayering.py**: Defines the `InterleavingAndLayering` class, determining interleaving patterns and layer counts for windings based on core area and wire gauges.
- **main.py**: Entry point, launching the main application window.
- **main_window.py**: Manages the main GUI, providing access to all analysis modules and displaying basic results.
- **MultiWinding.py**: Implements the multi-winding transformer design module, handling multiple secondary windings and plotting flux densities.
- **CoreLossCopperLossEst.py**: Computes and displays core and copper losses, visualizing loss breakdowns.
- **CoreSaturationAnalyzer.py**: Analyzes core saturation, calculating peak flux density and plotting the B-H curve.
- **MagneticFluxDensity.py**: Visualizes magnetic flux density waveforms for primary and secondary windings.
- **HFRA.py**: Performs harmonic and frequency response analysis, simulating mechanical and acoustic behavior.
- **LeakageAndStrayFieldEst.py**: Estimates leakage inductance and stray magnetic fields, plotting stray field strength.
- **WindingConfigOptimizer.py**: Optimizes winding configurations, minimizing leakage inductance through interleaving factor variations.
- **TempCorrectedResistanceModeling.py**: Models winding resistances and copper losses at ambient and operating temperatures, plotting loss comparisons.
- **ThermalSimulationFEM.py**: Conducts transient thermal simulations using FEM, animating temperature distributions.
- **VibrationAndAcousticNoiseModel.py**: Simulates vibrations and acoustic noise using FEM, visualizing displacement magnitudes.
- **plot.py**: Utility class for plotting magnetic flux density in the main window.

Each module interacts with the `Transformer` class, ensuring consistent calculations across the application.

## Physics and Mathematical Models

The tool employs physics-based models to simulate transformer behavior, using mathematical equations to describe electrical, magnetic, thermal, and mechanical phenomena. 

### Electrical Models
- **Turns Ratio**: Adjusts for voltage regulation:
  `turns_ratio = (Vp / Vs) * regulation_factor`
  where `regulation_factor = 1.05` (Power/Distribution) or `1.02` (others).
- **Current Calculation**:
  `Ip = total_power / (Vp * efficiency)`
  `Is = power / (Vs * efficiency)`
- **Winding Resistance (Temperature-Corrected)**:
  `R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))`
  where `resistivity_cu` = 1.68e-8 Ω·m, `alpha_cu` = 0.00393/°C, `T_ref` = 20°C.
- **Copper Losses**:
  `P_copper = I²R`
  Total copper loss sums primary and secondary contributions.

### Magnetic Models
- **Core Area (Empirical)**:
  `core_area = k * (power / frequency)^0.5`
  where `k` varies by transformer and core type.
- **Primary Turns**:
  `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`
  where `B_max` is the maximum flux density.
- **Secondary Turns**:
  `Ns = Np / turns_ratio`
- **Flux Density**:
  `B(t) = B_max * sin(2 * pi * freq * t)`
  Calculated for primary and secondary windings.
- **Mutual Inductance**:
  `Lp = Np^2 * mu * core_area / l`, `Ls = Ns^2 * mu * core_area / l`, `M = coupling_coeff * sqrt(Lp * Ls)`
  where `mu = mu_r * mu_0`, `mu_0 = 4 * pi * 1e-7 H/m`, `l = 0.1 m`.
- **Leakage Inductance**:
  `L_leak = leakage_factor * mu_0 * N^2 * window_area / bobbin_height`
  Adjusted by interleaving factor for winding configuration.
- **Stray Magnetic Field**:
  `H_stray = stray_field_factor * (Ip * Np) / (2 * pi * (d + core_dimension))`
  Calculated over distances d.

### Core Loss Models
- **Steinmetz Equation**:
  `P_core = k * freq^alpha * B_max^beta * core_volume`
  where `k`, `alpha`, `beta` are material-specific (e.g., Silicon Steel: k=0.1, alpha=1.6, beta=2.0).
- **Hysteresis Loss**:
  `P_hyst` = `k * freq^alpha * B_max^beta * core_volume`
- **Eddy Current Loss**:
  `P_eddy = 0.5 * k * freq^2 * B_max^2 * core_volume`, Total core loss = `P_hyst` + `P_eddy`.

### Thermal Models (FEM)
- **Heat Transfer Equation**:
  `rho * cp * dT/dt = div(k * grad(T)) + q`
  where `rho` is density, `cp` is specific heat, `k` is thermal conductivity, `q` is heat source.
- **Boundary Conditions**:
  `-k * grad(T) · n = h_total * (T - T_amb) + epsilon * sigma * (T^4 - T_amb^4)`
  where `h_total = h_conv + h_rad`, `h_rad = 4 * epsilon * sigma * T_ref^3`, `sigma = 5.67e-8 W/m^2·K^4`.
- **Heat Sources**:
  `q_core = core_loss / core_volume`, `q_winding = copper_loss / winding_volume`
- **Heat Flux**:
  `q = -k * grad(T)`

### Mechanical and Acoustic Models (FEM)
- **Vibration Equation**:
  `[M] * d^2u/dt^2 + [C] * du/dt + [K] * u = F`
  where `M` is mass matrix, `C` is damping matrix, `K` is stiffness matrix, F includes magnetostriction and Lorentz forces.
- **Magnetostriction Strain**:
  `epsilon_m` = `lambda_s * (B / B_s)^2`
  where `lambda_s` is magnetostriction coefficient, `B_s` = 2.0 T.
- **Magnetostriction Stress**:
  `sigma_m = E * epsilon_m`
  where `E` is Young’s modulus.
- **Electromagnetic Force**:
  `F_em = (B^2 * core_area) / (2 * mu_0 * depth)`
- **Lorentz Force**:
  `F_lorentz = I * B_leakage / depth`
- **Acoustic Power**:
  `P_acoustic = rho_air * c_air * S * v_rms^2`
  where `rho_air` = 1.225 kg/m^3, `c_air` = 343 m/s, `S` is surface area, `v_rms` is RMS velocity.
- **Sound Pressure Level**:
  `SPL = 20 * log10(p_rms / p_0)`
  where `p_rms = sqrt((P_acoustic * rho_air * c_air) / (4 * pi * r^2))`, `p_0` = 20e-6 Pa, `r` = 1 m.

## Algorithms Utilized

The tool employs several algorithms to perform calculations and simulations efficiently:

- **Core Sizing Algorithm**:
  - Uses an empirical formula to compute core area based on power and frequency.
  - Adjusts coefficients based on transformer and core type.
- **Turns and Wire Gauge Selection**:
  - Calculates turns ratios with a regulation factor.
  - Selects wire gauges from an AWG table using a current density of 2 A/mm², choosing the smallest AWG meeting the area requirement.
- **Interleaving and Layering Algorithm**:
  - Determines interleaving patterns (e.g., primary-secondary alternation) based on transformer type and winding counts.
  - Calculates layer counts using core area and wire diameters.
- **Finite Element Method (FEM) for Thermal Simulation**:
  - Discretizes the 2D domain into a 20x10 triangular mesh.
  - Assembles stiffness (`K`) and capacitance (`C`) matrices using element-wise contributions.
  - Applies boundary conditions for convection and radiation.
  - Solves the transient heat equation using implicit time-stepping and sparse linear solvers (spla.spsolve).
  - Computes heat fluxes from temperature gradients.
- **Finite Element Method (FEM) for Vibration Simulation**:
  - Uses a similar 20x10 triangular mesh.
  - Assembles mass (`M`), damping (`C`), and stiffness (`K`) matrices based on material properties (Young’s modulus, Poisson’s ratio, density).
  - Applies fixed boundary conditions at `y=0`.
  - Solves the harmonic vibration equation using a sparse least-squares solver (spla.lsqr) for complex displacements.
  - Computes acoustic noise from surface velocities.
- **Steinmetz Loss Calculation**:
  - Applies the Steinmetz equation to compute hysteresis and eddy current losses.
  - Uses material-specific parameters for accuracy.
- **Leakage Inductance Optimization**:
  - Iterates over interleaving factors (0.5 to 2.0) to compute leakage inductance.
  - Selects the configuration minimizing leakage inductance.
- **Plotting Algorithms**:
  - Uses Matplotlib for static plots (e.g., flux density, loss bar charts, B-H curves).
  - Implements animations for transient thermal simulations, updating temperature contours and hotspot positions.

These algorithms ensure accurate and efficient computation of transformer parameters, leveraging sparse matrix operations for FEM simulations and empirical models for rapid design iterations.
