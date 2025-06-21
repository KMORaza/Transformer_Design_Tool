# Transformer Design Tool

The Transformer Design Tool provides a suite of modules for analyzing and designing transformers, focusing on electrical, magnetic, thermal, mechanical, and acoustic performance.

## _Functioning_

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

## _Simulation Logic_

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

## _Design Logic_

The design logic focuses on creating a transformer that meets user-specified electrical and physical requirements while optimizing performance. 

- **Transformer Type and Core Selection**: Supports various transformer types (Power, Distribution, Isolation, Current, Potential, High-Frequency) and core configurations (Core-Type, Shell-Type). Parameters like core loss coefficient, maximum flux density, and coupling coefficient are adjusted accordingly.
- **Core Area Calculation**: Uses an empirical formula based on power, frequency, and transformer type to size the core, ensuring sufficient magnetic flux capacity.
- **Turns Ratio and Winding Design**: Calculates turns ratios with a regulation factor (1.05 for Power/Distribution, 1.02 otherwise) to account for voltage drops. Wire gauges are selected based on a current density of 2 A/mm².
- **Winding Configuration**: Supports Concentric or Interleaved configurations, affecting leakage inductance and stray fields. Interleaving patterns and layer counts are determined to optimize performance.
- **Loss Minimization**: Balances core and copper losses to achieve the specified efficiency, using material-specific Steinmetz parameters and temperature-corrected resistances.
- **Thermal Management**: Designs for safe operation by simulating temperature rise, ensuring hotspot temperatures remain within limits.
- **Mechanical Stability**: Minimizes vibrations and noise through material selection (e.g., Young’s modulus) and damping considerations.

The design integrates multiple constraints (e.g., saturation limits, fill factor, bobbin utilization) to produce a feasible transformer configuration.

## _Code Structure_

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

## _Physics and Mathematical Models_

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

## _Algorithms Utilized_

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

---

## ① Transformer Types

The Transformer Design Tool supports six transformer types, each with distinct applications and design considerations. The `Transformer` class in `transformer.py` centralizes calculations, adjusting parameters based on the selected type to model specific electrical, magnetic, thermal, and mechanical behaviors. 

- **Power Transformer**: Core loss coefficient = 0.02, max flux density = 1.5 T, coupling coefficient = 0.98
- **Distribution Transformer**: Core loss coefficient = 0.015, max flux density = 1.4 T, coupling coefficient = 0.97
- **Isolation Transformer**: Core loss coefficient = 0.01, max flux density = 1.3 T, coupling coefficient = 0.95
- **Current Transformer**: Core loss coefficient = 0.005, max flux density = 1.0 T, coupling coefficient = 0.99
- **Potential Transformer**: Core loss coefficient = 0.008, max flux density = 1.1 T, coupling coefficient = 0.98
- **High-Frequency Transformer**: Core loss coefficient = 0.03, max flux density = 0.8 T, coupling coefficient = 0.90

These parameters influence core sizing, loss calculations, and electromagnetic performance, ensuring type-specific accuracy.

### Power Transformers

**Description**: Used in power generation and transmission to step up/down voltages at high power levels (>100 kVA). Designed for high efficiency and minimal losses.

**Support**:
- **Core Sizing**: Larger core areas to handle high power, calculated using an empirical formula in `CoreSizing.py`.
- **Loss Calculations**: Emphasizes low core and copper losses, using higher core loss coefficient (0.02).
- **Turns Ratio**: Applies a 1.05 regulation factor to account for significant voltage drops, in `RatioAndWireGaugeCalc.py`.
- **Thermal and Vibration Analysis**: Robust thermal simulations in `ThermalSimulationFEM.py` and vibration analysis in `VibrationAndAcousticNoiseModel.py` to manage high power-induced stresses.

**Algorithms**:
- **Core Sizing**: Empirical formula adjusts core area based on high power and 50/60 Hz frequency.
- **Loss Optimization**: Iterates interleaving factors in `WindingConfigOptimizer.py` to minimize leakage inductance and losses.
- **FEM Simulations**: Uses sparse matrix solvers for thermal and vibration analyses, handling large meshes for large cores.

**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, where `k` adjusted for high power.
- **Turns**: `Np = (Vp * 10^8) / (4.44 * freq * core_area * B_max)`, `Ns = Np / (turns_ratio * 1.05)`
- **Core Loss (Steinmetz)**: `P_core = k * freq^1.6 * B_max^2 * core_volume`, `k = 0.1` for Silicon Steel.
- **Copper Loss**: `P_copper = I^2 * R`, with `R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))`, `resistivity_cu = 1.68e-8`, `alpha_cu = 0.00393`, `T_ref = 20°C`.

### Distribution Transformers

**Description**: Steps down voltages in distribution networks (e.g., 11 kV to 400V), typically 10-100 kVA. Focuses on reliability and efficiency at medium power levels.
**Support**:
- **Core Design**: Balances core size and cost, using a slightly smaller core area than Power transformers.
- **Loss Calculations**: Lower core loss coefficient (0.015) reflects efficient core materials.
- **Voltage Regulation**: Uses 1.05 regulation factor in `RatioAndWireGaugeCalc.py` for load variations.
- **Winding Optimization**: Supports Concentric and Interleaved configurations in `WindingConfigOptimizer.py` to reduce leakage inductance.
**Algorithms**:
- **Core Sizing**: Empirical formula in `CoreSizing.py` optimizes core area for medium power and frequency.
- **Wire Selection**: Selects AWG sizes based on 2 A/mm² current density, prioritizing cost-effective gauges.
- **Interleaving Optimization**: Iterates interleaving factors to minimize leakage inductance and improve efficiency.
**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, with `k` tuned for medium power.
- **Turns Ratio**: `turns_ratio = (Vp / Vs) * 1.05`
- **Core Loss**: `P_core = k * freq^1.6 * B_max^2 * core_volume`, `k = 0.1`, `B_max = 1.4 T`.
- **Leakage Inductance**: `L_leak = leakage_factor * mu_0 * N^2 * window_area / bobbin_height`, `leakage_factor = 0.05` (Concentric) or `0.02` (Interleaved), `mu_0 = 4 * pi * 10^-7`.

### Isolation Transformers

**Description**: Provides galvanic isolation for safety and noise reduction, often used in medical or sensitive electronic equipment. Operates at low to medium power.
**Support**:
- **Core Sizing**: Smaller cores due to lower power, calculated in `CoreSizing.py`.
- **Loss Calculations**: Lowest core loss coefficient (0.01) for high efficiency at low flux density (1.3 T).
- **Winding Design**: Emphasizes Interleaved configurations in `InterleavingAndLayering.py` to minimize leakage and stray fields.
- **Stray Field Analysis**: Detailed stray field calculations in `LeakageAndStrayFieldEst.py` to ensure low electromagnetic interference.
**Algorithms**:
- **Core Sizing**: Empirical formula reduces core area for low power.
- **Winding Optimization**: Prioritizes interleaving to reduce leakage inductance.
- **Stray Field Calculation**: Computes field strength over distances to find safe zones.
**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, `k` minimized for low power.
- **Turns**: `Np = (Vp * 10^8) / (4.44 * freq * core_area * 1.3)`, `Ns = Np / (turns_ratio * 1.02)`.
- **Stray Field**: `H_stray = stray_field_factor * (Ip * Np) / (2 * pi * (d + core_dimension))`, `stray_field_factor = 0.1` (Concentric) or `0.05` (Interleaved).
- **Mutual Inductance**: `M = coupling_coeff * sqrt(Lp * Ls)`, `coupling_coeff = 0.95`, `Lp = Np^2 * mu * core_area / l`, `l = 0.1 m`.

### Current Transformers

**Description**: Measures high currents by stepping down to a measurable level, used in metering and protection. Operates at low flux density to avoid saturation.
**Support**:
- **Core Sizing**: Small core areas due to low power, in `CoreSizing.py`.
- **Flux Density**: Lowest max flux density (1.0 T) to prevent saturation, enforced in `CoreSaturationAnalyzer.py`.
- **High Coupling**: Highest coupling coefficient (0.99) for accurate current measurement.
- **Turns Ratio**: Uses 1.02 regulation factor for precision, in `RatioAndWireGaugeCalc.py`.
**Algorithms**:
- **Core Sizing**: Empirical formula for minimal core area.
- **Saturation Check**: Compares peak flux density to saturation limit (`B_sat = 1.8 T` for Silicon Steel).
- **Wire Selection**: Selects fine gauges for secondary to handle low currents.
**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, `k` minimized.
- **Turns**: `Np = (Vp * 10^8) / (4.44 * freq * core_area * 1.0)`, `Ns = Np / (turns_ratio * 1.02)`.
- **Core Loss**: `P_core = k * freq^1.6 * B_max^2 * core_volume`, `k = 0.1`, `B_max = 1.0 T`.
- **Saturation**: `B_peak = Vp / (4.44 * freq * Np * core_area)`, compared to `B_sat = 1.8 T`.

### Potential Transformers

**Description**: Steps down high voltages for measurement, used in metering and protection. Designed for high accuracy and low losses.
**Support**:
- **Core Sizing**: Small cores for low power, in `CoreSizing.py`.
- **Loss Calculations**: Low core loss coefficient (0.008) and flux density (1.1 T) for precision.
- **Coupling**: High coupling coefficient (0.98) for accurate voltage measurement.
- **Winding Design**: Supports fine wire gauges in `RatioAndWireGaugeCalc.py` for low-current secondaries.
**Algorithms**:
- **Core Sizing**: Empirical formula for small core area.
- **Wire Selection**: Chooses AWG sizes for low secondary currents.
- **Loss Calculation**: Uses Steinmetz equation with low flux density.
**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, `k` small.
- **Turns**: `Np = (Vp * 10^8) / (4.44 * freq * core_area * 1.1)`, `Ns = Np / (turns_ratio * 1.02)`.
- **Core Loss**: `P_core = k * freq^1.6 * B_max^2 * core_volume`, `k = 0.1`, `B_max = 1.1 T`.
- **Mutual Inductance**: `M = 0.98 * sqrt(Lp * Ls)`, `Lp = Np^2 * mu * core_area / l`.

### High-Frequency Transformers

**Description**: Operates at high frequencies (kHz-MHz) in switch-mode power supplies and RF applications. Designed for low core losses and high efficiency at low flux density.
**Support**:
- **Core Material**: Prefers Ferrite (`B_sat = 0.5 T`, `mu_r = 2000`) in `transformer.py`.
- **Flux Density**: Lowest max flux density (0.8 T) to minimize core losses at high frequencies.
- **Loss Calculations**: Higher core loss coefficient (0.03) and Steinmetz parameters (`k=0.05`, `alpha=1.5`, `beta=2.5`) for Ferrite.
- **Winding Design**: Emphasizes Interleaved configurations in `WindingConfigOptimizer.py` to reduce leakage inductance.
- **HFRA**: Detailed harmonic analysis in `HFRA.py` for high-frequency response.
**Algorithms**:
- **Core Sizing**: Adjusts empirical formula for high frequency.
- **Loss Calculation**: Uses Ferrite-specific Steinmetz parameters.
- **Interleaving Optimization**: Minimizes leakage inductance for high-frequency operation.
- **Harmonic Analysis**: Simulates mechanical and acoustic responses at high frequencies.
**Physics and Mathematical Models**:
- **Core Area**: `core_area = k * (power / frequency)^0.5`, `k` adjusted for high frequency.
- **Turns**: `Np = (Vp * 10^8) / (4.44 * freq * core_area * 0.8)`, `Ns = Np / (turns_ratio * 1.02)`.
- **Core Loss**: `P_core = k * freq^1.5 * B_max^2.5 * core_volume`, `k = 0.05` (Ferrite).
- **Leakage Inductance**: `L_leak = 0.02 * mu_0 * N^2 * window_area / bobbin_height` (Interleaved).
- **Flux Density**: `B(t) = 0.8 * sin(2 * pi * freq * t)`.

---

## Core-Type and Shell-Type Configurations

The Transformer Design Tool supports two core configurations: **Core-Type** and **Shell-Type**, which are selectable via a dropdown menu in each module's GUI (e.g., Multi-Winding, Core Saturation, Thermal Simulation). These configurations influence magnetic flux paths, core sizing, winding arrangements, and loss calculations, impacting overall transformer performance. The key functionalities include:

- **Core-Type Configuration**: Features a single magnetic path where windings are wrapped around the core limbs. It is simpler to construct and typically used for high-voltage applications. The tool adjusts flux density and core area calculations to account for the single-path magnetic circuit.
- **Shell-Type Configuration**: Features a magnetic circuit where the core surrounds the windings, providing two parallel flux paths. This configuration is often used for low-voltage, high-current applications due to its ability to reduce leakage flux. The tool applies a flux density reduction factor and modifies interleaving patterns to optimize performance.

Both configurations are integrated into the `Transformer` class, which adjusts parameters like flux density, core area, and interleaving based on the selected core type. 

### Core Sizing Algorithm
- **Purpose**: Calculates the core cross-sectional area based on power, frequency, and core type.
- **Process**:
  - For **Core-Type**, uses a higher core constant to account for a single magnetic path, ensuring sufficient flux capacity.
  - For **Shell-Type**, applies a lower core constant due to dual flux paths, reducing the required core area.
  - Formula: `core_area = k * (power / frequency)^0.5`
    - `k = 0.01` for Core-Type, `k = 0.008` for Shell-Type (empirical values adjusted for transformer type).
- **Implementation**: Handled in `CoreSizing.calculate_core_area`, called by `Transformer.calculate`.

### Flux Density Adjustment
- **Purpose**: Adjusts maximum flux density to prevent saturation, accounting for core geometry.
- **Process**:
  - Core-Type: Uses full material-specific flux density (e.g., 1.8 T for Silicon Steel).
  - Shell-Type: Applies a reduction factor (`flux_density_factor = 0.9`) to account for distributed flux paths.
  - Ensures flux density does not exceed material saturation limits.
- **Implementation**: Applied in `Transformer.calculate` using `flux_density_factor`.

### Interleaving and Layering Algorithm
- **Purpose**: Determines winding interleaving patterns and layer counts, adjusted for core type.
- **Process**:
  - Core-Type: Assumes concentric winding placement around core limbs, with simpler interleaving patterns (e.g., primary-secondary alternation).
  - Shell-Type: Promotes interleaved configurations to minimize leakage flux due to the core enclosing the windings.
  - Uses core area and wire gauges to calculate layers: `layers = ceil(turns * wire_diameter^2 / window_area)`.
  - Shell-Type configurations prioritize more interleaving to reduce leakage inductance.
- **Implementation**: Handled in `InterleavingAndLayering.calculate`, called by `Transformer.calculate`.

### Loss Calculation Algorithm
- **Purpose**: Computes core and copper losses, adjusted for core type.
- **Process**:
  - Core-Type: Higher core losses due to concentrated flux paths.
  - Shell-Type: Lower core losses due to distributed flux, but potentially higher leakage inductance.
  - Uses Steinmetz equation for core losses and temperature-corrected resistance for copper losses.
- **Implementation**: Integrated in `Transformer.calculate`, with core loss coefficients adjusted by transformer type.

### FEM Simulation Algorithms (Thermal and Vibration)
- **Purpose**: Simulates thermal and mechanical behavior, accounting for core geometry.
- **Process**:
  - Core-Type: Models a rectangular core with windings around limbs, assuming uniform heat and force distribution.
  - Shell-Type: Models a core surrounding windings, adjusting mesh boundaries to reflect the enclosing structure.
  - Uses a 20x10 triangular mesh for FEM, with core region defined based on core_side = `sqrt(core_area)`.
- **Implementation**: In `ThermalSimulationFEM` and `VibrationAndAcousticNoiseModel`, with core dimensions derived from `Transformer.calculate`.

### Magnetic Models
- **Flux Density**:
  - Core-Type: Flux is concentrated in a single path, so `B_max = min(params["max_flux_density"], material["saturation_flux_density"])`.
  - Shell-Type: Flux is split across two paths, so `B_max = min(params["max_flux_density"] * 0.9, material["saturation_flux_density"])`.
  - Equation: `B(t) = B_max * sin(2 * pi * freq * t)`.
- **Primary Turns**:
  - `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`.
  - Core area and B_max are adjusted based on core type.
- **Mutual Inductance**:
  - `Lp = Np^2 * mu * core_area / l`, `Ls = Ns^2 * mu * core_area / l`, `M = coupling_coeff * sqrt(Lp * Ls)`.
  - Shell-Type uses a slightly lower coupling coefficient due to distributed flux paths.
- **Leakage Inductance**:
  - `L_leak = leakage_factor * mu_0 * N^2 * window_area / bobbin_height`.
  - Shell-Type has a lower leakage factor (0.02 for Interleaved) due to reduced leakage flux.
- **Stray Magnetic Field**:
  - `H_stray = stray_field_factor * (Ip * Np) / (2 * pi * (d + core_dimension))`.
  - Shell-Type has a higher stray_field_factor (0.5 for Interleaved) due to core geometry.

### Core Loss Models
- **Steinmetz Equation**:
  - `P_core = k * freq^alpha * B_max^beta * core_volume`.
  - Core-Type: Higher B_max increases losses.
  - Shell-Type: Lower B_max reduces losses.
  - Parameters (e.g., Silicon Steel: `k=0.1`, `alpha=1.6`, `beta=2.0`) are material-specific.
- **Hysteresis Loss**:
  - `P_hyst = k * freq^alpha * B_max^beta * core_volume`.
- **Eddy Current Loss**:
  - `P_eddy = 0.5 * k * freq^2 * B_max^2 * core_volume`.

### Electrical Models
- **Current Calculation**:
  - `Ip = total_power / (Vp * efficiency)`, `Is = power / (Vs * efficiency)`.
  - Identical for both configurations, but Shell-Type may have higher currents due to lower voltage applications.
- **Winding Resistance**:
  - `R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))`.
  - Length depends on core geometry (longer for Core-Type due to limb wrapping).

### Thermal Models (FEM)
- **Heat Transfer**:
  - `rho * cp * dT/dt = div(k * grad(T)) + q`.
  - Core-Type: Heat concentrates in core limbs, modeled with a rectangular core region.
  - Shell-Type: Heat distributes across a larger core area, modeled with an enclosing core geometry.
- **Heat Sources**:
  - `q_core = core_loss / core_volume`, `q_winding = copper_loss / winding_volume`.
  - Core volume is `core_area * 0.1` (0.1 m depth).
- **Boundary Conditions**:
  - `-k * grad(T) · n = h_total * (T - T_amb) + epsilon * sigma * (T^4 - T_amb^4)`.
  - Shell-Type has a larger boundary perimeter due to enclosing core.

### Mechanical and Acoustic Models (FEM)
- **Vibration Equation**:
  - `[M] * d^2u/dt^2 + [C] * du/dt + [K] * u = F`.
  - Core-Type: Forces concentrate on core limbs, leading to higher displacements.
  - Shell-Type: Forces distribute across the core, reducing peak displacements.
- **Magnetostriction Strain**:
  - `epsilon_m = lambda_s * (B / B_s)^2`.
  - Applied primarily in core region, adjusted for B_max differences.
- **Electromagnetic Force**:
  - `F_em = (B^2 * core_area) / (2 * mu_0 * depth)`.
  - Core-Type has higher F_em due to concentrated flux.
- **Lorentz Force**:
  - `F_lorentz = I * B_leakage / depth`.
  - Shell-Type has lower B_leakage, reducing F_lorentz.
- **Acoustic Power**:
  - `P_acoustic = rho_air * c_air * S * v_rms^2`.
  - Shell-Type has a larger surface area (S), potentially increasing noise.
- **Sound Pressure Level**:
  - `SPL = 20 * log10(p_rms / p_0)`, where `p_rms = sqrt((P_acoustic * rho_air * c_air) / (4 * pi * r^2))`.

### Implementation Details

- **Transformer.calculate**:
  - Central method in `transformer.py`, adjusts flux_density_factor (1.0 for Core-Type, 0.9 for Shell-Type) and calls helper classes.
  - Computes core area, turns, losses, inductances, and flux densities, tailoring results to core type.
- **CoreSizing.calculate_core_area**:
  - Uses different k values for Core-Type and Shell-Type, reflecting their magnetic path differences.
- **InterleavingAndLayering.calculate**:
  - Adjusts interleaving patterns, favoring more interleaving for Shell-Type to minimize leakage flux.
- **ThermalSimulationFEM**:
  - Models core geometry based on `core_side = sqrt(core_area)`, with Shell-Type having a larger effective core area in the mesh.
- **VibrationAndAcousticNoiseModel**:
  - Applies core-specific boundary conditions and forces, with Shell-Type distributing forces across a larger core area.

The core-type and shell-type configurations are seamlessly integrated, ensuring accurate modeling of magnetic, thermal, and mechanical behavior tailored to each geometry.

---

## Interleaving and Layering Configuration

The interleaving and layering configuration tools determine the optimal arrangement of primary and secondary windings in a transformer.

- **Interleaving Patterns**: Specifies how primary and secondary windings are alternated to reduce leakage inductance and improve coupling. Common patterns include simple alternation (e.g., primary-secondary-primary) or more complex sequences for multi-winding designs.
- **Layer Counts**: Calculates the number of layers for primary and secondary windings based on core geometry, wire gauges, and bobbin dimensions, ensuring efficient use of the winding window.
- **Optimization**: Evaluates different interleaving factors to minimize leakage inductance, as implemented in the `WindingConfigOptimizer` module.
- **Multi-Winding Support**: Handles multiple secondary windings in the `MultiWinding` module, assigning appropriate turns, gauges, and layers for each winding.

These tools are invoked when users specify transformer parameters (e.g., transformer type, core type, voltages, power, frequency) and winding configurations (Concentric or Interleaved). Results include interleaving patterns, layer counts, and associated electrical properties like leakage inductance and fill factor, displayed via text outputs and plots.

### Interleaving Pattern Determination (`InterleavingAndLayering.calculate`)
- **Input**: Transformer type, core type, primary turns (`Np`), total secondary turns (`Ns_total`), primary wire gauge, secondary wire gauge, core area.
- **Steps**:
  1. Determine the turns ratio: `Np / Ns_total`.
  2. Select interleaving pattern based on transformer type and turns ratio:
     - For High-Frequency transformers, use a pattern with frequent alternation (e.g., `P-S-P-S`).
     - For Power/Distribution, use a balanced pattern (e.g., `P-S-P`).
     - For others, use a simpler pattern (e.g., `P-S`).
  3. Adjust pattern based on core type:
     - Core-Type: Symmetrical alternation.
     - Shell-Type: Asymmetrical to account for flux distribution.
  4. Return the interleaving pattern as a string (e.g., `P-S-P`).
- **Output**: Interleaving pattern string.

### Layer Count Calculation (`InterleavingAndLayering.calculate`)
- **Input**: Same as above, plus core area for window size estimation.
- **Steps**:
  1. Estimate bobbin window area: `window_area = sqrt(core_area) * 0.1`.
  2. Calculate wire diameters from AWG gauges using a predefined table (e.g., AWG 10 = 5.261 mm²).
  3. Compute turns per layer:
     - Primary: `turns_per_layer_p = window_width / wire_diameter_p`.
     - Secondary: `turns_per_layer_s = window_width / wire_diameter_s`.
  4. Calculate required layers:
     - Primary: `primary_layers = ceil(Np / turns_per_layer_p)`.
     - Secondary: `secondary_layers = ceil(Ns_total / turns_per_layer_s)`.
  5. Adjust for interleaving:
     - If Interleaved, increase layer count to accommodate alternation (e.g., `primary_layers *= 2` for `P-S-P`).
- **Output**: Primary and secondary layer counts.

### Leakage Inductance Optimization (`WindingConfigOptimizer.calculate_optimization`)
- **Input**: Transformer type, core type, winding configuration, primary/secondary voltages, frequency, power, efficiency, bobbin dimensions, insulation thickness, primary/secondary layers.
- **Steps**:
  1. Validate inputs (e.g., positive values, efficiency between 0 and 1).
  2. Iterate over interleaving factors (`0.5` to `2.0` in steps of `0.1`).
  3. For each factor:
     - Adjust leakage factor: `leakage_factor_adj = leakage_factor * interleaving_factor`.
     - Call `Transformer.calculate` to compute leakage inductance.
     - Store primary leakage inductance (in µH).
  4. Plot leakage inductance vs. interleaving factor.
  5. Select configuration with minimum leakage inductance.
- **Output**: Optimized interleaving pattern, layer counts, leakage inductance, fill factor, and bobbin utilization.

### Multi-Winding Layer Assignment (`MultiWinding.calculate_design`)
- **Input**: Transformer type, core type, primary voltage, frequency, total power, efficiency, secondary voltages, and powers.
- **Steps**:
  1. Validate secondary winding data (positive voltages/powers, sum of powers ≤ total power).
  2. Call `Transformer.calculate` to obtain results for each secondary winding.
  3. Assign layers for each secondary winding using `InterleavingAndLayering.calculate`.
  4. Generate interleaving pattern to include all windings (e.g., `P-S1-S2-P`).
  5. Compute fill factor and window utilization.
- **Output**: Interleaving pattern, primary/secondary layers, and electrical properties for each winding.

### Winding Geometry
- **Window Area**:
  `window_area = bobbin_width * bobbin_height`
  Represents the available space for windings.
- **Insulation Area**:
  `insulation_area = insulation_thickness * (primary_layers + sum(secondary_layers))`
  Accounts for space occupied by insulation.
- **Winding Area**:
  `winding_area = window_area - insulation_area`
  Ensures positive area with a fallback value of `0.001 m²` if negative.
- **Fill Factor**:
  `fill_factor = (Np * wire_diameter_p^2 * primary_layers + sum(Ns * wire_diameter_s^2 * secondary_layers)) / winding_area`
  Measures winding efficiency in utilizing the window area.

### Leakage Inductance
- **Leakage Factor Adjustment**:
  `leakage_factor_adj = leakage_factor * interleaving_factor`
  where `leakage_factor = 0.05` (Concentric) or `0.02` (Interleaved).
- **Interleaving Factor**:
  `interleaving_factor = 1.0` (Concentric) or `(primary_layers + sum(secondary_layers)) / max(primary_layers, max(secondary_layers))` (Interleaved).
- **Primary Leakage Inductance**:
  `L_leak_p = leakage_factor_adj * mu_0 * Np^2 * window_area / bobbin_height`
  where `mu_0 = 4 * pi * 1e-7 H/m`.
- **Secondary Leakage Inductance**:
  `L_leak_s = leakage_factor_adj * mu_0 * Ns^2 * window_area / bobbin_height`
  Calculated for each secondary winding.

### Wire Gauge and Turns
- **Wire Diameter**:
  Derived from AWG table (e.g., AWG 10 = 5.261 mm², diameter = `sqrt(4 * area / pi)`).
- **Turns per Layer**:
  `turns_per_layer = window_width / wire_diameter`
  where `window_width` is approximated from core area or user input.
- **Layer Count**:
  `layers = ceil(turns / turns_per_layer)`
  Ensures all turns are accommodated.

### Magnetic Coupling
- **Mutual Inductance**:
  `Lp = Np^2 * mu * core_area / l`
  `Ls = Ns^2 * mu * core_area / l`
  `M = coupling_coeff * sqrt(Lp * Ls)`
  where `mu = mu_r * mu_0`, `l = 0.1 m`, and `coupling_coeff` varies by transformer type (e.g., 0.98 for Power).

### Equations

- Window Area:
  `window_area = bobbin_width * bobbin_height`
- Insulation Area:
  `insulation_area = insulation_thickness * (primary_layers + sum(secondary_layers))`
- Winding Area:
  `winding_area = window_area - insulation_area`
- Fill Factor:
  `fill_factor = (Np * wire_diameter_p^2 * primary_layers + sum(Ns * wire_diameter_s^2 * secondary_layers)) / winding_area`
- Interleaving Factor (Interleaved):
  `interleaving_factor = (primary_layers + sum(secondary_layers)) / max(primary_layers, max(secondary_layers))`
- Leakage Factor Adjustment:
  `leakage_factor_adj = leakage_factor * interleaving_factor`
- Primary Leakage Inductance:
  `L_leak_p = leakage_factor_adj * mu_0 * Np^2 * window_area / bobbin_height`
- Secondary Leakage Inductance:
  `L_leak_s = leakage_factor_adj * mu_0 * Ns^2 * window_area / bobbin_height`
- Primary Inductance:
  `Lp = Np^2 * mu * core_area / l`
- Secondary Inductance:
  `Ls = Ns^2 * mu * core_area / l`
- Mutual Inductance:
  `M = coupling_coeff * sqrt(Lp * Ls)`
- Turns per Layer:
  `turns_per_layer = window_width / wire_diameter`
- Layer Count:
  `layers = ceil(turns / turns_per_layer)`

### Integration in Transformer Calculations

The `Transformer` class orchestrates interleaving and layering calculations by:

1. Calling `InterleavingAndLayering.calculate` to obtain default interleaving patterns and layer counts.
2. Allowing user-specified layer counts (e.g., in `WindingConfigOptimizer`) to override defaults.
3. Computing leakage inductance and fill factor based on the interleaving factor and winding geometry.
4. Passing results to modules like `MultiWinding` and `WindingConfigOptimizer` for display and further optimization.

The tools ensure that winding configurations are physically feasible (e.g., positive winding area, realistic fill factor) and electrically optimal (e.g., minimized leakage inductance), supporting a wide range of transformer designs.

---

## Multi-Winding Feature in Transformer Design Tool

The Multi-Winding feature supports the design of transformers with one primary winding and multiple secondary windings, catering to applications requiring multiple output voltages (e.g., power supplies with diverse voltage rails). 

- **User Input**: Users specify transformer type, core configuration, primary voltage, frequency, total power, efficiency, and details for each secondary winding (voltage and power).
- **Design Calculations**: Computes core area, primary and secondary turns, wire gauges, interleaving patterns, layer counts, mutual inductances, core losses, and efficiency.
- **Validation**: Ensures secondary winding powers do not exceed total power and all inputs are valid.
- **Visualization**: Displays magnetic flux density waveforms for each secondary winding over time.
- **Output**: Provides detailed results, including core area, turns, wire gauges, layers, interleaving pattern, core loss, efficiency, and mutual inductances.

The feature integrates with the `Transformer` class to perform calculations, ensuring consistency with other modules.

### Algorithms

- **Input Validation Algorithm**:
  - Checks that primary voltage, frequency, total power, efficiency, and secondary voltages/powers are positive.
  - Ensures the sum of secondary powers does not exceed the total power.
  - Validates that all secondary winding fields are filled.
- **Core Sizing Algorithm**:
  - Uses the `CoreSizing` class to calculate core area based on total power and frequency.
  - Applies an empirical formula adjusted for transformer type and core configuration.
- **Turns and Wire Gauge Calculation**:
  - Iterates over each secondary winding to compute turns ratios using the `RatioAndWireGaugeCalc` class.
  - Calculates primary wire gauge based on total power and secondary wire gauges for each winding.
  - Selects wire gauges from an AWG table using a current density of 2 A/mm².
- **Turns Calculation**:
  - Computes primary turns using the transformer equation.
  - Calculates secondary turns for each winding based on turns ratios.
- **Interleaving and Layering**:
  - Uses the `InterleavingAndLayering` class to determine interleaving patterns and layer counts.
  - Bases calculations on transformer type, core area, primary turns, and sum of secondary turns.
- **Flux Density Plotting**:
  - Generates time-varying flux density waveforms for each secondary winding.
  - Plots waveforms using Matplotlib, with each secondary winding represented by a distinct line.

### Electrical Models
- **Turns Ratio**:
  - Adjusts for voltage regulation to account for losses:
    `turns_ratio = (Vp / Vs) * regulation_factor`
    where `regulation_factor = 1.05` for Power/Distribution transformers, `1.02` otherwise.
- **Current Calculation**:
  - Primary current:
    `Ip = total_power / (Vp * efficiency)`
  - Secondary current for each winding:
    `Is = power / Vs`
- **Winding Resistance**:
  - Uses temperature-corrected resistance for copper losses:
    `R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))`
    where `resistivity_cu = 1.68e-8` Ω·m, `alpha_cu = 0.00393/°C`, `T_ref = 20°C`.
- **Copper Losses**:
  - Calculated as:
    `P_copper = I^2 * R`
    Summed across primary and secondary windings.

### Magnetic Models
- **Core Area**:
  - Empirical formula from `CoreSizing`:
    `core_area = k * (total_power / frequency)^0.5`
    where `k` varies by transformer and core type.
- **Primary Turns**:
  - Derived from the transformer equation:
    `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`
    where `B_max` is the maximum flux density, adjusted by core type and material.
- **Secondary Turns**:
  - For each secondary winding:
    `Ns = Np / turns_ratio`
- **Flux Density**:
  - Time-varying flux density for primary:
    `B(t) = B_max * sin(2 * pi * freq * t)`
  - For each secondary winding:
    `B_sec(t) = B_max * sin(2 * pi * freq * t) * (Vs / Vp)`
- **Mutual Inductance**:
  - Primary inductance:
    `Lp = Np^2 * mu * core_area / l`
  - Secondary inductance for each winding:
    `Ls = Ns^2 * mu * core_area / l`
  - Mutual inductance:
    `M = coupling_coeff * sqrt(Lp * Ls)`
    where `mu = mu_r * mu_0`, `mu_0 = 4 * pi * 1e-7` H/m, `l = 0.1` m, and `coupling_coeff` depends on transformer type.
- **Core Losses**:
  - Uses the Steinmetz equation:
    `P_core = k * freq^alpha * B_max^beta * core_volume`
    where `k`, `alpha`, `beta` are material-specific (e.g., Silicon Steel: `k=0.1`, `alpha=1.6`, `beta=2.0`).
  - Splits into hysteresis and eddy current losses:
    `P_hyst = k * freq^alpha * B_max^beta * core_volume`
    `P_eddy = 0.5 * k * freq^2 * B_max^2 * core_volume`

---

## Magnetic Flux Density Visualization in Transformer Design Tool

The module, implemented in `MagneticFluxDensity.py` and supported by `transformer.py` and `plot.py`, visualizes the magnetic flux density waveforms for primary and secondary windings over time. 

### Functioning

The Magnetic Flux Density Visualization module allows users to analyze the time-varying magnetic flux density in a transformer's primary and secondary windings.

- **User Inputs**: Users specify transformer parameters such as transformer type (e.g., Power, Distribution), core configuration (Core-Type, Shell-Type), primary voltage, secondary voltage(s), frequency, power, efficiency, and time range for visualization.
- **Calculation**: The module computes magnetic flux density waveforms using the `Transformer` class, which calculates flux density based on input parameters and transformer characteristics.
- **Visualization**: Displays waveforms of magnetic flux density (in Tesla) versus time (in milliseconds) for primary and secondary windings using Matplotlib, integrated into a PyQt5 dialog window.
- **Output**: Provides numerical results, including maximum flux density, in a text output panel alongside the plot.

The module is accessible from the main application window (`main_window.py`) and uses a consistent dark-themed interface for user interaction.

### Key Physics Concepts
- **Magnetic Flux Density**: Represents the strength of the magnetic field within the transformer core, measured in Tesla (T). It varies sinusoidally with the applied AC voltage.
- **Faraday’s Law**: Relates the induced voltage to the rate of change of magnetic flux, used to determine the number of turns and flux density.
- **Core Saturation**: The maximum flux density (`B_max`) is constrained by the core material’s saturation limit to prevent saturation.
- **Transformer Equation**: Links voltage, frequency, core area, and turns to flux density.

### Mathematical Models and Equations
1. **Primary Turns Calculation**:
   The number of primary turns (`Np`) is calculated using the transformer equation to ensure the core operates within safe flux density limits:
   `Np = (Vp * 10^8) / (4.44 * freq * core_area * B_max)`
   where:
   - `Vp` is the primary voltage (V),
   - `freq` is the frequency (Hz),
   - `core_area` is the core cross-sectional area (m²),
   - `B_max` is the maximum allowable flux density (T), determined by transformer type and core material (e.g., 1.5 T for Power transformers with Silicon Steel).
2. **Secondary Turns Calculation**:
   Secondary turns (`Ns`) are determined based on the turns ratio, adjusted for voltage regulation:
   `turns_ratio = (Vp / Vs) * regulation_factor`
   `Ns = Np / turns_ratio`
   where:
   - `Vs` is the secondary voltage (V),
   - `regulation_factor` is 1.05 for Power/Distribution transformers, 1.02 for others.
3. **Magnetic Flux Density**:
   The flux density `B(t)` varies sinusoidally over time, proportional to the applied voltage:
   `B(t) = B_max * sin(2 * pi * freq * t)`
   For the primary winding, `B_max` is directly used. For secondary windings, the flux density is scaled by the voltage ratio:
   `B_secondary(t) = B_max * (Vs / Vp) * sin(2 * pi * freq * t)`
   where:
   - `t` is time (s),
   - `pi` is 3.14159.
4. **Core Area**:
   The core area is calculated empirically by the `CoreSizing` class:
   `core_area = k * (power / frequency)^0.5`
   where `k` is a constant dependent on transformer and core type, and `power` is the total power (VA).
5. **Material Parameters**:
   The maximum flux density (`B_max`) is constrained by the core material’s saturation flux density:
   - Silicon Steel: `saturation_flux_density = 1.8 T`
   - Ferrite: `saturation_flux_density = 0.5 T`
   `B_max` is the minimum of the transformer type’s limit (e.g., 1.5 T for Power) and the material’s saturation flux density, adjusted by a core type factor (1.0 for Core-Type, 0.9 for Shell-Type).

### Algorithms Utilized

The module employs algorithms for computation and visualization, leveraging the `Transformer` class for core calculations and Matplotlib for plotting. The key algorithms are:

1. **Flux Density Calculation Algorithm**:
   - **Input Processing**: Retrieves user inputs (transformer type, core type, voltages, frequency, power, efficiency) from the dialog window.
   - **Validation**: Ensures all inputs are positive and efficiency is between 0 and 1.
   - **Core Sizing**: Calls `CoreSizing.calculate_core_area` to compute the core area based on total power and frequency.
   - **Turns Calculation**: Uses the `RatioAndWireGaugeCalc` class to compute turns ratios and the transformer equation to determine primary and secondary turns.
   - **Flux Density Computation**:
     - Generates a time array `t = [i / (freq * time_points) for i in range(time_points)]` with `time_points = 1000` for one cycle.
     - Computes primary flux density: `primary_flux = [B_max * sin(2 * pi * freq * ti) for ti in t]`.
     - Computes secondary flux densities for each secondary winding: `secondary_flux = [B_max * (Vs / Vp) * sin(2 * pi * freq * ti) for ti in t]`.
     - Stores results in a list: `flux_densities = [primary_flux, secondary_flux_1, ...]`.
   - **Output**: Returns flux densities arrays and time points to the dialog for visualization.
2. **Visualization Algorithm**:
   - **Plot Initialization**: Creates a Matplotlib figure and axis within a `FigureCanvasQTAgg` widget in the dialog (`MagneticFluxDensityWindow`).
   - **Data Plotting**:
     - Clears the existing plot.
     - Iterates over secondary flux density arrays (excluding primary, as per `MultiWindingPlot` logic).
     - Plots each waveform using `ax.plot(flux_densities[i], label="Winding {i+1}")` with time converted to milliseconds.
     - Sets plot title (“Magnetic Flux Density (Multiple Windings)”), x-label (“Time (ms)”), y-label (“Flux Density (T)”), and adds a grid and legend.
   - **Rendering**: Calls `canvas.draw()` to update the plot.
   - **Error Handling**: Catches exceptions during plotting, logging errors to the console.
3. **Main Window Plotting Algorithm** (via `plot.py`):
   - **Purpose**: Visualizes flux density in the main window’s plot panel.
   - **Process**:
     - Initializes a Matplotlib figure in the `Plot` class.
     - Computes flux density using the `Transformer` class (similar to above).
     - Plots primary and secondary flux densities waveforms with distinct colors.
     - Updates axis labels, title, and legend, then redraws the canvas.
   - **Integration**: Triggered by the main window’s calculate button, using inputs from the main interface.

### Implementation Details

- **Module**: `MagneticFluxDensity.py` defines the `MagneticFluxDensityWindow` class, a PyQt5 dialog with input fields, a calculate button, and a plot widget (`FluxDensityPlot`). It interfaces with the `Transformer` class to compute flux densities.
- **Core Calculations**: The `Transformer.calculate` method in `transformer.py` handles flux density computation, returning a list of flux density arrays (`flux_densities`) and time points (`time_points`).
- **Plotting**:
  - `FluxDensityPlot` in `MagneticFluxDensity.py` plots secondary winding flux densities waveforms.
  - `MultiWindingPlot` in `MultiWinding.py` reuses similar logic for multi-winding transformers.
  - `Plot` in `plot.py` supports the main window’s visualization.
- **Dependencies**: Relies on `CoreSizing`, `RatioAndWireGaugeCalc`, and `InterleavingAndLayering` for supporting calculations (core area, turns ratios, etc.).

---

## Core Saturation Analysis

The Core Saturation Analysis feature of the Transformer Design Tool evaluates the magnetic core's saturation behavior, determining whether the transformer operates within safe flux density limits. Implemented in `CoreSaturationAnalyzer.py` and supported by the `Transformer` class in `transformer.py`, this module calculates peak flux density, saturation status, and safety margin, visualizing the B-H curve to aid in transformer design. 

### Functioning

The Core Saturation Analysis feature allows users to input transformer parameters and assess the magnetic core's saturation risk.

- **User Inputs**: Users specify transformer type (e.g., Power, Distribution, Isolation), core type (Core-Type, Shell-Type), core material (Silicon Steel, Ferrite), winding configuration (Concentric, Interleaved), primary and secondary voltages, frequency, power, efficiency, and operating temperatures.
- **Calculations**: The module computes:
  - Core area based on power and frequency.
  - Primary and secondary turns using voltage ratios.
  - Peak magnetic flux density (B_max).
  - Saturation status by comparing B_max to the material’s saturation flux density (B_sat).
  - Safety margin as a percentage of the difference between B_sat and B_max.
  - Magnetic field strength (H) corresponding to the flux density.
- **Outputs**: Results are displayed in a text area, including core saturation status, peak flux density, safety margin, and saturation flux density. A plot of the B-H curve visualizes the core’s magnetic behavior, highlighting the operating point relative to saturation.
- **Integration**: The feature interfaces with the `Transformer` class, which performs core calculations, ensuring consistency with other modules (e.g., Multi-Winding, Core Loss Estimation).

The analysis helps engineers ensure the transformer avoids saturation, which could lead to nonlinear behavior, increased losses, and distorted output.

### Magnetic Flux Density Calculation
The magnetic flux density in the core is determined using the transformer equation, which relates voltage, turns, frequency, core area, and flux density. The peak flux density is critical for assessing saturation.

- **Primary Turns**:
  `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`
  where:
  - `Np`: Number of primary turns.
  - `Vp`: Primary voltage (V).
  - `freq`: Operating frequency (Hz).
  - `B_max`: Maximum allowable flux density (T), adjusted for transformer type and core area.
  - `core_area`: Core cross-sectional area (m²).
- **Peak Flux Density**:
  Rearranging the transformer equation:
  `B_max = (Vp * 1e8) / (4.44 * freq * core_area * Np)`
  In practice, `B_max` is constrained by the transformer type’s maximum flux density and the core material’s saturation flux density (`B_sat`).

### Core Saturation Assessment
Saturation occurs when the core’s flux density exceeds the material’s saturation flux density (`B_sat`), causing the core to lose its magnetic permeability, leading to nonlinear behavior. The module compares `B_max` to `B_sat` to determine saturation status.

- **Saturation Flux Density**:
  Defined by core material:
  - Silicon Steel: `B_sat = 1.8 T`
  - Ferrite: `B_sat = 0.5 T`
- **Saturation Status**:
  - If `B_max < B_sat`, the core is unsaturated (safe operation).
  - If `B_max >= B_sat`, the core is saturated (risk of nonlinearity).
- **Safety Margin**:
  `safety_margin = ((B_sat - B_max) / B_sat) * 100`
  Expressed as a percentage, indicating how close the core is to saturation.

### Magnetic Field Strength (H)
The magnetic field strength is calculated to construct the B-H curve, which illustrates the core’s magnetization behavior.

- **Peak Magnetic Field Strength**:
  `H_peak = (Np * Ip) / l`
  where:
  - `Ip`: Primary current (A), calculated as `Ip = total_power / (Vp * efficiency)`.
  - `l`: Effective magnetic path length (m), assumed as 0.1 m.
- **H Variation**:
  `H(t) = H_peak * sin(2 * pi * freq * t)`
  where `t` is time over one cycle.
- **B-H Curve**:
  The B-H curve is approximated using:
  `B(t) = B_max * sin(2 * pi * freq * t)`
  The curve is plotted with `H(t)` on the x-axis and `B(t)` on the y-axis, showing the core’s magnetic response up to `B_max`.

### Core Area Calculation
The core area is computed using an empirical formula implemented in the `CoreSizing` class, ensuring sufficient magnetic flux capacity.

- **Core Area**:
  `core_area = k * (power / frequency)^0.5`
  where `k` is a coefficient dependent on transformer type (e.g., 0.01 for Power, 0.008 for High-Frequency) and core type (Core-Type or Shell-Type).

### Material Properties
The core material defines key magnetic parameters:
- **Silicon Steel**:
  - Saturation flux density: `B_sat = 1.8 T`
  - Relative permeability: `mu_r = 4000`
- **Ferrite**:
  - Saturation flux density: `B_sat = 0.5 T`
  - Relative permeability: `mu_r = 2000`

The effective permeability is:
`mu = mu_r * mu_0`
where `mu_0 = 4 * pi * 1e-7 H/m`.

### Transformer Type Constraints
Each transformer type imposes a maximum flux density limit, adjusted by core type (Shell-Type reduces flux density by a factor of 0.9):

- **Type-Specific Parameters**:
  - Power: `B_max_limit = 1.5 T`, coupling_coeff = 0.98
  - Distribution: `B_max_limit = 1.4 T`, coupling_coeff = 0.97
  - Isolation: `B_max_limit = 1.3 T`, coupling_coeff = 0.95
  - Current: `B_max_limit = 1.0 T`, coupling_coeff = 0.99
  - Potential: `B_max_limit = 1.1 T`, coupling_coeff = 0.98
  - High-Frequency: `B_max_limit = 0.8 T`, coupling_coeff = 0.90

The effective `B_max` is the minimum of the type-specific limit (adjusted for core type) and `B_sat`.

### Algorithms Utilized

The Core Saturation Analysis feature employs the following algorithms to perform calculations efficiently:

- **Core Sizing Algorithm**:
  - Implemented in `CoreSizing.calculate_core_area`.
  - Input: Transformer type, core type, total power, frequency.
  - Process: Applies the empirical formula `core_area = k * (power / frequency)^0.5`, selecting `k` based on transformer and core type.
  - Output: Core area (m²).
- **Turns Calculation Algorithm**:
  - Implemented in `Transformer.calculate`.
  - Input: Primary voltage (`Vp`), frequency (`freq`), core area, maximum flux density (`B_max`).
  - Process: Computes primary turns using `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`. Secondary turns are derived from turns ratios.
  - Output: Primary and secondary turns.
- **Flux Density Calculation Algorithm**:
  - Input: Primary voltage, frequency, core area, primary turns.
  - Process: Calculates peak flux density using the transformer equation. Ensures `B_max` does not exceed material or type-specific limits.
  - Output: Peak flux density (`B_max`).
- **Saturation Assessment Algorithm**:
  - Input: Peak flux density (`B_max`), material saturation flux density (`B_sat`).
  - Process:
    1. Compares `B_max` to `B_sat`.
    2. If `B_max >= B_sat`, flags saturation.
    3. Computes safety margin: `((B_sat - B_max) / B_sat) * 100`.
  - Output: Saturation status, safety margin.
- **Magnetic Field Strength Algorithm**:
  - Input: Primary turns, primary current, magnetic path length.
  - Process: Computes peak `H_peak = (Np * Ip) / l` and time-varying `H(t)` over one cycle.
  - Output: `H_peak`, `H(t)` for B-H curve.
- **B-H Curve Plotting Algorithm**:
  - Implemented in `CoreSaturationAnalyzer.BHCurvePlot.update_plot`.
  - Input: `B_max`, `H_peak`, frequency, time points.
  - Process:
    1. Generates time points over one cycle (`t = [i / (freq * 1000) for i in range(1000)]`).
    2. Computes `B(t) = B_max * sin(2 * pi * freq * t)` and `H(t) = H_peak * sin(2 * pi * freq * t)`.
    3. Plots `B(t)` vs. `H(t)` using Matplotlib, with annotations for saturation point.
  - Output: B-H curve visualization.

---

## Leakage and Stray Field Estimation 

The Leakage and Stray Field Estimation feature analyzes leakage inductance and stray magnetic fields in transformers, critical for assessing efficiency, electromagnetic interference, and design optimization. Implemented in `LeakageAndStrayFieldEst.py` and supported by calculations in `transformer.py`, this feature uses user inputs to compute and visualize results, leveraging physics-based models and numerical algorithms. 

### Functionality

The Leakage and Stray Field Estimation feature enables users to evaluate the leakage inductance of primary and secondary windings and the stray magnetic field strength around the transformer. 

- **Input Parameters**: Users specify transformer type (e.g., Power, Distribution, Isolation), core configuration (Core-Type, Shell-Type), winding configuration (Concentric, Interleaved), primary voltage, secondary voltage, frequency, power, efficiency, bobbin width, bobbin height, insulation thickness, and primary and secondary layer counts.
- **Calculations**:
  - Computes leakage inductance for the primary and each secondary winding, accounting for winding configuration and interleaving effects.
  - Estimates stray magnetic field strength as a function of distance from the transformer core.
- **Outputs**:
  - Displays leakage inductance values in microhenries (µH) for primary and secondary windings.
  - Reports the safe distance where stray field strength falls below 0.1 A/m.
  - Visualizes stray field strength versus distance using a plot generated by Matplotlib.
- **Integration**: Relies on the `Transformer` class in `transformer.py` to perform core calculations, ensuring consistency with other modules like core sizing and winding optimization.

The feature is accessed via a dialog window, where users input parameters, trigger calculations with a "Calculate" button, and view results in a text output area and a graphical plot.

### Leakage Inductance Model

- **Leakage Inductance Calculation**:
  Leakage inductance is proportional to the number of turns squared, core permeability, and window area, inversely proportional to bobbin height, and adjusted by a leakage factor and interleaving effects.
  `L_leak = leakage_factor * interleaving_factor * mu_0 * N^2 * window_area / bobbin_height`
  where:
  - `L_leak` is the leakage inductance (H).
  - `leakage_factor` is a configuration-specific factor (0.05 for Concentric, 0.02 for Interleaved).
  - `interleaving_factor` adjusts for winding interleaving, defined as:
    `interleaving_factor = 1.0` for Concentric, or `(primary_layers + sum(secondary_layers)) / max(primary_layers, max(secondary_layers))` for Interleaved.
  - `mu_0 = 4 * pi * 10^-7` H/m is the permeability of free space.
  - `N` is the number of turns (primary or secondary).
  - `window_area = bobbin_width * bobbin_height` (m²).
  - `bobbin_height` is the bobbin height (m).
- **Primary and Secondary Leakage**:
  - Primary leakage inductance uses primary turns (`Np`).
  - Secondary leakage inductance is calculated for each secondary winding using its turns (`Ns`).
  - The model assumes uniform flux distribution within the window area, modified by interleaving to reduce leakage in Interleaved configurations.
- **Window Area Adjustment**:
  The effective window area accounts for insulation:
  `winding_area = window_area - insulation_area`
  where `insulation_area = insulation_thickness * (primary_layers + sum(secondary_layers))`.
  If `winding_area` is negative, a default value of 0.001 m² is used to avoid singularities.

### Stray Magnetic Field Model
Stray magnetic fields result from magnetic flux escaping the core, posing risks for electromagnetic interference. The model approximates the field as originating from the primary winding current and turns, decaying with distance.

- **Stray Field Strength**:
  The magnetic field strength (`H`) at a distance `d` from the core is modeled using Ampere’s law, assuming a cylindrical field distribution around the core:
  `H_stray = stray_field_factor * (Ip * Np) / (2 * pi * (d + core_dimension))`
  where:
  - `H_stray` is the magnetic field strength (A/m).
  - `stray_field_factor` is configuration-specific (0.1 for Concentric, 0.5 for Interleaved, reflecting increased stray flux in Interleaved designs).
  - `Ip = total_power / (Vp * efficiency)` is the primary current (A).
  - `Np` is the number of primary turns.
  - `d` is the distance from the core (m), evaluated from 0.001 to 1.000 m in 0.001 m increments.
  - `core_dimension = sqrt(core_area)` approximates the core’s effective radius (m).
- **Safe Distance**:
  The safe distance is the smallest distance where `H_stray < 0.1` A/m, determined by iterating through calculated field strengths:
  `safe_distance = min(d for d, H in zip(distances, stray_fields) if H < 0.1, default=1.0)`

### Supporting Calculations
- **Core Area**:
  The core area is computed using the `CoreSizing` class:
  `core_area = k * (power / frequency)^0.5`
  where `k` is an empirical constant based on transformer and core type.
- **Primary Turns**:
  Determined using the transformer equation:
  `Np = (Vp * 10^8) / (4.44 * freq * core_area * B_max)`
  where `B_max` is the maximum flux density, limited by transformer type and core material (e.g., 1.5 T for Power, 0.8 T for High-Frequency).
- **Fill Factor**:
  The fill factor assesses winding efficiency within the bobbin:
  `fill_factor = (Np * wire_diameter^2 * primary_layers + sum(Ns * wire_diameter^2 * secondary_layers)) / winding_area`
  where `wire_diameter = 0.001` m (assumed for simplicity).
- **Material and Configuration Parameters**:
  - Winding configuration affects `leakage_factor` and `stray_field_factor`.
  - Core material (Silicon Steel, Ferrite) influences permeability (`mu_r`) and saturation flux density, used in turns calculations.
  - Transformer type adjusts `B_max` and coupling coefficients.

### Algorithms Utilized

The feature employs straightforward numerical algorithms to compute leakage inductance and stray fields, integrated within the `Transformer` class’s `calculate` method and visualized in `LeakageAndStrayFieldEst.py`. 

- **Leakage Inductance Calculation**:
  1. Retrieve inputs: transformer type, core type, winding configuration, voltages, frequency, power, efficiency, bobbin dimensions, insulation thickness, and layer counts.
  2. Validate inputs: Ensure voltages, frequency, power, efficiency, and dimensions are positive and within realistic ranges.
  3. Compute core area using the `CoreSizing` class.
  4. Calculate primary turns (`Np`) using the transformer equation.
  5. Determine secondary turns (`Ns`) via turns ratios from `RatioAndWireGaugeCalc`.
  6. Compute window area and insulation area to find effective `winding_area`.
  7. Calculate interleaving factor based on winding configuration and layer counts.
  8. Compute primary leakage inductance using `L_leak` equation with `Np`.
  9. Compute secondary leakage inductance for each secondary winding using `Ns`.
  10. Convert inductances to microhenries (µH) for output.
- **Stray Field Calculation**:
  1. Define distance array: `distances = [0.001, 0.002, ..., 1.000]` m (1000 points).
  2. Compute primary current (`Ip`) and core dimension (`sqrt(core_area)`).
  3. For each distance `d`, calculate stray field strength (`H_stray`) using the provided equation.
  4. Store results in a list (`stray_fields`).
  5. Find safe distance by iterating through `stray_fields` to identify the first distance where `H_stray < 0.1` A/m.
- **Result Visualization**:
  1. Update text output with leakage inductances, safe distance, and transformer parameters.
  2. Generate a plot of `stray_fields` versus `distances` using Matplotlib, with labeled axes and a grid.
- **Integration with Transformer Class**:
  - The `Transformer.calculate` method returns a dictionary containing `leakage_inductance_primary`, `leakage_inductance_secondary`, `stray_fields`, `distances`, and `safe_distance`.
  - The `LeakageAndStrayFieldEst` dialog extracts these results for display and plotting.

---

## Winding Configuration Optimization 

The Winding Configuration Optimization feature optimizes transformer winding parameters to minimize leakage inductance, enhancing performance. Implemented in `WindingConfigOptimizer.py`, this module uses the `Transformer` class to perform calculations and displays results via text output and plots. 

### Functioning

The Winding Configuration Optimization feature allows users to design a transformer with a single secondary winding, optimizing its winding configuration to reduce leakage inductance.

- **Input Parameters**: Users specify transformer type (e.g., Power, Distribution), core configuration (Core-Type, Shell-Type), winding configuration (Concentric, Interleaved), primary and secondary voltages, frequency, power, efficiency, bobbin dimensions (width, height), insulation thickness, and primary and secondary layer counts.
- **Calculations**: The module computes core area, turns ratios, wire gauges, winding resistances, core and copper losses, mutual inductance, leakage inductance, fill factor, and bobbin utilization. It optimizes the interleaving factor to minimize leakage inductance.
- **Output**: Displays transformer type, core type, winding configuration, core area, primary and secondary leakage inductances, fill factor, bobbin utilization, interleaving pattern, and layer counts. A plot shows leakage inductance versus interleaving factor.
- **Visualization**: A Matplotlib plot illustrates how leakage inductance varies with interleaving factor, aiding users in understanding optimization results.

The feature integrates with the `Transformer` class, ensuring consistent calculations across the software.

### Algorithms

The optimization process employs algorithms to compute transformer parameters and minimize leakage inductance. 

- **Input Validation Algorithm**:
  - Ensures all inputs (voltages, frequency, power, efficiency, bobbin dimensions, insulation thickness, layer counts) are positive, with efficiency between 0 and 1.
  - Validates numeric inputs, raising errors for invalid entries.
- **Core Sizing Algorithm**:
  - Uses the `CoreSizing` class to calculate core area based on power, frequency, transformer type, and core type.
  - Employs an empirical formula: `core_area = k * (power / frequency)^0.5`, where `k` varies by transformer and core type.
- **Turns and Wire Gauge Calculation Algorithm**:
  - Uses the `RatioAndWireGaugeCalc` class to compute the turns ratio with a regulation factor (1.05 for Power/Distribution, 1.02 otherwise): `turns_ratio = (Vp / Vs) * regulation_factor`.
  - Calculates primary current: `Ip = power / (Vp * efficiency)` and secondary current: `Is = power / Vs`.
  - Selects wire gauges from an AWG table based on a current density of 2 A/mm², choosing the smallest AWG meeting the required area: `wire_area = current / current_density`.
- **Winding Parameter Calculation Algorithm**:
  - Computes primary turns: `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`, where `B_max` is the maximum flux density.
  - Determines secondary turns: `Ns = Np / turns_ratio`.
  - Uses user-specified primary and secondary layer counts or defaults from the `InterleavingAndLayering` class.
- **Leakage Inductance Optimization Algorithm**:
  - Iterates over interleaving factors from 0.5 to 2.0 in steps of 0.1.
  - For each factor, recalculates leakage inductance using the `Transformer` class.
  - Stores primary leakage inductance values for plotting.
  - Outputs the configuration with the recommended interleaving pattern.
- **Fill Factor and Bobbin Utilization Algorithm**:
  - Calculates window area: `window_area = bobbin_width * bobbin_height`.
  - Computes insulation area: `insulation_area = insulation_thickness * (primary_layers + secondary_layers)`.
  - Determines winding area: `winding_area = window_area - insulation_area`.
  - Calculates fill factor: `fill_factor = (Np * wire_diameter^2 * primary_layers + Ns * wire_diameter^2 * secondary_layers) / winding_area`.
- **Plotting Algorithm**:
  - Uses Matplotlib to plot leakage inductance (in µH) versus interleaving factor.
  - Draws a line plot with red scatter points for clarity.
  - Sets axis limits to 90%–110% of interleaving factor range and 0–110% of maximum leakage inductance.

### Electrical Models
- **Turns Ratio**:
  - Accounts for voltage regulation:
    `turns_ratio = (Vp / Vs) * regulation_factor`
    where `regulation_factor` is 1.05 for Power/Distribution transformers, 1.02 otherwise.
- **Current Calculation**:
  - Primary current:
    `Ip = power / (Vp * efficiency)`
  - Secondary current:
    `Is = power / Vs`
- **Winding Resistance**:
  - Temperature-corrected resistance:
    `R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))`
    where `resistivity_cu = 1.68e-8` Ω·m, `alpha_cu = 0.00393/°C`, `T_ref = 20°C`, `length` is winding length, `area` is wire cross-sectional area.
- **Copper Losses**:
  - Power loss:
    `P_copper = I^2 * R`
    Summed for primary and secondary windings.

### Magnetic Models
- **Core Area**:
  - Empirical sizing:
    `core_area = k * (power / frequency)^0.5`
    where `k` is a type-specific constant.
- **Primary Turns**:
  - Based on transformer equation:
    `Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)`
    where `B_max` is the maximum flux density, adjusted by core type (1.0 for Core-Type, 0.9 for Shell-Type) and material saturation limits.
- **Secondary Turns**:
  - Derived from turns ratio:
    `Ns = Np / turns_ratio`
- **Mutual Inductance**:
  - Primary inductance:
    `Lp = Np^2 * mu * core_area / l`
  - Secondary inductance:
    `Ls = Ns^2 * mu * core_area / l`
  - Mutual inductance:
    `M = coupling_coeff * sqrt(Lp * Ls)`
    where `mu = mu_r * mu_0`, `mu_0 = 4 * pi * 1e-7` H/m, `l = 0.1` m, `coupling_coeff` is type-specific (e.g., 0.98 for Power).
- **Leakage Inductance**:
  - Primary leakage inductance:
    `L_leak_primary = leakage_factor * interleaving_factor * mu_0 * Np^2 * window_area / bobbin_height`
  - Secondary leakage inductance:
    `L_leak_secondary = leakage_factor * interleaving_factor * mu_0 * Ns^2 * window_area / bobbin_height`
    where `leakage_factor` is 0.05 (Concentric) or 0.02 (Interleaved), `interleaving_factor = 1.0` (Concentric) or `(primary_layers + secondary_layers) / max(primary_layers, secondary_layers)` (Interleaved).
- **Fill Factor**:
  - Ratio of conductor area to available window area:
    `fill_factor = (Np * wire_diameter^2 * primary_layers + Ns * wire_diameter^2 * secondary_layers) / winding_area`
    where `wire_diameter` is derived from wire gauge.

### Core Loss Models
- **Steinmetz Equation**:
  - Core loss:
    `P_core = k * freq^alpha * B_max^beta * core_volume`
    where `k`, `alpha`, `beta` are material-specific (e.g., Silicon Steel: `k=0.1`, `alpha=1.6`, `beta=2.0`), `core_volume = core_area * 0.1` m.
- **Hysteresis Loss**:
  - `P_hyst = k * freq^alpha * B_max^beta * core_volume`
- **Eddy Current Loss**:
  - `P_eddy = 0.5 * k * freq^2 * B_max^2 * core_volume`
  - Total core loss: `P_core = P_hyst + P_eddy`.

---

## Core Loss and Copper Loss Estimation

The Core Loss and Copper Loss Estimation implemented in `CoreLossCopperLossEst.py` and supported by `transformer.py`, calculates the core losses (hysteresis and eddy current) and copper losses in a transformer based on user inputs.

### Functionality

The Core Loss and Copper Loss Estimation feature enables users to analyze the power losses in a transformer, which are critical for determining efficiency and thermal performance.

- Accepts inputs such as transformer type (e.g., Power, Distribution), core type (Core-Type, Shell-Type), core material (Silicon Steel, Ferrite), winding configuration (Concentric, Interleaved), primary and secondary voltages, frequency, power, efficiency, ambient temperature, and operating temperature.
- Computes core losses, comprising hysteresis and eddy current losses, based on the transformer's magnetic properties and operating conditions.
- Calculates copper losses in the primary and secondary windings, accounting for temperature-corrected resistances.
- Outputs numerical results, including individual loss components and total losses, displayed in a text area.
- Visualizes the loss breakdown using a bar chart, showing hysteresis, eddy current, and copper losses for comparison.

The feature interfaces with the `Transformer` class in `transformer.py` to perform calculations, ensuring consistency with other modules of the software.

### Core Loss Model
Core losses arise from the magnetic behavior of the transformer core under alternating flux. They are divided into hysteresis and eddy current losses, calculated using the Steinmetz equation, which relates loss to frequency and flux density.

- **Hysteresis Loss**: Occurs due to the energy required to realign magnetic domains in the core material.
  ```
  P_hyst = k * freq^alpha * B_max^beta * core_volume
  ```
  where:
  - `P_hyst` is the hysteresis loss (W),
  - `k` is the Steinmetz coefficient (material-specific, e.g., 0.1 for Silicon Steel),
  - `freq` is the operating frequency (Hz),
  - `alpha` is the frequency exponent (e.g., 1.6 for Silicon Steel),
  - `B_max` is the maximum flux density (T),
  - `beta` is the flux density exponent (e.g., 2.0 for Silicon Steel),
  - `core_volume` is the core volume (m³), approximated as `core_area * 0.1`.

- **Eddy Current Loss**: Results from induced currents in the core material due to changing magnetic flux.
  ```
  P_eddy = 0.5 * k * freq^2 * B_max^2 * core_volume
  ```
  where parameters are as above, with the frequency squared term reflecting the quadratic dependence of eddy currents.

- **Total Core Loss**: Sum of hysteresis and eddy current losses.
  ```
  P_core = P_hyst + P_eddy
  ```

### Copper Loss Model
Copper losses occur due to resistive heating in the transformer windings, influenced by current and temperature-dependent resistance.

- **Winding Resistance (Temperature-Corrected)**: Resistance varies with temperature due to the temperature coefficient of copper.
  ```
  R = resistivity_cu * length / area * (1 + alpha_cu * (T - T_ref))
  ```
  where:
  - `R` is the winding resistance (Ω),
  - `resistivity_cu = 1.68e-8` Ω·m is the resistivity of copper,
  - `length` is the winding length (m), approximated as `turns * 0.1`,
  - `area` is the wire cross-sectional area (m²), derived from wire diameter (default 0.001 m),
  - `alpha_cu = 0.00393` /°C is the temperature coefficient of copper,
  - `T` is the temperature (°C, ambient or operating),
  - `T_ref = 20` °C is the reference temperature.

- **Current Calculation**:
  - Primary current:
    ```
    Ip = total_power / (Vp * efficiency)
    ```
  - Secondary current (for each secondary winding):
    ```
    Is = power / (Vs * efficiency)
    ```
  where `Vp` and `Vs` are primary and secondary voltages (V), `total_power` is the sum of secondary powers (VA), and `power` is the individual secondary power (VA).

- **Copper Losses**:
  - Primary copper loss:
    ```
    P_copper_primary = Ip^2 * R_primary
    ```
  - Secondary copper loss (summed over all secondary windings):
    ```
    P_copper_secondary = sum(Is^2 * R_secondary)
    ```
  - Total copper loss:
    ```
    P_copper_total = P_copper_primary + P_copper_secondary
    ```
  where `R_primary` and `R_secondary` are calculated at the operating temperature (default 80°C) unless specified otherwise.

### Magnetic Flux Density
The maximum flux density (`B_max`) is determined to ensure the core operates below saturation, adjusted by core type and material properties.
```
B_max = min(type_params[transformer_type]["max_flux_density"] * flux_density_factor, material_params[core_material]["saturation_flux_density"])
```
where:
- `flux_density_factor = 1.0` for Core-Type or `0.9` for Shell-Type,
- `type_params` provides transformer-specific flux density limits (e.g., 1.5 T for Power),
- `material_params` provides material-specific saturation flux density (e.g., 1.8 T for Silicon Steel).

### Core Area
The core area is calculated using an empirical formula in `CoreSizing.py`:
```
core_area = k * (power / frequency)^0.5
```
where `k` is a coefficient dependent on transformer and core type, and `power` is the total power (VA).

### Core Loss Calculation 
1. **Input Processing**:
   - Retrieve transformer type, core type, core material, frequency, and total power.
   - Select material-specific Steinmetz parameters (`k`, `alpha`, `beta`) and saturation flux density from `material_params`.
   - Select transformer-specific maximum flux density from `type_params`.
2. **Core Area Calculation**:
   - Use the `CoreSizing` class to compute `core_area` based on `power` and `frequency`.
3. **Flux Density Determination**:
   - Calculate `B_max` by taking the minimum of the transformer type’s flux density limit (adjusted by `flux_density_factor`) and the material’s saturation flux density.
4. **Core Volume Estimation**:
   - Approximate `core_volume` as `core_area * 0.1` (assuming a core depth of 0.1 m).
5. **Hysteresis Loss**:
   - Compute `P_hyst` using the Steinmetz equation with `freq`, `B_max`, and `core_volume`.
6. **Eddy Current Loss**:
   - Compute `P_eddy` using the modified Steinmetz equation with `freq^2` and `B_max^2`.
7. **Total Core Loss**:
   - Sum `P_hyst` and `P_eddy` to obtain `P_core`.

### Copper Loss Calculation 
1. **Input Processing**:
   - Retrieve primary voltage (`Vp`), secondary voltages (`Vs_list`), total power, individual secondary powers (`power_list`), efficiency, ambient temperature, and operating temperature.
2. **Current Calculation**:
   - Compute primary current (`Ip`) using `total_power / (Vp * efficiency)`.
   - Compute secondary currents (`Is_list`) for each secondary winding using `power / (Vs * efficiency)`.
3. **Turns Calculation**:
   - Calculate primary turns (`Np`) using:
     ```
     Np = (Vp * 1e8) / (4.44 * freq * core_area * B_max)
     ```
   - Calculate secondary turns (`Ns`) for each secondary using `Np / turns_ratio`, where `turns_ratio` is computed by `RatioAndWireGaugeCalc`.
4. **Wire Parameters**:
   - Use a default wire diameter (0.001 m) to compute wire area (`pi * (diameter / 2)^2`).
   - Estimate winding lengths as `turns * 0.1` (assuming 0.1 m per turn).
5. **Resistance Calculation**:
   - Compute primary and secondary resistances at operating temperature using the temperature-corrected resistance formula.
   - Use `resistivity_cu`, `alpha_cu`, and `T_ref` for accuracy.
6. **Copper Loss Calculation**:
   - Compute primary copper loss as `Ip^2 * R_primary`.
   - Compute secondary copper loss as the sum of `Is^2 * R_secondary` for each secondary winding.
   - Sum primary and secondary losses to obtain total copper loss.

### Plotting Algorithm
1. **Data Collection**:
   - Gather `hysteresis_loss`, `eddy_current_loss`, and `total_copper_loss` from the `Transformer` calculation results.
2. **Bar Chart Generation**:
   - Use Matplotlib to create a bar chart with three categories: Hysteresis Loss, Eddy Current Loss, and Copper Loss.
   - Assign distinct colors to each bar (e.g., red, green (‘g’), and blue).
   - Label each bar with loss values (W) above the bar for clarity.
   - Set the y-axis to represent power loss (W) and include a grid lines for reference.
3. **Plot Configuration**:
   - Set the plot title to "Loss Breakdown".
   - Configure the y-axis label as "Power Loss (W).
   - Ensure the plot canvas is cleared and redrawn to reflect the new data.

### Integration with Transformer Class

The feature relies on the `Transformer` class’s `calculate` method, which orchestrates the loss calculations and returns a dictionary containing:
- `core_loss`: Total core loss (W),
- `hysteresis_loss`: Hysteresis loss (W),
- `eddy_loss_current_loss`: Eddy loss (W),
- `total_loss_copper_loss`: Total copper loss (W),
- `copper_loss_primary_loss`: Primary copper loss (W),
- `copper_loss_copper_secondary_loss`: Secondary copper loss (W).

The `Transformer` class integrates helper classes:
- `CoreSizing` for core area calculations.
- `RatioAndWireGaugeCalc` for current and turns ratio calculations.

---

## Thermal Simulation (FEM) Feature in Transformer Design Tool

The Thermal Simulation enables transient thermal analysis of transformers, modeling temperature distribution, heat transfer mechanisms, and hotspot locations. Implemented in `ThermalSimulationFEM.py`, this module simulates heat generation from core and copper losses, providing insights into thermal performance under specified operating conditions. 

### Functioning

The Thermal Simulation (FEM) feature models the transient temperature distribution across a 2D transformer cross-section, accounting for heat generation from core and copper losses. Users input transformer parameters (type, core configuration, material, voltages, frequency, power, efficiency), thermal properties (core and winding thermal conductivities, emissivity, convection coefficient, airflow velocity), and simulation settings (ambient temperature, operating temperature, simulation time). 

- Calculates core and copper losses using the `Transformer` class.
- Performs a transient FEM simulation to compute temperature distribution over time.
- Identifies hotspot locations and calculates heat losses via conduction, convection, and radiation.
- Visualizes results through an animated contour plot of temperature distribution, with an optional heat flux vector display.

Results include maximum and average temperatures, hotspot coordinates, and heat loss breakdowns, displayed in a text output and animated plot.

### Simulation

1. **Input Collection**: Gathers user inputs for transformer parameters (e.g., primary voltage, frequency, power), thermal properties (e.g., core thermal conductivity, convection coefficient), and simulation duration.
2. **Input Validation**: Ensures inputs are valid (e.g., positive voltages, power, and thermal conductivities; emissivity between 0 and 1; simulation time positive).
3. **Loss Calculation**: Uses the `Transformer` class to compute core and copper losses based on input parameters, providing heat source terms for the simulation.
4. **FEM Setup**: Defines a 2D domain (core and windings) with a 20x10 triangular mesh, assigning material properties (thermal conductivity, density, specific heat) to core and winding regions.
5. **Heat Transfer Simulation**: Solves the transient heat equation using FEM, incorporating conduction, convection, and radiation boundary conditions.
6. **Result Processing**: Extracts maximum and average temperatures, hotspot positions, and heat flux for each time step, storing results for approximately 50 frames.
7. **Visualization**: Animates temperature distribution over time, optionally displaying heat flux vectors, and updates text output with final results.

The simulation runs for the specified duration, saving results at intervals to balance computational efficiency and visualization detail.

### Heat Transfer Equation
The transient heat equation governs temperature distribution:
```
rho * cp * dT/dt = div(k * grad(T)) + q
```
- `rho`: Material density (kg/m³; 7650 for silicon steel core, 8960 for copper windings).
- `cp`: Specific heat capacity (J/kg·K; 450 for core, 385 for windings).
- `T`: Temperature (K).
- `k`: Thermal conductivity (W/m·K; user-specified, e.g., 30 for core, 400 for windings).
- `q`: Heat source (W/m³), from core and copper losses.

### Heat Sources
Heat generation is derived from transformer losses:
```
q_core = core_loss / core_volume
q_winding = copper_loss / winding_volume
```
- `core_loss`: Total core loss (W), from hysteresis and eddy currents.
- `copper_loss`: Total copper loss (W), from primary and secondary windings.
- `core_volume`: Core volume (m³), calculated as `core_area * 0.1` (0.1 m depth).
- `winding_volume`: Winding volume (m³), calculated as `(domain_width * domain_height - core_area) * 0.1`.

### Boundary Conditions
Heat transfer at boundaries includes convection and radiation:
```
-k * grad(T) · n = h_total * (T - T_amb) + epsilon * sigma * (T^4 - T_amb^4)
```
- `n`: Normal vector to the boundary.
- `h_total`: Total heat transfer coefficient (W/m²·K), sum of convection (`h_conv`) and radiation (`h_rad`).
- `h_conv`: Convection coefficient (W/m²·K; user-specified, e.g., 10).
- `h_rad`: Radiation coefficient, approximated as:
  ```
  h_rad = 4 * epsilon * sigma * T_ref^3
  ```
  where `epsilon` is emissivity (0-1), `sigma` = 5.67e-8 W/m²·K^4, `T_ref` = (T_amb + T_operating)/2 (K).
- `T_amb`: Ambient temperature (K).

### Heat Flux
Heat flux is computed from temperature gradients:
```
q = -k * grad(T)
```
- `q`: Heat flux vector (W/m²), with components in x and y directions.
- `grad(T)`: Temperature gradient, calculated element-wise in FEM.

### Heat Losses
Total heat loss comprises conduction, convection, and radiation:
- **Conduction Loss**:
  ```
  conduction_loss = sum(heat_source) * (domain_width * domain_height * 0.1)
  ```
- **Convection Loss**:
  ```
  convection_loss = sum(h_conv * length_per_node * (T[n] - T_amb)) for boundary nodes n
  ```
- **Radiation Loss**:
  ```
  radiation_loss = sum(epsilon * sigma * length_per_node * (T[n]^4 - T_amb^4)) for boundary nodes n
  ```
- `length_per_node`: Boundary length per node, calculated as `perimeter / num_boundary_nodes`.

### FEM Mesh Generation
- Creates a 2D domain with core (square, side = sqrt(core_area)) and windings (10 mm thick annulus).
- Discretizes into a 20x10 grid, forming triangular elements (two triangles per grid cell).
- Nodes are defined at grid points, with coordinates `(x, y)`.

### Material Property Assignment
- Assigns thermal conductivity (`k`), density (`rho`), and specific heat (`cp`) to nodes.
- Core region: Nodes within `[winding_thickness, core_side + winding_thickness]` (x) and `[0, core_side]` (y) use core properties (e.g., k=30 W/m·K, rho=7650 kg/m³, cp=450 J/kg·K).
- Winding region: Remaining nodes use copper properties (e.g., k=400 W/m·K, rho=8960 kg/m³, cp=385 J/kg·K).

### FEM System Assembly
- Formulates the system: `[C] * dT/dt + [K] * T = F`.
- **Stiffness Matrix (K)**: For each triangular element:
  - Computes area: `area = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))`.
  - Calculates shape function derivatives: `b = [y1-y2, y2-y0, y0-y1]/(2*area)`, `c = [x2-x1, x0-x2, x1-x0]/(2*area)`.
  - Element stiffness: `Ke = k_elem * area * (outer(b, b) + outer(c, c))`.
  - Assembles into global `K` matrix.
- **Capacitance Matrix (C)**: Element capacitance: `Ce = rho_elem * cp_elem * area / 3 * eye(3)`.
- **Force Vector (F)**: Element heat source: `Fe = q_elem * area / 3 * ones(3)`, where `q_elem` is the average heat source over element nodes.

### Boundary Conditions
- Identifies boundary nodes (x=0, x=domain_width, y=0, y=domain_height).
- Applies convection and radiation: `K[n, n] += h_total * length_per_node`, `F[n] += h_total * T_amb * length_per_node`.

### Transient Simulation
- Uses a fixed time step: `dt = 10 s`.
- Computes time steps: `num_steps = t_max / dt`.
- Saves results every `num_steps // 50` steps for ~50 frames.
- Solves: `[C + dt * K] * T_new = C * T_old + dt * F` using sparse solver (`spla.spsolve`).
- Initial condition: `T = T_amb` (Kelvin) at all nodes.

### Heat Flux Calculation
- Computes temperature gradients per element: `grad_x_elem = dot(b, Te)`, `grad_y_elem = dot(c, Te)`.
- Averages gradients over nodes, calculates heat flux: `q = -k * [grad_x, grad_y]`.

### Animation
- Uses Matplotlib’s `FuncAnimation` to animate temperature contours over time.
- Updates contour plot, hotspot marker, and time text for each saved frame.
- Optionally plots heat flux vectors using `quiver`.


---

## Vibration and Acoustic Noise Modeling

The Vibration and Acoustic Noise Modeling feature of the Transformer Design Tool analyzes the mechanical vibrations and acoustic emissions of a transformer under operational conditions. This module, implemented in `VibrationAndAcousticNoiseModel.py`, uses finite element method (FEM) simulations to model the effects of magnetic forces, including magnetostriction and Lorentz forces, on the transformer's core and windings. 

### Functionality

The vibration and acoustic noise modeling feature performs the following tasks:

- **Input Processing**: Accepts user inputs for transformer parameters (type, core configuration, material, voltages, frequency, power, efficiency), mechanical properties (Young’s modulus for core and windings, magnetostriction coefficient, damping ratio), and environmental conditions (ambient and operating temperatures).
- **Vibration Simulation**: Computes mechanical displacements and velocities in the transformer’s core and windings due to magnetic forces using a 2D FEM model.
- **Acoustic Analysis**: Estimates the sound pressure level at a reference distance (1 meter) based on surface vibration velocities.
- **Result Visualization**: Displays a contour plot of vibration displacement magnitude (in micrometers) across the transformer domain, alongside numerical results for maximum displacement, velocity, and SPL.
- **Output Reporting**: Provides detailed results, including transformer type, core material, maximum vibration displacement, maximum vibration velocity, SPL, and dominant frequencies (twice and four times the input frequency).

The module integrates with the `Transformer` class to obtain electrical and magnetic parameters, ensuring consistency with other analyses in the tool.

### Simulation

1. **Input Validation**: Checks that all inputs are valid (e.g., positive voltages, efficiency between 0 and 1, positive mechanical properties, damping ratio between 0 and 1, temperatures between -50°C and 200°C).
2. **Transformer Parameter Calculation**: Uses the `Transformer` class to compute core area, maximum flux density, primary current, and leakage flux density based on user inputs.
3. **FEM Domain Setup**:
   - Defines a 2D domain representing the transformer core (square) and surrounding windings.
   - Discretizes the domain into a 20x10 grid of nodes, forming triangular elements.
   - Assigns material properties (Young’s modulus, density, Poisson’s ratio) to core and winding regions.
4. **Force Calculation**:
   - Computes magnetostriction-induced stresses in the core based on flux density and magnetostriction coefficient.
   - Calculates electromagnetic forces in the core and Lorentz forces in the windings due to leakage flux and current.
5. **Vibration Simulation**:
   - Assembles mass, damping, and stiffness matrices for the FEM model.
   - Applies fixed boundary conditions at the bottom edge (y=0).
   - Solves the harmonic vibration equation at twice the operating frequency to obtain displacement and velocity fields.
6. **Acoustic Noise Calculation**:
   - Computes acoustic power from the RMS surface velocity.
   - Estimates SPL at 1 meter using spherical wave propagation.
7. **Visualization and Output**:
   - Generates a contour plot of displacement magnitude (in micrometers).
   - Outputs numerical results, including maximum displacement, velocity, SPL, and dominant frequencies.

### Mechanical Vibration Model

The vibration behavior is modeled using the dynamic equation for a mechanical system:

`[M] * d^2 u / dt^2 + [C] * du / dt + [K] * u = F`
- **[M]**: Mass matrix, representing the mass distribution of the core and windings.
- **[C]**: Damping matrix, accounting for energy dissipation.
- **[K]**: Stiffness matrix, representing the elastic properties of the materials.
- **u**: Displacement vector (x and y components for each node).
- **F**: Force vector, including magnetostriction and Lorentz forces.

The system is solved in the frequency domain at the harmonic frequency `omega = 2 * pi * 2 * freq`, where `freq` is the input electrical frequency (e.g., 50 Hz).

### Magnetostriction
Magnetostriction in the core causes strain due to magnetic flux density:

`strain_m = lambda_s * (B / B_s))^2`

- `lambda_s`: Magnetostriction coefficient (e.g., 4e-6).
- `B`: Maximum flux density (from transformer calculations).
- `B_s`: Saturation flux density (2.0 T).

The resulting stress is converted to force:

`stress_m = E * strain_m`

- `E`: Young’s modulus of the core (e.g., 200e9 Pa).

This stress is applied as a volumetric force in the core region:

`F_m = stress_m / depth`

- `depth`: Assumed depth of the 2D model (0.1 m).

### Electromagnetic Force
An electromagnetic force acts on the core due to magnetic flux:

`F_em = (B^2 * core_area) / (2 * mu_0 * depth)`

- `core_area`: Core cross-sectional area (from transformer calculations).
- `mu_0`: Permeability of free space (4 * pi * 1e-7 H/m).

### Lorentz Force
In the windings, Lorentz forces arise from the interaction of primary current and leakage flux:

`F_lorentz = I_primary * B_leakage / depth`

- `I_primary`: Primary current (power / (Vp * efficiency)).
- `B_leakage`: Leakage flux density (default 0.01 T).

### Material Properties
- **Young’s Modulus**:
  - Core: `E_core` (e.g., 200e9 Pa for silicon steel).
  - Windings: `E_winding` (e.g., 110e9 Pa for copper).
- **Density**:
  - Core: 7650 kg/m³ (silicon steel).
  - Windings: 8960 kg/m³ (copper).
- **Poisson’s Ratio**:
  - Core: 0.3.
  - Windings: 0.34.

### Damping
Rayleigh damping is used:

`[C] = 2 * zeta * omega * [M]`

- `zeta`: Damping ratio (e.g., 0.02).
- `omega`: Angular frequency (2 * pi * 2 * freq).

### Acoustic Noise Model

Acoustic noise is estimated from surface vibrations:

1. **Acoustic Power**:
   `P_acoustic = rho_air * c_air * S * v_rms^2`
   - `rho_air`: Air density (1.225 kg/m³).
   - `c_air`: Speed of sound in air (343 m/s).
   - `S`: Surface area of the transformer (2 * (domain_width + domain_height) * depth).
   - `v_rms`: Root-mean-square velocity, computed as `sqrt(mean(v_x^2 + v_y^2))`, where `v_x` and `v_y` are nodal velocities.

2. **Sound Pressure**:
   `p_rms = sqrt((P_acoustic * rho_air * c_air) / (4 * pi * r^2))`
   - `r`: Reference distance (1 m).

3. **Sound Pressure Level**:
   `SPL = 20 * log10(p_rms / p_0)`
   - `p_0`: Reference pressure (20e-6 Pa).

The dominant frequencies of the noise are twice and four times the input frequency (e.g., 100 Hz and 200 Hz for 50 Hz input), due to the magnetic forces oscillating at twice the electrical frequency.

### Finite Element Method (FEM) 
- **Mesh Generation**:
  - Creates a 20x10 grid of nodes over the 2D domain (core + windings).
  - Forms triangular elements by connecting nodes, resulting in two triangles per grid cell.
- **Matrix Assembly**:
  - **Mass Matrix ([M])**: Computed element-wise using material density and element area, assuming a uniform distribution:
    `Me = (rho_elem * area * depth / 3) * I_6`
    where `rho_elem` is the average density, `area` is the triangle area, `depth` is 0.1 m, and `I_6` is a 6x6 identity matrix.
  - **Stiffness Matrix ([K])**: Assembled using the material’s Young’s modulus and Poisson’s ratio:
    `Ke = area * depth * (B_elem^T * D * B_elem)`
    where `D` is the constitutive matrix, and `B_elem` is the strain-displacement matrix derived from nodal coordinates.
  - **Damping Matrix ([C])**: Constructed using Rayleigh damping proportional to the mass matrix.
  - **Force Vector (F)**: Aggregates magnetostriction, electromagnetic, and Lorentz forces, distributed evenly across element nodes.
- **Boundary Conditions**:
  - Applies fixed constraints at nodes where `y = 0` by adding a large penalty (1e12) to the diagonal of `[K]` and setting corresponding forces to zero.
- **Harmonic Solution**:
  - Solves the complex system:
    `A * u = F_vec`
    where `A = -omega^2 * [M] + i * omega * [C] + [K]`.
  - Uses a sparse least-squares solver (`scipy.sparse.linalg.lsqr`) to compute real and imaginary displacements (`u_x`, `u_y`) and velocities (`v_x`, `v_y`).
  - Checks for non-finite solutions, defaulting to zero displacements if numerical issues arise.

### Acoustic Noise Calculation
- **Velocity RMS**:
  - Computes the RMS velocity from nodal velocities:
    `v_rms = sqrt(mean(v_x^2 + v_y^2))`
- **Acoustic Power and SPL**:
  - Calculates acoustic power using the surface area and RMS velocity.
  - Converts to sound pressure and SPL using spherical wave propagation formulas.
- **Frequency Identification**:
  - Reports dominant frequencies as `2 * freq` and `4 * freq`, reflecting the harmonic nature of magnetic forces.

### Visualization Algorithm
- **Contour Plot**:
  - Computes displacement magnitude:
    `displacement_magnitude = sqrt(u_x^2 + u_y^2)`
  - Converts to micrometers for visualization.
  - Uses Matplotlib’s `contourf` to plot a color-coded map with 50 levels, ensuring non-zero and non-NaN values.
  - Adds a colorbar with a label indicating displacement in micrometers.

### Debugging and Logging
- **Matrix Conditioning**: Estimates the condition number of the system matrix for small systems to monitor numerical stability.
- **Force and Displacement Logging**: Prints computed forces, displacement ranges, and solver status for debugging.
- **Fallback Visualization**: Uses a small non-zero displacement (1e-9 m) if results are invalid to prevent plotting failures.

---

## Harmonic and Frequency Response Analysis (HFRA)

The Harmonic and Frequency Response Analysis (HFRA) implemented in `HFRA.py`, simulates the mechanical and acoustic behavior of a transformer under harmonic excitation. It calculates displacement, resonant frequencies, and sound pressure levels, providing insights into vibration and noise characteristics. 

### Functionality

The HFRA module analyzes the transformer's response to harmonic forces induced by alternating magnetic fields. It models the transformer as a 2D structure with a core and windings, subjected to magnetostriction and Lorentz forces. Key functionalities include:

- **Input Parameters**: Accepts transformer type (e.g., Power, Distribution), core configuration (Core-Type, Shell-Type), core material (Silicon Steel, Ferrite), winding configuration (Concentric, Interleaved), primary and secondary voltages, frequency, power, efficiency, core and winding Young’s moduli, magnetostriction coefficient, and damping ratio.
- **Simulation Outputs**: Computes maximum vibration displacement, velocity, sound pressure level (SPL) at 1 meter, and dominant frequencies (2x and 4x the input frequency). Displays results in a text output and visualizes displacement magnitude across the transformer domain.
- **Visualization**: Plots a contour map of vibration displacement magnitude (in micrometers) using Matplotlib, showing spatial distribution of vibrations.

The module interfaces with the `Transformer` class to retrieve core area, flux density, and currents, then performs a finite element method (FEM) simulation to model vibrations and estimate acoustic noise.

### Simulation

1. **Input Collection**: Gathers user inputs from the GUI, including electrical parameters (voltages, frequency, power, efficiency) and mechanical properties (Young’s moduli, magnetostriction coefficient, damping ratio).
2. **Input Validation**: Ensures all inputs are positive and within realistic ranges (e.g., efficiency between 0 and 1, damping ratio between 0 and 1, temperatures between -50°C and 200°C).
3. **Transformer Parameter Calculation**: Calls the `Transformer.calculate` method to obtain core area, maximum flux density (`B`), primary current (`I_primary`), and leakage flux density (`B_leakage`).
4. **Vibration Simulation**: Uses FEM to solve the dynamic equation for mechanical displacement under harmonic forces, considering magnetostriction in the core and Lorentz forces in the windings.
5. **Acoustic Noise Calculation**: Estimates sound pressure level based on surface velocities derived from the vibration simulation.
6. **Result Display**: Outputs numerical results (displacement, velocity, SPL, frequencies) in a text field and visualizes displacement magnitude in a contour plot.

### Mechanical Vibration Model
The transformer is modeled as a 2D elastic structure with a finite element mesh. The governing equation for harmonic vibration is:

`[M] * d^2u/dt^2 + [C] * du/dt + [K] * u = F`

- `[M]`: Mass matrix, based on material densities.
- `[C]`: Damping matrix, using Rayleigh damping.
- `[K]`: Stiffness matrix, based on Young’s modulus and Poisson’s ratio.
- `u`: Displacement vector (x and y components).
- `F`: Force vector, including magnetostriction and Lorentz forces.

For harmonic excitation at angular frequency `omega = 2 * pi * (2 * freq)`, the equation is solved in the frequency domain:

`(-omega^2 * [M] + i * omega * [C] + [K]) * u = F`

where `i` is the imaginary unit.

### Magnetostriction Force
Magnetostriction in the core causes strain due to magnetic field variations:

`epsilon_m = lambda_s * (B / B_s)^2`

- `lambda_s`: Magnetostriction coefficient (e.g., 4e-6).
- `B`: Maximum flux density (from `Transformer` class).
- `B_s`: Saturation flux density (2.0 T).

The resulting stress is:

`sigma_m = E_core * epsilon_m`

- `E_core`: Core Young’s modulus (e.g., 200e9 Pa).

The force per unit volume in the core (x-direction) is:

`F_m = sigma_m / depth`

where `depth = 0.1 m`.

### Electromagnetic Force
The core experiences an electromagnetic force due to the magnetic field:

`F_em = (B^2 * core_area) / (2 * mu_0 * depth)`

- `core_area`: Core cross-sectional area (from `Transformer` class).
- `mu_0 = 4 * pi * 1e-7 H/m`: Permeability of free space.

### Lorentz Force
In the windings, the Lorentz force arises from the interaction of current and leakage flux:

`F_lorentz = I_primary * B_leakage / depth`

- `I_primary = power / (Vp * efficiency)`: Primary current.
- `B_leakage`: Leakage flux density (default 0.01 T).

### Damping
Rayleigh damping is used:

`[C] = 2 * zeta * omega * [M]`

- `zeta`: Damping ratio (e.g., 0.02).
- `omega = 2 * pi * (2 * freq)`: Angular frequency of the first harmonic.

### Acoustic Noise Model
Acoustic power is calculated from surface velocities:

`P_acoustic = rho_air * c_air * S * v_rms^2`

- `rho_air = 1.225 kg/m^3`: Air density.
- `c_air = 343 m/s`: Speed of sound in air.
- `S = 2 * (domain_width + domain_height) * 0.1`: Surface area of the transformer.
- `v_rms`: Root-mean-square velocity, derived from displacement solution.

Sound pressure level at 1 meter is:

`p_rms = sqrt((P_acoustic * rho_air * c_air) / (4 * pi * r^2))`

`SPL = 20 * log10(p_rms / p_0)`

- `r = 1.0 m`: Distance from the transformer.
- `p_0 = 20e-6 Pa`: Reference pressure.

### Material Properties
- Core: Density = 7650 kg/m^3 (Silicon Steel), Poisson’s ratio = 0.3.
- Windings: Density = 8960 kg/m^3 (Copper), Poisson’s ratio = 0.34.

### Algorithms Utilized

- **Finite Element Method (FEM)**:
  - **Mesh Generation**: Creates a 20x10 triangular mesh over a 2D domain (core side length = `sqrt(core_area)`, winding thickness = 0.01 m).
  - **Element Assembly**:
    - Computes element-wise mass (`M`), stiffness (`K`), and force (`F`) matrices.
    - Stiffness matrix uses Young’s modulus (`E`), Poisson’s ratio (`nu`), and element geometry:
      `D = (E / (1 - nu^2)) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu)/2]]`
      `K_e = area * depth * (B^T * D * B)`
      where `B` is the strain-displacement matrix derived from element node coordinates.
    - Mass matrix: `M_e = (rho * area * depth / 3) * I`, where `rho` is density.
    - Force vector: Averages magnetostriction, electromagnetic, and Lorentz forces over elements.
  - **Boundary Conditions**: Fixes nodes at `y=0` by adding large penalties (1e12) to diagonal entries of `K`.
  - **System Solution**: Solves the harmonic equation using a sparse least-squares solver (`spla.lsqr`):
    `A = -omega^2 * M + i * omega * C + K`
    `u = lsqr(A, F_vec)`
    Extracts real (displacement) and imaginary (velocity) components: `u_x`, `u_y`, `v_x = omega * imag(u_x)`, `v_y = omega * imag(u_y)`.

- **Acoustic Noise Calculation**:
  - Computes RMS velocity from displacement solution.
  - Calculates acoustic power and sound pressure level using the equations above.
  - Reports dominant frequencies as 2x and 4x the input frequency, reflecting harmonic content.

- **Visualization**:
  - Computes displacement magnitude: `sqrt(u_x^2 + u_y^2)`.
  - Generates a contour plot using Matplotlib’s `contourf`, scaling displacement to micrometers.
  - Handles edge cases (e.g., zero or NaN displacements) by setting a minimum visualization value (1e-9 m).
