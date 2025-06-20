import math
from CoreSizing import CoreSizing
from RatioAndWireGaugeCalc import RatioAndWireGaugeCalc
from InterleavingAndLayering import InterleavingAndLayering

class Transformer:
    def calculate(self, transformer_type, core_type, Vp, Vs_list, freq, power_list, efficiency, core_material="Silicon Steel", winding_config="Concentric", bobbin_width=0.02, bobbin_height=0.03, insulation_thickness=0.0005, primary_layers=None, secondary_layers=None, ambient_temp=20.0, operating_temp=80.0):
        # Input validation
        if not all(Vs > 0 for Vs in Vs_list) or Vp <= 0 or freq <= 0 or not all(power > 0 for power in power_list) or efficiency <= 0 or efficiency > 1:
            raise ValueError("All voltages, frequency, and power must be positive; efficiency must be between 0 and 1.")
        if ambient_temp < -50 or operating_temp < -50 or operating_temp > 200:
            raise ValueError("Temperatures must be between -50°C and 200°C.")
        
        # Transformer type-specific parameters
        type_params = {
            "Power": {"core_loss_coeff": 0.02, "max_flux_density": 1.5, "coupling_coeff": 0.98},
            "Distribution": {"core_loss_coeff": 0.015, "max_flux_density": 1.4, "coupling_coeff": 0.97},
            "Isolation": {"core_loss_coeff": 0.01, "max_flux_density": 1.3, "coupling_coeff": 0.95},
            "Current": {"core_loss_coeff": 0.005, "max_flux_density": 1.0, "coupling_coeff": 0.99},
            "Potential": {"core_loss_coeff": 0.008, "max_flux_density": 1.1, "coupling_coeff": 0.98},
            "High-Frequency": {"core_loss_coeff": 0.03, "max_flux_density": 0.8, "coupling_coeff": 0.90}
        }
        
        # Core material parameters
        material_params = {
            "Silicon Steel": {
                "saturation_flux_density": 1.8,
                "relative_permeability": 4000,
                "steinmetz_k": 0.1,
                "steinmetz_alpha": 1.6,
                "steinmetz_beta": 2.0
            },
            "Ferrite": {
                "saturation_flux_density": 0.5,
                "relative_permeability": 2000,
                "steinmetz_k": 0.05,
                "steinmetz_alpha": 1.5,
                "steinmetz_beta": 2.5
            }
        }
        
        # Winding configuration parameters
        winding_params = {
            "Concentric": {"leakage_factor": 0.05, "stray_field_factor": 0.1},
            "Interleaved": {"leakage_factor": 0.02, "stray_field_factor": 0.5}
        }
        
        # Core type adjustments
        flux_density_factor = 1.0 if core_type == "Core-Type" else 0.9
        
        params = type_params.get(transformer_type, type_params["Power"])
        material = material_params.get(core_material, material_params["Silicon Steel"])
        winding = winding_params.get(winding_config, winding_params["Concentric"])
        core_loss_coeff = params["core_loss_coeff"]
        max_flux_density = min(params["max_flux_density"] * flux_density_factor, material["saturation_flux_density"])
        coupling_coeff = params["coupling_coeff"]
        leakage_factor = winding["leakage_factor"]
        stray_field_factor = winding["stray_field_factor"]
        mu_r = material["relative_permeability"]
        mu_0 = 4 * math.pi * 1e-7
        mu = mu_r * mu_0
        steinmetz_k = material["steinmetz_k"]
        steinmetz_alpha = material["steinmetz_alpha"]
        steinmetz_beta = material["steinmetz_beta"]
        
        # Calculate total power for core sizing
        total_power = sum(power_list)
        
        # Calculate core area using CoreSizing
        core_sizing = CoreSizing()
        core_area = core_sizing.calculate_core_area(transformer_type, core_type, total_power, freq)
        
        # Calculate turns ratio and wire gauges for each secondary
        ratio_calc = RatioAndWireGaugeCalc()
        turns_ratios = []
        secondary_wire_gauges = []
        for Vs, power in zip(Vs_list, power_list):
            turns_ratio, _, secondary_wire_gauge = ratio_calc.calculate(Vp, Vs, power, efficiency, transformer_type)
            turns_ratios.append(turns_ratio)
            secondary_wire_gauges.append(secondary_wire_gauge)
        
        # Use total power for primary wire gauge calculation
        primary_wire_gauge = ratio_calc.calculate(Vp, Vs_list[0], total_power, efficiency, transformer_type)[1]
        
        # Calculate number of turns using transformer equation
        primary_turns = (Vp * 1e8) / (4.44 * freq * core_area * max_flux_density)
        secondary_turns = [primary_turns / tr for tr in turns_ratios]
        
        # Calculate interleaving and layering
        interleaving = InterleavingAndLayering()
        interleaving_pattern, default_primary_layers, default_secondary_layers = interleaving.calculate(
            transformer_type, core_type, primary_turns, sum(secondary_turns), 
            primary_wire_gauge, secondary_wire_gauges[0], core_area
        )
        
        # Override default layers if specified
        primary_layers = primary_layers if primary_layers is not None else default_primary_layers
        secondary_layers_list = secondary_layers if secondary_layers is not None else [default_secondary_layers] * len(Vs_list)
        
        # Calculate mutual inductance for each secondary
        l = 0.1
        Lp = primary_turns ** 2 * mu * core_area / l
        mutual_inductances = []
        for Ns in secondary_turns:
            Ls = Ns ** 2 * mu * core_area / l
            M = coupling_coeff * math.sqrt(Lp * Ls)
            mutual_inductances.append(M)
        
        # Calculate currents
        Ip = total_power / (Vp * efficiency) if Vp * efficiency != 0 else 0.0
        Is_list = [power / (Vs * efficiency) if Vs * efficiency != 0 else 0.0 for Vs, power in zip(Vs_list, power_list)]
        
        # Calculate winding resistances (temperature-corrected)
        resistivity_cu = 1.68e-8
        alpha_cu = 0.00393
        T_ref = 20.0
        wire_diameter = 0.001
        wire_area = math.pi * (wire_diameter / 2) ** 2
        winding_length_primary = primary_turns * 0.1
        winding_length_secondary = [Ns * 0.1 for Ns in secondary_turns]
        
        # Resistances at ambient temperature
        R_primary_ambient = resistivity_cu * winding_length_primary / wire_area * (1 + alpha_cu * (ambient_temp - T_ref))
        R_secondary_ambient = [resistivity_cu * length / wire_area * (1 + alpha_cu * (ambient_temp - T_ref)) for length in winding_length_secondary]
        
        # Resistances at operating temperature
        R_primary_operating = resistivity_cu * winding_length_primary / wire_area * (1 + alpha_cu * (operating_temp - T_ref))
        R_secondary_operating = [resistivity_cu * length / wire_area * (1 + alpha_cu * (operating_temp - T_ref)) for length in winding_length_secondary]
        
        # Calculate copper losses
        copper_loss_primary_ambient = Ip ** 2 * R_primary_ambient
        copper_loss_secondary_ambient = sum(Is ** 2 * R for Is, R in zip(Is_list, R_secondary_ambient))
        total_copper_loss_ambient = copper_loss_primary_ambient + copper_loss_secondary_ambient
        
        copper_loss_primary_operating = Ip ** 2 * R_primary_operating
        copper_loss_secondary_operating = sum(Is ** 2 * R for Is, R in zip(Is_list, R_secondary_operating))
        total_copper_loss_operating = copper_loss_primary_operating + copper_loss_secondary_operating
        
        # Calculate core losses using Steinmetz equation
        core_volume = core_area * 0.1
        hysteresis_loss = steinmetz_k * freq ** steinmetz_alpha * max_flux_density ** steinmetz_beta * core_volume
        eddy_current_loss = 0.5 * steinmetz_k * (freq ** 2) * (max_flux_density ** 2) * core_volume
        total_core_loss = hysteresis_loss + eddy_current_loss
        
        # Use operating temperature for default copper losses
        core_loss = total_core_loss
        total_copper_loss = total_copper_loss_operating
        copper_loss_primary = copper_loss_primary_operating
        copper_loss_secondary = copper_loss_secondary_operating
        
        # Calculate flux density variation
        flux_densities = []
        time_points = 1000
        t = [i / (freq * time_points) for i in range(time_points)]
        primary_flux = [max_flux_density * math.sin(2 * math.pi * freq * ti) for ti in t]
        flux_densities.append(primary_flux)
        for i in range(len(Vs_list)):
            secondary_flux = [max_flux_density * math.sin(2 * math.pi * freq * ti) * (Vs_list[i] / Vp)
                             for ti in t]
            flux_densities.append(secondary_flux)
        
        # Calculate magnetic field strength (H)
        H_peak = (primary_turns * Ip) / l if l != 0 else 0.0
        H = [H_peak * math.sin(2 * math.pi * freq * ti) for ti in t]
        
        # Calculate leakage inductance
        window_area = bobbin_width * bobbin_height
        insulation_area = insulation_thickness * (primary_layers + sum(secondary_layers_list))
        winding_area = window_area - insulation_area if window_area > insulation_area else 0.001
        fill_factor = (primary_turns * wire_diameter**2 * primary_layers + sum(Ns * wire_diameter**2 * sl for Ns, sl in zip(secondary_turns, secondary_layers_list))) / winding_area
        interleaving_factor = 1.0 if winding_config == "Concentric" else (primary_layers + sum(secondary_layers_list)) / max(primary_layers, max(secondary_layers_list))
        leakage_factor_adj = leakage_factor * interleaving_factor
        leakage_inductance_primary = leakage_factor_adj * mu_0 * primary_turns ** 2 * window_area / bobbin_height
        leakage_inductance_secondary = []
        for Ns, sl in zip(secondary_turns, secondary_layers_list):
            L_leak_sec = leakage_factor_adj * mu_0 * Ns ** 2 * window_area / bobbin_height
            leakage_inductance_secondary.append(L_leak_sec)
        
        # Calculate stray fields
        distances = [i * 0.001 for i in range(1, 1001)]
        stray_fields = []
        core_dimension = math.sqrt(core_area)
        for d in distances:
            H_stray = stray_field_factor * (Ip * primary_turns) / (2 * math.pi * (d + core_dimension))
            stray_fields.append(H_stray)
        
        safe_distance = next((d for d, H in zip(distances, stray_fields) if H < 0.1), 1.0)
        
        return {
            'turns_ratios': turns_ratios,
            'primary_turns': primary_turns,
            'secondary_turns': secondary_turns,
            'core_loss': core_loss,
            'hysteresis_loss': hysteresis_loss,
            'eddy_current_loss': eddy_current_loss,
            'copper_loss_primary': copper_loss_primary,
            'copper_loss_secondary': copper_loss_secondary,
            'total_copper_loss': total_copper_loss,
            'copper_loss_primary_ambient': copper_loss_primary_ambient,
            'copper_loss_secondary_ambient': copper_loss_secondary_ambient,
            'total_copper_loss_ambient': total_copper_loss_ambient,
            'R_primary_ambient': R_primary_ambient,
            'R_secondary_ambient': R_secondary_ambient,
            'R_primary_operating': R_primary_operating,
            'R_secondary_operating': R_secondary_operating,
            'efficiency': efficiency,
            'flux_densities': flux_densities,
            'time_points': t,
            'core_area': core_area,
            'primary_wire_gauge': primary_wire_gauge,
            'secondary_wire_gauges': secondary_wire_gauges,
            'interleaving_pattern': interleaving_pattern,
            'primary_layers': primary_layers,
            'secondary_layers': secondary_layers_list,
            'mutual_inductances': mutual_inductances,
            'max_flux_density': max_flux_density,
            'H_field': H,
            'H_peak': H_peak,
            'saturation_flux_density': material["saturation_flux_density"],
            'leakage_inductance_primary': leakage_inductance_primary,
            'leakage_inductance_secondary': leakage_inductance_secondary,
            'stray_fields': stray_fields,
            'distances': distances,
            'safe_distance': safe_distance,
            'fill_factor': fill_factor,
            'window_area': window_area
        }