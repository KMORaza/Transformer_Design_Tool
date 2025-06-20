import math
from CoreSizing import CoreSizing
from RatioAndWireGaugeCalc import RatioAndWireGaugeCalc
from InterleavingAndLayering import InterleavingAndLayering

class Transformer:
    def calculate(self, transformer_type, core_type, Vp, Vs_list, freq, power_list, efficiency):
        # Transformer type-specific parameters
        type_params = {
            "Power": {"core_loss_coeff": 0.02, "max_flux_density": 1.5, "coupling_coeff": 0.98},
            "Distribution": {"core_loss_coeff": 0.015, "max_flux_density": 1.4, "coupling_coeff": 0.97},
            "Isolation": {"core_loss_coeff": 0.01, "max_flux_density": 1.3, "coupling_coeff": 0.95},
            "Current": {"core_loss_coeff": 0.005, "max_flux_density": 1.0, "coupling_coeff": 0.99},
            "Potential": {"core_loss_coeff": 0.008, "max_flux_density": 1.1, "coupling_coeff": 0.98},
            "High-Frequency": {"core_loss_coeff": 0.03, "max_flux_density": 0.8, "coupling_coeff": 0.90}
        }
        
        # Core type adjustments
        flux_density_factor = 1.0 if core_type == "Core-Type" else 0.9  # Shell-Type has lower flux density
        
        params = type_params.get(transformer_type, type_params["Power"])
        core_loss_coeff = params["core_loss_coeff"]
        max_flux_density = params["max_flux_density"] * flux_density_factor
        coupling_coeff = params["coupling_coeff"]
        
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
        interleaving_pattern, primary_layers, secondary_layers = interleaving.calculate(
            transformer_type, core_type, primary_turns, sum(secondary_turns), 
            primary_wire_gauge, secondary_wire_gauges[0], core_area
        )
        
        # Calculate mutual inductance for each secondary
        # L = N^2 * mu * A / l (simplified inductance model)
        mu = 4 * math.pi * 1e-7  # Permeability of free space
        l = 0.1  # Effective magnetic path length (m, assumed)
        Lp = primary_turns ** 2 * mu * core_area / l
        mutual_inductances = []
        for Ns in secondary_turns:
            Ls = Ns ** 2 * mu * core_area / l
            M = coupling_coeff * math.sqrt(Lp * Ls)
            mutual_inductances.append(M)
        
        # Core loss estimation (proportional to total power)
        core_loss = total_power * (1 - efficiency) + (total_power * core_loss_coeff)
        
        # Calculate flux density variation for each secondary
        flux_densities = []
        for i in range(len(Vs_list)):
            flux_density = [max_flux_density * math.sin(2 * math.pi * freq * t / 1000) 
                           for t in range(1000)]
            flux_densities.append(flux_density)
        
        return {
            'turns_ratios': turns_ratios,
            'primary_turns': primary_turns,
            'secondary_turns': secondary_turns,
            'core_loss': core_loss,
            'efficiency': efficiency,
            'flux_densities': flux_densities,
            'core_area': core_area,
            'primary_wire_gauge': primary_wire_gauge,
            'secondary_wire_gauges': secondary_wire_gauges,
            'interleaving_pattern': interleaving_pattern,
            'primary_layers': primary_layers,
            'secondary_layers': [secondary_layers] * len(Vs_list),  # Simplified: same layers for all secondaries
            'mutual_inductances': mutual_inductances
        }