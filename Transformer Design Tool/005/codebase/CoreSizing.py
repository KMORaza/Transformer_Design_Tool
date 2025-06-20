import math

class CoreSizing:
    def calculate_core_area(self, transformer_type, core_type, power, frequency):
        # Transformer type-specific core sizing constants
        type_params = {
            "Power": {"k": 0.0012},
            "Distribution": {"k": 0.0011},
            "Isolation": {"k": 0.0010},
            "Current": {"k": 0.0008},
            "Potential": {"k": 0.0009},
            "High-Frequency": {"k": 0.0007}
        }
        
        # Core type adjustment
        core_factor = 1.0 if core_type == "Core-Type" else 1.2  # Shell-Type requires larger core area
        
        # Get base constant for transformer type
        params = type_params.get(transformer_type, type_params["Power"])
        k = params["k"] * core_factor
        
        # Calculate core area (m^2) using empirical formula: A = k * (VA^0.5 / f)
        core_area = k * (math.sqrt(power) / frequency)
        
        return core_area