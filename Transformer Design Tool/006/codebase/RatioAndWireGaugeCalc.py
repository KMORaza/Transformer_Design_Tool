import math

class RatioAndWireGaugeCalc:
    def calculate(self, Vp, Vs, power, efficiency, transformer_type):
        # Calculate currents
        Ip = power / (Vp * efficiency)  # Primary current
        Is = power / Vs  # Secondary current
        
        # Adjust turns ratio for voltage regulation (simplified)
        regulation_factor = 1.05 if transformer_type in ["Power", "Distribution"] else 1.02
        turns_ratio = (Vp / Vs) * regulation_factor
        
        # Wire gauge calculation based on current density (2 A/mm^2 typical)
        current_density = 2.0  # A/mm^2
        # AWG table: {AWG: area in mm^2}
        awg_table = {
            10: 5.261, 12: 3.309, 14: 2.081, 16: 1.309, 18: 0.823,
            20: 0.518, 22: 0.326, 24: 0.205, 26: 0.129, 28: 0.081
        }
        
        # Calculate required wire area (mm^2)
        primary_wire_area = Ip / current_density
        secondary_wire_area = Is / current_density
        
        # Select AWG size (smallest AWG with area >= required)
        primary_wire_gauge = min((awg for awg, area in awg_table.items() if area >= primary_wire_area), default=10)
        secondary_wire_gauge = min((awg for awg, area in awg_table.items() if area >= secondary_wire_area), default=10)
        
        return turns_ratio, primary_wire_gauge, secondary_wire_gauge