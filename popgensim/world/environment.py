from dataclasses import dataclass
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
import numpy as np

@dataclass
class EnvironmentConfig:
    """All environmental parameters in one place"""
    # Energy parameters
    
    free_energy_rate: float = 0.04
    holding_cost: float = 0.04
    
    # solar_temperature: float = 6000  # K
    # surface_temperature: float = 288  # K  
    # heat_sink_temperature: float = 3  # K (cosmic background)
    
    # # Spatial parameters
    grid_size: int = 100
    diffusion_rate: float = 0.1
    
    # # Temporal parameters
    # day_length: int = 10  # timesteps
    # season_length: int = 100  # timesteps

class Environment:
    """Manages physical environment and resources"""
    def __init__(self, config: EnvironmentConfig):
        self.config = config
        # self.grid = np.zeros((config.grid_size, config.grid_size))
        # self.temperature_map = np.full((config.grid_size, config.grid_size), 
                                    #   config.surface_temperature)
        self.time = 0
        
    # def calculate_free_energy_flux(self, x: int, y: int) -> float:
    #     """Free energy available at a location"""
    #     # Solar input varies with time (day/night, seasons)
    #     solar_phase = np.sin(2 * np.pi * self.time / self.config.day_length)
    #     seasonal_phase = np.sin(2 * np.pi * self.time / self.config.season_length)
        
    #     T_hot = self.config.solar_temperature * max(0, solar_phase)
    #     T_cold = self.temperature_map[x, y]
    #     T_sink = self.config.heat_sink_temperature
        
    #     # Energy flux from gradient
    #     flux_in = (T_hot**4 - T_cold**4) * 1e-8  # Simplified Stefan-Boltzmann
    #     flux_out = (T_cold**4 - T_sink**4) * 1e-8
        
    #     return max(0, flux_in - flux_out)
    
    def update(self):
        """Update environmental conditions"""
        self.time += 1
        # Could add: weather, climate change, etc.
        
    # def __getattr__(self, attr):
    #     if hasattr(self.config, attr):
    #         return getattr(self.config, attr)
