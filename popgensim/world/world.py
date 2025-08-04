
from dataclasses import dataclass
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
import numpy as np


from ..population.population import Population, PopulationConfig
from ..genetics.genes import Gene
from ..genetics.genome import Genome
from .environment import Environment, EnvironmentConfig
from .energetics import EnergyCalculator

class World:
    """Top-level container managing the entire simulation"""
    def __init__(self, env_config: EnvironmentConfig):
        self.environment = Environment(env_config)
        self.populations: Dict[str, Population] = {}
        self.generation = 0
        
        self.energetics = EnergyCalculator()
        
        # Global tracking
        self.extinct_species: List[str] = []
        # self.invasion_queue: List[PopulationConfig] = []
        
    def add_population(self, pop_config: PopulationConfig, genome: Genome) -> None:
        """Add a new population/species"""
        # Factory pattern for creating populations
        population = self._create_population(pop_config, genome)
        self.populations[pop_config.species_name] = population
    
    def _create_population(self, config: PopulationConfig, genome: Genome) -> Population:
        """Factory method - could be overridden for different population types"""
        # Import here to avoid circular imports
        return Population(config, genome)
    
    def _calculate_total_energy(self):
        #TBI
        return 1.0
    
    def timestep(self):
        """Coordinate all populations and environment"""
        self.generation += 1
        
        # 1. Update environment
        self.environment.update()
        
        # 2. Calculate resource availability
        total_energy_flux = self._calculate_total_energy()
        
        # 3. Populations compete for resources
        # energy_allocation = self._allocate_energy(total_energy_flux)
        
        # 4. Each population evolves
        for species_name, population in list(self.populations.items()):
            # population.available_energy = energy_allocation.get(species_name, 0)
            # population.timestep(self.environment)
            population.timestep()
            
            # Check extinction
            if len(population.individuals) == 0:
                self.extinct_species.append(species_name)
                del self.populations[species_name]
        
        # 5. Check for speciation events
        # self._check_speciation()
        
        # 6. Process invasions
        # self._process_invasions()
    
    # def _allocate_energy(self, total_flux: float) -> Dict[str, float]:
    #     """Divide energy among populations based on competition"""
    #     if not self.populations:
    #         return {}
        
    #     # Simple version - proportional to biomass
    #     allocations = {}
    #     total_demand = sum(pop.get_energy_demand() 
    #                       for pop in self.populations.values())
        
    #     if total_demand == 0:
    #         return {name: total_flux / len(self.populations) 
    #                for name in self.populations}
        
    #     for name, pop in self.populations.items():
    #         fraction = pop.get_energy_demand() / total_demand
    #         allocations[name] = total_flux * fraction
        
    #     return allocations
    
    def _check_speciation(self):
        """Detect if populations have diverged enough to speciate"""
        pass
        # for name, population in self.populations.items():
        #     # Simple version - check genetic distance
        #     if self._population_has_diverged(population):
        #         # Create new species
        #         new_config = self._create_species_config(population)
        #         self.invasion_queue.append(new_config)

# Builder pattern for complex initialization
class WorldBuilder:
    """Fluent interface for building complex worlds"""
    def __init__(self):
        self.env_config = EnvironmentConfig()
        self.population_configs = []
    
    # def with_temperature_gradient(self, hot: float, cold: float) -> 'WorldBuilder':
    #     self.env_config.solar_temperature = hot
    #     self.env_config.surface_temperature = cold
    #     return self
    
    def with_grid_size(self, size: int) -> 'WorldBuilder':
        self.env_config.grid_size = size
        return self
    
    def add_species(self, config: PopulationConfig, genome: Genome) -> 'WorldBuilder':
        self.population_configs.append((config, genome))
        return self
    
    def build(self) -> World:
        world = World(self.env_config)
        for pop_config in self.population_configs:
            world.add_population(*pop_config)
        return world
    
    
    
def main():
    # Usage example
    def create_simple_ecosystem():
        genome = Genome()
        for n in range(10):
            genome.add_gene(Gene(f"gene{n}"))
        
        world = (WorldBuilder()
            # .with_temperature_gradient(hot=6000, cold=288)
            # .with_grid_size(50)
            .add_species(PopulationConfig(
                species_name="phototrophs",
                initial_size=100,
            ),
            genome)
            .add_species(PopulationConfig(
                species_name="chemotrophs", 
                initial_size=50,
            ),
            genome)
            .build()
        )
        return world

if __name__=="__main__":
    main()