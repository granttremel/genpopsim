from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Tuple

class StrategicGeneType(Enum):
    # Energy acquisition
    METABOLIC_HARVESTER = "harvester"      # Increases energy capture
    METABOLIC_SPECIALIST = "specialist"     # Better in specific conditions
    
    # Energy allocation
    REPRODUCTION_TIMING = "repro_timing"    # When to reproduce (energy threshold)
    OFFSPRING_INVESTMENT = "offspring_inv"   # Energy per offspring
    OFFSPRING_NUMBER = "offspring_num"       # How many offspring to produce
    
    # Complexity management
    GENE_DUPLICATOR = "duplicator"          # Chance to duplicate genes
    GENE_SILENCER = "silencer"              # Can turn off costly genes
    COMPLEXITY_TOLERANCE = "complexity_tol"  # Reduces complexity penalty
    
    # Energy consumption
    ENERGY_CONSUMER = "consumer"            # Uses energy for some benefit
    STORAGE_PROTEIN = "storage"             # Stores energy for later
    DEFENSE_SYSTEM = "defense"              # Costs energy but provides protection
    SENSOR = "sensor"                       # Costs energy but provides information
    
    # Social/cooperative
    PUBLIC_GOOD = "public_good"             # Produces benefit for all
    CHEATER = "cheater"                     # Exploits public goods
    KIN_RECOGNITION = "kin_recognition"     # Preferentially helps relatives

@dataclass
class StrategicAllele:
    """Allele that encodes behavioral strategy"""
    id: int
    gene_name: str
    value: float  # Strategic parameter value
    energy_cost: float = 0.0  # Ongoing energy cost
    energy_production: float = 0.0  # Energy production rate
    
class StrategicIndividual(ThermodynamicIndividual):
    """Individual with evolvable strategies"""
    
    def get_strategy(self, strategy_type: str) -> float:
        """Get strategic parameter from genome"""
        # Search chromosomes for strategy genes
        for chromosome in self.chromosomes:
            for allele in chromosome.alleles:
                if hasattr(allele, 'gene_name') and strategy_type in allele.gene_name:
                    return allele.value
        return 0.5  # Default strategy
    
    @property
    def reproduction_threshold(self) -> float:
        """Energy level at which to reproduce"""
        base_threshold = 1.0
        timing_gene = self.get_strategy("repro_timing")
        # Gene value 0-1 maps to threshold 0.5-2.0
        return base_threshold * (0.5 + 1.5 * timing_gene)
    
    @property
    def offspring_investment_ratio(self) -> float:
        """Fraction of energy to give each offspring"""
        investment_gene = self.get_strategy("offspring_inv")
        # Maps to 0.1-0.9 of available energy
        return 0.1 + 0.8 * investment_gene
    
    @property
    def offspring_number_strategy(self) -> int:
        """How many offspring to produce at once"""
        number_gene = self.get_strategy("offspring_num")
        # Maps to 1-5 offspring
        return max(1, int(1 + 4 * number_gene))
    
    def calculate_metabolic_rate(self) -> float:
        """Total energy production minus consumption"""
        production = 0.0
        consumption = 0.0
        
        for chromosome in self.chromosomes:
            for allele in chromosome.alleles:
                if hasattr(allele, 'energy_production'):
                    production += allele.energy_production
                if hasattr(allele, 'energy_cost'):
                    consumption += allele.energy_cost
        
        # Base metabolic efficiency
        base_rate = self.metabolic_efficiency
        
        # Modified by producer/consumer genes
        return base_rate + production - consumption
    
    def should_duplicate_genes(self) -> bool:
        """Decide whether to duplicate genes during reproduction"""
        dup_rate = self.get_strategy("duplicator")
        return np.random.random() < (dup_rate * 0.01)  # Max 1% chance
    
    def complexity_cost_modifier(self) -> float:
        """How well organism handles complexity"""
        tolerance = self.get_strategy("complexity_tol")
        # High tolerance reduces complexity penalty
        return 1.0 - (0.5 * tolerance)

class EnergyBalancedGene(Gene):
    """Gene with explicit energy costs/benefits"""
    def __init__(self, name: str, gene_type: GeneType, 
                 energy_cost: float = 0.0,
                 energy_production: float = 0.0,
                 complexity_cost: float = 0.01):
        super().__init__(name, gene_type)
        self.energy_cost = energy_cost
        self.energy_production = energy_production
        self.complexity_cost = complexity_cost

# Example: Create a metabolic pathway with trade-offs
def create_photosynthesis_pathway() -> List[EnergyBalancedGene]:
    """Energy-producing pathway with multiple steps"""
    genes = []
    
    # Light harvesting complex - cheap but essential
    genes.append(EnergyBalancedGene(
        "light_harvester",
        GeneType.PATHWAY,
        energy_cost=0.05,
        energy_production=0.3,
        complexity_cost=0.02
    ))
    
    # Electron transport - expensive but amplifies production
    genes.append(EnergyBalancedGene(
        "electron_transport",
        GeneType.PATHWAY,
        energy_cost=0.1,
        energy_production=0.5,
        complexity_cost=0.03
    ))
    
    # Carbon fixation - very expensive, high payoff
    genes.append(EnergyBalancedGene(
        "carbon_fixation",
        GeneType.PATHWAY,
        energy_cost=0.2,
        energy_production=0.8,
        complexity_cost=0.05
    ))
    
    return genes

def create_defense_system() -> List[EnergyBalancedGene]:
    """Costly defense that might provide advantage"""
    return [
        EnergyBalancedGene(
            "toxin_production",
            GeneType.MODULAR,
            energy_cost=0.15,  # Expensive!
            energy_production=0.0,  # No direct benefit
            complexity_cost=0.02
        ),
        EnergyBalancedGene(
            "toxin_resistance",
            GeneType.MODULAR,
            energy_cost=0.05,  # Cheaper than producing
            energy_production=0.0,
            complexity_cost=0.01
        )
    ]

# Population-level effects
class StrategicPopulation(ThermodynamicPopulation):
    """Population where strategies can evolve"""
    
    def calculate_public_goods(self) -> float:
        """Sum of public goods produced by population"""
        total_public_goods = 0.0
        for ind in self.individuals:
            public_good_level = ind.get_strategy("public_good")
            if public_good_level > 0:
                # Producer pays cost
                ind.energy -= public_good_level * 0.1
                total_public_goods += public_good_level
        
        return total_public_goods
    
    def distribute_public_goods(self, total_goods: float):
        """All individuals benefit from public goods"""
        if len(self.individuals) == 0:
            return
            
        per_capita_benefit = total_goods / len(self.individuals)
        
        for ind in self.individuals:
            # Cheaters get extra benefit
            cheater_level = ind.get_strategy("cheater")
            benefit = per_capita_benefit * (1 + cheater_level)
            ind.energy += benefit