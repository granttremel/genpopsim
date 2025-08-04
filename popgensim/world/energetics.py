import numpy as np
from typing import List, Dict, Tuple
# from . import Allele, Chromosome, Individual, Population
from ..core.allele import Allele
from ..core.chromosome import Chromosome
from ..population.individual import Individual
from ..genetics.genes import Gene
from ..genetics.genome import Genome
from ..world.environment import Environment

class EnergyCalculator:
    """Calculate fitness based on gene network"""
    
    def __init__(self):
        
        pass
    
    def calculate_energy_change(self, individual:Individual, genome: Genome, environment:Environment) -> float:
        chromosomes = individual.chromosomes
        gain = 0.0
        loss = 0.0
        for c in chromosomes:
            expression = genome.calculate_expression(c.alleles)
            for g,a in c.alleles.items():
                if hasattr(a.params,"benefit"):
                    gain += a.expression * a.params.benefit
                    
                if hasattr(a.params,"cost"):
                    loss += a.expression * a.params.cost
                
        loss += individual.energy * environment.config.holding_cost
    
        return gain, loss

    def _get_alleles_for_gene(self, gene_name: str, chromosomes: List[Chromosome]) -> List[Allele]:
        """Get all alleles for a gene across chromosomes"""
        alleles = []
        for chromosome in chromosomes:
            allele = chromosome.get_allele(gene_name)
            if allele:
                alleles.append(allele)
        return alleles