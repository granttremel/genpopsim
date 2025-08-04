

import numpy as np
import random
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
from dataclasses import dataclass, field
from enum import Enum
import scipy.stats as stats

from ..core.allele import Allele
from ..core.chromosome import Chromosome

class Individual:
    """An individual in the population"""
    def __init__(self, name:str, chromosomes: List[Chromosome], ploidy: int = 2, energy:float=0.0):
        self.name = name
        self.chromosomes = chromosomes
        self.genes = list(self.chromosomes[0].alleles.keys())
        self.ploidy = ploidy
        self.age = 0
        self.fitness_calculator = None
        self.energy:float = energy
    
    @property
    def offspring_ratio(self)->float:
        if "offspring_ratio" in self.genes:
            return sum([c.alleles["offspring_ratio"].params.ratio for c in self.chromosomes])/2
        else:
            return 0.25
        
    @property
    def offspring_number(self)->int:
        if "offspring_number" in self.genes:
            return sum([c.alleles["offspring_number"].params.number for c in self.chromosomes])//2
        else:
            return 3
    
    @property
    def metabolic_efficiency(self) -> float:
        # base_efficiency = self._calculate_fitness()
        base_efficiency = 1.0
        return base_efficiency
    
    @property
    def maintenance_cost(self) -> float:
        """Energy cost to maintain organization (fight entropy)"""
        # Could depend on:
        # - Genome complexity (more genes = more maintenance)
        # - Age (older = higher maintenance)
        # - Size/ploidy (more cells = more cost)
        base_cost = 0.1
        age_factor = 1 + (self.age * 0.01)  # Aging increases entropy
        complexity_factor = 1 + (len(self.chromosomes) * 0.01)
        return base_cost * age_factor * complexity_factor
    
    def harvest_energy(self, available_free_energy: float, competition_factor: float = 1.0) -> float:
        """Extract free energy from environment"""
        # Energy harvested depends on:
        # - Metabolic efficiency (from genes)
        # - Available free energy in environment
        # - Competition (more individuals = less per capita)
        harvested = self.metabolic_efficiency * available_free_energy * competition_factor
        self.energy += harvested
        return harvested
    
    def pay_maintenance(self) -> bool:
        """Pay entropy tax - return False if can't afford it (death)"""
        cost = self.maintenance_cost
        if self.energy >= cost:
            self.energy -= cost
            self.age += 1
            return True
        else:
            return False  # Death by energy depletion
    
    def can_reproduce(self, reproduction_cost: float) -> bool:
        """Check if enough energy to reproduce"""
        return self.energy >= reproduction_cost
    
    def mutate_allele(self, allele: Allele, 
                    mutation_rate: float = 0.001, 
                    effect_size: float = 0.1
                    ) -> Allele:
        """
        Mutate an allele with realistic fitness effects.
        
        Args:
            allele: The allele to potentially mutate
            mutation_rate: Probability of mutation per allele per generation
            beneficial_prob: Probability that a mutation is beneficial (vs deleterious)
            effect_size: Standard deviation of mutational effects
        
        Returns:
            Either the original allele or a new mutated allele
        """
        if random.random() > mutation_rate:
            return allele  # No mutation
        else:
            return allele
        
        benefit_change = np.random.normal(0, effect_size)
        new_benefit = np.clip(allele.benefit + benefit_change,0.0,1.0)
        
        cost_change = np.random.normal(0, effect_size)
        new_cost = np.clip(allele.cost + cost_change,0.0,1.0)
            
        
        return Allele(id=-1,gene_name=allele.gene_name, benefit=new_benefit, cost = new_cost)
    
    def get_gamete(self, energy_cost:float) -> List[Chromosome]:
        """
        For sexual reproduction: return one chromosome from each homologous pair.
        For now, assuming free recombination (no linkage).
        """
        if energy_cost > self.energy:
            print('individual cant pay energy cost!')
            return []
        self.energy -= energy_cost
        if self.ploidy == 1:
            # Haploid: return copy of chromosome
            return [self.chromosomes[0].copy()]
        elif self.ploidy == 2:
            # Diploid: randomly select one from each homologous pair
            newchrom = Chromosome()
            # gamete = []
            for g in self.genes:
                # Randomly choose maternal or paternal
                chosen = random.choice(range(self.ploidy))
                chosen_allele = self.chromosomes[chosen].alleles[g]
                mutated_allele = self.mutate_allele(chosen_allele)
                newchrom.add_allele(mutated_allele)
            return [newchrom]
        else:
            raise NotImplementedError(f"Ploidy {self.ploidy} not yet supported")
        
    def __str__(self):
        return f"Individual(name={self.name},age={self.age},energy={self.energy:0.3f})"
