

import numpy as np
import random
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
from dataclasses import dataclass, field
from enum import Enum
import scipy.stats as stats

class ReproductionMode(Enum):
    ASEXUAL = "asexual"
    SEXUAL = "sexual"

@dataclass
class Allele:
    """Single allele with associated fitness value"""
    id: int
    gene_name: str
    fitness: float  # Between 0 and 1
    _frequency: int=0
    discovered: int=0
    
    def __str__(self):
        return f"Allele(id={self.id},fitness={self.fitness:0.3f},{self._frequency}x)"

    def __repr__(self):
        return str(self)
    
    def __hash__(self):
        return hash((self.gene_name, self.id))


class Chromosome:
    """Container for alleles"""
    def __init__(self, genes: Iterable[str]):
        self.genes: List[str] =  list(genes)
        self.alleles: Dict[str,Allele] = {}
        self.n_genes = len(self.genes)
    
    def copy(self):
        """Create a copy of the chromosome"""
        return Chromosome.from_alleles(self.alleles.values())
    
    @classmethod
    def from_alleles(cls, alleles:Iterable[Allele]):
        chr = cls([a.gene_name for a in alleles])
        
        for a in alleles:
            chr.add_allele(a)
            
        return chr
        
    def add_allele(self, allele:Allele):
        if allele.gene_name in self.genes:
            self.alleles[allele.gene_name] = allele
            allele._frequency += 1
        else:
            print('gene not in chromosome!')
    
    def get_allele(self, gene_name:str):
        
        return self.alleles.get(gene_name)
    
    def done(self):
        b = all([self.alleles[g] for g in self.alleles])
        return b
    
    def __len__(self):
        return len(self.genes)
    
    def __str__(self):
        content = []
        for gene in self.genes:
            content.append(f"{gene}:{self.alleles[gene]}")
        contentstr = "\n".join(content)
        return f"Chromosome({contentstr})"
    
    def __repr__(self):
        return str(self)

class Individual:
    """An individual in the population"""
    def __init__(self, name:str, chromosomes: List[Chromosome], ploidy: int = 2):
        self.name = name
        self.chromosomes = chromosomes
        self.genes = list(self.chromosomes[0].genes)
        self.ploidy = ploidy
        self.age = 0
        self.fitness_calculator = None
        self._fitness:float = -1.0
    
    @property
    def fitness(self) -> float:
        return self._fitness
    
    def mutate_allele(self, allele: Allele, 
                    mutation_rate: float = 0.001, 
                    beneficial_prob: float = 0.5, 
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
        
        # Mutation occurs
        # self.allele_counter += 1
        
        # Determine if beneficial or deleterious
        if random.random() < beneficial_prob:
            # Beneficial mutation - increase fitness
            fitness_change = abs(np.random.normal(0, effect_size))
            new_fitness = min(1.0, allele.fitness + fitness_change)
        else:
            # Deleterious mutation - decrease fitness
            fitness_change = -abs(np.random.normal(0, effect_size))
            new_fitness = max(0.0, allele.fitness + fitness_change)
        
        return Allele(id=-1,gene_name=allele.gene_name, fitness=new_fitness)
    
    def get_gamete(self) -> List[Chromosome]:
        """
        For sexual reproduction: return one chromosome from each homologous pair.
        For now, assuming free recombination (no linkage).
        """
        if self.ploidy == 1:
            # Haploid: return copy of chromosome
            return [self.chromosomes[0].copy()]
        elif self.ploidy == 2:
            # Diploid: randomly select one from each homologous pair
            newchrom = Chromosome(self.chromosomes[0].genes)
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
        return f"Individual(name={self.name},age={self.age},fitness={self.fitness:0.3f})"
