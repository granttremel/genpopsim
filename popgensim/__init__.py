import numpy as np
import random
from typing import List, Dict, Tuple, Optional, Set, Union
from dataclasses import dataclass, field
from enum import Enum
import scipy.stats as stats

from .plotting import Tracker
from .genetics import Gene, Genome

class ReproductionMode(Enum):
    ASEXUAL = "asexual"
    SEXUAL = "sexual"

@dataclass
class Allele:
    """Single allele with associated fitness value"""
    id: int
    gene: str
    fitness: float  # Between 0 and 1
    _frequency: Optional[int]=0
    
    def __str__(self):
        return f"Allele(id={self.id},fitness={self.fitness:0.3f},{self._frequency}x)"

    def __repr__(self):
        return str(self)


class Chromosome:
    """Container for alleles"""
    def __init__(self, genes: List[str]):
        self.genes =  {g:None for g in genes}
        self.n_genes = len(self.genes)
    
    def copy(self):
        """Create a copy of the chromosome"""
        return Chromosome([Allele(a.id, a.fitness) for a in self.alleles])
    
    def add_allele(self, allele:Allele):
        if allele.gene in self.genes:
            self.genes[allele.gene] = allele
            allele._frequency += 1
        else:
            print('gene not in chromosome!')
    
    def done(self):
        b = all([self.genes[g] for g in self.genes])
        return b
    
    def __getitem__(self, ind):
        if isinstance(ind, int):
            for c, a in self.genes.items():
                if a.id == ind:
                    return a
        elif isinstance(ind, str):
            return self.genes.get(ind, None)
    
    def __len__(self):
        return len(self.genes)
    
    def __str__(self):
        content = []
        for gene in self.genes:
            content.append(f"{gene}:{self.genes[gene]}")
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
        self._fitness = None
    
    @property
    def fitness(self) -> float:
        """Calculate fitness using epistatic interactions"""
        if self._fitness is None:
            if self.fitness_calculator:
                self._fitness = self.fitness_calculator.calculate_fitness(self)
            else:
                # Fallback to additive
                all_fitness_values = []
                for chromosome in self.chromosomes:
                    all_fitness_values.extend([a.fitness for a in chromosome.genes.values()])
                self._fitness = np.mean(all_fitness_values) if all_fitness_values else 0.0
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
        
        return Allele(id=-1,gene=allele.gene, fitness=new_fitness)
    
    def get_gamete(self) -> List[Chromosome]:
        """
        For sexual reproduction: return one chromosome from each homologous pair.
        For now, assuming free recombination (no linkage).
        """
        if self.ploidy == 1:
            # Haploid: return copy of chromosome
            return self.chromosomes[0].copy()
        elif self.ploidy == 2:
            # Diploid: randomly select one from each homologous pair
            newchrom = Chromosome(self.chromosomes[0].genes)
            # gamete = []
            for g in self.genes:
                # Randomly choose maternal or paternal
                chosen = random.choice(range(self.ploidy))
                chosen_allele = self.chromosomes[chosen][g]
                mutated_allele = self.mutate_allele(chosen_allele)
                newchrom.add_allele(mutated_allele)
            return newchrom
        else:
            raise NotImplementedError(f"Ploidy {self.ploidy} not yet supported")
        
    def __str__(self):
        return f"Individual(name={self.name},age={self.age},fitness={self.fitness:0.3f})"

class Population:
    """Population of individuals with configurable genetics"""
    def __init__(self, 
            size: int,
          #  genes: List[str],
            genome:Genome,
            ploidy: int = 2,
            reproduction_mode: ReproductionMode = ReproductionMode.SEXUAL,
            parity:int=2,
            reproductive_rate = 0.25,
            mortality_rate = 0.1,
            allele_fitness_range: Tuple[float, float] = (0.3, 0.7),
            capacity: int = 0,
            lifespan = 20):
        
        self.size = size
        # self.n_genes_per_chromosome = len(genes)
        self.genome = genome
        self.ploidy = ploidy
        self.gene_pool = {}
        
        self.reproduction_mode = reproduction_mode
        self.parity = parity
        self.reproductive_rate = reproductive_rate
        self.reproductive_fitness_strength = 0.1
        self._mortality_rate = mortality_rate
        
        self.individuals: List[Individual] = []
        if not capacity:
            self.capacity = 5*size
        else:
            self.capacity = capacity
        self.generation = 0
        self.n_ever = 0
        self.allele_counter = 0
        self.lifespan = lifespan
        
        self.tracker = Tracker(self)
        
        # Initialize population
        self._initialize_population(2, allele_fitness_range)
    
    def _create_random_allele(self, gene_name, fitness_range: Tuple[float, float]) -> Allele:
        """Create a random allele with fitness in specified range"""
        # self.allele_counter += 1
        na = len(self.gene_pool[gene_name])
        fitness = np.random.uniform(*fitness_range)
        new_allele = Allele(id=na,gene=gene_name,fitness=fitness)
        self.gene_pool[gene_name].append(new_allele)
        return new_allele
    
    def _initialize_population(self, alleles_per_gene:int, fitness_range: Tuple[float, float]):
        """Create initial population with random alleles"""
        
        for g in self.genome.genes:
            if not g in self.gene_pool:
                self.gene_pool[g] = []
            for na in range(alleles_per_gene):
                _ = self._create_random_allele(g, fitness_range)
                
        for ni in range(self.size):
            chromosomes = []
            name = f"indi{ni}"
            # Create ploidy sets of chromosomes
            for _ in range(self.ploidy):
                newchr = Chromosome(self.genome.genes.keys())
                for g in self.genome.genes:
                    a = random.choice(self.gene_pool[g])
                    newchr.add_allele(a)
                chromosomes.append(newchr)
            self.individuals.append(Individual(name, chromosomes, self.ploidy))

    def reproduce_single(self, parent1: Individual, parent2: Optional[Individual] = None) -> Individual:
        """Create offspring from one or two parents"""
        if self.reproduction_mode == ReproductionMode.ASEXUAL:
            # Simple cloning with potential for mutation (to be added later)
            offspring_chromosomes = [c.copy() for c in parent1.chromosomes]
            newname = f"Indi{self.n_ever}"
            self.n_ever += 1
            newind = Individual(newname, offspring_chromosomes, self.ploidy)
            self.register_new_individual(newind)
            return newind
        
        elif self.reproduction_mode == ReproductionMode.SEXUAL:
            if parent2 is None:
                raise ValueError("Sexual reproduction requires two parents")
            # print(f'reproducing with {parent1} and {parent2}')
            
            # Get gametes from each parent
            gamete1 = parent1.get_gamete()
            gamete2 = parent2.get_gamete()
            
            # Combine gametes
            offspring_chromosomes = [gamete1, gamete2]
            newname = f"Indi{self.n_ever}"
            self.n_ever += 1
            newind = Individual(newname, offspring_chromosomes, self.ploidy)
            self.register_new_individual(newind)
            return newind
    
    def reproduce(self, n_offspring, parent1:Individual, parent2: Optional[Individual] = None):
        offs = list()
        for n in range(n_offspring):
            offs.append(self.reproduce_single(parent1, parent2 = parent2))
        # print(f'{parent1} and {parent2} produced {n_offspring} offspring!')
        return offs
    
    def register_new_individual(self, new_indi):
        for c in new_indi.chromosomes:
            for g in c.genes:
                if c[g].id < 0:
                    #new allele!
                    newid = len(self.gene_pool[g])
                    c[g].id = newid
                    self.gene_pool[g].append(c[g])
        
        pass
    
    def mortality_rate(self, individual: Individual) -> float:
        """Calculate death rate based on fitness (lower fitness = higher death rate)"""
        # Simple inverse relationship for now
        # Fitness of 1.0 -> death rate near 0
        # Fitness of 0.0 -> death rate of 1.0
        if callable(self._mortality_rate):
            mr = self._mortality_rate(self.generation)
        else:
            mr = self._mortality_rate
        return mr * (1.0 - individual.fitness)*len(self.individuals)/self.capacity
        # return -np.log(individual.fitness)*len(self.individuals)/self.capacity
    
    def timestep(self):
        """Execute one generation: reproduction and mortality"""
        self.generation += 1
        
        n_repros = int(len(self.individuals) * self.reproductive_rate)
        for nr in range(n_repros):
            # Reproduction
            if self.reproduction_mode == ReproductionMode.SEXUAL:
                # Random mating
                if len(self.individuals) >= 2:
                    parent1, parent2 = random.sample(self.individuals, 2)
                    pair_fit = parent1.fitness + parent2.fitness
                    n_offspring = int(stats.poisson.rvs(self.parity)*pair_fit)
                    n_offspring = min(n_offspring, 2*self.parity)
                    offspring = self.reproduce(n_offspring, parent1, parent2)
                    self.individuals.extend(offspring)
            else:
                # Asexual reproduction
                if self.individuals:
                    parent = random.choice(self.individuals)
                    n_offspring = stats.poisson.rvs(self.parity)
                    offspring = self.reproduce(n_offspring, parent)
                    self.individuals.append(offspring)
        
        # Mortality (Poisson process approximation)
        survivors = []
        for individual in self.individuals:
            if individual.age > self.lifespan:
                # print(f'{individual} died of old age!')
                continue
            death_prob = self.mortality_rate(individual)  # Scale factor
            if random.random() > death_prob:
                individual.age+=1
                survivors.append(individual)
            else:
                pass
                # print(f'{individual} succumbed to natural causes!')
        self.individuals = survivors
        self.tracker.update()
        
        
    def get_allele_frequencies(self) -> dict:
        """Calculate frequency of each allele in the population"""
        allele_counts = {}
        total_alleles = 0
        
        for individual in self.individuals:
            for chromosome in individual.chromosomes:
                for allele in chromosome.alleles:
                    allele_counts[allele.id] = allele_counts.get(allele.id, 0) + 1
                    total_alleles += 1
        
        if total_alleles == 0:
            return {}
        
        return {aid: count/total_alleles for aid, count in allele_counts.items()}
    
    def get_fitness_distribution(self) -> List[float]:
        """Get list of all individual fitness values"""
        return [ind.fitness for ind in self.individuals]
    
    
if __name__ == "__main__":
    
    max_gen = 200
    
    # mr = lambda gen:1+np.sin(gen/20)/2
    mr = 0.9
    
    # Create a population
    pop = Population(
        size=100,
        capacity=500,
        # n_chromosomes=5,
        genes = [f'gene{n}' for n in range(10)],
        ploidy=2,
        parity=5,
        reproduction_mode=ReproductionMode.SEXUAL,
        reproductive_rate=0.1,
        mortality_rate=mr,
        allele_fitness_range=(0.1, 0.9)
    )
    
    print(f"Initial population size: {len(pop.individuals)}")
    print(f"Initial mean fitness: {np.mean(pop.get_fitness_distribution()):.3f}")
    
    # nindis = list()
    # meanfit = list()
    
    # Run for a few generations
    for _ in range(max_gen):
        pop.timestep()
        # nindis.append(len(pop.individuals))
        # meanfit.append(np.mean([i.fitness for i in pop.individuals]))
        print(f"Generation {pop.generation}: size={len(pop.individuals)}, "
              f"mean fitness={np.mean(pop.get_fitness_distribution()):.3f}")
    
    pop.tracker.plot("age","./testage.png")
    pop.tracker.plot("number_individuals","./testnumindis.png")
    pop.tracker.plot("fitness","./testfitness.png")
    pop.tracker.plot("allele_frequencies","./testallele_freq.png")
    pop.tracker.plot("fitness_vs_frequency","./testfitfreqs.png")
    