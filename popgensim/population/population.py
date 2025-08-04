
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
from typing import TYPE_CHECKING
import random
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from itertools import chain
import os
from dataclasses import dataclass
from enum import Enum

from ..core.allele import Allele
from ..core.chromosome import Chromosome
from .individual import Individual
from ..genetics.genes import Gene
from ..genetics.genome import Genome
from ..world.energetics import EnergyCalculator
from ..world.environment import Environment

# if TYPE_CHECKING:
#     from ..world.environment import Environment
# from .energetics import EnergyCalculator

class ReproductionMode(Enum):
    ASEXUAL = "asexual"
    SEXUAL = "sexual"

@dataclass
class PopulationConfig:
    """Configuration for a population/species"""
    species_name: str
    initial_size: int
    
    # Genetic parameters
    # n_genes: int = 10
    ploidy: int = 2
    mutation_rate: float = 0.001
    
    # Reproduction parameters
    reproduction_mode: ReproductionMode = ReproductionMode.SEXUAL
    parity: float = 2.0
    reproductive_rate: float = 0.02
    
    
    # Metabolic parameters
    base_metabolic_rate: float = 0.1
    reproduction_cost: float = 0.5
    lifespan: int = 20
    
    # Niche parameters TO BE IMPLEMENTED
    # optimal_temperature: float = 288
    # temperature_tolerance: float = 20
    # preferred_energy_source: str = "photosynthesis"

class Population:
    """Population of individuals with configurable genetics"""
    def __init__(self, config:PopulationConfig, genome:Genome):
        
        self.size = config.initial_size
        self.genome = genome
        self.ploidy = config.ploidy
        self.gene_pool:Dict[str,List[Allele]] = {}
        
        self.reproduction_mode = config.reproduction_mode
        self.parity = config.parity
        self.reproductive_rate = config.reproductive_rate
        self.reproductive_fitness_strength = 0.1
        
        self.individuals: List[Individual] = []
            
        self.generation = 0
        self.n_ever = 0
        self.allele_counter = 0
        self.lifespan = config.lifespan
        
        self.energy_calc = EnergyCalculator()
        self.tracker = Tracker(self)
        
        
    
    def _initialize_population(self, expression_range: Tuple[float, float], energy_range:Tuple[float,float]):
        """Create initial population with random alleles"""
        
        cs = [self.genome._create_random_chromosome(expression_range=expression_range) for _ in range(self.ploidy)]
        
        for g in self.genome.genes:
            if not g in self.gene_pool:
                self.gene_pool[g] = []
            for c in cs:
                self.gene_pool[g].append(c.alleles[g])
                
        for ni in range(self.size):
            chromosomes = []
            name = f"indi{ni}"
            nrg = random.uniform(*energy_range)
            for _ in range(self.ploidy):
                newchr = Chromosome()
                for g in self.genome.genes:
                    newchr.add_allele(self.gene_pool[g][random.choice(range(self.ploidy))])
                chromosomes.append(newchr)
            indi = Individual(name, chromosomes, self.ploidy, energy=nrg)
            self.register_new_individual(indi)

    def reproduce_single(self, parent1: Individual, parent2: Optional[Individual] = None) -> Individual:
        p1_energy = parent1.offspring_ratio*parent1.energy/parent1.offspring_number
        """Create offspring from one or two parents"""
        if self.reproduction_mode == ReproductionMode.ASEXUAL:
            # Simple cloning with potential for mutation (to be added later)
            offspring_chromosomes = [c.copy() for c in parent1.chromosomes]
            newname = f"Indi{self.n_ever}"
            newind = Individual(newname, offspring_chromosomes, self.ploidy)
            self.register_new_individual(newind)
            return newind
        
        elif self.reproduction_mode == ReproductionMode.SEXUAL:
            if parent2 is None:
                raise ValueError("Sexual reproduction requires two parents")
            p2_energy = parent2.offspring_ratio*parent2.energy/parent2.offspring_number
            
            # Get gametes from each parent
            gamete1 = parent1.get_gamete(p1_energy)
            gamete2 = parent2.get_gamete(p2_energy)
            
            if not gamete1 or not gamete2:
                print('Parents {parent1} and {parent2} failed to reproduce!')
            
            # Combine gametes
            offspring_chromosomes = gamete1 + gamete2
            newname = f"Indi{self.n_ever}"
            newind = Individual(newname, offspring_chromosomes, self.ploidy, p1_energy+p2_energy)
            self.register_new_individual(newind)
            return newind
        
        else:
            return None
    
    def reproduce(self, parent1:Individual, parent2: Optional[Individual] = None) -> Iterable[Individual]:
        n_offspring = (parent1.offspring_number + parent2.offspring_number)//2
        offs = list()
        for n in range(n_offspring):
            offs.append(self.reproduce_single(parent1, parent2 = parent2))
        return offs
    
    def register_new_individual(self, new_indi:Individual):
        for c in new_indi.chromosomes:
            for g,a in c.alleles.items():
                if c.alleles[g].id < 0:
                    #new allele!
                    newid = len(self.gene_pool[g])
                    a.id = newid
                    # a.discovered = self.generation
                    self.gene_pool[g].append(a)
        # self.fitcalc.calculate_reward(new_indi, self.genome)
        self.individuals.append(new_indi)
        self.n_ever+=1
    
    def timestep(self, environment:Environment):
    # def timestep(self):
        
        """Execute one generation: reproduction and mortality"""
        self.generation += 1
        
        costs = dict()
        for indi in self.individuals:
            gain, cost = self.energy_calc.calculate_energy_change(indi, self.genome, environment)
            indi.energy += gain
            costs[indi.name] = cost
        
        n_repros = int(len(self.individuals) * self.reproductive_rate)
        for nr in range(n_repros):
            # Reproduction
            if self.reproduction_mode == ReproductionMode.SEXUAL:
                # Random mating
                if len(self.individuals) >= 2:
                    parent1, parent2 = random.sample(self.individuals, 2)
                    # pair_fit = parent1.energy + parent2.energy
                    offspring = self.reproduce(parent1, parent2)
                    self.individuals.extend(offspring)
            else:
                # Asexual reproduction
                if self.individuals:
                    parent = random.choice(self.individuals)
                    offspring = self.reproduce(parent)
                    self.individuals.extend(offspring)
        
        survivors = []
        for indi in self.individuals:
            indi.energy -= costs.get(indi.name,0.0)
            if indi.age > self.lifespan:
                # print(f"{indi.name} perished of old age!")
                continue
            if indi.energy < 0:
                print(f"{indi.name} perished of no energy!")
                # print(f"{indi.energy}, {gain}, {cost}")
                continue
            indi.age+=1
            survivors.append(indi)
        self.individuals = survivors
        self.tracker.update()
        
        
    def get_allele_frequencies(self) -> dict:
        """Calculate frequency of each allele in the population"""
        allele_counts = {}
        total_alleles = 0
        
        for individual in self.individuals:
            for chromosome in individual.chromosomes:
                for allele in chromosome.alleles.values():
                    allele_counts[allele.id] = allele_counts.get(allele.id, 0) + 1
                    total_alleles += 1
        
        if total_alleles == 0:
            return {}
        
        return {aid: count/total_alleles for aid, count in allele_counts.items()}
    
    def get_fitness_distribution(self) -> List[float]:
        """Get list of all individual fitness values"""
        return [ind.energy for ind in self.individuals]
    
    def __repr__(self):
        
        outstrs = []
        outstrs.append(f'generation={self.generation}')
        outstrs.append(f'size={len(self.individuals)}')
        outstr = ','.join(outstrs)
        
        return f"Population({outstr})"

class Tracker:
    def __init__(self, population:Population, outdir:str = ""):
        self.population = population
        self.number_generations: int = 0
        self.number_individuals: List[int] = []
        self.sum_energies = []
        self.allele_frequencies: Dict[str,Dict[int,List[float]]] = {} #indexed by gene id, timestep
        self.ages: List[List[int]] = []
        if not outdir:
            self.outdir = "./out"
    
    def update(self):
        indis = self.population.individuals
        self.number_generations += 1
        self.number_individuals.append(len(indis))
        
        counter = {}
        ages=[]
        energies = []
        for ind in indis:
            ages.append(ind.age)
            energies.append(ind.energy)
            for c in ind.chromosomes:
                for g,a in c.alleles.items():
                    if not g in counter:
                        counter[g] = {}
                    if not a.id in counter[g]:
                        counter[g][a.id] = 0
                    counter[g][a.id] += 1
        
        for g in counter:
            for a in counter[g]:
                if not g in self.allele_frequencies:
                    self.allele_frequencies[g] = {}
                if not a in self.allele_frequencies[g]:
                    self.allele_frequencies[g][a] = []
                self.allele_frequencies[g][a].append(counter[g][a])
                    
        self.ages.append(ages)
        self.sum_energies.append(sum(energies))
    
    def plot(self, quantity:str, outfile:str, **kwargs):
        
        if not outfile.endswith(".png"):
            outfile = outfile + ".png"
        outpath = os.path.join(self.outdir, outfile)
        
        if quantity == "number_individuals":
            f,ax = plt.subplots()
            ax.plot(self.number_individuals)
            
            ax.set_title("population size over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("number of individuals")
            f.savefig(outpath)
            
        elif quantity == "energy":
            f,ax = plt.subplots()
            ax.plot(self.sum_energies)
            
            ax.set_title("energy over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("sum of energy")
            f.savefig(outpath)
            
        elif quantity == "age":
            meanages = [np.mean(age) for age in self.ages]
            f,ax = plt.subplots()
            ax.plot(meanages)
            
            ax.set_title("ages over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("mean age")
            f.savefig(outpath)
            
        elif quantity == "allele_frequencies":
            
            genes = kwargs.get('genes',[])
            # rel = kwargs.get('relative',False)
                
            for g in self.allele_frequencies:
                if genes and g not in genes:
                    continue
                
                # if rel:
                    # scale = [sum([self.allele_frequencies[g][a][n] for n in ])]
                
                f,ax = plt.subplots()
                lg = []
                _outf = outpath.format(gene=g)    
                for a in self.allele_frequencies[g]:
                    allele = self.population.gene_pool[g][a]
                    timeseries = self.allele_frequencies[g][a]
                    if len(timeseries) < 10:
                        continue
                    if max(timeseries) < 10:
                        continue
                    
                    dom = range(0,len(timeseries))
                    ax.plot(dom, timeseries)
                    lg.append(str(self.population.gene_pool[g][a]))
                ax.legend(lg)
                f.savefig(_outf)
                
        elif quantity == "fitness_vs_frequency":
            
            afits = [a.benefit for a in chain(*self.population.gene_pool.values())]
            afreqs = [a._frequency for a in chain(*self.population.gene_pool.values())]

            f,ax = plt.subplots()
            ax.scatter(afits, afreqs)
            ax.set_title("fitness vs frequency among alleles")
            ax.set_xlabel("fitness")
            ax.set_ylabel("frequency")
            f.savefig(outpath)
        
if __name__ == "__main__":
    
    max_gen = 200
    
    # mr = lambda gen:1+np.sin(gen/20)/2
    mr = 0.9
    
    genome = Genome()
    for n in range(10):
        genome.add_gene(Gene(f"gene{n}"))
    
    popconfig = PopulationConfig(
        species_name="species1",
        initial_size=100,
        ploidy=2,
    )
    
    # Create a population
    pop = Population(popconfig, genome)
    
    print(f"Initial population size: {len(pop.individuals)}")
    print(f"Initial mean fitness: {np.mean(pop.get_fitness_distribution()):.3f}")

    # env = Environment()

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
    



