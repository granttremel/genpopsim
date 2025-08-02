
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
import random
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from itertools import chain
import os

from .core import Individual, Chromosome, Allele, ReproductionMode
from .genetics import Genome, Gene, GeneType
from .epistasis import FitnessCalculator


class Population:
    """Population of individuals with configurable genetics"""
    def __init__(self, 
            size: int,
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
        self.genome = genome
        self.ploidy = ploidy
        self.gene_pool:Dict[str,List[Allele]] = {}
        
        self.fitcalc = FitnessCalculator()
        
        self.reproduction_mode = reproduction_mode
        self.parity = parity
        self.reproductive_rate = reproductive_rate
        self.reproductive_fitness_strength = 0.1
        self._mortality_rate = mortality_rate
        
        self.individuals: List[Individual] = []
        # self.perished: List[Individual] = []
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
        na = len(self.gene_pool[gene_name])
        fitness = np.random.uniform(*fitness_range)
        new_allele = Allele(id=na,gene_name=gene_name,fitness=fitness,discovered=self.generation)
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
            indi = Individual(name, chromosomes, self.ploidy)
            self.register_new_individual(indi)
            # self.individuals.append()

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
            offspring_chromosomes = gamete1 + gamete2
            newname = f"Indi{self.n_ever}"
            self.n_ever += 1
            newind = Individual(newname, offspring_chromosomes, self.ploidy)
            self.register_new_individual(newind)
            return newind
        
        else:
            return None
    
    def reproduce(self, n_offspring, parent1:Individual, parent2: Optional[Individual] = None) -> Iterable[Individual]:
        offs = list()
        for n in range(n_offspring):
            offs.append(self.reproduce_single(parent1, parent2 = parent2))
        # print(f'{parent1} and {parent2} produced {n_offspring} offspring!')
        return offs
    
    def register_new_individual(self, new_indi:Individual):
        for c in new_indi.chromosomes:
            for g,a in c.alleles.items():
                if c.alleles[g].id < 0:
                    #new allele!
                    newid = len(self.gene_pool[g])
                    a.id = newid
                    a.discovered = self.generation
                    self.gene_pool[g].append(a)
        self.fitcalc.calculate_fitness(new_indi, self.genome)
        self.individuals.append(new_indi)
    
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
                    self.individuals.extend(offspring)
        
        # Mortality (Poisson process approximation)
        survivors = []
        for individual in self.individuals:
            if individual.age > self.lifespan:
                continue
            death_prob = self.mortality_rate(individual)  # Scale factor
            if random.random() > death_prob:
                individual.age+=1
                survivors.append(individual)
            else:
                pass
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
        return [ind.fitness for ind in self.individuals]

class Tracker:
    def __init__(self, population:Population, outdir:str = ""):
        self.population = population
        self.number_generations: int = 0
        self.number_individuals: List[int] = []
        self.mean_fitnesses = []
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
        fitnesses = []
        for ind in indis:
            ages.append(ind.age)
            fitnesses.append(ind.fitness)
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
        self.mean_fitnesses.append(np.mean(fitnesses))
    
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
            
        elif quantity == "fitness":
            f,ax = plt.subplots()
            ax.plot(self.mean_fitnesses)
            
            ax.set_title("fitness over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("mean fitness")
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
                    
                    dom = range(self.population.gene_pool[g][a].discovered,self.population.gene_pool[g][a].discovered+len(timeseries))
                    ax.plot(dom, timeseries)
                    lg.append(str(self.population.gene_pool[g][a]))
                ax.legend(lg)
                f.savefig(_outf)
                
        elif quantity == "fitness_vs_frequency":
            
            afits = [a.fitness for a in chain(*self.population.gene_pool.values())]
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
    
    # Create a population
    pop = Population(
        size=100,
        capacity=500,
        genome = genome,
        ploidy=2,
        parity=5,
        reproduction_mode=ReproductionMode.SEXUAL,
        reproductive_rate=0.1,
        mortality_rate=mr,
        allele_fitness_range=(0.1, 0.9)
    )
    
    print(f"Initial population size: {len(pop.individuals)}")
    print(f"Initial mean fitness: {np.mean(pop.get_fitness_distribution()):.3f}")

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
    



