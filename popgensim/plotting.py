

from typing import List, Dict, Tuple, Optional, Set, Union
from itertools import chain
import matplotlib.pyplot as plt
import numpy as np

class Tracker:
    def __init__(self, population):
        self.population = population
        self.number_generations: int = 0
        self.number_individuals: List[int] = []
        self.mean_fitnesses = []
        self.allele_frequencies: Dict[Tuple[str,int],List[float]] = {} #indexed by gene id, timestep
        self.ages: List[List[int]] = []
    
    def update(self):
        
        self.number_generations += 1
        self.number_individuals.append(len(self.population.individuals))
        
        counter = {}
        ages=[]
        fitnesses = []
        for ind in self.population.individuals:
            ages.append(ind.age)
            fitnesses.append(ind.fitness)
            for c in ind.chromosomes:
                for g in c.genes:
                    a = c.genes[g]
                    aind = a.id
                    ct = counter.get((g,aind),0)+1
                    counter[(g,aind)] = ct
        
        for g in counter:
            afs = self.allele_frequencies.get(g,[]) 
            afs.append(counter[g])
            self.allele_frequencies[g] = afs
            
        self.ages.append(ages)
        self.mean_fitnesses.append(np.mean(fitnesses))
        
    def plot(self, quantity, outfile, **kwargs):
        
        if quantity == "number_individuals":
            f,ax = plt.subplots()
            ax.plot(self.number_individuals)
            
            ax.set_title("population size over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("number of individuals")
            
            f.savefig(outfile)
        elif quantity == "fitness":
            f,ax = plt.subplots()
            ax.plot(self.mean_fitnesses)
            
            ax.set_title("fitness over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("mean fitness")
            
            f.savefig(outfile)
        elif quantity == "age":
            meanages = [np.mean(age) for age in self.ages]
            f,ax = plt.subplots()
            ax.plot(meanages)
            
            ax.set_title("ages over generation")
            ax.set_xlabel("generation")
            ax.set_ylabel("mean age")
            
            f.savefig(outfile)
        elif quantity == "allele_frequencies":
            
            genes = kwargs.get('genes',[])
                
            afs = {}
            
            for g,a in self.allele_frequencies:
                if genes and g not in genes:
                    continue
                timeseries = self.allele_frequencies[(g,a)]
                if len(timeseries) < 10:
                    continue
                if not g in afs:
                    afs[g] = []
                afs[g].append(timeseries)
                
            for g in afs:
                _outf = outfile.format(gene=g)    
                lg = []
                f,ax = plt.subplots()
                for i,ts in enumerate(afs[g]):
                    nts = len(ts)
                    dom = range(self.number_generations - nts,self.number_generations)
                    ax.plot(dom, ts)
                    lg.append(str(self.population.gene_pool[g][i]))
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
            
            f.savefig(outfile)
        