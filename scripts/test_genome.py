
from popgensim.genetics import GeneType,Gene,Genome
from popgensim import Allele,Chromosome,Individual,Population,ReproductionMode
from popgensim.epistasis import EpistaticFitnessCalculator
import numpy as np
import matplotlib.pyplot as plt

genome = Genome()

for i in range(3):
    gene = Gene(f"enzyme_{i}", GeneType.PATHWAY, 
                pathway_id="glycolysis", pathway_position=i)
    genome.add_gene(gene)



# Add regulatory gene
reg_gene = Gene("glycolysis_regulator", GeneType.REGULATORY)
genome.add_gene(reg_gene)
for i in range(3):
    genome.add_regulatory_relationship("glycolysis_regulator", f"enzyme{i}")


genome.add_gene(Gene("proteinA", GeneType.PAIRWISE))
genome.add_gene(Gene("proteinB", GeneType.PAIRWISE))
genome.add_interaction("proteinA", "proteinB", "synergistic", strength=2.0)

genome.add_gene(Gene("proteinC",GeneType.INDEPENDENT))

print(genome)


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

nindis = list()
meanfit = list()

# Run for a few generations
for _ in range(max_gen):
    pop.timestep()
    nindis.append(len(pop.individuals))
    meanfit.append(np.mean([i.fitness for i in pop.individuals]))
    print(f"Generation {pop.generation}: size={len(pop.individuals)}, "
            f"mean fitness={np.mean(pop.get_fitness_distribution()):.3f}")
    
f,ax = plt.subplots()
ax.plot(nindis)
f.savefig('./pop.png')

f,ax = plt.subplots()
ax.plot(meanfit)
f.savefig('./meanfit.png')

meanages = [np.mean(age) for age in pop.tracker.ages]
f,ax = plt.subplots()
ax.plot(meanages)
f.savefig('./meanages.png')
