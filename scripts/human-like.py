
from popgensim.genetics import GeneType,Gene,Genome
from popgensim import Allele,Chromosome,Individual,Population,ReproductionMode
from popgensim.world.energetics import EnergyCalculator
import numpy as np
import matplotlib.pyplot as plt

genome = Genome()

for i in range(3):
    gene = Gene(f"enzyme{i}", GeneType.PATHWAY, 
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

genome.add_gene(Gene("proteinD", GeneType.PAIRWISE))
genome.add_gene(Gene("proteinE", GeneType.PAIRWISE))
genome.add_interaction("proteinD", "proteinE", "synergistic", strength=2.0)

print(genome)

max_gen = 100

# mr = lambda gen:1+np.sin(gen/20)/2
mr = 0.1

# Create a population
pop = Population(
    size=100,
    genome=genome,
    capacity=200,
    ploidy=2,
    parity=5,
    reproduction_mode=ReproductionMode.SEXUAL,
    reproductive_rate=0.015,
    mortality_rate=mr,
    allele_fitness_range=(0.1, 0.9),
    lifespan = 65
)

print(f"Initial population size: {len(pop.individuals)}")
print(f"Initial mean fitness: {np.mean(pop.get_fitness_distribution()):.3f}")

# Run for a few generations
for _ in range(max_gen):
    pop.timestep()
    print(f"Generation {pop.generation}: size={len(pop.individuals)}, "
            f"mean fitness={np.mean(pop.get_fitness_distribution()):.3f}")
    
pop.tracker.plot("fitness","human-like_fitness")
pop.tracker.plot("age","human-like_age")
pop.tracker.plot("allele_frequencies","human-like_alleles_{gene}")

