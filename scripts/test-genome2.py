
from popgensim import *
from popgensim.genetics.genes import *
from popgensim.population.population import *
import random

genome = Genome()

benefit_range = (0.8,1.0)
cost_range = (0.1,0.3)
gain_range = (0.9,1.1)

genes = [GenericGene(f'gene{n}',
                    benefit_range=benefit_range,
                    cost_range=cost_range) 
        for n in range(7)
]

regulatory = [RegulatoryGene(f'reggene{n}',
                            gain_range=gain_range,
                            cost_range=cost_range) 
            for n in range(3)
]

offgene = OffspringNumberGene()
genome.add_gene(offgene)

for g in genes:
    genome.add_gene(g)

#make random pathways
npaths = 3
path_lens = 3

for np in range(npaths):
    gs = random.sample(genes, path_lens)
    genome.add_pathway(f'pathway{np}', *gs)
    
# make random regulatory network
nregs = 10

for nr in range(nregs):
    reg = random.sample(regulatory, 1)[0]
    tgt = random.sample(genes,1)[0]
    genome.add_regulatory_relationship(reg, tgt)
print(genome)

popconfig = PopulationConfig("species1", 100, lifespan = 40)
pop = Population(popconfig, genome)

pop._initialize_population((0.9,1.1),(0.8,1.0))

econf = EnvironmentConfig(holding_cost = 0.2)
env = Environment(econf)

print(pop)
ngens = 100
for n in range(ngens):
    pop.timestep(env)
    print(pop)
    if len(pop.individuals) < 1:
        break
# print(pop)
if pop.individuals:
    print(pop.individuals[-1].chromosomes)

pop.tracker.plot("energy","pop_energy")
