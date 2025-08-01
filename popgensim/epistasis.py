import numpy as np
from typing import List, Dict, Tuple
from . import Allele, Chromosome, Individual, Population

class EpistaticFitnessCalculator:
    """Calculate fitness with various types of epistatic interactions"""
    
    def __init__(self, n_genes: int, epistasis_type: str = 'pairwise'):
        self.n_genes = n_genes
        self.epistasis_type = epistasis_type
        
        # For pairwise epistasis: which genes interact
        self.interaction_pairs = self._generate_interaction_pairs()
        
        # For pathway epistasis: gene order in pathway
        self.pathway_order = list(range(n_genes))
        np.random.shuffle(self.pathway_order)
        
        # For modular epistasis: assign genes to modules
        self.gene_modules = self._assign_modules()
    
    def _generate_interaction_pairs(self) -> List[Tuple[int, int]]:
        """Generate random pairs of interacting genes"""
        pairs = []
        # Each gene interacts with 1-3 others
        for i in range(self.n_genes):
            n_interactions = np.random.randint(1, 4)
            partners = np.random.choice(
                [j for j in range(self.n_genes) if j != i],
                size=min(n_interactions, self.n_genes - 1),
                replace=False
            )
            for j in partners:
                if (i, j) not in pairs and (j, i) not in pairs:
                    pairs.append((i, j))
        return pairs
    
    def _assign_modules(self, n_modules: int = 3) -> Dict[int, int]:
        """Assign each gene to a functional module"""
        assignments = {}
        for i in range(self.n_genes):
            assignments[i] = i % n_modules
        return assignments
    
    def calculate_fitness(self, individual) -> float:
        """Calculate fitness with epistatic interactions"""
        # First, get all alleles organized by gene position
        alleles_by_gene = self._organize_alleles_by_gene(individual)
        
        if self.epistasis_type == 'additive':
            # No epistasis - simple mean
            return np.mean([a.fitness for alleles in alleles_by_gene.values() 
                           for a in alleles])
        
        elif self.epistasis_type == 'pairwise':
            # Fitness modified by gene-gene interactions
            base_fitness = []
            interaction_effects = []
            
            # Base fitness from individual genes
            for gene_id, alleles in alleles_by_gene.items():
                base_fitness.extend([a.fitness for a in alleles])
            
            # Interaction effects
            for i, j in self.interaction_pairs:
                if i in alleles_by_gene and j in alleles_by_gene:
                    # Synergistic: high + high = bonus
                    # Antagonistic: high + high = penalty
                    fitness_i = np.mean([a.fitness for a in alleles_by_gene[i]])
                    fitness_j = np.mean([a.fitness for a in alleles_by_gene[j]])
                    
                    # Synergistic epistasis example
                    if fitness_i > 0.7 and fitness_j > 0.7:
                        interaction_effects.append(0.2)  # Bonus
                    # Antagonistic epistasis example
                    elif fitness_i > 0.8 and fitness_j < 0.3:
                        interaction_effects.append(-0.3)  # Incompatibility
                    else:
                        interaction_effects.append(0)
            
            base = np.mean(base_fitness)
            interactions = np.mean(interaction_effects) if interaction_effects else 0
            return max(0, min(1, base + interactions))
        
        elif self.epistasis_type == 'pathway':
            # Sequential pathway - weakest link matters most
            pathway_fitnesses = []
            
            for gene_id in self.pathway_order:
                if gene_id in alleles_by_gene:
                    gene_fitness = np.mean([a.fitness for a in alleles_by_gene[gene_id]])
                    pathway_fitnesses.append(gene_fitness)
            
            if not pathway_fitnesses:
                return 0.5
            
            # Fitness limited by weakest step, but not entirely
            # Geometric mean gives more weight to weak links than arithmetic mean
            return np.power(np.prod(pathway_fitnesses), 1/len(pathway_fitnesses))
        
        elif self.epistasis_type == 'modular':
            # Different modules contribute independently
            module_fitnesses = {}
            
            for gene_id, alleles in alleles_by_gene.items():
                module = self.gene_modules.get(gene_id, 0)
                if module not in module_fitnesses:
                    module_fitnesses[module] = []
                module_fitnesses[module].extend([a.fitness for a in alleles])
            
            # Each module's fitness is the mean of its genes
            # Overall fitness is product (or mean) of module fitnesses
            module_scores = []
            for module, fitnesses in module_fitnesses.items():
                module_scores.append(np.mean(fitnesses))
            
            # Multiplicative across modules
            return np.power(np.prod(module_scores), 1/len(module_scores))
        
        elif self.epistasis_type == 'regulatory':
            # Some genes act as regulators that modify others
            # First 20% of genes are regulators
            n_regulators = max(1, self.n_genes // 5)
            
            regulator_fitness = []
            regulated_fitness = []
            
            for gene_id, alleles in alleles_by_gene.items():
                gene_fitness = np.mean([a.fitness for a in alleles])
                if gene_id < n_regulators:
                    regulator_fitness.append(gene_fitness)
                else:
                    regulated_fitness.append(gene_fitness)
            
            if not regulator_fitness or not regulated_fitness:
                return 0.5
            
            # Regulators act as multipliers on regulated genes
            regulator_effect = np.mean(regulator_fitness)
            base_fitness = np.mean(regulated_fitness)
            
            # Regulator effect: good regulators (>0.5) enhance, bad ones suppress
            modifier = 0.5 + regulator_effect  # Range 0.5-1.5
            return max(0, min(1, base_fitness * modifier))
    
    # def _organize_alleles_by_gene(self, individual) -> Dict[int, List]:
    #     """Organize alleles by their gene position"""
    #     # Assuming genes are at same position across chromosomes
    #     alleles_by_gene = {}
        
    #     for chrom_idx, chromosome in enumerate(individual.chromosomes):
    #         for allele_idx, allele in enumerate(chromosome.alleles):
    #             gene_id = allele_idx  # Gene position
    #             if gene_id not in alleles_by_gene:
    #                 alleles_by_gene[gene_id] = []
    #             alleles_by_gene[gene_id].append(allele)
        
    #     return alleles_by_gene

# Modified Individual class to use epistatic fitness
# class EpistaticIndividual(Individual):
#     def __init__(self, chromosomes: List[Chromosome], ploidy: int = 2, 
#                  fitness_calculator=None):
#         super().__init__(chromosomes, ploidy)
#         self.fitness_calculator = fitness_calculator
    
#     @property
#     def fitness(self) -> float:
#         """Calculate fitness using epistatic interactions"""
#         if self._fitness is None:
#             if self.fitness_calculator:
#                 self._fitness = self.fitness_calculator.calculate_fitness(self)
#             else:
#                 # Fallback to additive
#                 all_fitness_values = []
#                 for chromosome in self.chromosomes:
#                     all_fitness_values.extend([a.fitness for a in chromosome.alleles])
#                 self._fitness = np.mean(all_fitness_values) if all_fitness_values else 0.0
#         return self._fitness