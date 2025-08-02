import numpy as np
from typing import List, Dict, Tuple
# from . import Allele, Chromosome, Individual, Population
from .core import Individual, Chromosome, Allele
from .genetics import Genome, Gene, GeneType

class FitnessCalculator:
    """Calculate fitness based on gene network"""
    
    def calculate_fitness(self, individual:Individual, genome: Genome) -> float:
        """Calculate epistatic fitness based on gene interactions"""
        fitness_components = {}
        chromosomes = individual.chromosomes
        # Process each gene
        for gene_name, gene in genome.genes.items():
            alleles = self._get_alleles_for_gene(gene_name, chromosomes)
            if not alleles:
                continue
            
            base_fitness = np.mean([a.fitness for a in alleles])
            
            if gene.gene_type == GeneType.INDEPENDENT:
                fitness_components[gene_name] = base_fitness
            
            elif gene.gene_type == GeneType.PAIRWISE:
                # Calculate with interactions
                interaction_effect = 0
                for partner_name in gene.interaction_partners:
                    partner_alleles = self._get_alleles_for_gene(partner_name, chromosomes)
                    if partner_alleles:
                        partner_fitness = np.mean([a.fitness for a in partner_alleles])
                        
                        # Apply interaction rules
                        if partner_name in gene.synergistic_pairs:
                            if base_fitness > 0.6 and partner_fitness > 0.6:
                                interaction_effect += 0.2
                        elif partner_name in gene.antagonistic_pairs:
                            if base_fitness > 0.7 and partner_fitness > 0.7:
                                interaction_effect -= 0.3
                        
                        # General interaction based on strength
                        strength = gene.interaction_strengths.get(partner_name, 1.0)
                        interaction_effect += (partner_fitness - 0.5) * strength * 0.1
                
                fitness_components[gene_name] = max(0, min(1, base_fitness + interaction_effect))
            
            elif gene.gene_type == GeneType.REGULATORY:
                # Regulators affect their targets
                regulatory_power = base_fitness
                for target_name in gene.regulates:
                    # Store regulatory effect to apply to targets
                    fitness_components[f"{gene_name}_regulates_{target_name}"] = regulatory_power
            
            elif gene.gene_type == GeneType.REGULATED:
                # Base fitness modified by regulators
                regulatory_modifier = 1.0
                for regulator_name in gene.regulated_by:
                    reg_key = f"{regulator_name}_regulates_{gene_name}"
                    if reg_key in fitness_components:
                        # Regulators act as multipliers
                        reg_effect = fitness_components[reg_key]
                        regulatory_modifier *= (0.5 + reg_effect)  # Range 0.5-1.5
                
                fitness_components[gene_name] = base_fitness * regulatory_modifier
            
            else:
                fitness_components[gene_name] = base_fitness
        
        # Calculate pathway fitness (weakest link)
        for pathway_id, gene_names in genome.pathways.items():
            pathway_fitnesses = []
            for gene_name in gene_names:
                if gene_name and gene_name in fitness_components:
                    pathway_fitnesses.append(fitness_components[gene_name])
            
            if pathway_fitnesses:
                # Geometric mean emphasizes weak links
                pathway_fitness = np.power(np.prod(pathway_fitnesses), 1/len(pathway_fitnesses))
                fitness_components[f"pathway_{pathway_id}"] = pathway_fitness
        
        # Calculate module fitness
        for module_id, gene_names in genome.modules.items():
            module_fitnesses = []
            for gene_name in gene_names:
                if gene_name in fitness_components:
                    module_fitnesses.append(fitness_components[gene_name])
            
            if module_fitnesses:
                module_fitness = np.mean(module_fitnesses)
                fitness_components[f"module_{module_id}"] = module_fitness
        
        # Combine all fitness components
        all_fitnesses = [f for f in fitness_components.values() if isinstance(f, (int, float))]
        fitness = np.mean(all_fitnesses) if all_fitnesses else 0.5
        individual._fitness = fitness
        return fitness
    
    def _get_alleles_for_gene(self, gene_name: str, chromosomes: List[Chromosome]) -> List[Allele]:
        """Get all alleles for a gene across chromosomes"""
        alleles = []
        for chromosome in chromosomes:
            allele = chromosome.get_allele(gene_name)
            if allele:
                alleles.append(allele)
        return alleles