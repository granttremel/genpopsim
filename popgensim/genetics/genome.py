from typing import List, Dict, Tuple, Optional, Set, Union, Iterable
import numpy as np

from ..core.allele import Allele
from ..core.chromosome import Chromosome
from .genes import Gene, GenericGene, RegulatoryGene


class Genome:
    """Population-level gene interaction network"""
    def __init__(self):
        self.genes: Dict[str, Gene] = {}  # {gene_name: Gene}
        self.pathways: Dict[str, List[str]] = {}  # {pathway_id: [gene_names in order]}
        # self.interactions: Dict[Tuple[str,str],float] = {}
        self.regulators: Dict[str,Set[str]] = {}
        self.regulated: Dict[str,Set[str]] = {}
        
        self.rmax = 5
    
    def add_gene(self, gene: Gene):
        """Add a gene to the genome"""
        self.genes[gene.name] = gene
    
    def add_pathway(self, pathway_id:str, *genes):
        self.pathways[pathway_id] = [str(g) for g in genes]
    
    # def add_interaction(self, gene1_name: str, gene2_name: str, strength: float = 1.0):
    #     """Add pairwise interaction between genes"""
    #     self.interactions[(gene1_name,gene2_name)] = strength
    
    def add_regulatory_relationship(self, regulator:RegulatoryGene, target: Gene|str):
        """Add regulatory relationship"""
        if not regulator in self.genes:
            self.add_gene(regulator)
        if not str(regulator) in self.regulators:
            self.regulators[str(regulator)] = set()
        self.regulators[str(regulator)].add(str(target))
            
    def calculate_expression(self, alleles:Dict[str,Allele]):
        expression = {g:a.expression for g,a in alleles.items()}
        expression = self._apply_regulatory_effects(expression, alleles)
        expression = self._apply_pathway_effects(expression)
        for g, a in alleles.items():
            a.expression = expression[g]
        
        return expression
            
    def _apply_regulatory_effects(self, expression: Dict[str, float], 
                                  alleles: Dict[str, 'Allele']) -> Dict[str, float]:
        """Apply regulatory gene effects"""
        # Calculate in dependency order to handle regulatory cascades
        calculated = set()
        max_iterations = len(expression) * 2  # Prevent infinite loops
        
        for _ in range(max_iterations):
            changed = False
            
            for gene_name, gene in self.genes.items():
                if gene_name in calculated:
                    continue
                regby = self.regulated.get(gene_name,[])
                # Check if all regulators have been calculated
                if all(reg in calculated or reg not in expression 
                       for reg in regby):
                    
                    if regby and gene_name in expression:
                        # Apply regulatory modifications
                        reg_modifier = 1.0
                        
                        for regulator_name in regby:
                            if regulator_name in expression:
                                reg_expr = expression[regulator_name]
                                # Regulatory effect: genetically encoded gain and expression
                                reg_effect = gene.gain * reg_expr
                                reg_modifier *= reg_effect
                        
                        # Apply bounds: 1/r to r where r=10
                        expression[gene_name] *= np.clip(reg_modifier, 1/self.rmax, self.rmax)
                        
                        # Check for silencing
                        # if expression[gene_name] < 0.05:
                        #     expression[gene_name] = 0  # Silenced
                    
                    calculated.add(gene_name)
                    changed = True
            
            if not changed:
                break
        
        return expression
    
    def _apply_pathway_effects(self, expression: Dict[str, float]) -> Dict[str, float]:
        """Apply pathway bottleneck effects"""
        for pathway in self.pathways.values():
            bottleneck_factor = 1.0
            
            for gene_name in pathway:
                if gene_name in expression:
                    # Current gene limited by upstream bottleneck
                    expression[gene_name] *= bottleneck_factor
                    
                    gene_expr = expression[gene_name]
                    if gene_expr == 0:  # Silenced gene blocks pathway
                        bottleneck_factor = 0
                    else:
                        bottleneck_factor = min(bottleneck_factor, gene_expr)
        
        return expression
    
    def _create_random_chromosome(self,expression_range: Tuple[float,float]) -> Chromosome:
        
        chrom = Chromosome()
        for gn,g in self.genes.items():
            a = g.random_allele(expression_range=expression_range)
            chrom.add_allele(a)
        return chrom
    
    def __str__(self):
        
        outstrs = []
        
        allgenes = ','.join(self.genes.keys())
        outstrs.append(f"genes=[{allgenes}]")
        if self.pathways:
            allpways = ','.join(self.pathways.keys())
            outstrs.append(f"pathways=[{allpways}]")
        if self.regulators:
            regstrs = []
            for r in self.regulators:
                regs = self.regulators[r]
                regstr = ','.join(regs)
                if len(regs) > 1:
                    regstr = f'({regstr})'
                regstrs.append(f"{r}->{regstr}")
            regstr = ','.join(regstrs)
            outstrs.append(f'regulatory=[{regstr}]')
        
        outstr = ','.join(outstrs)
        
        return f"Genome({outstr})"       

