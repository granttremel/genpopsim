
from typing import List, Dict, Tuple, Optional, Set, Union
from dataclasses import dataclass, field
from enum import Enum

class GeneType(Enum):
    INDEPENDENT = "independent"          # No interactions
    PAIRWISE = "pairwise"          # Interacts with specific partners
    PATHWAY = "pathway"            # Part of sequential pathway
    MODULAR = "modular"            # Part of functional module
    REGULATORY = "regulatory"      # Regulates other genes
    REGULATED = "regulated"        # Is regulated by other genes

@dataclass
class Gene:
    """Gene with interaction information"""
    name: str
    gene_type: GeneType
    
    # Interaction partners (for pairwise)
    interaction_partners: Set[str] = field(default_factory=set)
    
    # Pathway information
    pathway_id: Optional[str] = None
    pathway_position: Optional[int] = None
    
    # Module information
    module_id: Optional[str] = None
    
    # Regulatory relationships
    regulates: Set[str] = field(default_factory=set)      # Genes this regulates
    regulated_by: Set[str] = field(default_factory=set)   # Genes that regulate this
    
    # Interaction strength/type parameters
    interaction_strengths: Dict[str, float] = field(default_factory=dict)
    synergistic_pairs: Set[str] = field(default_factory=set)
    antagonistic_pairs: Set[str] = field(default_factory=set)
        

class Genome:
    """Population-level gene interaction network"""
    def __init__(self):
        self.genes: Dict[str, Gene] = {}  # {gene_name: Gene}
        self.pathways: Dict[str, List[str]] = {}  # {pathway_id: [gene_names in order]}
        self.modules: Dict[str, Set[str]] = {}    # {module_id: {gene_names}}
    
    def add_gene(self, gene: Gene):
        """Add a gene to the genome"""
        self.genes[gene.name] = gene
        
        # Update pathway registry
        if gene.pathway_id:
            if gene.pathway_id not in self.pathways:
                self.pathways[gene.pathway_id] = []
            # Insert at correct position or append
            if gene.pathway_position is not None:
                # Ensure list is long enough
                while len(self.pathways[gene.pathway_id]) <= gene.pathway_position:
                    self.pathways[gene.pathway_id].append(None)
                self.pathways[gene.pathway_id][gene.pathway_position] = gene.name
            else:
                self.pathways[gene.pathway_id].append(gene.name)
        
        # Update module registry
        if gene.module_id:
            if gene.module_id not in self.modules:
                self.modules[gene.module_id] = set()
            self.modules[gene.module_id].add(gene.name)
    
    def add_interaction(self, gene1_name: str, gene2_name: str, 
                       interaction_type: str = "neutral", strength: float = 1.0):
        """Add pairwise interaction between genes"""
        if gene1_name in self.genes and gene2_name in self.genes:
            gene1 = self.genes[gene1_name]
            gene2 = self.genes[gene2_name]
            
            gene1.interaction_partners.add(gene2_name)
            gene2.interaction_partners.add(gene1_name)
            
            gene1.interaction_strengths[gene2_name] = strength
            gene2.interaction_strengths[gene1_name] = strength
            
            if interaction_type == "synergistic":
                gene1.synergistic_pairs.add(gene2_name)
                gene2.synergistic_pairs.add(gene1_name)
            elif interaction_type == "antagonistic":
                gene1.antagonistic_pairs.add(gene2_name)
                gene2.antagonistic_pairs.add(gene1_name)
    
    def add_regulatory_relationship(self, regulator_name: str, target_name: str):
        """Add regulatory relationship"""
        if regulator_name in self.genes and target_name in self.genes:
            self.genes[regulator_name].regulates.add(target_name)
            self.genes[target_name].regulated_by.add(regulator_name)
            
    def __str__(self):
        
        outstrs = []
        
        allgenes = ','.join(self.genes.keys())
        outstrs.append(f"genes=[{allgenes}]")
        if self.pathways:
            allpways = ','.join(self.pathways.keys())
            outstrs.append(f"pathways=[{allpways}]")
        if self.modules:
            mods = ','.join(self.modules)
            outstrs.append(f"modules=[{mods}]")
        
        outstr = ','.join(outstrs)
        
        return f"Genome({outstr})"       
