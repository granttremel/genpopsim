

from typing import List, Dict, Tuple, Optional, Set, Union, Iterable
from.allele import Allele

class Chromosome:
    """Container for alleles"""
    def __init__(self):
        self.alleles: Dict[str,Allele] = {}
        self.n_genes = len(self.alleles)
    
    def copy(self):
        """Create a copy of the chromosome"""
        return Chromosome.from_alleles(self.alleles.values())
    
    @classmethod
    def from_alleles(cls, alleles:Iterable[Allele]):
        chr = cls()
        
        for a in alleles:
            chr.add_allele(a)
            
        return chr
        
    def add_allele(self, allele:Allele):
        self.alleles[allele.gene_name] = allele
        allele._frequency += 1
    
    def get_allele(self, gene_name:str):
        
        return self.alleles.get(gene_name)
    
    def done(self):
        b = all([self.alleles[g] for g in self.alleles])
        return b
    
    def __len__(self):
        return len(self.alleles)
    
    def __str__(self):
        content = []
        for gene in self.alleles:
            content.append(f"{gene}:{self.alleles[gene]}")
        contentstr = "\n".join(content)
        return f"Chromosome({contentstr})"
    
    def __repr__(self):
        return str(self)
