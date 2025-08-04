
from dataclasses import dataclass, field
from typing import List, Dict, Tuple,Iterable, Optional, Set, Union
from types import SimpleNamespace

@dataclass
class Allele:
    
    id: int
    gene_name: str  
    expression: float=1.0
    params: SimpleNamespace = field(default_factory=SimpleNamespace)
    _frequency:int = 0
    
    @classmethod
    def create(cls, id: int, gene_name: str, expression: float=1.0, **kwargs):
        """Factory method"""
        return cls(id, gene_name, expression, SimpleNamespace(**kwargs))
    
    def __str__(self):
        # return f"Allele(id={self.id},benefit={self.benefit:0.3f},cost={self.cost:0.3f},{self._frequency}x)"
        outstrs = []
        outstrs.append(f"id={self.id}")
        for p,v in self.params.__dict__.items():
            outstrs.append(f"{p}={v:0.3g}")
        outstr = ','.join(outstrs)
        return f"Allele({outstr})"

    def __repr__(self):
        return str(self)
    
    def __hash__(self):
        return hash((self.gene_name, self.id))

