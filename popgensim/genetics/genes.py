
from typing import List, Dict, Tuple, Optional, Set, Union, Iterable
from typing import TYPE_CHECKING
import numpy as np
from dataclasses import dataclass, field
from enum import Enum
import random
from abc import ABC, abstractmethod

from .mutations import MutationStrategy
from ..core.allele import Allele

class Gene(ABC):
    """Gene with interaction information"""
    def __init__(self, name: str, **ranges):
        self.name = name
        self.ranges:Dict[str,Tuple[float,float]]={}
        self.allele_counter = 0
    
    @abstractmethod
    def create_allele(self, expression:float) -> Allele:
        """Factory method to create an allele of this gene"""
        pass
    
    @abstractmethod
    def random_allele(self, expression_range: Tuple[float,float]) -> Allele:
        """Create a random allele for initial population"""
        pass
    
    # @abstractmethod
    # def get_mutation_strategy(self) -> MutationStrategy:
    #     """Return appropriate mutation strategy for this gene type"""
    #     pass
    
    @abstractmethod
    def mutate_allele(self, allele:Allele)->Allele:
        
        pass
    
    def _next_allele_id(self) -> int:
        """Generate unique allele ID"""
        self.allele_counter += 1
        return self.allele_counter
    
    def __str__(self):
        return self.name
    
    def __repr__(self):
        
        return f"{type(self).__name__}({self.name})"

class GenericGene(Gene):
    
    params=['benefit','cost']
    def __init__(self, name, benefit_range=(0.1,0.9), cost_range=(0.1,0.9)):
        super().__init__(name)
        self.ranges['benefit'] = benefit_range
        self.ranges['cost'] = cost_range
    
    def create_allele(self, expression:float=1.0, **kwargs) -> Allele:
        """Factory method to create an allele of this gene"""
        benefit = kwargs.get('benefit')
        cost = kwargs.get('cost')
        return Allele.create(self._next_allele_id(),self.name,expression,  benefit=benefit, cost=cost)
    
    def random_allele(self, expression_range: Tuple[float,float],) -> Allele:
        """Create a random allele for initial population"""
        
        expression = random.uniform(*expression_range)
        benefit = random.uniform(*self.ranges["benefit"])
        cost = random.uniform(*self.ranges["cost"])
        
        return self.create_allele(expression, benefit=benefit, cost=cost)
    
    def mutate_allele(self, allele:Allele)->Allele:
        return allele
    
    
    def make_regulatory(self):
        return RegulatoryGene(self.name)
    
    
class RegulatoryGene(Gene):
    
    params=['gain','cost']
    def __init__(self, name, gain_range=(-1,1), cost_range=(0.1,0.9)):
        super().__init__(name)
        self.ranges["gain"]=gain_range
        self.ranges["cost"]=cost_range
        
    def create_allele(self, expression:float=1.0, **kwargs) -> Allele:
        """Factory method to create an allele of this gene"""
        gain=kwargs.get("gain")
        cost=kwargs.get("cost")
        return Allele.create(self._next_allele_id(), self.name, expression, gain=gain, cost=cost)
    
    def random_allele(self, expression_range: Tuple[float,float],) -> Allele:
        """Create a random allele for initial population"""
        
        expression = random.uniform(*expression_range)
        gain = random.uniform(*self.ranges["gain"])
        cost = random.uniform(*self.ranges["cost"])
        return self.create_allele(expression, gain=gain, cost=cost)
    
    def mutate_allele(self, allele:Allele)->Allele:
        return allele
    
class OffspringNumberGene(Gene):
    """Gene controlling how many offspring to produce"""
    def __init__(self, number_range=(1,5)):
        super().__init__("offspring_number")
        self.ranges["number"] = number_range
    
    def create_allele(self, expression=1.0, **kwargs) -> Allele:
        number = int(kwargs.get('number'))
        return Allele.create(self._next_allele_id(), self.name,expression,number=number)
    
    def random_allele(self, expression_range: Tuple[float,float],) -> Allele:
        """Create a random allele for initial population"""
        
        expression = random.uniform(*expression_range)
        number = random.choice(range(*self.ranges["number"]))
        return self.create_allele(expression, number=number)
    
    def mutate_allele(self, allele:Allele)->Allele:
        return allele
    