
from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np

from ..core.allele import Allele

class MutationStrategy(ABC):
    @abstractmethod
    def mutate(self, allele: Allele) -> Allele:
        pass

class ContinuousMutation(MutationStrategy):
    """For continuous traits like metabolic rate"""
    def __init__(self, std_dev: float = 0.1, bounds: Tuple[float, float] = (0, 1)):
        self.std_dev = std_dev
        self.bounds = bounds
    
    def mutate(self, allele: Allele) -> Allele:
        new_value = np.random.normal(allele.value, self.std_dev)
        new_value = np.clip(new_value, *self.bounds)
        return Allele(
            id=allele.id + 10000,  # New ID
            gene_name=allele.gene_name,
            value=new_value,
            parameters=allele.parameters.copy()
        )

class DiscreteMutation(MutationStrategy):
    """For discrete traits like offspring number"""
    def __init__(self, max_change: int = 1, bounds: Tuple[int, int] = (1, 10)):
        self.max_change = max_change
        self.bounds = bounds
    
    def mutate(self, allele: Allele) -> Allele:
        change = np.random.randint(-self.max_change, self.max_change + 1)
        new_value = int(allele.value) + change
        new_value = np.clip(new_value, *self.bounds)
        return Allele(
            id=allele.id + 10000,
            gene_name=allele.gene_name,
            value=float(new_value),
            parameters=allele.parameters.copy()
        )

class RegulatoryMutation(MutationStrategy):
    
    def __init__(self, std_dev: float = 0.1, bounds: Tuple[float, float] = (0, 1)):
        self.std_dev = std_dev
        self.bounds = bounds
    
    def mutate(self, allele: Allele) -> Allele:
        new_value = np.random.normal(allele.value, self.std_dev)
        new_value = np.clip(new_value, *self.bounds)
        return Allele(
            id=allele.id + 10000,  # New ID
            gene_name=allele.gene_name,
            value=new_value,
            parameters=allele.parameters.copy()
        )