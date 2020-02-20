import numpy as np
from coffea.processor.dataframe import LazyDataFrame
from pdb import set_trace
from copy import copy

class MaskedLazyDataFrame(LazyDataFrame):
    def   __init__(self):
        self.mask_ = np.array([])
    
    @classmethod
    def from_lazy(self, df):
        ret = MaskedLazyDataFrame()
        ret.__dict__.update(
          copy(df.__dict__)
        )
        ret.mask_ = np.ones(df.size).astype(bool)
        return ret
    
    def __getitem__(self, k):
        if isinstance(k, np.ndarray):
            ret = MaskedLazyDataFrame.from_lazy(self)
            ret.mask_ = self.mask_.copy()
            ret.mask_[ret.mask_] = k
            return ret
        else:
            return super().__getitem__(k)[self.mask_]
    
    def __setitem__(self, k, val):
        if self.mask_.all(): # the mask is not active
            super().__setitem__(k, val)
        else:
            raise ValueError('You cannot add fields to a MaskedLazyDataFrame after masking')
    
    @property
    def size(self):
        return self.mask_.sum()
