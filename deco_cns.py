'''
Created on 13 Aug 2016

@author: vke
'''

import sys
import copy

class deco_cns_ctx:

    def __init__(self):
        self.load = 0
    
    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __str__(self):
        return "CNS %u%%" % (round(self.load * 100))
    
    def update_constant(self, time, p_amb):
        pass
    
    def update_asc_desc(self, time, p_begin, p_end):
        pass
    
    def decoinfo_update(self, p_amb, di):
        
        return di
