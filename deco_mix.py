'''
Created on 13 Aug 2016

@author: vke
'''

import sys
import copy

class deco_mix:
    #float f_o2;
    #float f_he;
    #bool enabled;
    #float switch_depth;

    def __init__(
        self,
        f_o2 = 0.21,
        f_he = 0.00,
        enabled = False,
        switch_depth = 0):
        self.f_o2 = f_o2
        self.f_he = f_he
        self.enabled = enabled
        self.switch_depth = switch_depth
        self.volume = 0
    
    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __str__(self):
        if (self.f_he == 0.0):
            if (self.f_o2 == 0.21):
                return "AIR"
            else:
                return "NX%u" % (round(self.f_o2 * 100))
        else:
            return "TX%u/%u" % (round(self.f_o2 * 100), round(self.f_he * 100))

    def get_enabled(self):
        return self.enabled

    def get_fhe(self):
        return self.f_he

    def get_fn2(self):
        return 1.0 - self.f_he - self.f_o2

    def get_fo2(self):
        return self.f_o2

    def get_switch_depth(self):
        return self.switch_depth

    def set_enabled(self, enabled):
        self.enabled = enabled

    def set_fhe(self, fhe):
        self.f_he = fhe

    def set_fo2(self, fo2):
        self.f_o2 = fo2

    def set_switch_depth(self, switch_depth):
        self.switch_depth = switch_depth

class deco_mix_ctx:
    #constants
    DECO_MIX_AIR = 0
    DECO_NUM_MIXES = 12
    mixes = [ ]

    def __init__(self):
        for imix in range(0, deco_mix_ctx.DECO_NUM_MIXES):
            self.mixes.append(deco_mix())
        self.first_mix = deco_mix_ctx.DECO_MIX_AIR
        self.active_mix = deco_mix_ctx.DECO_MIX_AIR

    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        newone.mixes = copy.deepcopy(self.mixes, memo)
        return newone

    def __str__(self):
        result = ""
        for imix in range(0, deco_mix_ctx.DECO_NUM_MIXES):
            result += "Mix %u: %s @ %u m %u l\n" % (imix, self.get_mix(imix), self.get_mix(imix).get_switch_depth(), self.get_mix(imix).volume)
        return result

    def get_active_mix(self):
        return self.active_mix

    def get_best_mix(self, depth):
        best_mix = self.get_active_mix()
        best_depth = sys.maxsize
        for i in range(0, deco_mix_ctx.DECO_NUM_MIXES):
            mix = self.get_mix(i)
            if (mix.get_enabled() == True and
                mix.get_switch_depth() >= depth and
                mix.get_switch_depth() < best_depth):
                best_depth = mix.get_switch_depth()
                best_mix = i
        return best_mix

    def get_first_mix(self):
        return self.first_mix

    def get_mix(self, imix):
        return self.mixes[imix]

    def set_active_mix(self, mix):
        self.active_mix = mix

    def set_first_mix(self, mix):
        self.first_mix = mix

    def set_mix(self, imix, mix):
        self.mixes[imix] = mix

    def switch_mix(self, imix):
        self.active_mix = imix
    
    def clear_mix_volumes(self):
        for i in range(0, deco_mix_ctx.DECO_NUM_MIXES):
            self.mixes[i].volume = 0

