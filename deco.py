'''
Created on 8 Jan 2015

@author: vke
'''

import numpy
import copy
import deco_mix

class deco_settings:
    def __init__(
                 self,
                 gf_low = 0.5,
                 gf_high = 0.8,
                 nofly_pressure = 0.6,
                 ascent_rate = 10,
                 stop_size = 3,
                 last_stop = 3):
        self.gf_low = gf_low
        self.gf_high = gf_high
        self.nofly_pressure = nofly_pressure
        self.ascent_rate = ascent_rate
        self.stop_size = stop_size
        self.last_stop = last_stop

    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __str__(self):
        return "gf: %u/%u\nnofly_pressure: %f\nascent_rate: %u\nstop_size: %u\nlast_stop: %u\n" % (
                                round(self.gf_low * 100),
                                round(self.gf_high * 100),
                                self.nofly_pressure,
                                self.ascent_rate,
                                self.stop_size,
                                self.last_stop)

class deco_info:
    N_STOPS = 200
    def __init__(self, 
            dtype = numpy.float32):
        self.dtype = dtype
        self.stops = numpy.zeros(deco_info.N_STOPS, dtype=[('mix', int), ('time', dtype), ('depth', dtype)])
        self.lead_tissue = 0
        self.lead_tissue_gf = 0
        self.allowed_gf = 0.0
        self.gf_low_depth = 0.0
        self.nostop_time = 0
        self.time_to_surface = 0
        self.cns_load = 0.0
        self.otu_load = 0.0
        self.mix_ctx = deco_mix.deco_mix_ctx()

    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        newone.stops = copy.deepcopy(self.stops, memo)
        newone.mix_ctx = copy.deepcopy(self.mix_ctx, memo)
        return newone

    def __str__(self):
        result = ""
        for stop in self.stops:
            if (stop['depth'] == 0):
                break
            result += "%s @ %u m for %u min.\n" % (
                                self.mix_ctx.get_mix(stop['mix']),
                                stop['depth'],
                                stop['time'])
        result += "Total deco time: %u min.\n" % (self.time_to_surface)
        return result

class deco_zhl16:
    # Constants
    N_TISSUES = 16
    P_H2O_BUHLMANN = 0.0627
    P_H2O_NAVY = 0.0567
    P_H2O_SCHREINER = 0.0493
    DECO_ZHL16_MAX_NDL_TIME = 999
    p_h2o = P_H2O_NAVY
    MAX_STOP_TIME = 999
    halftime_n2 = None
    a_n2 = None
    b_n2 = None
    halftime_he = None
    a_he = None
    b_he = None
    
    settings = None
    mix_ctx = None

    def __init__(
            self,
            settings = deco_settings(),
            mix_ctx = deco_mix.deco_mix_ctx(),
            p_surface = 1.0,
            f_n2_init = 0.79,
            f_he_init = 0.0,
            p_h2o = P_H2O_NAVY,
            salinity = 1.0,
            dtype = numpy.float32):
        
        # inputs
        self.settings = settings
        self.mix_ctx = mix_ctx
        self.p_surface = p_surface
        self.salinity = salinity
        self.dtype = dtype
        self.ascend_speed = 10.0
        self.imprecision_tolerance = 0.001
        
        # outputs
        self.p_gf_low = 0.0
        
        self.p_tissues_n2 = numpy.full(deco_zhl16.N_TISSUES, (p_surface - p_h2o) * f_n2_init, dtype=dtype)
        self.p_tissues_he = numpy.full(deco_zhl16.N_TISSUES, (p_surface - p_h2o) * f_he_init, dtype=dtype)
        self.supersaturation_n2 = numpy.full(deco_zhl16.N_TISSUES, 0, dtype=dtype)
        self.supersaturation_he = numpy.full(deco_zhl16.N_TISSUES, 0, dtype=dtype)
        self.a = numpy.full(deco_zhl16.N_TISSUES, 0, dtype=dtype)
        self.b = numpy.full(deco_zhl16.N_TISSUES, 0, dtype=dtype)
        
    def __copy__(self):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def __deepcopy__(self, memo):
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        newone.settings = copy.deepcopy(self.settings, memo)
        newone.mix_ctx = copy.deepcopy(self.mix_ctx, memo)
        newone.p_tissues_n2 = copy.deepcopy(self.p_tissues_n2, memo)
        newone.p_tissues_he = copy.deepcopy(self.p_tissues_he, memo)
        newone.supersaturation_he = copy.deepcopy(self.supersaturation_he, memo)
        newone.supersaturation_n2 = copy.deepcopy(self.supersaturation_n2, memo)
        newone.a = copy.deepcopy(self.a, memo)
        newone.b = copy.deepcopy(self.b, memo)
        return newone
    
    def __str__(self):
        return "settings:\n%smixes:\n%s\ntissues:\nn2: %s\nhe: %s\n" % (self.settings, self.mix_ctx, self.p_tissues_n2, self.p_tissues_he)
    
    def update_settings(self, settings):
        self.settings = settings
    
    def update_mix_ctx(self, mix_ctx):
        self.mix_ctx = mix_ctx
    
    def update_constant(self, time, p_amb):
        # Get acive mix data
        active_mix = self.mix_ctx.get_mix(self.mix_ctx.get_active_mix())
        p_n2 = (p_amb - self.p_h2o) * active_mix.get_fn2()
        p_he = (p_amb - self.p_h2o) * active_mix.get_fhe()

        # Update integrated supersaturation
        ss_n2 = (
            time *
            (-self.p_tissues_n2 + p_n2) +
            (
                (
                    numpy.power(2, -(time / self.halftime_n2)) *
                    (-1 + numpy.power(2, time / self.halftime_n2)) *
                    self.halftime_n2 *
                    (self.p_tissues_n2 - p_n2)
                ) /
                (
                    numpy.log(2)
                )
            )
        )
        self.supersaturation_n2 += numpy.where(ss_n2 > 0, ss_n2, 0)
        ss_he = (
            time *
            (-self.p_tissues_he + p_he) +
            (
                (
                    numpy.power(2, -(time / self.halftime_he)) *
                    (-1 + numpy.power(2, time / self.halftime_he)) *
                    self.halftime_he *
                    (self.p_tissues_he - p_he)
                ) /
                (
                    numpy.log(2)
                )
            )
        )
        self.supersaturation_he += numpy.where(ss_he > 0, ss_he, 0)
        
        # Update tissue loadings.
        self.p_tissues_n2 += numpy.multiply((p_n2 - self.p_tissues_n2),
            (1 - numpy.power(2, -1 * time / self.halftime_n2)))
        self.p_tissues_he += numpy.multiply((p_he - self.p_tissues_he),
            (1 - numpy.power(2, -1 * time / self.halftime_he)))

        # Update weighted a and b values.
        self.a = numpy.divide(self.a_n2 * self.p_tissues_n2 + self.a_he * self.p_tissues_he,
            self.p_tissues_n2 + self.p_tissues_he,
            where=numpy.not_equal(self.p_tissues_n2 + self.p_tissues_he, 0))
        self.b = numpy.divide(self.b_n2 * self.p_tissues_n2 + self.b_he * self.p_tissues_he,
            self.p_tissues_n2 + self.p_tissues_he,
            where=numpy.not_equal(self.p_tissues_n2 + self.p_tissues_he, 0))
    
    def update_asc_desc(self, time, p_begin, p_end):
        # Get acive mix data
        active_mix = self.mix_ctx.get_mix(self.mix_ctx.get_active_mix())
        
        p_begin_alv_n2 = (p_begin - self.p_h2o) * active_mix.get_fn2()
        p_end_alv_n2 = (p_begin - self.p_h2o) * active_mix.get_fn2()
        p_begin_alv_he = (p_end - self.p_h2o) * active_mix.get_fhe()
        p_end_alv_he = (p_end - self.p_h2o) * active_mix.get_fn2()
        R_n2 = ((p_end - self.p_h2o) * active_mix.get_fn2() - p_begin_alv_n2) / time
        R_he = ((p_end - self.p_h2o) * active_mix.get_fhe() - p_begin_alv_he) / time

        # Update integrated supersaturation
        ss_n2 = (
            (time * time * (p_begin_alv_n2 - p_end_alv_he + time * R_n2)) / (2 * time) +
            (
                numpy.exp(-time * (numpy.log(2) / self.halftime_n2)) *
                (
                    -(numpy.log(2) / self.halftime_n2) * self.p_tissues_n2 +
                    (numpy.log(2) / self.halftime_n2) * p_begin_alv_n2 -
                    R_n2 +
                    numpy.exp(time * (numpy.log(2) / self.halftime_n2)) *
                    (
                        R_n2 + (numpy.log(2) / self.halftime_n2) *
                        (self.p_tissues_n2 - p_begin_alv_n2 - time * R_n2)
                    )
                )
            ) /
            (
                (numpy.log(2) / self.halftime_n2) * (numpy.log(2) / self.halftime_n2)
            )
        )
        self.supersaturation_n2 += numpy.where(ss_n2 > 0, ss_n2, 0)
        ss_he = (
            (time * time * (p_begin_alv_he - p_end_alv_he + time * R_he)) / (2 * time) +
            (
                numpy.exp(-time * (numpy.log(2) / self.halftime_he)) *
                (
                    -(numpy.log(2) / self.halftime_he) * self.p_tissues_he +
                    (numpy.log(2) / self.halftime_he) * p_begin_alv_he -
                    R_he +
                    numpy.exp(time * (numpy.log(2) / self.halftime_he)) *
                    (
                        R_he + (numpy.log(2) / self.halftime_he) *
                        (self.p_tissues_he - p_begin_alv_he - time * R_he)
                    )
                )
            ) /
            (
                (numpy.log(2) / self.halftime_he) * (numpy.log(2) / self.halftime_he)
            )
        )
        self.supersaturation_he += numpy.where(ss_he > 0, ss_he, 0)

        # Update tissues gas loadings
        self.p_tissues_n2 = p_begin_alv_n2 + R_n2 * (time - 1 / (numpy.log(2) / self.halftime_n2)) \
            - (p_begin_alv_n2 - self.p_tissues_n2 - (R_n2 / (numpy.log(2) / self.halftime_n2))) \
            * numpy.exp(-(numpy.log(2) / self.halftime_n2) * time)
        self.p_tissues_he = p_begin_alv_he + R_he * (time - 1 / (numpy.log(2) / self.halftime_he)) \
            - (p_begin_alv_he - self.p_tissues_he - (R_he / (numpy.log(2) / self.halftime_he))) \
            * numpy.exp(-(numpy.log(2) / self.halftime_he) * time)

        # Update weighted a and b values.
        self.a = numpy.divide(self.a_n2 * self.p_tissues_n2 + self.a_he * self.p_tissues_he,
            self.p_tissues_n2 + self.p_tissues_he,
            where=numpy.not_equal(self.p_tissues_n2 + self.p_tissues_he, 0))
        self.b = numpy.divide(self.b_n2 * self.p_tissues_n2 + self.b_he * self.p_tissues_he,
            self.p_tissues_n2 + self.p_tissues_he,
            where=numpy.not_equal(self.p_tissues_n2 + self.p_tissues_he, 0))
    
    def ceiling_get_pressure(self, gf):
        return ((self.p_tissues_n2 + self.p_tissues_he - self.a * gf) / (gf / self.b - gf + 1)).max()
    
    def depth_to_pressure(self, depth):
        return self.p_surface + depth * 1.0 * 1.0 * 0.1
    
    def pressure_to_depth(self, p_amb):
        if (p_amb < self.p_surface):
            return 0
        return (p_amb - self.p_surface) / 1.0 / 1.0 / 0.1
    
    def gf_current_get(self, p_amb):
        gf = ((self.p_tissues_n2 + self.p_tissues_he - p_amb) / ((p_amb / self.b + self.a) - p_amb)).max()
        if (gf < 0):
            return 0
        return gf

    def gf_allowed_get(self, p_amb):
        if (p_amb <= self.p_surface):
            return self.settings.gf_high
        
        if (p_amb >= self.p_gf_low):
            return self.settings.gf_low
        
        return (p_amb - self.p_surface) \
            * ((self.settings.gf_high - self.settings.gf_low) \
            / (0 - (self.p_gf_low - self.p_surface))) \
            + self.settings.gf_high
    
    def lead_tissue_get(self, p_amb):
        max_gf = 0
        lead_tissue = 0
        
        for i in range(0, self.N_TISSUES):
            this_gf = (self.p_tissues_n2[i] + self.p_tissues_he[i] - p_amb) \
                / ((p_amb / self.b[i] + self.a[i]) - p_amb)
            
            if this_gf > max_gf:
                max_gf = this_gf
                lead_tissue = i
            
        return lead_tissue

    def ndt_get(self, p_amb):
        ndt = 0

        # Create a temporary context copy.
        temp_ctx = copy.deepcopy(self)
        
        while (ndt < deco_zhl16.DECO_ZHL16_MAX_NDL_TIME):
            # Simulate ascent to surface.
            asc_ctx = copy.deepcopy(temp_ctx)
            
            # Only simulate ascent if there is a pressure difference.
            if (p_amb > self.p_surface):
                asc_ctx.update_asc_desc(
                time = asc_ctx.pressure_to_depth(p_amb) / asc_ctx.settings.ascent_rate,
                p_begin = p_amb,
                p_end = asc_ctx.p_surface)
            
            # Check ceiling against simulated ascent.
            if (asc_ctx.ceiling_get_pressure(asc_ctx.settings.gf_high) >= asc_ctx.p_surface):
                break
            
            # Advance one minute.
            temp_ctx.update_constant(
                time = 1.0,
                p_amb = p_amb)
            
            ndt += 1
        
        return ndt

    def decotime_update(self, p_amb, di):
        # Create a temporary context copy.
        temp_ctx = copy.deepcopy(self)
        
        di.stops = numpy.zeros(di.N_STOPS, dtype=[('mix', int), ('depth', self.dtype), ('time', self.dtype)])
        di.time_to_surface = 0
        stop_time = 0
        
        # Determine first ceiling pressure.
        p_ceiling = temp_ctx.ceiling_get_pressure(self.settings.gf_low)
        
        # Determine lowest depth where no stop is required (worst-case).
        # This is where we will start looking for required stops.
        stop_depth = numpy.ceil(self.pressure_to_depth(p_ceiling))
        
        # First stop candidate is the first higher depth on the 3m raster.
        if (stop_depth % self.settings.stop_size):
            stop_depth = stop_depth - (stop_depth % self.settings.stop_size) + self.settings.stop_size
        
        # Raise pressure of GF_LOW to the first stop depth.
        if (self.depth_to_pressure(stop_depth) > self.p_gf_low + self.imprecision_tolerance):
            self.p_gf_low = self.depth_to_pressure(stop_depth)
        # Lower pressure of GF_LOW in case the diver is still below deco zone.
        #elif (self.depth_to_pressure(stop_depth) < self.p_gf_low - 0.001 and
        #    p_amb > self.depth_to_pressure(stop_depth)):
        #        self.p_gf_low = self.depth_to_pressure(stop_depth)
        
        # If we are deeper than the first stop candidate, simulate ascent accordingly.
        if (p_amb > self.depth_to_pressure(stop_depth) + self.imprecision_tolerance):
            
            time = (self.pressure_to_depth(p_amb) - stop_depth) / self.ascend_speed
            
            temp_ctx.update_asc_desc(
                time = time,
                p_begin = p_amb,
                p_end = self.depth_to_pressure(stop_depth))
            
            # Calculate remaining time to fill up minute.
            remainder = numpy.ceil(time) - time
            
            # Stay at first stop depth for remaining time.
            if (remainder > 1 / 60):
                temp_ctx.update_constant(
                    time = remainder,
                    p_amb = self.depth_to_pressure(stop_depth))
                stop_time += 1
        
        i_stop = 0
        
        while (stop_depth >= self.settings.last_stop and stop_time < deco_zhl16.MAX_STOP_TIME):
            # Switch to better mix if available.
            best_mix = self.mix_ctx.get_best_mix(stop_depth)
            if (best_mix != temp_ctx.mix_ctx.get_active_mix()):
                temp_ctx.mix_ctx.set_active_mix(best_mix)

            if (stop_depth == self.settings.last_stop):
                clearance_depth = 0
            else:
                clearance_depth = stop_depth - self.settings.stop_size
            
            # Update current gf.
            gf = self.gf_allowed_get(self.depth_to_pressure(clearance_depth))
            
            # Check if a stop is required using the current gf.
            p_ceiling = temp_ctx.ceiling_get_pressure(gf)
            
            if (p_ceiling < self.depth_to_pressure(clearance_depth)):
                # No (more) stopping required.
                if (stop_time == 0):
                    # No stop required so do the minimum stop of one minute
                    # but take ascent time from last stop into account. */
                    temp_ctx.update_constant(
                        time = 1 - self.settings.stop_size / self.settings.ascent_rate,
                        p_amb = self.depth_to_pressure(stop_depth))
                    
                    stop_time += 1
            
                # Commit stop to the deco info result.
                if (stop_time > 0):
                    if (i_stop < di.N_STOPS):
                        di.stops[i_stop]['mix'] = temp_ctx.mix_ctx.get_active_mix()
                        di.stops[i_stop]['time'] = stop_time
                        di.stops[i_stop]['depth'] = stop_depth
                    i_stop += 1
                    di.time_to_surface += stop_time
                
                # Simulate ascent to next stop.
                temp_ctx.update_asc_desc(
                    time = (stop_depth - clearance_depth) / self.ascend_speed,
                    p_begin = self.depth_to_pressure(stop_depth),
                    p_end = self.depth_to_pressure(clearance_depth))
                
                # Advance stop depth.
                stop_depth -= self.settings.stop_size

                # Reset local stop data.
                stop_time = 0
                
            else:
                # Stop is required, so do it.
                if (stop_time == 0):
                    # This is the first waiting segment on this stop
                    # so subtract the ascent time from last stop.
                    temp_ctx.update_constant(
                        time = 1 - self.settings.stop_size / self.settings.ascent_rate,
                        p_amb = self.depth_to_pressure(stop_depth))
                else:
                    # This is not the first waiting segment on this stop
                    # so do the full minute.
                    temp_ctx.update_constant(
                        time = 1,
                        p_amb = self.depth_to_pressure(stop_depth))
                    
                stop_time += 1
        
        # Finally add ascent time up to first stop to the TTS value.
        if (p_amb > self.depth_to_pressure(di.stops[0]['depth'])):
            # Calculate time for ascent.
            time = (self.pressure_to_depth(p_amb) - di.stops[0]['depth']) \
                / self.ascend_speed
            
            # Add initial ascent from ambient to first stop.deepcopy
            di.time_to_surface += numpy.ceil(time)
            
        return di

    def decoinfo_update(self, p_amb):
        di = deco_info()
        di.nostop_time = self.ndt_get(p_amb)
        if (di.nostop_time == 0):
            # Calculate deco plan.
            di = self.decotime_update(p_amb, di)
            # Update gf and gf_low depth.
            di.allowed_gf = self.gf_allowed_get(p_amb)
            di.gf_low_depth = self.pressure_to_depth(self.p_gf_low)
        else:
            di.time_to_surface = 0
            # Update gf and gf_low depth. Note that ndt is calculated
            # against gf_high at surface pressure only.
            di.allowed_gf = self.settings.gf_high
            di.gf_low_depth = 0.0
        
        # Update lead tissue number and gf. */
        di.lead_tissue = self.lead_tissue_get(p_amb)
        di.lead_tissue_gf = self.gf_current_get(p_amb)
        
        return di

class deco_zhl16b(deco_zhl16):
    # Constants
    halftime_n2 = numpy.array(
        [5.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0,
         109.0, 146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0], dtype=numpy.float32)
    a_n2 = numpy.array(
        [1.1696, 1.0, 0.8618, 0.7562, 0.6667, 0.56, 0.4947, 0.45, 
         0.4187, 0.3798, 0.3497, 0.3223, 0.2850, 0.2737, 0.2523, 0.2327], dtype=numpy.float32)
    b_n2 = numpy.array(
        [0.5578, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.891,
         0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653], dtype=numpy.float32)
    halftime_he = numpy.array(
        [1.88, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11,
         41.2, 55.19, 70.69, 90.34, 115.29, 147.42, 188.24, 240.03], dtype=numpy.float32)
    a_he = numpy.array(
        [1.6189, 1.383, 1.1919, 1.0458, 0.922, 0.8205, 0.7305, 0.6502,
         0.595, 0.5545, 0.5333, 0.5189, 0.5181, 0.5176, 0.5172, 0.5119], dtype=numpy.float32)
    b_he = numpy.array(
        [0.477, 0.5747, 0.6527, 0.7223, 0.7582, 0.7957, 0.8279, 0.8553,
         0.8757, 0.8903, 0.8997, 0.9073, 0.9122, 0.9171, 0.9217, 0.9267], dtype=numpy.float32)

class deco_zhl16c(deco_zhl16):
    # Constants
    halftime_n2 = numpy.array(
        [4.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0,
         109.0, 146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0], dtype=numpy.float32)
    a_n2 = numpy.array(
        [1.2599, 1.0, 0.8618, 0.7562, 0.62, 0.5043, 0.441, 0.4, 
         0.375, 0.35, 0.3295, 0.3065, 0.2835, 0.261, 0.248, 0.2327], dtype=numpy.float32)
    b_n2 = numpy.array(
        [0.505, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.891,
         0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653], dtype=numpy.float32)
    halftime_he = numpy.array(
        [1.51, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11,
         41.2, 55.2, 70.7, 90.2, 115.2, 147.2, 188.0, 240.0], dtype=numpy.float32)
    a_he = numpy.array(
        [1.7424, 1.383, 1.1919, 1.0458, 0.922, 0.8205, 0.7305, 0.6502,
         0.595, 0.5545, 0.5333, 0.5189, 0.5181, 0.5176, 0.5172, 0.5119], dtype=numpy.float32)
    b_he = numpy.array(
        [0.4245, 0.5747, 0.6527, 0.7223, 0.7582, 0.7957, 0.8279, 0.8553,
         0.8757, 0.8903, 0.8997, 0.9073, 0.9122, 0.9171, 0.9217, 0.9267], dtype=numpy.float32)

if __name__ == '__main__':
    settings = deco_settings(last_stop = 6)
    mixes = deco_mix.deco_mix_ctx()
    mixes.set_mix(1, deco_mix.deco_mix(f_o2 = 0.21, f_he = 0.35, enabled = False))
    mixes.set_mix(2, deco_mix.deco_mix(f_o2 = 0.50, f_he = 0.25, enabled = True, switch_depth = 21))
    mixes.set_mix(3, deco_mix.deco_mix(f_o2 = 1, f_he = 0, enabled = True, switch_depth = 6))
    mixes.set_active_mix(1)
    deco_ctx = deco_zhl16b(settings, mixes)
    
    deco_ctx.update_asc_desc(5, 1.0, 6)
    deco_ctx.update_constant(45, 6)

    di = deco_ctx.decoinfo_update(6)

    print (deco_ctx)
    print (di)
    

