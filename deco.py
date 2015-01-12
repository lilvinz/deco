'''
Created on 8 Jan 2015

@author: vke
'''

import numpy
import copy

class DecoZHL16:
    # Constants
    n_tissues = 16
    p_h2o_buhlmann = 0.0627
    p_h2o_navy = 0.0567
    p_h2o_schreiner = 0.0493

    def __init__(
            self,
            p_surface = 1.0,
            f_n2_init = 0.79,
            f_he_init = 0.0,
            gf_low = 0.1,
            gf_high = 0.9,
            p_h2o = p_h2o_navy,
            salinity = 1.0,
            dtype = numpy.float32):
        
        # inputs
        self.p_surface = p_surface
        self.gf_low = gf_low
        self.gf_high = gf_high
        self.p_h2o = p_h2o
        self.salinity = salinity
        self.dtype = dtype
        
        self.n_stops = 200
        self.ascend_speed = 10.0
        self.last_stop = 3
        self.stop_size = 3
        
        # outputs
        self.p_gf_low = 0.0
        
        self.p_tissues_n2 = numpy.full(DecoZHL16.n_tissues, (p_surface - p_h2o) * f_n2_init, dtype=dtype)
        self.p_tissues_he = numpy.full(DecoZHL16.n_tissues, (p_surface - p_h2o) * f_he_init, dtype=dtype)
        self.a = numpy.full(DecoZHL16.n_tissues, 0, dtype=dtype)
        self.b = numpy.full(DecoZHL16.n_tissues, 0, dtype=dtype)
        
        self.stops = numpy.zeros(self.n_stops, dtype=[('time', dtype), ('depth', dtype)])
        self.time_to_surface = 0
        
    def __str__(self):
        return "%s(%r)" % (self.p_tissues_n2, self.p_tissues_he)
    
    def update_constant(self, time, p_ambient, f_n2, f_he):
        # Update tissue loadings.
        self.p_tissues_n2 += ((p_ambient - self.p_h2o) * f_n2 - self.p_tissues_n2) \
            * (1 - numpy.power(2, -1 * time / self.halftime_n2))
        self.p_tissues_he += ((p_ambient - self.p_h2o) * f_he - self.p_tissues_he) \
            * (1 - numpy.power(2, -1 * time / self.halftime_he))
        
        # Update weighted a and b values.
        if (self.p_tissues_n2.all() + self.p_tissues_he.all()) != 0:
            self.a = (self.a_n2 * self.p_tissues_n2 + self.a_he * self.p_tissues_he) \
                / (self.p_tissues_n2 + self.p_tissues_he);
            self.b = (self.b_n2 * self.p_tissues_n2 + self.b_he * self.p_tissues_he) \
                / (self.p_tissues_n2 + self.p_tissues_he);
            
    def update_asc_desc(self, time, p_start, p_end, f_n2, f_he):
        p_start_alv_n2 = self.dtype((p_start - self.p_h2o) * f_n2)
        p_start_alv_he = self.dtype((p_start - self.p_h2o) * f_he)
        R_n2 = self.dtype(((p_end - self.p_h2o) * f_n2 - p_start_alv_n2) / time)
        R_he = self.dtype(((p_end - self.p_h2o) * f_he - p_start_alv_he) / time)
        
        for i in range(0, self.n_tissues):
            # Update tissues gas loadings
            k = self.dtype(numpy.log(2) / self.halftime_n2[i])
            self.p_tissues_n2[i] = p_start_alv_n2 + R_n2 * (time - 1 / k) \
                - (p_start_alv_n2 - self.p_tissues_n2[i] - (R_n2 / k)) \
                * numpy.exp(-k * time)
            
            k = self.dtype(numpy.log(2) / self.halftime_he[i])
            self.p_tissues_he[i] = p_start_alv_he + R_he * (time - 1 / k) \
                - (p_start_alv_he - self.p_tissues_he[i] - (R_he / k)) \
                * numpy.exp(-k * time)
            
            # Update the weighted a and b values
            if (self.p_tissues_n2[i] + self.p_tissues_he[i]) != 0:
                self.a[i] = (self.a_n2[i] * self.p_tissues_n2[i] + self.a_he[i] * self.p_tissues_he[i]) \
                    / (self.p_tissues_n2[i] + self.p_tissues_he[i])
                self.b[i] = (self.b_n2[i] * self.p_tissues_n2[i] + self.b_he[i] * self.p_tissues_he[i]) \
                    / (self.p_tissues_n2[i] + self.p_tissues_he[i])
            
    def ceiling_get_pressure(self, gf):
        #return (self.p_tissues_n2 + self.p_tissues_he - self.a * gf) / (gf / self.b - gf + 1).max()
        ceiling = 0
        for i in range(0, self.n_tissues):
            this_ceiling = self.dtype((self.p_tissues_n2[i] + self.p_tissues_he[i] - self.a[i] * gf)
                / (gf / self.b[i] - gf + 1))
                
            if this_ceiling > ceiling:
                ceiling = this_ceiling
                
        return ceiling;
    
    def depth_to_pressure(self, depth):
        return self.p_surface + depth * 1.0 * 1.0 * 0.1
    
    def pressure_to_depth(self, p_ambient):
        if (p_ambient < self.p_surface):
            return 0
        return (p_ambient - self.p_surface) / 1.0 / 1.0 / 0.1
    
    def gf_current_get(self, p_ambient):
        #return (self.p_tissues_n2 + self.p_tissues_he - p_ambient) / (((p_ambient / self.b + self.a) - p_ambient)).max()
        
        max_gf = 0
        
        for i in range(0, self.n_tissues):
            this_gf = self.dtype((self.p_tissues_n2[i] + self.p_tissues_he[i] - p_ambient)
                / ((p_ambient / self.b[i] + self.a[i]) - p_ambient))
                
            if this_gf > max_gf:
                max_gf = this_gf
                
        return max_gf;
    
    def gf_allowed_get(self, p_ambient):
        if (p_ambient <= self.p_surface):
            return self.gf_high
        
        if (p_ambient >= self.p_gf_low):
            return self.gf_low
        
        return self.dtype((p_ambient - self.p_surface)
            * ((self.gf_high - self.gf_low)
            / (0 - (self.p_gf_low - self.p_surface)))
            + self.gf_high)
    
    def lead_tissue_get(self, p_ambient):
        max_gf = 0
        lead_tissue = 0
        
        for i in range(0, self.n_tissues):
            this_gf = (self.p_tissues_n2[i] + self.p_tissues_he[i] - p_ambient) \
                / ((p_ambient / self.b[i] + self.a[i]) - p_ambient)
                
            if this_gf > max_gf:
                max_gf = this_gf
                lead_tissue = i
                
        return lead_tissue;
    
    def decotime_update(self, p_ambient, f_n2, f_he):
        # Create a temporary context copy.
        temp_ctx = copy.deepcopy(self)
        
        self.stops = numpy.zeros(self.n_stops, dtype=[('depth', self.dtype), ('time', self.dtype)])
        self.time_to_surface = 0
        stop_time = 0
        
        # Determine first ceiling pressure.
        _p_ceiling = self.dtype(temp_ctx.ceiling_get_pressure(self.gf_low))
        
        # Determine lowest depth where no stop is required (worst-case).
        # This is where we will start looking for required stops.
        stop_depth = numpy.ceil(self.pressure_to_depth(_p_ceiling))
        
        # First stop candidate is the first higher depth on the 3m raster.
        if stop_depth % self.stop_size:
            stop_depth = stop_depth - (stop_depth % self.stop_size) + self.stop_size
        
        # Raise pressure of GF_LOW to the first stop depth.
        if (self.depth_to_pressure(stop_depth) > self.p_gf_low + 0.001):
            self.p_gf_low = self.depth_to_pressure(stop_depth)
        # Lower pressure of GF_LOW in case the diver is still below deco zone.
        #elif (self.depth_to_pressure(stop_depth) < self.p_gf_low - 0.001 and
        #    p_ambient > self.depth_to_pressure(stop_depth)):
        #        self.p_gf_low = self.depth_to_pressure(stop_depth)
        
        # If we are deeper than the first stop candidate, simulate ascent accordingly.
        if (p_ambient > self.depth_to_pressure(stop_depth)):
            
            time = (self.pressure_to_depth(p_ambient) - stop_depth) \
                / self.ascend_speed
            
            temp_ctx.update_asc_desc(
                time = time,
                p_start = p_ambient,
                p_end = self.depth_to_pressure(stop_depth),
                f_n2 = f_n2,
                f_he = f_he)
            
            # Calculate remaining time to fill up minute.
            remainder = numpy.ceil(time) - time
            
            # Stay at first stop depth for remaining time.
            if (remainder > 0):
                temp_ctx.update_constant(
                    time = remainder,
                    p_ambient = self.depth_to_pressure(stop_depth),
                    f_n2 = f_n2,
                    f_he = f_he)
                ++stop_time
                
        i_stop = 0
        
        while (stop_depth >= self.last_stop):
            clearance_depth = stop_depth - self.stop_size
            
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
                        1 - self.stop_size / self.ascend_speed,
                        self.depth_to_pressure(stop_depth),
                        f_n2 = f_n2,
                        f_he = f_he)
                    
                    stop_time += 1
            
                # Commit stop to the deco info result.
                if (stop_time > 0):
                    if (i_stop < self.n_stops):
                        self.stops[i_stop]['time'] = stop_time
                        self.stops[i_stop]['depth'] = stop_depth
                    i_stop += 1
                    self.time_to_surface += stop_time
                
                # Simulate ascent to next stop.
                temp_ctx.update_asc_desc(
                    (stop_depth - clearance_depth) / self.ascend_speed,
                    self.depth_to_pressure(stop_depth),
                    self.depth_to_pressure(clearance_depth),
                    f_n2 = f_n2,
                    f_he = f_he)
                
                # Advance stop depth.
                stop_depth -= self.stop_size

                # Reset local stop data.
                stop_time = 0
                
            else:
                # Stop is required, so do it.
                if (stop_time == 0):
                    # This is the first waiting segment on this stop
                    # so subtract the ascent time from last stop.
                    temp_ctx.update_constant(
                        1.0 - self.stop_size / self.ascend_speed,
                        self.depth_to_pressure(stop_depth),
                        f_n2 = f_n2,
                        f_he = f_he)
                else:
                    # This is not the first waiting segment on this stop
                    # so do the full minute.
                    temp_ctx.update_constant(
                        1.0,
                        self.depth_to_pressure(stop_depth),
                        f_n2 = f_n2,
                        f_he = f_he)
                    
                stop_time += 1
        
        # Finally add ascent time up to first stop to the TTS value.
        if (p_ambient > self.depth_to_pressure(self.stops[0]['depth'])):
            # Calculate time for ascent.
            time = \
                (self.pressure_to_depth(p_ambient) - self.stops[0]['depth']) \
                / self.ascend_speed
            
            # Add initial ascent from ambient to first stop.
            self.time_to_surface += int(time + 0.999)

class DecoZHL16B(DecoZHL16):
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
