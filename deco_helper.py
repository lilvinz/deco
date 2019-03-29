import copy
import numpy
import deco


class profile:
    i = 0
    ctx = None
    di = None
    def __init__(self, runtime, ctx, maxentries):
        self.runtime = runtime
        self.ctx = copy.deepcopy(ctx)
        self.data = numpy.zeros(
            maxentries,
            dtype=[('time', ctx.dtype), ('depth', ctx.dtype), ('p_amb_n2', ctx.dtype, 16), ('p_amb_he', ctx.dtype, 16), ('p_ceiling', ctx.dtype), ('gf_allowed', ctx.dtype),
                   ('gf_current', ctx.dtype), ('p_tissues_n2', ctx.dtype, 16), ('p_tissues_he', ctx.dtype, 16),
                   ('supersaturation_n2', ctx.dtype, 16), ('supersaturation_he', ctx.dtype, 16)])
        self.ctx.mix_ctx.clear_mix_volumes()

    def update(self, time, depth_begin, depth_end, sac = 20):
        if (depth_begin == depth_end):
            self.ctx.update_constant(
                time = time,
                p_amb = self.ctx.depth_to_pressure(depth_end))
            self.ctx.mix_ctx.mixes[self.ctx.mix_ctx.active_mix].volume += \
                sac * time * self.ctx.depth_to_pressure(depth_end)
        else:
            self.ctx.update_asc_desc(
                time = time,
                p_begin = self.ctx.depth_to_pressure(depth_begin),
                p_end = self.ctx.depth_to_pressure(depth_end))
            self.ctx.mix_ctx.mixes[self.ctx.mix_ctx.active_mix].volume += \
                sac * time * \
                self.ctx.depth_to_pressure(depth_begin + (depth_begin - depth_end) / 2)

        self.di = deco.deco_info()
        self.di = self.ctx.decotime_update(
            p_amb = self.ctx.depth_to_pressure(depth_end),
            di = self.di)

        self.runtime += time
        self.data[self.i]['time'] = self.runtime
        self.data[self.i]['depth'] = depth_end
        self.data[self.i]['p_amb_n2'] = depth_end / 10 * self.ctx.mix_ctx.mixes[self.ctx.mix_ctx.active_mix].get_fn2()
        self.data[self.i]['p_amb_he'] = depth_end / 10 * self.ctx.mix_ctx.mixes[self.ctx.mix_ctx.active_mix].get_fhe()
        numpy.copyto(self.data[self.i]['p_tissues_n2'], self.ctx.p_tissues_n2)
        numpy.copyto(self.data[self.i]['p_tissues_he'], self.ctx.p_tissues_he)
        numpy.copyto(self.data[self.i]['supersaturation_n2'], self.ctx.supersaturation_n2)
        numpy.copyto(self.data[self.i]['supersaturation_he'], self.ctx.supersaturation_he)
        self.data[self.i]['gf_current'] = self.ctx.gf_current_get(self.ctx.depth_to_pressure(depth_end))
        self.data[self.i]['gf_allowed'] = self.ctx.gf_allowed_get(self.ctx.depth_to_pressure(depth_end))
        self.data[self.i]['p_ceiling'] = self.ctx.ceiling_get_pressure(self.data[self.i]['gf_allowed'])
        self.i += 1
        
    def strip(self):
        self.data = self.data[:self.i]
