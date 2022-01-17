# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:36:22 2020

@author: feierabend
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:35:32 2020

@author: feierabend
"""

import numpy as np
cimport cython

DTYPE = np.float64

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def calc_heat_transfer(double[:] wall_temp, double[:] fluid_temp, 
                       double[:] capacity_rate, double[:] heat_coeff,
                       int flow_direction):

    assert tuple(capacity_rate.shape) == tuple(wall_temp.shape)
    assert tuple(heat_coeff.shape) == tuple(wall_temp.shape)  
    
    cdef Py_ssize_t elements = wall_temp.shape[0]
    cdef Py_ssize_t nodes = fluid_temp.shape[0]
    heat = np.zeros(elements, dtype=DTYPE)
    temp_result = np.zeros(nodes, dtype=DTYPE)
    cdef double[:] temp_result_view = temp_result
    cdef double[:] heat_view = heat
    cdef double fluid_avg = 0.0
    cdef double fluid_out_old = 5e5
    cdef double fluid_in = 0.0
    cdef double fluid_out = 0.0
    cdef double q = 0.0
    cdef double error = 1e3
    cdef double delta_temp = 0.0
    cdef int iter = 0
    cdef int itermax = 10
    cdef fluid_temp_avg = 0.0
    cdef int i
    temp_result_view[:] = fluid_temp
    if flow_direction == -1:
        for i in reversed(range(elements)):
            fluid_avg = (fluid_temp[i + 1] + fluid_temp[i]) * 0.5
            fluid_out_old = 5e5
            error = 1e3
            iter = 0
            itermax = 10
            while error > 1e-5 and iter <= itermax:
                fluid_in = fluid_temp[i + 1]
                
                if wall_temp[i] == fluid_in:
                    delta_temp = wall_temp[i] - fluid_avg
                else:
                    temp_diff_ratio = (wall_temp[i] - fluid_avg) \
                        / (wall_temp[i] - fluid_in)
                    if temp_diff_ratio > 0.0:
                        delta_temp = wall_temp[i] - fluid_avg
                    else:
                        delta_temp = wall_temp[i] - fluid_in
                    
                q = heat_coeff[i] * delta_temp        
                fluid_out = fluid_in + q/capacity_rate[i]
                if fluid_in < wall_temp[i]:
                    fluid_out = min(wall_temp[i] - 1e-4, fluid_out)
                else:
                    fluid_out = max(wall_temp[i] + 1e-4, fluid_out)
                fluid_avg = (fluid_in + fluid_out) * 0.5
                error = abs(fluid_out_old - fluid_out)/fluid_out
                fluid_out_old = fluid_out
                iter += 1             
            temp_result_view[i] = fluid_out
            heat_view[i] = heat_coeff[i] * (wall_temp[i] - fluid_avg)
    else:         
        for i in range(elements):
            fluid_avg = (fluid_temp[i + 1] + fluid_temp[i]) * 0.5
            fluid_out_old = 5e5
            error = 1e3
            iter = 0
            itermax = 10
            while error > 1e-5 and iter <= itermax:
                fluid_in = fluid_temp[i]
                if wall_temp[i] == fluid_in:
                    delta_temp = wall_temp[i] - fluid_avg
                else:
                    temp_diff_ratio = (wall_temp[i] - fluid_avg) \
                        / (wall_temp[i] - fluid_in)
                    if temp_diff_ratio > 0.0:
                        delta_temp = wall_temp[i] - fluid_avg
                    else:
                        delta_temp = wall_temp[i] - fluid_in
                        
                delta_temp = wall_temp[i] - fluid_avg
                q = heat_coeff[i] * delta_temp        
                fluid_out = fluid_in + q/capacity_rate[i]
                if fluid_in < wall_temp[i]:
                    fluid_out = min(wall_temp[i] - 1e-4, fluid_out)
                else:
                    fluid_out = max(wall_temp[i] + 1e-4, fluid_out)
                fluid_avg = (fluid_in + fluid_out) * 0.5
                error = abs(fluid_out_old - fluid_out)/fluid_out
                fluid_out_old = fluid_out
                iter += 1
            temp_result_view[i + 1] = fluid_out
            heat_view[i] = heat_coeff[i] * (wall_temp[i] - fluid_avg)
    return temp_result, heat