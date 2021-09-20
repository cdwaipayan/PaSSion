import argparse
import numpy as np
import vegas
import math

def angular_partition_function(x):    
    return math.exp( lamor*(x**2 - 1.0) )

# Parameter to change
nitn   = 100
neval  = 5000000
lamor  = 900.0

print( 'nitn = %s  neval = %s \n' % (nitn, neval) )
print( 'lamor = %s \n' % (lamor) )

integ = vegas.Integrator([[0, 1]])
result = integ(angular_partition_function, nitn=nitn, neval=neval)

print(result.summary())
print('result = %s    Q = %.2f' % (result, result.Q))
print('Orientational Free Energy per Particle = %s' % -math.log(result.mean) )
