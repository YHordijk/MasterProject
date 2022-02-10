import numpy as np
import matplotlib.pyplot as plt


def random_data(m, n, noise=.1):
	#will generate m random datapoints of n variables

	#first create alpha
	a = 2*np.random.rand(n+1,1)-1 #column

	x = np.random.rand(m, n+1)
	x[:,0] = 1
	y = x@a + (2*np.random.rand(m,1)-1)*noise

	return {'a':a, 'x':x, 'y':y}

