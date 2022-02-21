import numpy as np
import matplotlib.pyplot as plt


def random_data(m, n, noise=.1):
	#will generate m random datapoints of n variables

	#first create alpha
	a = 2*np.random.rand(n,1)-1 #column

	x = np.random.rand(m, n)
	y = x@a + (2*np.random.rand(m,1)-1)*noise

	return {'a':a, 'x':x, 'y':y}


def non_linear_data(n=2, N=1000):
	d = 2
	c = 1
	X = np.random.rand(n, N)
	K = np.zeros((N,N))
	print(X)
	for i, x in enumerate(X.T):
		for j, y in enumerate(X.T):
			K[i,j] = (x@y.T+c)**d
	return K

