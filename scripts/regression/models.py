import numpy as np 
import matplotlib.pyplot as plt
import data_gen
import sklearn



def linear_regression(data):
	#derivation: 
	x = data['x']
	y = data['y']
	#perform regression
	a = np.linalg.inv(x.T@x) @ x.T @ y
	#calculate correlation
	m = y.size
	e = y - x@a
	s = y - np.sqrt(y.T@y)/m
	R2 = 1 - (e.T@e)/(s.T@s)
	return {'a':a, 'R2':R2}


def ridge_regression(data):
	...



data = data_gen.random_data(600, 12)
a = linear_regression(data)

print(a['a'])
print(data['a'])
print(a['R2'])