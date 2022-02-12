import numpy as np 
import matplotlib.pyplot as plt
import data_gen


def correlation(y, y_pred):
	m = y.size 
	ymean = np.mean(y)
	vary = np.sum((y-ymean)**2)
	vary_pred = np.sum((y-y_pred)**2)

	return 1 - vary_pred/vary


def linear_regression(data):
	#derivation: 
	x = data['x']
	y = data['y']
	#perform regression
	a = np.linalg.inv(x.T@x) @ x.T @ y
	return {'a':a, 'R2':correlation(y, x@a)}


def ridge_regression(data, lam):
	#derivation: 
	x = data['x']
	y = data['y']
	#perform regression
	a = np.linalg.inv(x.T@x + lam) @ x.T @ y
	return {'a':a, 'R2':correlation(y, x@a)}




data = data_gen.random_data(600, 1, noise=0.1)
plt.scatter(data['x'][:,1:], data['y'])
for l in [1000,100,10,1, .1, .01, .001, .0001, .00001]:
	a = ridge_regression(data,l)
	plt.scatter(data['x'][:,1:], data['x']@a['a'])
# plt.plot(data['x'])


# 
plt.show()
