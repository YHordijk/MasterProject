import numpy as np 
import matplotlib.pyplot as plt
import data_gen
import sklearn.kernel_ridge
import itertools


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



class FeatureMap:
	def __call__(self, x):
		NotImplemented


class LinearFM(FeatureMap):
	def __call__(self, x):
		return x

class PolynomialFM(FeatureMap):
	def __init__(self, n):
		self.n = n

	def __call__(self, x):
		return np.array(list(itertools.combinations_with_replacement(x, self.n))).reshape(-1,1)


class QuadraticFM(FeatureMap):
	def __call__(self, x):
		return np.array(list(itertools.combinations_with_replacement(x, 2))).reshape(-1,1)




class KRR:
	def __init__(self, kernel='polynomial', **kwargs):
		self.kernel_type = kernel
		if kernel == 'polynomial':
			self.hyperparameters = {'alpha':	kwargs.get('alpha',  1),
									'degree': 	kwargs.get('degree', 1)}

	# def train(self, X, Y):
	# 	'''
	# 	X is an array where each row is a datapoint. The number of columns is te number of independent variables
	# 	Y is an array with observed outcomes corresponding to the data in X
	# 	'''

	# 	N = Y.size
	# 	m = X.shape[1]

	# def get_kernel(self):
		


# help(sklearn.kernel_ridge)

# data = data_gen.random_data(600, 1, noise=0.1)
# # plt.scatter(data['x'][:,1:], data['y'])
# # for l in [1000,100,10,1, .1, .01, .001, .0001, .00001]:
# # 	a = ridge_regression(data,l)
# # 	plt.scatter(data['x'][:,1:], data['x']@a['a'])

# # plt.show()

# model = KRR()
# model.train(data['x'], data['y'])


U = [1, 2, 3, 4]
V = [5, 6, 7, 8]

fm = PolynomialFM(5)

fu = fm(U)
fv = fm(V)

print(fu, fv)
print(fu.T @ fv)