import numpy as np 
import matplotlib.pyplot as plt
# import data_gen
# import sklearn.kernel_ridge
import itertools


def correlation(y, y_pred):
	m = y.size 
	ymean = np.mean(y)
	vary = np.sum((y-ymean)**2)
	vary_pred = np.sum((y-y_pred)**2)

	return 1 - vary_pred/vary


def linear_regression(x, y, add_ones=True):
	if add_ones:
		x = np.hstack((np.ones((x.shape[0],1)), x))
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


class MLModel:
	def __init__(self):
		self.a = np.array([])
		self.description = 'Template Machine Learning Model'
		self.trained = False

	def __repr__(self):
		s = '' 
		s += self.description + '\n'
		s += f'This model has{" not"*(not self.trained)} been trained'
		return s

	def train(self, Xtrain, Ytrain):
		'''
		Method that takes training X (m x n) and Y (m x 1) arrays
		and trains the model to set the parameter "a"
		'''
		self.trained = True
		NotImplemented

	def evaluate(self, Xtest, Ytest):
		'''
		Method that takes test X (m x n) and Y (m x 1) arrays
		and evaluates the model, that is predict the Ytest using
		trained parameters. Then it will calculate the R2 value
		'''
		NotImplemented

	def predict(self, X):
		'''
		Method that takes X (m x n) array and 
		predicts the Y-values corresponding to these datapoints
		Does not calculate R2
		'''
		NotImplemented


class LR(MLModel):
	def __init__(self):
		super().__init__()
		self.description = 'Linear Regression Model'

	def train(self, Xtrain, Ytrain, add_ones=True):
		if add_ones:
			Xtrain = np.hstack((np.ones((Xtrain.shape[0],1)), Xtrain))
		#perform regression
		self.a = np.linalg.inv(Xtrain.T@Xtrain) @ Xtrain.T @ Ytrain
		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def evaluate(self, Xtest, Ytest, add_ones=True):
		if add_ones:
			Xtest = np.hstack((np.ones((Xtest.shape[0],1)), Xtest))
		Ypred = Xtest@self.a
		R2 = correlation(Ytest, Ypred)
		return R2

	def predict(self, X, add_ones=True):
		if add_ones:
			X = np.hstack((np.ones((X.shape[0],1)), X))
		return X@self.a


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



class KRR(MLModel):
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
		

