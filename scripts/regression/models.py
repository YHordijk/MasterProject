import numpy as np 
import matplotlib.pyplot as plt
# import data_gen
# import sklearn.kernel_ridge
import itertools
import kernels


def correlation(y, y_pred):
	m = y.size 
	ymean = np.mean(y)
	vary = np.sum((y-ymean)**2)
	vary_pred = np.sum((y-y_pred)**2)

	return 1 - vary_pred/vary


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


class RR(MLModel):
	def __init__(self):
		super().__init__()
		self.lam = None
		self.description = f'Ridge Regression Model'

	def train(self, Xtrain, Ytrain, lam=1, add_ones=True):
		if add_ones:
			Xtrain = np.hstack((np.ones((Xtrain.shape[0],1)), Xtrain))
		n = Xtrain.shape[1]
		#perform regression
		self.a = np.linalg.inv(Xtrain.T@Xtrain + np.eye(n)*lam) @ Xtrain.T @ Ytrain
		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X, add_ones=True):
		if add_ones:
			X = np.hstack((np.ones((X.shape[0],1)), X))
		return X@self.a


class KRR(MLModel):
	def __init__(self, kernel):
		'''
		kernel: instance of one of the subclasses of kernels.Kernel
		'''
		super().__init__()
		self.lam = None
		self.kernel = kernel
		self.description = f'Kernel Ridge Regression Model [{kernel}]'

	def train(self, Xtrain, Ytrain, lam=0, add_ones=True, lambda_convergence=1e-3):
		if add_ones:
			Xtrain = np.hstack((np.ones((Xtrain.shape[0],1)), Xtrain))
		m = Xtrain.shape[0]
		#perform regression
		K = self.kernel(Xtrain)
		self.a = Xtrain.T @ np.linalg.inv(np.eye(m)*lam + K) @ Ytrain
		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

		
if __name__ == '__main__':
	# m = KRR(kernels.Polynomial(2))
	m = RR()
	# print(m)
	X = np.random.rand(50, 5)
	# print(X)
	correct_a = np.array([3, 1, 0, .2, 1.5, 0])
	Y = X@correct_a[1:] + correct_a[0]
	# plt.plot(correct_a, linewidth=5)
	# print(Y)
	for lam in np.linspace(0,1,10):
		m.train(X, Y, lam)
		print(m.a)
		# Ypred = m.predict(X)
		# plt.plot(m.a)
		plt.scatter(lam,m.R2train)

	# plt.scatter(Y, Ypred)
	plt.show()


