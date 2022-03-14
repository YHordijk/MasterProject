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


def add_ones(X):
	return np.hstack((np.ones((X.shape[0],1)), X))


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

	def train(self, Xtrain, Ytrain):
		#perform regression
		self.a = np.linalg.inv(Xtrain.T@Xtrain) @ Xtrain.T @ Ytrain
		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def evaluate(self, Xtest, Ytest):
		Ypred = Xtest@self.a
		R2 = correlation(Ytest, Ypred)
		return R2

	def predict(self, X):
		return X@self.a


class RR(MLModel):
	def __init__(self):
		super().__init__()
		self.lam = None
		self.description = f'Ridge Regression Model'

	def train(self, Xtrain, Ytrain, lam=1):
		n = Xtrain.shape[1]
		# G = kernels.Gram(Xtrain)
		# w = np.linalg.inv(G + np.eye(n)*lam) @ Ytrain
		# self.a = Xtrain.T @ w
		self.a = np.linalg.inv(Xtrain.T@Xtrain + np.eye(n)*lam) @ Xtrain.T @ Ytrain

		# Xtrain.T@Xtrain @ self.a + np.eye(n)*lam @ self.a = Xtrain.T @ Ytrain
		# np.eye(n)*lam @ self.a = Xtrain.T @ Ytrain - Xtrain.T@Xtrain @ self.a
		# np.eye(n)*lam @ self.a = Xtrain.T @ (Ytrain - Xtrain @ self.a)
		# self.a = np.linalg.inv(np.eye(n)*lam) @ Xtrain.T @ (Ytrain - Xtrain @ self.a)

		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X):
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

	def train(self, Xtrain, Ytrain, lam=1, lambda_convergence=1e-3):
		m = Xtrain.shape[0]
		K = self.kernel(Xtrain)
		self.w = np.linalg.inv(K + np.eye(m)*lam) @ Ytrain
		self.a = Xtrain.T @ self.w
		self.R2train = correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X):
		K = self.kernel(X)
		Ypred = K.T @ self.w
		return Ypred

		
if __name__ == '__main__':
	# m = KRR(kernels.Polynomial(3))
	model = RR()

	print(model)
	X = np.random.rand(100, 300)
	X = add_ones(X)
	# a0 = np.array([3, 1, 0, .2, 1.5, 0])
	a0 = np.random.rand(X.shape[1])*2
	zeroidx = np.random.choice(np.arange(a0.size), replace=False, size=150)
	a0[zeroidx] = 0
	Y = X@a0

	cmap = plt.get_cmap('hot')
	

	R2s = []
	lams = np.linspace(0.00001,100,100)
	# plt.plot(a0, '--k', linewidth=4)
	for i, lam in enumerate(lams):
		model.train(X, Y, lam=lam)
		R2s.append(model.R2train)
		plt.plot(np.abs(model.a-a0), color=cmap(i/100))


	plt.show()
	plt.plot(lams, R2s)
	plt.show()