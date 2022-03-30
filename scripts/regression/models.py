import numpy as np 
import matplotlib.pyplot as plt
# import data_gen
# import sklearn.kernel_ridge
import itertools
try:
	import regression.kernels as kernels
	import regression.MLutils as MLutils
except:
	import kernels
	import MLutils



class MLModel:
	def __init__(self):
		self.a = np.array([])
		self.description = 'Template Machine Learning Model'
		self.trained = False

	def __repr__(self):
		s = '' 
		s += self.description
		# s += f'This model has{" not"*(not self.trained)} been trained'
		return s

	def train(self, Xtrain, Ytrain):
		'''
		Method that takes training X (m x n) and Y (m x 1) arrays
		and trains the model to set the parameter "a"
		'''
		self.trained = True
		NotImplemented

	def multi_train(self, X, Y, f, n):
		'''
		Method that trains multiple times and then averages the coefficients
		'''
		NotImplemented

	def evaluate(self, Xtest, Ytest):
		'''
		Method that takes test X (m x n) and Y (m x 1) arrays
		and evaluates the model, that is predict the Ytest using
		trained parameters. Then it will calculate the R2 value
		'''
		Ypred = self.predict(Xtest)
		R2 = MLutils.correlation(Ytest, Ypred)
		return R2

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
		self.R2train = MLutils.correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X):
		return X@self.a


class RR(MLModel):
	def __init__(self):
		super().__init__()
		self.lam = None
		self.description = f'Ridge Regression Model'

	def train(self, Xtrain, Ytrain, lam=1):
		n, m = Xtrain.shape
		#standard way
		# self.a = np.linalg.inv(Xtrain.T@Xtrain + np.eye(n)*lam) @ Xtrain.T @ Ytrain
		#KRR way
		w = np.linalg.inv(np.eye(n)*lam + kernels.Gram(Xtrain)) @ Y
		self.a = Xtrain.T @ w
		self.R2train = MLutils.correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def multi_train(self, X, Y, f=0.5, N=50, lam=1, silent=False):
		assert type(N) is int, f'N must be of type int, not {type(N)}'
		assert N > 0, 'N must be larger than 0'
		assert type(lam) in [int, float], f'lam must be numeric type, not {type(lam)}'
		assert lam > 0, f'lam must be larger than 0'
		assert X.shape[0] == Y.shape[0], f'Dimensions of X ({X.shape[0]}) and Y ({Y.shape[0]}) do not match'
		assert 0 < f < 1, f'f must be between 0 and 1'

		if not silent:
			print(f'Training {N} times with split={f*100}% and lambda={lam}')

		n, m = X.shape

		R2s_train = []
		R2s_test = []
		coeffs = []
		for i in range(N):
			(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, f)
			self.train(Xtrain, Ytrain, lam)
			R2s_train.append(self.R2train)
			R2s_test.append(self.evaluate(Xtest, Ytest))
			coeffs.append(self.a)

		R2s_train_mean = np.mean(R2s_train)
		R2s_train_std  = np.std(R2s_train)

		R2s_test_mean  = np.mean(R2s_test)
		R2s_test_std   = np.std(R2s_test)

		if not silent:
			print(f'Training complete!')
			print(f'\tR2 (training) = {R2s_train_mean: .5f} ± {R2s_train_std: .5f}')
			print(f'\tR2 (test)     = {R2s_test_mean: .5f} ± {R2s_test_std: .5f}')

		self.R2train = R2s_train_mean
		self.R2train_std = R2s_train_std
		self.R2test = R2s_test_mean
		self.R2test_std = R2s_test_std
		self.a = np.mean(np.array(coeffs), axis=0)


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
		self.R2train = MLutils.correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X):
		K = self.kernel(X)
		Ypred = K.T @ self.w
		return Ypred


class GPR(MLModel):
	def __init__(self, kernel):
		'''
		kernel: instance of one of the subclasses of kernels.Kernel
		'''
		super().__init__()
		self.lam = None
		self.kernel = kernel
		self.description = f'Gaussian Process Regression Model [{kernel}]'

	def train(self, Xtrain, Ytrain, noise=0):
		'''
		Training a GPR is simply precomputing the kernel
		'''
		print(Xtrain)
		plt.imshow(self.kernel(Xtrain))
		plt.show()
		self.invcov = np.linalg.inv(self.kernel(Xtrain) + np.eye(Xtrain.shape[0]) * (noise + np.finfo(float).eps)) #(K(X,X) + noise*I)^{-1}
		self.Xtrain = Xtrain
		self.Ytrain = Ytrain

	def predict(self, X, return_var=False):
		'''
		Posterior has mean mu = K(X,X*) (K(X,X) + noise*I)^{-1} Y
		and covariance matrix sigma = K(X*,X*) - K(X,X*) (K(X,X) + noise*I)^{-1} K(X*, X)
		and the diagonal of sigma is the variance vector
		'''
		#different covariance matrices
		Kns = self.kernel(self.Xtrain, X) #K(X,X*)
		Kss = self.kernel(X, X) #K(X*,X*)
		Ksn = self.kernel(X, self.Xtrain) #K(X*,X)

		mean = Kns.T @ self.invcov @ self.Ytrain
		if return_var:
			cov = Kss - Kns.T @ self.invcov @ Ksn.T
			var = np.diag(cov).reshape(-1,1)
			return mean, var
		return mean

	def evaluate(self, Xtest, Ytest):
		'''
		Method that takes test X (m x n) and Y (m x 1) arrays
		and evaluates the model, that is predict the Ytest using
		trained parameters. Then it will calculate the R2 value
		'''
		Ypred = self.predict(Xtest)[0]
		R2 = MLutils.correlation(Ytest, Ypred)
		return R2


		
if __name__ == '__main__':
	model = KRR(kernels.Polynomial(3))
	# model = RR()

	print(model)
	X = np.random.rand(100, 300)
	X = MLutils.add_ones(X)
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