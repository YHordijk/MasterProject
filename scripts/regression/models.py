import numpy as np 
import matplotlib.pyplot as plt
# import data_gen
# import sklearn.kernel_ridge
import itertools
import utility
try:
	import regression.kernels as kernels
	import regression.MLutils as MLutils
except:
	import kernels
	import MLutils

import sklearn.linear_model as sklm
import sklearn.kernel_ridge as skkr
import sklearn.model_selection as skms
import sklearn.gaussian_process as skgs
import skopt, skopt.space

np.seterr(under='ignore')



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

	def multi_train(self, X, Y, f=0.5, N=50, lam=1, silent=False):
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
		R2 = MLutils.correlation(Ytest, Ypred)**2
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
		w = np.linalg.inv(np.eye(n)*lam + kernels.Gram(Xtrain)) @ Ytrain
		self.a = Xtrain.T @ w
		self.R2train = MLutils.correlation(Ytrain, Xtrain@self.a)**2
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
			# plt.plot(self.a, alpha=1/75, c='k')
		# plt.ylabel('Coefficient')
		# plt.xlabel('Feature index')
		

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

		# plt.plot(self.a, c='r')

		# plt.show()

		MLutils.detect_coefficient_type(coeffs)


	def predict(self, X):
		return X@self.a


class SK_RR(RR):
	def __init__(self):
		super().__init__()
		self.description = f'SKLearn Ridge Regression Model'

	def train(self, Xtrain, Ytrain, lam=1):
		self.model = sklm.Ridge(alpha=lam)
		self.model.fit(Xtrain, Ytrain)
		self.R2train = self.model.score(Xtrain, Ytrain)
		self.trained = True

	def predict(self, X):
		return self.model.predict(X)

	def optimize(self, X, Y, iters=100):
		def on_step(r):
			utility.loading_bar(len(r['x_iters']), iters)

		self.model = sklm.Ridge()
		space = {'alpha':skopt.space.Real(0,10)}

		search = skopt.BayesSearchCV(self.model, space, n_iter=iters, cv=3, verbose=0, scoring='r2')
		print(f'Optimizing model {self.model}')
		search.fit(X, Y, callback=on_step)
		print(f'found best parameters with score {search.best_score_}:')
		for p, v in search.best_params_.items():
			print(p, v)
		self.model = sklm.Ridge(**search.best_params_)
		self.model.fit(X,Y)


class KRR(MLModel):
	def __init__(self, kernel):
		'''
		kernel: instance of one of the subclasses of kernels.Kernel
		'''
		super().__init__()
		self.lam = None
		self.kernel = kernel
		self.description = f'Kernel Ridge Regression Model [{kernel}]'

	def train(self, Xtrain, Ytrain, lam=1):
		m = Xtrain.shape[0]
		K = self.kernel(Xtrain)
		self.w = np.linalg.inv(K + np.eye(m)*lam) @ Ytrain
		self.a = Xtrain.T @ self.w
		self.R2train = MLutils.correlation(Ytrain, Xtrain@self.a)
		self.trained = True

	def predict(self, X):
		# K = self.kernel(X)
		# Ypred = K.T @ self.w
		Ypred = X@self.a
		return Ypred


class SK_KRR(KRR):
	def __init__(self, kernel, order=None, gamma=None, coef0=None):
		'''
		kernel: one of [‘additive_chi2’, ‘chi2’, ‘linear’, ‘poly’, ‘polynomial’, ‘rbf’, ‘laplacian’, ‘sigmoid’, ‘cosine’]
		'''
		self.lam = None
		self.kernel = kernel
		self._order = order
		self._gamma = gamma
		self._coef0 = coef0
		self.description = f'SKLearn Kernel Ridge Regression Model [{kernel}]'

	def train(self, Xtrain, Ytrain, lam=1):
		self.model = skkr.KernelRidge(kernel=self.kernel, alpha=lam, gamma=self._gamma, degree=self._order, coef0=self._coef0)
		self.model.fit(Xtrain, Ytrain)
		self.R2train = self.model.score(Xtrain, Ytrain)
		self.trained = True

	def predict(self, X):
		return self.model.predict(X)

	def cross_validate(self, X, Y, cv=3):
		self.model = skkr.KernelRidge(kernel=self.kernel)
		cv_res = skms.cross_validate(self.model,X,Y,cv=cv,return_estimator=True)
		smax = 0
		for model, score in zip(cv_res['estimator'], cv_res['test_score']):
			print(model.get_params())
			if score > smax:
				smax = score
				best_model = model
		self.model = best_model


	def optimize(self, X, Y, iters=100):
		def on_step(r):
			utility.loading_bar(len(r['x_iters']), iters)

		self.model = skkr.KernelRidge(kernel=self.kernel)
		if self.kernel in ['rbf', 'laplacian']:
			space = {'alpha': skopt.space.Real(0, 5),
					 'gamma': skopt.space.Real(1e-10, 10)}
		elif self.kernel == 'poly':
			space = {'alpha': skopt.space.Real(0, 2),
					 'gamma': skopt.space.Real(0, 2),
					 'coef0': skopt.space.Real(0, 10),
					 'degree': skopt.space.Integer(2,4)}
		elif self.kernel == 'sigmoid':
			space = {'alpha': skopt.space.Real(0, 2),
					 'gamma': skopt.space.Real(-5, 5),
					 'coef0': skopt.space.Real(-5, 5)}

		search = skopt.BayesSearchCV(self.model, space, n_iter=iters, cv=3, verbose=0, scoring='r2')
		print(f'Optimizing model {self.model} using Bayesian Search')
		search.fit(X, Y, callback=on_step)
		print(f'found best parameters with score {search.best_score_}:')
		for p, v in search.best_params_.items():
			print(p, v)
		self.model = skkr.KernelRidge(kernel=self.kernel, **search.best_params_)
		self.model.fit(X,Y)


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


class SK_GPR(GPR):
	def __init__(self, kernel):
		'''
		kernel: instance of one of the subclasses of kernels.Kernel
		'''
		super().__init__(kernel)
		self.description = f'SK Gaussian Process Regression Model [{kernel}]'

	def train(self, Xtrain, Ytrain, noise=1e-10):
		self.model = skgp.GaussianProcessRegressor(alpha=noise)



		
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