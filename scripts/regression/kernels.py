import numpy as np 
import matplotlib.pyplot as plt
import scipy.spatial.distance as scipy_dist

def Gram(X, X2=None):
	'''
	The Gram matrix is a matrix where each element Gij = <Xi, Xj>
	Also the inner product matrix
	'''
	if X2 is None:
		return X@X.T
	else:
		return X@X2.T


def distance2(X, X2=None):
	'''
	The distance matrix D between all points in X
	Written in terms of the Gram matrix:
	Dij = Gii + Gjj - 2Gij
	'''

	if X2 is None:
		G = Gram(X)
		d = np.diag(G).reshape(-1,1)
		return d + d.T - 2*G
	else:
		return scipy_dist.cdist(X, X2)**2


class Kernel:
	def __init__(self):
		NotImplemented

	def __call__(self, X):
		NotImplemented

	def __repr__(self):
		return f'{type(self).__name__} Kernel' 


class Linear(Kernel):
	def __call__(self, X):
		return Gram(X)


class Polynomial(Kernel):
	def __init__(self, order):
		assert type(order) is int, f'order must be integer, not {type(order)}'
		assert order >= 1, f'order must be >=1'

		self.order = order

	def __call__(self, X, X2=None):
		'''
		For each each element i, j in K we calculate (<Xi, Xj> + 1)^d.
		In numpy vectorized form this becomes (XX^T + 1)^d, where XX^T
		is the Gram matrix and d is the order of the polynomial kernel
		'''
		if X2 is None:
			return (X@X.T + 1)**self.order
		else:
			return (X@X2.T + 1)**self.order

	def __repr__(self):
		return super().__repr__() + f' (d = {self.order})'


class Gaussian(Kernel):
	def __init__(self, sigma2):
		assert type(sigma2) in [float, int], f'sigma2 must be a floating point value, not {type(sigma2)}'
		assert sigma2 != 0, 'sigma2 cannot be zero'

		self.sigma2 = sigma2

	def __call__(self, X, X2=None):
		'''
		For each element i, j in K we calculate exp(-||Xi - Xj||^2/(2sigma^2))
		In numpy vector notation this becomes exp(D2/(2sigma^2)
		D is the distance matrix here (Gii + Gjj -2Gij)
		'''
		D2 = distance2(X, X2)
		K = self.normalization() * np.exp(-D2/(2*self.sigma2))
		return K

	def normalization(self):
		return (2*np.pi*self.sigma2)**-.5

	def __repr__(self):
		return super().__repr__() + f' (s^2 = {self.sigma2})'



if __name__ == '__main__':
	X = np.linspace(-1, 1, 600).reshape(-1,1)
	kernel = Gaussian(1)
	K = kernel(X)

	plt.title(kernel)
	plt.imshow(K, origin='lower', extent=(X.min(),X.max(),X.min(),X.max()))
	plt.ylabel(r'$X_i$')
	plt.xlabel(r'$X_j$')
	plt.colorbar()
	plt.show()