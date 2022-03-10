import numpy as np 
import matplotlib.pyplot as plt



def Gram(X):
	'''
	The Gram matrix is a matrix where each element Gij = <Xi, Xj>
	Also the inner product matrix
	'''
	return X@X.T

def distance(X):
	'''
	The distance matrix D between all points in X
	Written in terms of the Gram matrix:
	Dij = Gii + Gjj - 2Gij
	'''
	G = Gram(X)
	d = np.diag(G).reshape(-1,1)
	return d + d.T - 2*G


class Kernel:
	def __init__(self):
		NotImplemented

	def __call__(self, X):
		NotImplemented

	def __repr__(self):
		return f'{type(self).__name__} Kernel' 


class Polynomial(Kernel):
	def __init__(self, order):
		assert type(order) is int, f'order must be integer, not {type(order)}'
		assert order >= 1, f'order must be >=1'
		self.order = order

	def __call__(self, X):
		'''
		For each each element i, j in K we calculate (<Xi, Xj> + 1)^d.
		In numpy vectorized form this becomes (XX^T + 1)^d, where XX^T
		is the Gram matrix and d is the order of the polynomial kernel
		'''
		K = (X@X.T + 1)**self.order
		return K

	def __repr__(self):
		return super().__repr__() + f' (d = {self.order})'

class RadialBasisFunction(Kernel):
	def __init__(self, sigma2):
		assert type(sigma2) in [float, int], f'sigma2 must be a floating point value, not {type(sigma2)}'
		assert sigma2 != 0, 'sigma2 cannot be zero'
		self.sigma2 = sigma2

	def __call__(self, X):
		'''
		For each element i, j in K we calculate exp(-||Xi - Xj||^2/(2sigma^2))
		In numpy vector notation this becomes exp(D/(2sigma^2)
		D is the distance matrix here (Gii + Gjj -2Gij)
		'''
		D = distance(X)
		K = np.exp(-D/(2*self.sigma2))
		return K

	def __repr__(self):
		return super().__repr__() + f' (s^2 = {self.sigma2})'






if __name__ == '__main__':
	X = np.linspace(-1, 1, 600).reshape(-1,1)
	Kp = Polynomial(2)
	K = Kp(X)
	plt.imshow(K, origin='low')
	plt.show()