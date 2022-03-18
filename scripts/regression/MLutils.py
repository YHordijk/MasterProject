import numpy as np
import matplotlib.pyplot as plt

def correlation(y, y_pred):
	m = y.size 
	ymean = np.mean(y)
	vary = np.sum((y-ymean)**2)
	vary_pred = np.sum((y-y_pred)**2)
	return 1 - vary_pred/vary


def add_ones(X):
	return np.hstack((np.ones((X.shape[0],1)), X))


def normalize(X):
	'''
	Normalizes all data in X by column
	'''
	return (X - np.min(X, axis=0))/np.max(X, axis=0)


def split_data(X, Y, f, randomize=True):
	'''
	Function that will split data into a training and a test set
	sizes of datasets are controlled by f \in [0,1]
	training: f*100 %
	test: (1-f)*100 %
	'''

	m = X.shape[0]
	assert Y.shape[0] == m, f'X and Y have different number of datapoints X: {X.shape[0]}, Y: {Y.shape[0]}'
	assert 0 <= f <= 1, 'f must be between 0 and 1'

	if randomize:
		idx = np.arange(m)
		np.random.shuffle(idx)
		X = X[idx]
		Y = Y[idx]

	trainM = int(f*m)
	testM  = int((1-f)*m)
	if trainM + testM < m:
		testM += 1

	return (X[:trainM,:], Y[:trainM,:]), (X[trainM:,:], Y[trainM:,:])



def detect_coefficient_type(coeffs, plot=True):
	'''
	This function will determine for every coefficient its type
	type 1 - coefficient consistently close to zero
	type 2 - coefficient consistently large and either positive or negative
	type 3 - coefficient consistently large but not consistently positive or negative

	Takes m x n matrix with m coefficient sets and n coefficients in each set
		list of arrays or list of lists also allowed
	returns vector with n elements denoting the type (1, 2 or 3)
	'''

	coeffs = np.array(coeffs)
	# print(coeffs.shape)
	m, n = coeffs.shape[0], coeffs.shape[1]
	coeffs =coeffs.reshape(m,n)
	mean = coeffs.mean(axis=0)
	std = coeffs.std(axis=0)
	fraction_under_zero = np.count_nonzero(coeffs < 0, axis=0)/n
	print(fraction_under_zero)

	if plot:
		plt.xlabel('Feature index')
		plt.ylabel('Coefficient')
		plt.plot(coeffs.T, alpha=10/m, c='k')
		plt.plot(mean,c='r')
		plt.plot(std,c='b')
		plt.show()


