import numpy as np
import scipy.stats as scipy_stats
import matplotlib.pyplot as plt

def correlation(y, y_pred):
	return scipy_stats.pearsonr(y.flatten(), y_pred.flatten())[0]


def RSS(y, y_pred):
	''' 
	Calculates the residual sum of squares
	'''
	return np.sum((y-y_pred)**2)/y.size

def add_ones(X):
	return np.hstack((np.ones((X.shape[0],1)), X))

def prepare_data(X, normalize=True, center=True):
	Xmin = np.min(X, axis=0)
	Xmax = np.max(X, axis=0)
	if normalize:
		X = (X - Xmin)/(Xmax - Xmin)

	mu   = np.mean(X, axis=0)
	if center:
		X = X - mu

	norm_data = {'normalized':normalize, 'centered':center, 'min': Xmin, 'max':Xmax, 'mu':mu}
	return X, norm_data

def unprepare_data(X, norm_data):
	if norm_data['centered']:
		X = X + norm_data['mu']
	if norm_data['normalized']:
		X = X * (norm_data['max'] - norm_data['min']) + norm_data['min']

	return X




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
		plt.plot(coeffs.T, alpha=10/m, c='k', label='features')
		plt.plot(mean,c='r', label='feature mean')
		plt.plot(std,c='b',label='feature std')
		plt.show()


