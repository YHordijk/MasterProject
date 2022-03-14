import numpy as np


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


