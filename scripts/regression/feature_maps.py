import numpy as np
import itertools
import matplotlib.pyplot as plt
try:
	import models, MLutils
except:
	import regression.models as models
	import regression.MLutils as MLutils



class FeatureMap:
	def __init__(self):
		NotImplemented

	def __call__(self, X):
		NotImplemented

	def __repr__(self):
		return f'{type(self).__name__} Feature Map' 


class Polynomial(FeatureMap):
	def __init__(self, order):
		assert type(order) is int, f'order must be integer, not {type(order)}'
		assert order >= 1, f'order must be >=1'

		self.order = order

	def __call__(self, X):
		features = []
		for x in X:
			feature = []
			for fs in itertools.combinations_with_replacement(x, self.order):
				p = 1
				for f in fs:
					p *= f
				feature.append(p)
			features.append(feature)
		return np.array(features)

	def __repr__(self):
		return f'{type(self).__name__} Feature Map (d = {self.order})'

	def get_feature_labels(self, labels):
		feature_labels = []
		for fs in itertools.combinations_with_replacement(labels, self.order):
			feature_labels.append(' * '.join(fs))
		return feature_labels

	def single(self, x):
		feature = []
		for fs in itertools.combinations_with_replacement(x, self.order):
			p = 1
			for f in fs:
				p *= f
			feature.append(p)
		return feature


class CombinedFeatures:
	def __init__(self, featuremaps):
		self.featuremaps = featuremaps

	def __call__(self, X):
		features = []
		for x in X:
			feature = []
			for fm in self.featuremaps:
				feature += fm.single(x)
			features.append(feature)
		return np.array(features)

	def __repr__(self):
		fmstrings = []
		for fm in self.featuremaps:
			if type(fm) is Polynomial:
				fmstrings.append(f'Polynomial({fm.order})')
		return f'{type(self).__name__} [{" + ".join(fmstrings)}]'

	def get_feature_labels(self, labels):
		feature_labels = []
		for fm in self.featuremaps:
			feature_labels += (fm.get_feature_labels(labels))
		return feature_labels

	

if __name__ == '__main__':
	cfm = CombinedFeatures([Polynomial(1), Polynomial(2), Polynomial(3)])
	m = 10000
	X = np.random.normal(0, 2, m).reshape(-1,1)
	Y = 1/3*X**3 + X**2
	plt.scatter(X, Y)
	plt.show()
	X = cfm(X) #featurize the input
	X = MLutils.add_ones(X)
	(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, 0.8)

	ML = models.LR()
	ML.train(Xtrain, Ytrain)

	plt.scatter(ML.predict(Xtrain), Ytrain, c='b')
	plt.scatter(ML.predict(Xtest), Ytest, c='r')

	print(ML.a)


	plt.show()

