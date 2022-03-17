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
		'''	
		Generate feature space given X. Iterate through X (by row)
		and get the combinations with replacement up to order self.order.
		Put these into list of lists and then return as new array.
		'''
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
		'''
		Generate the feature labels. We do this in the same way as
		generating feature spaces. Given a list of labels (as strings)
		we generate the combinations with replacement up to order self.order.
		Then count for each combination how many times each label occurs
		and then contract it, i.e. use powers (x*y^2) instead of full name (x*y*y)
		'''
		feature_labels = []
		for fs in itertools.combinations_with_replacement(labels, self.order):
			#fs is one combination
			label_dict = {}
			#track how many of each label in fs
			for f in fs:
				if not f in label_dict:
					label_dict[f] = 1
				else:
					label_dict[f] += 1
			#generate the label strings
			label_strings = []
			for l, p in label_dict.items():
				#if l is the empty string skip it (used for baseline)
				if l != '':
					if p > 1:
						#use power notation if l occurs more than once
						label_strings.append(f'{l}^{p}')
					else:
						label_strings.append(f'{l}')
			feature_labels.append('*'.join(label_strings))
		return feature_labels

	def single(self, x):
		'''
		Feature space for a single datapoint, used with CombinedFeatures
		as the iteration through X occurs in the parent.
		'''
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
	cfm = Polynomial(3) #feature map
	ML = models.LR() #ML model

	#generate some data
	m = 500
	X = np.random.normal(0, 3, m).reshape(-1,1)
	Y = 1/3*X**3 + X**2 + np.random.normal(0, 6, m).reshape(-1,1) + .5
	
	#generate the feature space
	phi = MLutils.add_ones(X)
	phi = cfm(phi) 

	#train the model
	ML.train(phi, Y)

	#plot prediction
	Xpred = np.linspace(-10,10,50).reshape(-1,1)
	phipred = MLutils.add_ones(Xpred)
	phipred = cfm(phipred)
	Ypred = ML.predict(phipred)

	plt.plot(Xpred, Ypred, '--k', label=f'Prediction, R2={ML.R2train: .5f}')
	plt.scatter(X, Y, alpha=0.25, label='Input')
	plt.title(f'{ML.description} using {cfm}')
	plt.xlabel(r'$x$')
	plt.ylabel(r'$y$')
	plt.gca().set_xlim(-10, 10)
	plt.gca().set_ylim(-150, 200)
	plt.legend()
	
	#print the predicted formula
	features = cfm.get_feature_labels(['','x'])
	s = []
	for a, f in zip(ML.a, features):
		s.append(f'{float(a):.5f}{f}')
	print('y =', ' + '.join(s))

	plt.show()