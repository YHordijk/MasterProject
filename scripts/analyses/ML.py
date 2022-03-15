import ML_params as MLparams
import regression.models as MLmodels
import regression.MLutils as MLutils
import regression.feature_maps as feature_maps
import matplotlib.pyplot as plt
import numpy as np
from utility import hartree2kcalmol as h2k
from utility import hartree2eV as h2e


def main():
	ML = MLmodels.RR()
	#get all data:
	parameters = MLparams.define_parameters()
	X = MLparams.get_all_columns()
	Y = MLparams.get_column('Eact')

	#remove outliers
	oi = MLparams.outlier_idx(X) + MLparams.outlier_idx(Y)
	X = np.delete(X, oi, axis=0)
	Y = np.delete(Y, oi, axis=0)
	print(f'Removed {len(oi)} outliers')

	pl = feature_maps.Polynomial
	cfm = feature_maps.CombinedFeatures([pl(1), pl(2), pl(3)])
	X = cfm(X)
	# X = MLutils.normalize(X)
	X = MLutils.add_ones(X)
	print(f'There are {X.shape[1]} features per datapoint using {cfm}')

	#split into train and test sets
	train, test = MLutils.split_data(X, Y, .8)
	Xtrain, Ytrain = train
	Xtest, Ytest = test

	#train the model
	ML.train(Xtrain,Ytrain, .1)

	#evaluate the model
	R2test = ML.evaluate(Xtest, Ytest)
	Ypredtrain = ML.predict(Xtrain)
	Ypredtest = ML.predict(Xtest)

	print(f'Training with {Xtrain.shape[0]} datapoints: R2 = {ML.R2train}')
	print(f'Testing with {Xtest.shape[0]} datapoints: R2 = {R2test}')
	plt.title(rf'Prediction of $\Delta E^‡$ using {ML.description}' + f'\nusing {cfm}')
	plt.scatter(Ytrain*h2k(1), Ypredtrain*h2k(1), color='blue', label=rf'Training, $r^2$={ML.R2train:.3f}')
	plt.scatter(Ytest*h2k(1), Ypredtest*h2k(1), color='red', label=rf'Testing, $r^2$={R2test:.3f}')

	minX = min(Ytrain.min(), Ytest.min())*h2k(1)
	maxX = max(Ytrain.max(), Ytest.max())*h2k(1)
	minY = min(Ypredtrain.min(), Ypredtest.min())*h2k(1)
	maxY = max(Ypredtrain.max(), Ypredtest.max())*h2k(1)

	mi = min(minX, minY) *1.1
	ma = max(maxX, maxY) *1.1

	plt.plot((mi, ma), (mi, ma), '--k', linewidth=2)
	plt.xlabel(r'$\Delta E^‡$ (kcal/mol)')
	plt.ylabel(r'$\Delta E^‡_{pred}$ (kcal/mol)')

	plt.gca().set_xlim([mi, ma])
	plt.gca().set_ylim([mi, ma])

	plt.legend()
	plt.show()

	sortidx = np.argsort(np.abs(ML.a.flatten()))
	labels = ['baseline'] + cfm.get_feature_labels(parameters)
	print('  Weight | Feature')
	for i in sortidx[::-1]:
		print(f'{float(ML.a[i]): .5f} | {labels[i]}')


	# R2s_training = []
	# R2s_evaluate = []
	# lams = np.linspace(0.01, .2, 30)
	# for lam in lams:
	# 	R2_train = []
	# 	R2_eval = []
	# 	for _ in range(15):
	# 		(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .8)
	# 		#train the model
	# 		ML.train(Xtrain,Ytrain, lam)
	# 		R2_train.append(ML.R2train)
	# 		R2_eval.append(ML.evaluate(Xtest, Ytest))
	# 	R2s_training.append(R2_train)
	# 	R2s_evaluate.append(R2_eval)


	# plt.plot(lams, np.mean(R2s_training, axis=1), c='b', label='Training')
	# plt.plot(lams, np.mean(R2s_evaluate, axis=1), c='r', label='Test')
	# plt.scatter([[l for _ in range(15)] for l in lams], R2s_training, c='b', alpha=.25)
	# plt.scatter([[l for _ in range(15)] for l in lams], R2s_evaluate, c='r', alpha=.25)
	# plt.ylabel(r'$r^2$')
	# plt.xlabel(r'$\lambda$')
	# plt.legend()
	# plt.show()


	# R2s_training = []
	# R2s_evaluate = []
	# fs = np.linspace(0.1, .90, 30)
	# for f in fs:
	# 	R2_train = []
	# 	R2_eval = []
	# 	for _ in range(15):
	# 		(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, f)
	# 		#train the model
	# 		ML.train(Xtrain,Ytrain, 1)
	# 		R2_train.append(ML.R2train)
	# 		R2_eval.append(ML.evaluate(Xtest, Ytest))
	# 	R2s_training.append(R2_train)
	# 	R2s_evaluate.append(R2_eval)


	# plt.plot(fs, np.mean(R2s_training, axis=1), c='b', label='Training')
	# plt.plot(fs, np.mean(R2s_evaluate, axis=1), c='r', label='Test')
	# plt.scatter([[l for _ in range(15)] for l in fs], R2s_training, c='b', alpha=.25)
	# plt.scatter([[l for _ in range(15)] for l in fs], R2s_evaluate, c='r', alpha=.25)
	# plt.ylabel(r'$r^2$')
	# plt.xlabel(r'$f$')
	# plt.legend()
	# plt.gca().set_ylim((0,1))
	# plt.show()

	ML.multi_train(X, Y, f=0.5, N=50, lam=0.1)
	print(ML.evaluate(X, Y))