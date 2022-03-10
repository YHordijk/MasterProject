import ML_params as MLparams
import regression.models as MLmodels
import matplotlib.pyplot as plt
import numpy as np
from utility import hartree2kcalmol as h2k
from utility import hartree2eV as h2e


def main():
	LR = MLmodels.LR()
	#get all data:
	parameters = MLparams.define_parameters()
	X = MLparams.get_all_columns()
	Y = MLparams.get_column('Eact')

	#remove outliers
	oi = MLparams.outlier_idx(X) + MLparams.outlier_idx(Y)
	X = np.delete(X, oi, axis=0)
	Y = np.delete(Y, oi, axis=0)
	print(f'Removed {len(oi)} outliers')

	#split into train and test sets
	train, test = MLparams.split_data(X, Y, .8)
	Xtrain, Ytrain = train
	Xtest, Ytest = test

	#train the model
	LR.train(Xtrain,Ytrain)

	#evaluate the model
	R2test = LR.evaluate(Xtest, Ytest)
	Ypredtrain = LR.predict(Xtrain)
	Ypredtest = LR.predict(Xtest)

	print(f'Training with {Xtrain.shape[0]} datapoints: R2={LR.R2train}')
	print(f'Testing with  {Xtest.shape[0]} datapoints: R2 = {R2test}')
	plt.title(r'Prediction of $\Delta E^‡$ using Linear Regression model')
	plt.scatter(Ytrain*h2k(1), Ypredtrain*h2k(1), color='blue', label=rf'Training, $r^2$={LR.R2train:.3f}')
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



	# X = MLparams.get_column('LUMO_energy')
	# Y = MLparams.get_column('Eact')
	# oi = oi + MLparams.outlier_idx(Y)
	# X = np.delete(X, oi, axis=0)
	# Y = np.delete(Y, oi, axis=0)

	# LR = MLmodels.linear_regression(X, Y)
	# print(LR)


	

	# print(Xtrain.shape)
	# print(Xtest.shape)
	# print(Xtest.shape[0]+Xtrain.shape[0])



	# for i, E in enumerate(['Eact', 'Gact']):
	# 	plt.subplot(1,2,i+1)
	# 	Y = MLparams.get_column(E)
	# 	X = MLparams.get_column('LUMO_energy')

	# 	# oi = MLparams.outlier_idx(Y, 3) + MLparams.outlier_idx(X, 3)
	# 	# X = np.delete(X, oi).reshape(-1,1) * h2e(1)
	# 	# Y = np.delete(Y, oi).reshape(-1,1) * h2k(1)

	# 	LR = MLmodels.linear_regression(X, Y)
	# 	Xpa = X.min() - (X.max()-X.min())*.1
	# 	Xpb = X.max() + (X.max()-X.min())*.1
	# 	Ypred = np.array([[1,Xpa], [1,Xpb]])@LR['a']

	# 	plt.scatter(X, Y, alpha=.5)
	# 	plt.plot([Xpa, Xpb], Ypred, '--r', label=rf'LR ($r^2$={LR["R2"]:.4f})')
	# 	if E == 'Eact':
	# 		plt.ylabel(r'$\Delta E_{act}$ (kcal/mol)')
	# 	else:
	# 		plt.ylabel(r'$\Delta G_{act}$ (kcal/mol)')
	# 	plt.xlabel(r'LUMO energy (eV)')
	# 	plt.legend()
	# plt.show()

