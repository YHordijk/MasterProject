import ML_params as MLparams
import regression.models as MLmodels
import regression.MLutils as MLutils
import regression.feature_maps as feature_maps
import regression.kernels as kernels
import regression.models as models
import matplotlib.pyplot as plt
import os
import numpy as np
from utility import hartree2kcalmol as h2k
from utility import hartree2eV as h2e
import moviepy.editor as mvp


def make_gif(out, frames, fps=14):
	# frames = [frames_folder + f for f in os.listdir(frames_folder)]
	clip = mvp.ImageSequenceClip(frames, fps=fps)
	clip.write_gif(out, fps=fps)


def main(model='GPR'):
	if model == 'RR':
		ML = MLmodels.RR()
		#get all data:
		parameters = MLparams.define_parameters()
		X = MLparams.get_all_columns()
		Y = MLparams.get_column('Eact')

		#remove outliers
		oi = MLparams.outlier_idx(X,2.5) + MLparams.outlier_idx(Y,2.5)
		X = np.delete(X, oi, axis=0)
		Y = np.delete(Y, oi, axis=0)
		print(f'Removed {len(oi)} outliers')

		cfm = feature_maps.Polynomial(3)
		X = MLutils.add_ones(X)
		X = cfm(X)
		X = MLutils.normalize(X)
		
		print(f'There are {X.shape[1]} features per datapoint using {cfm}')

		#split into train and test sets
		train, test = MLutils.split_data(X, Y, .8)
		Xtrain, Ytrain = train
		Xtest, Ytest = test

		#train the model
		ML.train(Xtrain,Ytrain, 1)

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

		


		R2s_training = []
		R2s_evaluate = []
		lams = np.linspace(0.01, .2, 30)
		for lam in lams:
			R2_train = []
			R2_eval = []
			for _ in range(15):
				(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .8)
				#train the model
				ML.train(Xtrain,Ytrain, lam)
				R2_train.append(ML.R2train)
				R2_eval.append(ML.evaluate(Xtest, Ytest))
			R2s_training.append(R2_train)
			R2s_evaluate.append(R2_eval)


		plt.plot(lams, np.mean(R2s_training, axis=1), c='b', label='Training')
		plt.plot(lams, np.mean(R2s_evaluate, axis=1), c='r', label='Test')
		plt.scatter([[l for _ in range(15)] for l in lams], R2s_training, c='b', alpha=.25)
		plt.scatter([[l for _ in range(15)] for l in lams], R2s_evaluate, c='r', alpha=.25)
		plt.ylabel(r'$r^2$')
		plt.xlabel(r'$\lambda$')
		plt.legend()
		plt.show()


		R2s_training = []
		R2s_evaluate = []
		fs = np.linspace(0.1, .90, 30)
		for f in fs:
			R2_train = []
			R2_eval = []
			for _ in range(15):
				(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, f)
				#train the model
				ML.train(Xtrain,Ytrain, 1)
				R2_train.append(ML.R2train)
				R2_eval.append(ML.evaluate(Xtest, Ytest))
			R2s_training.append(R2_train)
			R2s_evaluate.append(R2_eval)


		plt.plot(fs, np.mean(R2s_training, axis=1), c='b', label='Training')
		plt.plot(fs, np.mean(R2s_evaluate, axis=1), c='r', label='Test')
		plt.scatter([[l for _ in range(15)] for l in fs], R2s_training, c='b', alpha=.25)
		plt.scatter([[l for _ in range(15)] for l in fs], R2s_evaluate, c='r', alpha=.25)
		plt.ylabel(r'$r^2$')
		plt.xlabel(r'$f$')
		plt.legend()
		plt.gca().set_ylim((0,1))
		plt.show()

		ML.multi_train(X, Y, f=0.5, N=150, lam=1)
		sortidx = np.argsort(np.abs(ML.a.flatten()))
		labels = cfm.get_feature_labels([''] + parameters)
		print('  Weight | Feature')
		for i in sortidx[::-1]:
			print(f'{float(ML.a[i]): .5f} | {labels[i]}')


	elif model == 'GPR':
		parameters = MLparams.define_parameters()
		X = MLparams.get_all_columns()
		Y = MLparams.get_column('Eact')
		(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .7)

		print(Xtrain.shape, Ytrain.shape, Xtest.shape, Ytest.shape)

		vary_noise = False
		vary_sigma = True
		vary_both  = True

		out_path_idx = 1

		if vary_noise:
			out_path = f'analyses/ML_plots/GPR/vary_noise_{out_path_idx}'
			if not os.path.exists(out_path):
				os.mkdir(out_path)
			ns = np.linspace(-8, 2, 100)
			frames = []
			R2s = []
			for i, n in enumerate(ns):
				noise = 10**n
				kernel_sigma2 = 2.31e1
				kernel = kernels.Gaussian(kernel_sigma2)
				model = models.GPR(kernel)
				model.train(Xtrain, Ytrain, noise)
				
				plt.subplots(figsize=(10,5))
				plt.subplot(1,2,1)
				Ypred = model.predict(Xtrain)
				plt.scatter(Ytrain, Ypred, label='Train', alpha=.2, c='k')

				Ypred = model.predict(Xtest)
				plt.scatter(Ytest, Ypred, label='Test', c='r')
				R2s.append(MLutils.correlation(Ypred, Ytest)**2)
				
				plt.suptitle(model.__repr__() + f'\nNoise = {noise}')
				plt.xlabel(r'$Y$')
				plt.ylabel(r'$Y_{pred}$')
				plt.gca().set_xlim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
				plt.gca().set_ylim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
				plt.plot((Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), ( Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), 'k--')
				plt.legend()

				plt.subplot(1,2,2)
				plt.plot(ns[:i+1], R2s, c='b')
				plt.xlabel(r'log$_{10}(\sigma_n^2)$')
				plt.ylabel(r'$r^2$ correlation')
				plt.gca().set_xlim(ns.min(), ns.max())
				plt.gca().set_ylim(0,1)

				plt.savefig(os.path.join(out_path, f'vary_noise{i}.png'))
				plt.close()
				
				frames.append(os.path.join(out_path, f'vary_noise{i}.png'))

			make_gif(os.path.join(out_path, 'noise_anim.gif'), frames, fps=14)

		if vary_sigma:
			out_path = f'analyses/ML_plots/GPR/vary_sigma_{out_path_idx}'
			if not os.path.exists(out_path):
				os.mkdir(out_path)
			ss = np.linspace(-2, 4, 100)
			frames = []
			R2s = []
			for i, s in enumerate(ss):
				noise = 2.78e1
				kernel_sigma2 = float(10**s)
				kernel = kernels.Gaussian(kernel_sigma2)
				model = models.GPR(kernel)
				model.train(Xtrain, Ytrain, noise)

				plt.subplots(figsize=(10,5))
				plt.subplot(1,2,1)
				Ypred, var = model.predict(Xtrain)
				plt.scatter(Ytrain, Ypred, label='Train', alpha=.2, c='k')

				Ypred, var = model.predict(Xtest)
				plt.scatter(Ytest, Ypred, label='Test', c='r')
				R2s.append(MLutils.correlation(Ypred, Ytest)**2)
				
				plt.suptitle(model.__repr__() + f'\nNoise = {noise}')
				plt.xlabel(r'$Y$')
				plt.ylabel(r'$Y_{pred}$')
				plt.gca().set_xlim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
				plt.gca().set_ylim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
				plt.plot((Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), ( Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), 'k--')
				plt.legend()

				plt.subplot(1,2,2)
				plt.plot(ss[:i+1], R2s, c='b')
				plt.xlabel(r'log$_{10}(\sigma_n^2)$')
				plt.ylabel(r'$r^2$ correlation')
				plt.gca().set_xlim(ss.min(), ss.max())
				plt.gca().set_ylim(0,1)
				plt.savefig(os.path.join(out_path, f'vary_sigma{i}.png'))
				plt.close()
				
				frames.append(os.path.join(out_path, f'vary_sigma{i}.png'))

			make_gif(os.path.join(out_path, 'sigma_anim.gif'), frames, fps=14)

		if vary_both:
			out_path = f'analyses/ML_plots/GPR/vary_both_{out_path_idx}'
			if not os.path.exists(out_path):
				os.mkdir(out_path)
			ns = np.linspace(-10, 3, 100)
			ss = np.linspace(-1, 6, 100)
			R2s = np.zeros((ss.size, ns.size))
			R2best = (0,0)
			for i, n in enumerate(ns):
				for j, s in enumerate(ss):
					noise = 10**n
					kernel_sigma2 = float(10**s)
					kernel = kernels.Gaussian(kernel_sigma2)
					model = models.GPR(kernel)
					model.train(Xtrain, Ytrain, noise)
					Ypred, var = model.predict(Xtest)

					R2 = MLutils.correlation(Ypred, Ytest)**2
					if R2 >= R2s.max():
						print('new best (R2)', R2, noise, kernel_sigma2)
						R2best = (noise, kernel_sigma2)
						R2bestidx = (n, s)
					R2s[j,i] = R2

			extent = ns.min(), ns.max(), ss.min(), ss.max()
			plt.imshow(R2s, extent=extent, origin='lower', interpolation='antialiased', aspect='auto')
			plt.xlabel(r'log$_{10}(\sigma^2_n)$')
			plt.ylabel(r'log$_{10}(\sigma^2_K)$')
			plt.annotate(rf'Best fit, $r^2={R2s.max()}$' + '\n' + rf'noise$={R2best[0]:.5E}$, $\sigma^2={R2best[1]:.5E}$', 
					xy=R2bestidx, xytext=(0.3,.8), textcoords='figure fraction', arrowprops=dict(arrowstyle="->"), color='white',
					annotation_clip=False)
			plt.scatter(*R2bestidx, c='k')
			plt.savefig(os.path.join(out_path, 'heatmap_r2.png'))
			plt.close()

			noise = R2best[0]
			kernel_sigma2 = R2best[1]
			kernel = kernels.Gaussian(kernel_sigma2)
			model = models.GPR(kernel)
			model.train(Xtrain, Ytrain, noise)

			plt.figure(figsize=(9,9))

			Ypred, var = model.predict(Xtrain)
			plt.scatter(Ytrain, Ypred, label='Train', alpha=.2, c='k')
			Ypred, var = model.predict(Xtest)
			plt.scatter(Ytest, Ypred, label='Test', c='r')

			plt.suptitle(model.__repr__() + f'\nNoise = {noise}')
			plt.xlabel(r'$Y$')
			plt.ylabel(r'$Y_{pred}$')
			plt.gca().set_xlim(Ytest.min()-1, Ytest.max()+1)
			plt.gca().set_ylim(Ytest.min()-1, Ytest.max()+1)
			plt.plot((Ytest.max()+1, Ytest.min()-1), (Ytest.max()+1, Ytest.min()-1), 'k--')
			plt.legend()

			plt.savefig(os.path.join(out_path, f'best_fit_r2.png'))
			plt.close()