import ML_params as MLparams
import regression.models as MLmodels
import regression.MLutils as MLutils
import regression.feature_maps as feature_maps
import regression.kernels as kernels
import regression.models as models
import utility
import matplotlib.pyplot as plt
import os
import numpy as np
from utility import hartree2kcalmol as h2k
from utility import hartree2eV as h2e
import moviepy.editor as mvp

join = os.path.join

def make_gif(out, frames, fps=14):
	# frames = [frames_folder + f for f in os.listdir(frames_folder)]
	clip = mvp.ImageSequenceClip(frames, fps=fps)
	clip.write_gif(out, fps=fps)


def main(model='KRR', animate=True, out_path_idx=2, kernel=''):
	# if model == 'KRR':
	# 	# kernel = kernels.Polynomial(3)
	# 	# ML = MLmodels.KRR(kernel)
	# 	ML = MLmodels.SK_KRR('linear')
	# 	#get all data:
	# 	parameters = MLparams.define_parameters()
	# 	X = MLparams.get_all_columns()
	# 	Y = MLparams.get_column('Eact')

	# 	#remove outliers
	# 	oi = MLparams.outlier_idx(X, 2.5) + MLparams.outlier_idx(Y, 2.5)
	# 	X = np.delete(X, oi, axis=0)
	# 	Y = np.delete(Y, oi, axis=0)
	# 	print(f'Removed {len(oi)} outliers')
	# 	print(f'We have {X.shape[1]} features')

	# 	#split into train and test sets
	# 	train, test = MLutils.split_data(X, Y, .8)
	# 	Xtrain, Ytrain = train
	# 	Xtest, Ytest = test

	# 	#train the model
	# 	ML.train(Xtrain,Ytrain, .0000000001)

	# 	#evaluate the model
	# 	R2test = ML.evaluate(Xtest, Ytest)
	# 	Ypredtrain = ML.predict(Xtrain)
	# 	Ypredtest = ML.predict(Xtest)

	# 	print(f'Training with {Xtrain.shape[0]} datapoints: R2 = {ML.R2train}')
	# 	print(f'Testing with {Xtest.shape[0]} datapoints: R2 = {R2test}')
	# 	plt.title(rf'Prediction of $\Delta E^‡$ using {ML.description}')
	# 	plt.scatter(Ytrain*h2k(1), Ypredtrain*h2k(1), color='blue', label=rf'Training, $r^2$={ML.R2train:.3f}')
	# 	plt.scatter(Ytest*h2k(1), Ypredtest*h2k(1), color='red', label=rf'Testing, $r^2$={R2test:.3f}')

	# 	minX = min(Ytrain.min(), Ytest.min())*h2k(1)
	# 	maxX = max(Ytrain.max(), Ytest.max())*h2k(1)
	# 	minY = min(Ypredtrain.min(), Ypredtest.min())*h2k(1)
	# 	maxY = max(Ypredtrain.max(), Ypredtest.max())*h2k(1)

	# 	mi = min(minX, minY) *1.1
	# 	ma = max(maxX, maxY) *1.1

	# 	plt.plot((mi, ma), (mi, ma), '--k', linewidth=2)
	# 	plt.xlabel(r'$\Delta E^‡$ (kcal/mol)')
	# 	plt.ylabel(r'$\Delta E^‡_{pred}$ (kcal/mol)')

	# 	plt.gca().set_xlim([mi, ma])
	# 	plt.gca().set_ylim([mi, ma])

	# 	plt.legend()
	# 	plt.show()


	if model in ['RR', 'KRR']:
		if model == 'RR':
			ML = MLmodels.RR()
		elif model == 'KRR':
			# kernel = kernels.Polynomial(3)
			# ML = MLmodels.KRR(kernel)
			ML = MLmodels.SK_KRR(kernel)

			print(ML)

		#get all data:
		parameters = MLparams.define_parameters()
		X = MLparams.get_all_columns()
		Y = MLparams.get_column('Eact')
		
		#remove outliers
		oi = MLparams.outlier_idx(X,2.5) + MLparams.outlier_idx(Y,2.5)
		X = np.delete(X, oi, axis=0)
		Y = np.delete(Y, oi, axis=0)
		X = MLutils.normalize(X)
		Y = MLutils.normalize(Y)
		X, Xmu = MLutils.center(X)
		Y, Ymu = MLutils.center(Y)
		print(f'Removed {len(oi)} outliers')

		# cfm = feature_maps.Polynomial(3)
		# X = MLutils.add_ones(X)
		# X = cfm(X)
		# X = MLutils.normalize(X)
		
		# print(f'There are {X.shape[1]} features per datapoint using {cfm}')

		#split into train and test sets
		train, test = MLutils.split_data(X, Y, .8)
		Xtrain, Ytrain = train
		Xtest, Ytest = test

		ML.cross_validate(X, Y)
		ML.optimize(X,Y)


		if animate:
			out_path = f'analyses/ML_plots/{model}/vary_lambda_{out_path_idx}'
			if not os.path.exists(out_path):
				os.mkdir(out_path)

			ls = np.linspace(0, 1, 50)
			files = []
			for i, l in enumerate(ls):
				ML.train(Xtrain, Ytrain, l)
				#evaluate the model
				R2test = ML.evaluate(Xtest, Ytest)
				Ypredtrain = ML.predict(Xtrain)
				Ypredtest = ML.predict(Xtest)

				print(f'Training with {Xtrain.shape[0]} datapoints: R2 = {ML.R2train}')
				print(f'Testing with {Xtest.shape[0]} datapoints: R2 = {R2test}')
				plt.title(rf'Prediction of $\Delta E^‡$ using {ML.description}')
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
				plt.savefig(f'{out_path}/{i}.png')
				files.append(f'{out_path}/{i}.png')
				plt.close()
			make_gif(f'{out_path}/animate.gif', files)

		else:
			#train the model
			ML.train(Xtrain,Ytrain, 1)

			#evaluate the model
			R2test = ML.evaluate(Xtest, Ytest)
			Ypredtrain = ML.predict(Xtrain)
			Ypredtest = ML.predict(Xtest)

			print(f'Training with {Xtrain.shape[0]} datapoints: R2 = {ML.R2train}')
			print(f'Testing with {Xtest.shape[0]} datapoints: R2 = {R2test}')
			plt.title(rf'Prediction of $\Delta E^‡$ using {ML.description}')
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
			plt.savefig(f'{out_path}/vary_lam.png')
			plt.close()

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
		plt.savefig(f'{out_path}/vary_dist.png')
		plt.close()


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
		plt.savefig(f'{out_path}/coefficients_multi.png')
		plt.close()


		#vary both f and l:
		fs = np.linspace(0.1, .90, 60)
		lams = np.linspace(0.0, 10, 60)
		R2s = np.zeros((lams.size, fs.size)) -1 
		R2best = (0,0)
		for i, f in enumerate(fs):
			for j, lam in enumerate(lams):
				utility.loading_bar(i*fs.size + j, fs.size * lams.size-1, 50)
				(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, f)
				# noise = 10**n
				# kernel_sigma2 = float(10**s)
				# kernel = kernels.Gaussian(kernel_sigma2)
				# model = models.GPR(kernel)
				# model.train(Xtrain, Ytrain, noise)
				# Ypred = model.predict(Xtest)
				ML.train(Xtrain, Ytrain, lam)
				Ypred = ML.predict(Xtest)

				R2 = MLutils.correlation(Ypred, Ytest)**2
				if R2 >= R2s.max():
					print('new best (R2)', R2, f, lam)
					R2best = (f, lam)
					R2bestidx = (f, lam)
				if not np.isnan(R2):
					R2s[j,i] = R2
				else:
					R2s[j,i] = 0

		extent = fs.min(), fs.max(), lams.min(), lams.max()
		plt.imshow(R2s, extent=extent, origin='lower', interpolation='antialiased', aspect='auto')
		plt.xlabel(r'$f$')
		plt.ylabel(r'$\lambda$')
		plt.annotate(rf'Best fit, $r^2={R2s.max()}$' + '\n' + rf'f$={R2best[0]:.5E}$, $\lambda={R2best[1]:.5E}$', 
				xy=R2bestidx, xytext=(0.3,.8), textcoords='figure fraction', arrowprops=dict(arrowstyle="->"), color='white',
				annotation_clip=False)
		plt.scatter(*R2bestidx, c='k')
		plt.savefig(os.path.join(out_path, 'heatmap_r2.png'))
		plt.close()

		#vary kernel params
		if kernel in ['rbf', 'laplacian', 'poly', 'chi2', 'additive_chi2']:
			...

		# ML.multi_train(X, Y, f=0.5, N=150, lam=1)
		# sortidx = np.argsort(np.abs(ML.a.flatten()))
		# labels = cfm.get_feature_labels([''] + parameters)
		# print('  Weight | Feature')
		# for i in sortidx[::-1]:
		# 	print(f'{float(ML.a[i]): .5f} | {labels[i]}')


	elif model == 'GPR':
		parameters = MLparams.define_parameters()
		X = MLparams.get_all_columns()
		X, Xmu = MLutils.center(X)
		Y = MLparams.get_column('Eact')
		Y, Ymu = MLutils.center(Y)
		(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .7)

		vary_noise = False
		vary_sigma = False
		vary_both  = True

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
				# kernel = kernels.Gaussian(kernel_sigma2)
				kernel = kernels.Polynomial(3)
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
			ns = np.linspace(-6, 10, 100)
			ss = np.linspace(-10,10, 100)
			R2s = np.zeros((ss.size, ns.size)) -1 
			R2best = (0,0)
			for i, n in enumerate(ns):
				for j, s in enumerate(ss):
					utility.loading_bar(i*ns.size + j, ns.size * ss.size-1, 50)
					noise = 10**n
					kernel_sigma2 = float(10**s)
					kernel = kernels.Gaussian(kernel_sigma2)
					model = models.GPR(kernel)
					model.train(Xtrain, Ytrain, noise)
					Ypred = model.predict(Xtest)

					R2 = MLutils.correlation(Ypred, Ytest)**2
					if R2 >= R2s.max():
						print('new best (R2)', R2, noise, kernel_sigma2)
						R2best = (noise, kernel_sigma2)
						R2bestidx = (n, s)
					if not np.isnan(R2):
						R2s[j,i] = R2
					else:
						R2s[j,i] = 0

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

			Ypred = model.predict(Xtrain)
			plt.scatter(Ytrain, Ypred, label='Train', alpha=.2, c='k')
			Ypred = model.predict(Xtest)
			plt.scatter(Ytest, Ypred, label='Test', c='r')

			plt.suptitle(model.__repr__() + f'\nNoise = {noise}')
			plt.xlabel(r'$E^‡$ (kcal/mol)')
			plt.ylabel(r'$E^‡_{pred}$ (kcal/mol)')
			plt.gca().set_xlim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
			plt.gca().set_ylim(Ytest.min()-.1*Ytest.min(), Ytest.max()+.1*Ytest.max())
			plt.plot((Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), ( Ytest.max()+.1*Ytest.max(), Ytest.min()-.1*Ytest.min()), 'k--')
			plt.legend()

			plt.savefig(os.path.join(out_path, f'best_fit_r2.png'))
			plt.close()



def optimize(model='KRR', kernel='rbf', out_path_idx=0):
	def plot_prediction(ML, Xtrain, Ytrain, Xtest, Ytest):
		Ypred_train = ML.predict(Xtrain)
		Ypred_test = ML.predict(Xtest)

		Ypred_train = MLutils.unprepare_data(Ypred_train, Y_norm_data) * utility.hartree2kcalmol(1)
		Ypred_test =  MLutils.unprepare_data(Ypred_test, Y_norm_data)  * utility.hartree2kcalmol(1)
		Ytrain = 	  MLutils.unprepare_data(Ytrain, Y_norm_data)      * utility.hartree2kcalmol(1)
		Ytest = 	  MLutils.unprepare_data(Ytest, Y_norm_data)       * utility.hartree2kcalmol(1)

		plt.scatter(Ytrain, Ypred_train)
		plt.scatter(Ytest, Ypred_test)
		plt.xlabel(r'$E^‡$ (kcal/mol)')
		plt.ylabel(r'$E^‡_{pred}$ (kcal/mol)')
		plt.title(f'Prediction performance for {ML}')

		mi = min(Ytest.min()-.1*Ytest.min(), Ytrain.min()-.1*Ytrain.min())
		ma = max(Ytest.max()+.1*Ytest.max(), Ytrain.max()+.1*Ytrain.max())

		plt.gca().set_xlim(mi, ma)
		plt.gca().set_ylim(mi, ma)
		plt.plot((ma, mi), (ma, mi), 'k--')


	#get all data:
	parameters = MLparams.define_parameters()
	X = MLparams.get_all_columns()
	# X = MLutils.add_ones(X)
	Y = MLparams.get_column('Eact')
	
	#remove outliers
	oi = MLparams.outlier_idx(X,2.5) + MLparams.outlier_idx(Y,2.5)
	X = np.delete(X, oi, axis=0)
	Y = np.delete(Y, oi, axis=0)
	X, X_norm_data = MLutils.prepare_data(X, True, True)
	Y, Y_norm_data = MLutils.prepare_data(Y, True, True)
	print(f'Removed {len(oi)} outliers')

	

	if model == 'KRR':
		ML = MLmodels.SK_KRR(kernel)
		out_path = f'analyses/ML_plots/optimize/{model}_{kernel}_{out_path_idx}'

	if model == 'RR':
		fm = feature_maps.Polynomial(3)
		X = MLutils.add_ones(X)
		X = fm(X)
		ML = MLmodels.SK_RR()
		out_path = f'analyses/ML_plots/optimize/{model}_{out_path_idx}'

	if not os.path.exists(out_path):
		os.mkdir(out_path)

	(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .7)
	ML.optimize(Xtrain, Ytrain)
	plot_prediction(ML, Xtrain, Ytrain, Xtest, Ytest)
	plt.savefig(join(out_path, 'prediction'))
	plt.close()


		
