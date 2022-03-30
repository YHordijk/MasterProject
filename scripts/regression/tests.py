import models, feature_maps, kernels, MLutils
import matplotlib.pyplot as plt
import numpy as np
import moviepy.editor as mvp
import os


def make_gif(out, frames, fps=14):
	# frames = [frames_folder + f for f in os.listdir(frames_folder)]
	clip = mvp.ImageSequenceClip(frames, fps=fps)
	clip.write_gif(out, fps=fps)




target = 'all'
#tests
if target in ['all', 'GPR']:
	simple = True
	vary_noise = True
	vary_sigma = True
	vary_both = True
	multi_variate = False
	
	m = 10
	xmus = [-4,-1,4]
	xsigs = [.1, 1, .3]
	out_path_idx = 3
	# xmus = [0]
	# xsigs = [3]

	# Ytruth = lambda X: 1/300*X**3 + X**2/100 + .5
	Ytruth = lambda X: np.cos(X) + X**2/10

	#generate data
	X = np.random.multivariate_normal(xmus, np.eye(len(xsigs))*xsigs, m//len(xmus)).flatten().reshape(-1,1)
	Y = Ytruth(X) + np.random.normal(0, 6, X.size).reshape(-1,1)/75

	(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .7, randomize=True)


	if simple:
		out_path = f'test_plots/GPR/simple{out_path_idx}'
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		noise = 0
		kernel_sigma2 = 1

		kernel = kernels.Gaussian(kernel_sigma2)
		model = models.GPR(kernel)
		
		model.train(Xtrain, Ytrain, noise)
		R2 = model.evaluate(Xtest, Ytest)**2
		error = MLutils.RSS(Ytest, model.predict(Xtest))
		
		plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
		mean, var = model.predict(plotX)
		std = np.sqrt(var)
		
		plt.figure(figsize=(9,9))
		plt.fill_between(plotX.flatten(), (mean+2.58*std).flatten(), (mean-2.58*std).flatten(), color=(229/255,229/255,1), label='99% CI')
		plt.fill_between(plotX.flatten(), (mean+1.96*std).flatten(), (mean-1.96*std).flatten(), color=(205/255,205/255,1), label='95% CI')
		plt.fill_between(plotX.flatten(), (mean+1.15*std).flatten(), (mean-1.15*std).flatten(), color=(184/255,184/255,1), label='75% CI')
		plt.plot(plotX, mean, label=rf'$\mu_{{test}}, r^2={R2}$')
		Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
		plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
		plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
		plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
		xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
		plt.xlabel(xlabel)
		plt.ylabel(rf'$y$')
		plt.suptitle(model.__repr__() + f'\nNoise = {noise:.10f}')
		plt.legend()
		plt.gca().set_xlim(X.min()-1, X.max()+1)
		plt.gca().set_ylim(Y.min()-1, Y.max()+1)

		plt.savefig(os.path.join(out_path, f'default.png'))
		plt.close()

		noise = 0
		kernel_sigma2 = 1

		kernel = kernels.Gaussian(kernel_sigma2)
		model = models.GPR(kernel)
		
		model.train(Xtrain, Ytrain, noise)
		R2 = model.evaluate(Xtest, Ytest)**2
		error = MLutils.RSS(Ytest, model.predict(Xtest))
		
		plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
		mean, var = model.predict(plotX)
		std = np.sqrt(var)
		
		plt.figure(figsize=(9,9))
		Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
		plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
		plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
		plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
		xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
		plt.xlabel(xlabel)
		plt.ylabel(rf'$y$')
		plt.suptitle(model.__repr__() + f'\nNoise = {noise:.10f}')
		plt.legend()
		plt.gca().set_xlim(X.min()-1, X.max()+1)
		plt.gca().set_ylim(Y.min()-1, Y.max()+1)

		plt.savefig(os.path.join(out_path, f'default_nofit.png'))
		plt.close()


	if vary_noise:
		out_path = f'test_plots/GPR/vary_noise{out_path_idx}'
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		R2s = []
		errors = []
		ns = np.linspace(-10, 2, 100)
		frames = []
		for i, n in enumerate(ns):
			noise = 10**n
			kernel_sigma2 = 1

			kernel = kernels.Gaussian(kernel_sigma2)
			model = models.GPR(kernel)
			
			model.train(Xtrain, Ytrain, noise)
			R2 = model.evaluate(Xtest, Ytest)**2
			R2s.append(R2)
			error = MLutils.RSS(Ytest, model.predict(Xtest))
			errors.append(error)
			
			plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
			mean, var = model.predict(plotX)
			std = np.sqrt(var)
			
			fig, ax1 = plt.subplots(figsize=(16,9))
			plt.subplot(1,2,1)
			plt.fill_between(plotX.flatten(), (mean+2.58*std).flatten(), (mean-2.58*std).flatten(), color=(229/255,229/255,1), label='99% CI')
			plt.fill_between(plotX.flatten(), (mean+1.96*std).flatten(), (mean-1.96*std).flatten(), color=(205/255,205/255,1), label='95% CI')
			plt.fill_between(plotX.flatten(), (mean+1.15*std).flatten(), (mean-1.15*std).flatten(), color=(184/255,184/255,1), label='75% CI')
			plt.plot(plotX, mean, label=rf'$\mu_{{test}}, r^2={R2}$')
			Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
			plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
			plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
			plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
			xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
			plt.xlabel(xlabel)
			plt.ylabel(rf'$y$')
			plt.suptitle(model.__repr__() + f'\nNoise = {noise:.10f}')
			plt.legend()
			plt.gca().set_xlim(X.min()-1, X.max()+1)
			plt.gca().set_ylim(Y.min()-1, Y.max()+1)

			ax1 = plt.subplot(1,2,2)
			ax1.plot(ns[:i+1], R2s, c='b')
			ax1.tick_params('y', colors='b')
			ax2 = ax1.twinx()
			ax2.plot(ns[:i+1], errors, c='r')
			plt.xlabel(r'log$_{10}($noise$)$')
			ax1.set_ylabel(r'$r^2$ correlation')
			ax1.set_xlim(ns.min(), ns.max())
			ax1.set_ylim(0,1)
			ax2.set_xlim(ns.min(), ns.max())
			ax2.set_ylabel(r'Residual sum of squares')
			ax2.tick_params('y', colors='r')

			plt.savefig(os.path.join(out_path, f'{i}.png'))
			plt.close()
			frames.append(os.path.join(out_path, f'{i}.png'))

		make_gif(os.path.join(out_path, 'anim.gif'), frames, fps=14)

	if vary_sigma:
		out_path = f'test_plots/GPR/vary_sigma{out_path_idx}'
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		R2s = []
		errors = []
		ss = np.linspace(-2, 4, 100)
		frames = []
		for i, s in enumerate(ss):
			noise = 0.000001
			kernel_sigma2 = float(10**s)

			kernel = kernels.Gaussian(kernel_sigma2)
			model = models.GPR(kernel)
			
			model.train(Xtrain, Ytrain, noise)
			R2 = model.evaluate(Xtest, Ytest)**2
			R2s.append(R2)
			error = MLutils.RSS(Ytest, model.predict(Xtest))
			errors.append(error)
			
			plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
			mean, var = model.predict(plotX)
			std = np.sqrt(var)
			
			fig, ax1 = plt.subplots(figsize=(16,9))
			plt.subplot(1,2,1)
			plt.fill_between(plotX.flatten(), (mean+2.58*std).flatten(), (mean-2.58*std).flatten(), color=(229/255,229/255,1), label='99% CI')
			plt.fill_between(plotX.flatten(), (mean+1.96*std).flatten(), (mean-1.96*std).flatten(), color=(205/255,205/255,1), label='95% CI')
			plt.fill_between(plotX.flatten(), (mean+1.15*std).flatten(), (mean-1.15*std).flatten(), color=(184/255,184/255,1), label='75% CI')
			plt.plot(plotX, mean, label=rf'$\mu_{{test}}, r^2={R2}$')
			Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
			plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
			plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
			plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
			xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
			plt.xlabel(xlabel)
			plt.ylabel(rf'$y$')
			plt.suptitle(model.__repr__() + f'\nNoise = {noise:.10f}')
			plt.legend()
			plt.gca().set_xlim(X.min()-1, X.max()+1)
			plt.gca().set_ylim(Y.min()-1, Y.max()+1)

			ax1 = plt.subplot(1,2,2)
			ax1.plot(ss[:i+1], R2s, c='b')
			ax2 = ax1.twinx()
			ax2.plot(ss[:i+1], errors, c='r')
			plt.xlabel(r'log$_{10}($noise$)$')
			ax1.set_ylabel(r'$r^2$ correlation')
			ax1.set_xlim(ss.min(), ss.max())
			ax1.tick_params('y', colors='b')
			ax1.set_ylim(0,1)
			ax2.set_xlim(ss.min(), ss.max())
			ax2.set_ylabel(r'Residual sum of squares')
			ax2.tick_params('y', colors='r')

			plt.savefig(os.path.join(out_path, f'{i}.png'))
			plt.close()
			frames.append(os.path.join(out_path, f'{i}.png'))

		make_gif(os.path.join(out_path, 'anim.gif'), frames, fps=14)

	if vary_both:
		out_path = f'test_plots/GPR/vary_both{out_path_idx}'
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		fig = plt.figure(figsize=(9,9))
		ns = np.linspace(-14, 3, 10)
		ss = np.linspace(-4, 5, 10)
		R2s = np.zeros((ss.size, ns.size))
		errors = np.zeros((ss.size, ns.size)) + 10000000000
		for i, n in enumerate(ns):
			for j, s in enumerate(ss):
				noise = 10**n
				kernel_sigma2 = float(10**s)
				kernel = kernels.Gaussian(kernel_sigma2)
				model = models.GPR(kernel)
				
				model.train(Xtrain, Ytrain, noise)
				R2 = model.evaluate(Xtest, Ytest)**2
				error = MLutils.RSS(Ytest, model.predict(Xtest))
				if R2 >= R2s.max():
					print('new best (R2)', R2, noise, kernel_sigma2)
					R2best = (noise, kernel_sigma2)
					R2bestidx = (n, s)
				R2s[j,i] = R2
				if error <= errors.min():
					print('new best (error)', error, noise, kernel_sigma2)
					errorbest = (noise, kernel_sigma2)
					errorbestidx = (n, s)
				errors[j,i] = error

		extent = ns.min(), ns.max(), ss.min(), ss.max()
		plt.imshow(R2s, extent=extent, origin='lower', interpolation='antialiased', aspect='auto')
		plt.xlabel(r'log$_{10}(\sigma^2_n)$')
		plt.ylabel(r'log$_{10}(\sigma^2_K)$')
		plt.annotate(rf'Best fit, $r^2={R2s.max()}$' + '\n' + rf'noise$={R2best[0]}$, $\sigma^2={R2best[1]}$', 
				xy=R2bestidx, xytext=(0.5,0.7), textcoords='figure fraction', arrowprops=dict(arrowstyle="->"), color='white')
		plt.scatter(*R2bestidx, c='k')
		plt.savefig(os.path.join(out_path, 'heatmap_r2.png'))
		plt.close()

		extent = ns.min(), ns.max(), ss.min(), ss.max()
		plt.imshow(np.log10(errors), extent=extent, origin='lower', interpolation='antialiased', aspect='auto')
		plt.xlabel(r'log$_{10}(\sigma^2_n)$')
		plt.ylabel(r'log$_{10}(\sigma^2_K)$')
		plt.annotate(rf'Best fit, $error={errors.min()}$' + '\n' + rf'noise$={errorbest[0]}$, $\sigma^2={errorbest[1]}$', 
				xy=errorbestidx, xytext=(0.5,0.7), textcoords='figure fraction', arrowprops=dict(arrowstyle="->"), color='white')
		plt.scatter(*errorbestidx, c='k')
		plt.savefig(os.path.join(out_path, 'heatmap_error.png'))
		plt.close()


		kernel = kernels.Gaussian(R2best[1])
		model = models.GPR(kernel)
		
		model.train(Xtrain, Ytrain, R2best[0])
		R2 = model.evaluate(Xtest, Ytest)**2
		
		plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
		mean, var = model.predict(plotX)
		std = np.sqrt(var)
		
		plt.figure(figsize=(9,9))
		plt.fill_between(plotX.flatten(), (mean+2.58*std).flatten(), (mean-2.58*std).flatten(), color=(229/255,229/255,1), label='99% CI')
		plt.fill_between(plotX.flatten(), (mean+1.96*std).flatten(), (mean-1.96*std).flatten(), color=(205/255,205/255,1), label='95% CI')
		plt.fill_between(plotX.flatten(), (mean+1.15*std).flatten(), (mean-1.15*std).flatten(), color=(184/255,184/255,1), label='75% CI')
		plt.plot(plotX, mean, label=rf'$\mu_{{test}}, r^2={R2}$')
		Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
		plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
		plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
		plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
		xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
		plt.xlabel(xlabel)
		plt.ylabel(rf'$y$')
		plt.suptitle(model.__repr__() + f'\nNoise = {R2best[0]}')
		plt.legend()
		plt.gca().set_xlim(X.min()-1, X.max()+1)
		plt.gca().set_ylim(Y.min()-1, Y.max()+1)
		plt.savefig(os.path.join(out_path, 'best_fit_r2.png'))

		kernel = kernels.Gaussian(errorbest[1])
		model = models.GPR(kernel)
		
		model.train(Xtrain, Ytrain, errorbest[0])
		R2 = model.evaluate(Xtest, Ytest)**2
		
		plotX = np.linspace(X.min()-1, X.max()+1, 1000).reshape(-1,1)
		mean, var = model.predict(plotX)
		std = np.sqrt(var)
		
		plt.figure(figsize=(9,9))
		plt.fill_between(plotX.flatten(), (mean+2.58*std).flatten(), (mean-2.58*std).flatten(), color=(229/255,229/255,1), label='99% CI')
		plt.fill_between(plotX.flatten(), (mean+1.96*std).flatten(), (mean-1.96*std).flatten(), color=(205/255,205/255,1), label='95% CI')
		plt.fill_between(plotX.flatten(), (mean+1.15*std).flatten(), (mean-1.15*std).flatten(), color=(184/255,184/255,1), label='75% CI')
		plt.plot(plotX, mean, label=rf'$\mu_{{test}}, r^2={R2}$')
		Xtruth = np.linspace(X.min()-1, X.max()+1, 100)
		plt.plot(Xtruth, Ytruth(Xtruth), c='k', label='Ground Truth')
		plt.scatter(Xtrain, Ytrain, c='g', label=f'Train (N={Ytrain.size})')
		plt.scatter(Xtest, Ytest, c='r', label=f'Test (N={Ytest.size})')
		xlabel = rf'$x \sim $' + ' + '.join([rf'$\mathcal{{N}}({mu}, {sig})$' for mu, sig in zip(xmus, xsigs)])
		plt.xlabel(xlabel)
		plt.ylabel(rf'$y$')
		plt.suptitle(model.__repr__() + f'\nNoise = {errorbest[0]}')
		plt.legend()
		plt.gca().set_xlim(X.min()-1, X.max()+1)
		plt.gca().set_ylim(Y.min()-1, Y.max()+1)
		plt.savefig(os.path.join(out_path, 'best_fit_error.png'))


	if multi_variate:
		out_path = f'test_plots/GPR/multi_variate{out_path_idx}'
		vary_noise = True
		vary_sigma = True
		vary_both = True
		if not os.path.exists(out_path):
			os.mkdir(out_path)

		mu = [0, 1, 3, -1]
		sigma = np.random.rand(len(mu), len(mu)) *4 + 1
		X = np.random.multivariate_normal(mu, sigma, 400)
		Y = (X[:,0]**2/5 + 1/10*X[:,1] + np.exp(-X[:,2])/5 + np.random.normal(0, 6, X.shape[0])/40).reshape(-1,1)

		(Xtrain, Ytrain), (Xtest, Ytest) = MLutils.split_data(X, Y, .8)

		if vary_noise:
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
				Ypred, var = model.predict(Xtrain)
				plt.scatter(Ytrain, Ypred, label='Train', alpha=.2, c='k')

				Ypred, var = model.predict(Xtest)
				plt.scatter(Ytest, Ypred, label='Test', c='r')
				R2s.append(MLutils.correlation(Ypred, Ytest)**2)
				
				plt.suptitle(model.__repr__() + f'\nNoise = {noise}')
				plt.xlabel(r'$Y$')
				plt.ylabel(r'$Y_{pred}$')
				plt.gca().set_xlim(Ytest.min()-1, Ytest.max()+1)
				plt.gca().set_ylim(Ytest.min()-1, Ytest.max()+1)
				plt.plot((Ytest.max()+1, Ytest.min()-1), (Ytest.max()+1, Ytest.min()-1), 'k--')
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
			ss = np.linspace(-2, 4, 100)
			frames = []
			R2s = []
			for i, s in enumerate(ss):
				noise = 2.78e-3
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
				plt.gca().set_xlim(Ytest.min()-1, Ytest.max()+1)
				plt.gca().set_ylim(Ytest.min()-1, Ytest.max()+1)
				plt.plot((Ytest.max()+1, Ytest.min()-1), (Ytest.max()+1, Ytest.min()-1), 'k--')
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



			# Ypred, var = model.predict(Xtrain)
			# plt.scatter(Ytrain, Ypred)
			# plt.show()