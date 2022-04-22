from analyses import profiles, functional_comparison, draw, trends, ML
import matplotlib.pyplot as plt


mode = 'ML'
if mode == 'ML':
	# ML.main('KRR', out_path_idx=6, kernel='rbf')
	ML.optimize(model='KRR', kernel='laplacian', out_path_idx=4)

elif mode == 'view':
	for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']:
	# for Rcat in ['I2']:
		for R1 in ['OMe']:
			for R2 in ['o-FPh']:
				draw.draw('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':Rcat})

	# draw.draw('achiral_catalyst', {'R1':'NMe2', 'R2':'m-FPh', 'Rcat':'SnCl4'})

# elif mode == 'vibe_check':
