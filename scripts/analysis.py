from analyses import profiles, functional_comparison, draw, trends, ML
import matplotlib.pyplot as plt


mode = 'ML'
if mode == 'ML':
	# ML.main('KRR', out_path_idx=6, kernel='rbf')
	ML.optimize(model='KRR', kernel='rbf', out_path_idx=1)

elif mode == 'view':
	# for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']:
	# 	for R1 in ['Et']:
	# 		for R2 in ['p-FPh', 'o-FPh', 'm-FPh', 'Ph', 'tBu']:
	# 			draw.draw('achiral_catalyst', {'R1':R1, 'R2':'p-FPh', 'Rcat':Rcat})

	draw.draw('achiral_catalyst', {'R1':'NMe2', 'R2':'Ph', 'Rcat':'I2'})

# elif mode == 'vibe_check':
