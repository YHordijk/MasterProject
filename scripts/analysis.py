from analyses import profiles, functional_comparison, draw, trends, ML
import matplotlib.pyplot as plt


mode = 'ML'
if mode == 'ML':
	ML.main()

elif mode == 'view':
	# for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']:
	for Rcat in ['ZnCl2']:
		for R1 in ['Me', 'NH2', 'OMe']:
			draw.draw('achiral_catalyst', {'R1':R1, 'R2':'p-FPh', 'Rcat':Rcat})
		# draw.draw('achiral_catalyst', {'R1':R1, 'R2':'Ph', 'Rcat':Rcat})

	# draw.draw('achiral_catalyst', {'R1':'Br', 'R2':'Ph', 'Rcat':'AlF3'})
