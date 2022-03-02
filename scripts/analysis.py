from analyses import profiles, functional_comparison, draw, trends
import matplotlib.pyplot as plt

# functional_comparison.main()
# plt.show()

# plt.figure()
# profiles.achiral_catalyst(R1='F', R2='Ph')
# plt.figure()
# profiles.achiral_catalyst(R1='F', R2='tBu')
# plt.show()
# profiles.achiral_catalyst_Eact()
# plt.show()

for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']:
	print(Rcat)
	draw.draw('achiral_catalyst', {'R1':'OMe', 'R2':'tBu', 'Rcat':Rcat})
	draw.draw('achiral_catalyst', {'R1':'OMe', 'R2':'Ph', 'Rcat':Rcat})

# trends.Eact_EDA()