from analyses import profiles, functional_comparison, draw, trends
import matplotlib.pyplot as plt

functional_comparison.main()
plt.show()

plt.figure()
profiles.achiral_catalyst(R1='H', R2='Ph')
plt.figure()
profiles.achiral_catalyst(R1='H', R2='tBu')
plt.show()
# profiles.achiral_catalyst_Eact()
# plt.show()
# draw.draw('achiral_catalyst', {'R1':'F', 'R2':'tBu', 'Rcat':'SnCl4'})


trends.Eact_EDA()