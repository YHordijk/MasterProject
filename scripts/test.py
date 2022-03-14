import itertools



X = [
	['x1a', 'x1b', 'x1c'],
	['x2a', 'x2b', 'x2c'],
	['x3a', 'x3b', 'x3c'],
	['x4a', 'x4b', 'x4c'],
	]


for x in X:
	print(x)
	features = itertools.combinations_with_replacement(x, 3)
	for feature in features:
		print('*'.join(feature))




featurelist = []
for x in X:
	features = itertools.combinations_with_replacement(x, 2)
	feature_products = []
	for feature in features:
		p = 1
		for f in feature:
			p *= f

		feature_products.append(p)
	featurelist.append(feature_products)
print(featurelist)

