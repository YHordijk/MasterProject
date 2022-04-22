import os
import tensorflow as tf
from tensorflow.keras import layers, models
import ML_params 
import matplotlib.pyplot as plt



X, Y = ML_params.tf_dataset(.8)

model = models.Sequential([
	layers.Input(shape=(X.shape[1])),
	layers.Dense(128, activation='relu'),
	layers.Dense(128, activation='relu'),
	layers.Dense(2, activation='softmax')])

model.compile(optimizer='adam',
              loss='mean_squared_error',
              metrics=['mse'])

model.fit(X, Y, epochs=1500, validation_split=.2)
model.evaluate(X, Y)


Ypred = model.predict(X)


plt.scatter(Y[:,0], Ypred[:,1])
plt.show()