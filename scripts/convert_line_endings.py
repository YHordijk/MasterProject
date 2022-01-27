import os


for file in os.listdir():
	try:
		with open(file, 'U') as f:
			text = f.read()
		with open(file, 'w') as f:
			f.write(text)
	except: pass

