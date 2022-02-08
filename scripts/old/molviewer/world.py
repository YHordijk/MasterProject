import numpy as np
import glfw
from OpenGL.GL import *
from math import sin, cos

class World:
	'''
	Container class representing a world
	The objects must have a number of methods
	- draw
	'''

	def __init__(self, **kwargs):
		self.objs = kwargs.get('objects', [])

	def __iter__(self):
		self.iter_counter = 0
		return self

	def __next__(self):
		try:
			self.iter_counter += 1
			return self.objs[self.iter_counter-1]
		except IndexError: raise StopIteration

	def add_object(self, obj):
		self.objs.append(obj)

	def start(self):
		[o.start() for o in self.objs if hasattr(o, 'start')]



class cube:
	def __init__(self, **kwargs):
		self.position = kwargs.get('position', [0,0,0])
		self.orientation = kwargs.get('orientation', [0,0,0,0])
		self.edge_length = kwargs.get('edge_length', 1)


	def vertices(self):
		l = self.edge_length/2
		a, e = [-l, -l, -l], [-l, -l,  l]
		b, f = [-l,  l, -l], [-l,  l,  l]
		c, g = [ l, -l, -l], [ l, -l,  l]
		d, h = [ l,  l, -l], [ l,  l,  l]

		v = a+b + a+c + c+d + b+d + a+e + b+f + d+h + c+g + e+f + e+g + g+h + f+h

		return v

	def pre_draw(self):
		glLoadIdentity()
		glRotatef(*self.orientation)
		glTranslatef(*self.position)


	def draw(self):
		vertices = self.vertices()
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, vertices)
		self.pre_draw()
		glDrawArrays(GL_LINES, 0, len(vertices)//3)



class circle2D:
	def __init__(self, **kwargs):
		self.position = kwargs.get('position', [0,0,0])
		self.radius = kwargs.get('radius', 1)
		self.resolution = kwargs.get('resolution', 10)
		self.vertices = self.get_vertices()


	def get_vertices(self):
		v = []
		for i in range(self.resolution):
			v.append(cos(6.28318530718*i/self.resolution)/2)
			v.append(sin(6.28318530718*i/self.resolution)/2)
			v.append(0)
			if i > 0:
				v.append(cos(6.28318530718*i/self.resolution)/2)
				v.append(sin(6.28318530718*i/self.resolution)/2)
				v.append(0)
			if i == self.resolution-1:
				v.append(.5)
				v.append(0)
				v.append(0)
		return v

	def pre_draw(self):
		glLoadIdentity()
		glScale(self.radius, self.radius, 1)
		glTranslatef(*self.position)

	def draw(self):
		vertices = self.vertices
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, vertices)
		self.pre_draw()
		glDrawArrays(GL_LINES, 0, len(vertices)//3)


class rectangle2D:
	def __init__(self, **kwargs):
		self.position = kwargs.get('position', [0,0,0])
		self.orientation = kwargs.get('orientation', [0,0,0,0])
		self.edge_lengths = kwargs.get('edge_lengths', [1,1])


	def get_vertices(self):
		l1 = self.edge_lengths[0]/2
		l2 = self.edge_lengths[1]/2
		a = [-l1, -l2, 0]
		b = [ l1, -l2, 0]
		c = [ l1,  l2, 0]
		d = [-l1,  l2, 0]

		v = a+b + b+c + c+d + d+a

		return v

	def pre_draw(self):
		glLoadIdentity()
		ct = glfw.get_time()
		# glScale(sin(ct*2), sin(ct*2), sin(ct*2))
		# glRotatef(*self.orientation)
		glRotatef(self.orientation[0], 1,0,0)
		glRotatef(self.orientation[1], 0,1,0)
		glRotatef(self.orientation[2], 0,0,1)
		glRotatef(ct*60, 1,1,1)
		# glRotatef(0, ct,ct,0)
		glTranslatef(*self.position)


	def draw(self):
		vertices = self.get_vertices()
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, vertices)
		self.pre_draw()
		glDrawArrays(GL_LINES, 0, len(vertices)//3)



if __name__ == '__main__':
	world = World(objects=[1,2,3,4,5])

	for obj in world:
		print(obj)