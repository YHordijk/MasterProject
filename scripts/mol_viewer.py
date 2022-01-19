import molviewer.world as world
import molviewer.screen as screen
import molviewer.molecule as chemmolecule
from OpenGL.GL import *
import glfw
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
from math import cos, sin
import molviewer.data as data
import molviewer.text as text

'''
n  n-1 o

0  0   0
0  1   0
1  0   1
1  1   0
'''


class atom(world.circle2D):
	def __init__(self, **kwargs):
		self.orientation = kwargs.get('orientation', [0,0,0])
		super().__init__(**kwargs)
		self.position = kwargs.get('position', [0,0,0])
		self.element = kwargs.get('element', 1)
		self.molecule = kwargs.get('molecule')
		self.colors = (2*self.resolution) * data.ATOM_COLOURS[self.element]

	def get_vertices(self):
		v = []
		for i in range(self.resolution):
			v.append(self.radius*cos(6.28318530718*i/self.resolution)/2)
			v.append(self.radius*sin(6.28318530718*i/self.resolution)/2)
			v.append(0)
			if i > 0:
				v.append(self.radius*cos(6.28318530718*i/self.resolution)/2)
				v.append(self.radius*sin(6.28318530718*i/self.resolution)/2)
				v.append(0)
			if i == self.resolution-1:
				v.append(self.radius*.5)
				v.append(0)
				v.append(0)

		vs = np.array(v)
		vs = vs.reshape(vs.size//3, 3)
		
		ca, sa = cos(self.orientation[0]), sin(self.orientation[0])
		cb, sb = cos(self.orientation[1]), sin(self.orientation[1])
		cc, sc = cos(self.orientation[2]), sin(self.orientation[2])

		Ra = np.array([[ca, -sa, 0], [sa, ca,0],[0,0,1]])
		Rb = np.array([[cb,0,sb],[0,1,0],[-sb,0,cb]])
		Rc = np.array([[1,0,0],[0,cc,-sc],[0,sc,cc]])

		vs = (Ra@Rb@Rc@vs.T).T

		vs = vs + self.position

		return vs.flatten().tolist()
		


class bond(world.rectangle2D):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.positions = kwargs.get('positions')
		self.elements = kwargs.get('elements')
		self.molecule = kwargs.get('molecule')
		self.vertices = self.get_vertices()
		self.colors = self.get_colors()

	def get_vertices(self):
		# l1 = self.edge_lengths[0]/2
		# l2 = self.edge_lengths[1]/2
		# a = [-l1, -l2, 0]
		# b = [ l1, -l2, 0]
		# c = [ l1,  l2, 0]
		# d = [-l1,  l2, 0]

		# vs = np.array(a+b + b+c + c+d + d+a)
		# vs = vs.reshape(vs.size//3, 3)
		
		# ca, sa = cos(self.orientation[0]), sin(self.orientation[0])
		# cb, sb = cos(self.orientation[1]), sin(self.orientation[1])
		# cc, sc = cos(self.orientation[2]), sin(self.orientation[2])

		# Ra = np.array([[ca, -sa, 0], [sa, ca,0],[0,0,1]])
		# Rb = np.array([[cb,0,sb],[0,1,0],[-sb,0,cb]])
		# Rc = np.array([[1,0,0],[0,cc,-sc],[0,sc,cc]])

		# vs = (Ra@Rb@Rc@vs.T).T

		# vs = vs + self.position

		# return vs.flatten().tolist()

		# return self.positions[0].tolist() + self.positions[1].tolist()


		a, b = self.positions[0], self.positions[1]
		resolution = 10

		v = a.tolist() + b.tolist()
		# v.append()

		# for i in range(resolution):
		# 	v.append()

		return v


	def get_colors(self):
		return data.ATOM_COLOURS[self.elements[0]] + data.ATOM_COLOURS[self.elements[1]]



class MoleculeDrawer:
	def __init__(self, molecules, **kwargs):
		if type(molecules) is not list:
			if type(molecules) is str:
				self.molecules = [chemmolecule.load_molecule(molecules)]
			elif type(molecules) is dict:
				self.molecules = list(molecules.values())
			else:
				self.molecules = [molecules]

		else:
			self.molecules = molecules

		print(self.molecules)

		self.molidx = 0


		self.scale = kwargs.get('scale', 1/5)
		self.position = kwargs.get('position', [0,0,0])
		self.orientation = kwargs.get('orientation', [0,0,0,0])
		self.vertices = []
		self.colors = []
		for m in self.molecules:
			self.positions = m.positions * self.scale
			self.elements = m.atom_numbers
			self.radii = m.radii
			self.bonds = m.unique_bonds
			self.sub_objects = self.get_sub_objects()

			vertices = []
			for o in self.sub_objects:
				vertices = vertices + o.vertices

			vertices = vertices + [0,0,.02,0,0,0,0,.02,0,0,0,0,.02,0,0,0,0,0]
			self.vertices.append(vertices)

			colors = []
			for o in self.sub_objects:
				colors = colors + o.colors
			colors = (np.asarray(colors)/255).tolist()
			colors = colors + [0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0,0]
			self.colors.append(colors)


	def get_sub_objects(self):
		objs = [atom(position=p, radius=r*self.scale, element=e, resolution=50) for p, e, r in zip(self.positions, self.elements, self.radii)]
		for i in np.linspace(0,1,10):
			objs = objs + [atom(position=p, radius=r*self.scale, orientation=[0,np.pi*i,0], element=e, resolution=50) for p, e, r in zip(self.positions, self.elements, self.radii)]
			objs = objs + [atom(position=p, radius=r*self.scale, orientation=[0,0,np.pi*i], element=e, resolution=50) for p, e, r in zip(self.positions, self.elements, self.radii)]
		for b in self.bonds:
			p1, p2 = np.asarray(self.positions[b[0]]), np.asarray(self.positions[b[1]])
			e1, e2 = self.elements[b[0]], self.elements[b[1]]
			objs.append(bond(positions=[p1, p2], elements=[e1, e2]))
		return objs


	def mainloop(self):
		self.last_buttons = {'left':False, 'right':False, 'lmb':0, 'lmb_up':0, 'lmb_down':0}
		self.scr = screen.Screen(screen_width=1600, screen_height=900)
		wrld = world.World()
		self.wrld = wrld
		wrld.add_object(self)
		self.last_mouse_pos = [0,0]
		self.mouse_rx = 0
		self.mouse_ry = 0 
		self.mouse_d_last = [0,0]
		self.scr.mainloop(wrld)


	def pre_draw(self):
		
		glLoadIdentity()
		glRotatef(*self.orientation)
		mouse_pos = glfw.get_cursor_pos(self.scr.window)
		mouse_d = [mouse_pos[0] - self.last_mouse_pos[0], mouse_pos[1] - self.last_mouse_pos[1]]
		if glfw.get_mouse_button(self.scr.window, glfw.MOUSE_BUTTON_LEFT):
			self.mouse_rx += mouse_d[0]
			self.mouse_ry += mouse_d[1]

		glRotatef(self.mouse_ry/5, 0, 0, 1)
		glRotatef(self.mouse_rx/5, 1, 0, 0)
		glTranslatef(*self.position)
		self.last_mouse_pos = mouse_pos
		self.button_action()
		

	def start(self):
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, self.vertices[self.molidx])
		glEnableClientState(GL_COLOR_ARRAY)
		glColorPointer(3, GL_FLOAT, 0, self.colors[self.molidx])


	def draw(self):
		self.pre_draw()
		glDrawArrays(GL_LINES, 0, len(self.vertices[self.molidx])//3)
		# text.render_text(self.scr.window, self.molecules[self.molidx].name, 50., 50., 1, (255,255,255))

	def button_action(self):
		lmb = glfw.get_mouse_button(self.scr.window, glfw.MOUSE_BUTTON_LEFT)
		buttons = {'left': glfw.get_key(self.scr.window, glfw.KEY_LEFT),
				   'right': glfw.get_key(self.scr.window, glfw.KEY_RIGHT),
				   'lmb_down': int(lmb and not self.last_buttons['lmb']),
				   'lmb_up': int(self.last_buttons['lmb'] and not lmb), 
				   'lmb': lmb,}

		if buttons['left'] and not self.last_buttons['left']:
			self.molidx += 1
			self.molidx = self.molidx%len(self.molecules)
			self.start()
			print(f'Currently drawing: {self.molecules[self.molidx].name}')
		elif buttons['right'] and not self.last_buttons['right']:
			self.molidx -= 1
			self.molidx = self.molidx%len(self.molecules)
			self.start()
			print(f'Currently drawing: {self.molecules[self.molidx].name}')

		self.last_buttons = buttons

		



if __name__ == '__main__':
	moldrawer = MoleculeDrawer('water')
	moldrawer.mainloop()
