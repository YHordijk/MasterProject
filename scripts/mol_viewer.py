import molviewer.world as world
import molviewer.screen as screen
import molviewer.molecule as chemmolecule
from OpenGL.GL import *
import glfw
from OpenGL.GLU import *
import numpy as np
from math import cos, sin
import molviewer.data as data

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
	def __init__(self, molecule, **kwargs):
		if type(molecule) is str:
			self.molecule = chemmolecule.load_molecule(molecule)
		else:
			self.molecule = molecule


		self.scale = kwargs.get('scale', 1/5)
		self.position = kwargs.get('position', [0,0,0])
		self.positions = self.molecule.positions * self.scale
		self.orientation = kwargs.get('orientation', [0,0,0,0])
		self.elements = self.molecule.atom_numbers
		self.radii = self.molecule.radii
		self.bonds = self.molecule.unique_bonds
		self.sub_objects = self.get_sub_objects()

		self.vertices = []
		for o in self.sub_objects:
			self.vertices = self.vertices + o.vertices

		self.colors = []
		for o in self.sub_objects:
			self.colors = self.colors + o.colors
		self.colors = (np.asarray(self.colors)/255).tolist()


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
		self.scr = screen.Screen(screen_width=1600, screen_height=900)

		wrld = world.World()
		self.wrld = wrld
		wrld.add_object(self)
		self.last_mouse_pos = [0,0]
		self.mouse_rx = 0
		self.mouse_ry = 0 
		self.mouse_d_last = [0,0]
		self.scr.mainloop(wrld)


    def rotate_target(self, delta):
        right = glm.normalize(glm.cross(self.target - self.eye, self.up))
        M = glm.mat4(1)
        M = glm.translate(M, self.eye)
        M = glm.rotate(M, delta.y, right)
        M = glm.rotate(M, delta.x, self.up)
        M = glm.translate(M, -self.eye)
        self.target = glm.vec3(M * glm.vec4(self.target, 1.0))
        

	def pre_draw(self):
		glLoadIdentity()
		glRotatef(*self.orientation)
		mouse_pos = glfw.get_cursor_pos(self.scr.window)
		if glfw.get_mouse_button(self.scr.window, glfw.MOUSE_BUTTON_LEFT):
			dx = np.sqrt((mouse_pos[0] - self.last_mouse_pos[0])**2)
			dy = np.sqrt((mouse_pos[1] - self.last_mouse_pos[1])**2)
			self.mouse_rx += np.sqrt((mouse_pos[0] - self.last_mouse_pos[0])**2) - self.mouse_d_last[0]
			self.mouse_ry += np.sqrt((mouse_pos[1] - self.last_mouse_pos[1])**2) - self.mouse_d_last[1]
			self.mouse_d_last = [dx, dy]


		
		glRotatef(self.mouse_ry/5, 0, 0, 1)
		glRotatef(self.mouse_rx/5, 0, 1, 0)
		glTranslatef(*self.position)

		


	def start(self):
		glEnableClientState(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, self.vertices)
		glEnableClientState(GL_COLOR_ARRAY)
		glColorPointer(3, GL_FLOAT, 0, self.colors)


	def draw(self):
		self.pre_draw()
		glDrawArrays(GL_LINES, 0, len(self.vertices)//3)



if __name__ == '__main__':
	moldrawer = MoleculeDrawer('water')
	moldrawer.mainloop()
