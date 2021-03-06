import pygame as pg
import pygame.locals
import numpy as np
import random 
import molviewer2.data as data
import periodictable as pt
import skimage.draw
import matplotlib.pyplot as plt
from time import perf_counter
from math import cos, sin, pi, acos, pi


l2_norm = lambda u, v: np.linalg.norm(u-v)

# atom_img = './data/images/atom_default.png'


class Screen:
	def __init__(self, *args, **kwargs):
		self.size = kwargs.get('size', (500,300))
		self.background_color = kwargs.get('background_color', (0,0,0))
		print(self.background_color)
		self.main_display = pg.display.set_mode(self.size, pg.locals.HWSURFACE | pg.locals.DOUBLEBUF | pg.locals.RESIZABLE)
		self.molecule_surf = pg.surface.Surface(self.size, pg.locals.HWSURFACE | pg.locals.DOUBLEBUF | pg.locals.RESIZABLE | pg.locals.SRCALPHA)
		self.camera_position = [0,0]
		self.camera_orientation = [0,0,0]
		self.camera_z = 6
		self.time = 0
		self.screen_center = (self.size[0]/2, self.size[1]/2)
		self.draw_mode = 'normal'
		self.set_projection_plane()


	def prepare_atom_bonds_imgs(self, mol, pos=(.3,.3)):
		def gaussian(size, pos, m, one_dim=False):
			x, y = np.meshgrid(np.linspace(0,1,size[1]) - pos[0], np.linspace(0,1,size[0]) - pos[1])
			if one_dim:
				dst = np.sqrt(y*y)
			else:
				dst = np.sqrt(x*x+y*y)
			gauss = np.exp(-((dst)/m)**2)+0.1
			return gauss

		def generate_atom_img(n, size=500):
			r = size
			c = data.ATOM_COLOURS[n]
			atom_img = np.zeros((r,r))

			#draw circle
			rr, cc = skimage.draw.disk((int(r/2),int(r/2)),int(r/2))
			atom_img[rr, cc] = 1

			#add gaussians as highlight
			gauss1 = gaussian((r,r), pos, .25)
			gauss2 = 0.5 * gaussian((r,r), pos, .75)
			gauss = (gauss1+gauss2)/np.max(gauss1+gauss2) #double gaussian highlight

			#define surface to draw atom_img to
			surf = pg.surface.Surface((r,r))
			a = gauss*atom_img
			a = (a/a.max() * 255).astype(int) #normalize
			pg.pixelcopy.array_to_surface(surf, a + a*256 + a*256*256) #code for rgb in pygame internal
			atom_img = surf 

			atom_img.set_colorkey((0,0,0))
			atom_img.fill(c, special_flags=pg.locals.BLEND_RGB_MULT)

			return atom_img

		def generate_single_bond_img(n1, n2, size=20):
			size = (size,20)
			half_size = (size[0], size[1]//2)

			c1 = data.ATOM_COLOURS[n1]
			c2 = data.ATOM_COLOURS[n2]

			#generate bond_img
			bond_img = np.ones(half_size)

			gauss1 = gaussian(half_size, (.5/2,.5), .25, one_dim=True)
			gauss2 = 0.5 * gaussian(half_size, (.5/2,.5), .75, one_dim=True)
			gauss = (gauss1+gauss2)/np.max(gauss1+gauss2) #double gaussian highlight

			surf = pg.surface.Surface(size, pg.locals.SRCALPHA)
			surf_half = pg.surface.Surface(half_size, pg.locals.SRCALPHA)
			a = gauss*bond_img
			a = (a/a.max() * 255).astype(int) #normalize

			pg.pixelcopy.array_to_surface(surf_half, a + a*256 + a*256*256 + a*256*256*256)
			surf1 = surf_half.copy()
			surf2 = surf_half.copy()

			surf1.fill(c1, special_flags=pg.locals.BLEND_RGB_MULT)
			surf2.fill(c2, special_flags=pg.locals.BLEND_RGB_MULT)

			surf.blit(surf1, (0, 0))
			surf.blit(surf2, (0, size[1]//2))

			return surf

		def generate_double_bond_img(n1, n2, size=40, spacing=4):
			size = (size,20)
			sbn_size = ((size[0])//2, 20)
			sbn_img = generate_single_bond_img(n1, n2, size[0])

			surf = pg.surface.Surface(size, pg.locals.SRCALPHA)
			sbn_img = pg.transform.scale(sbn_img, sbn_size)

			surf.blit(sbn_img, (0, 0))
			surf.blit(sbn_img, (sbn_size[0], 0))

			return surf

		def generate_triple_bond_img(n1, n2, size=60, spacing=4):
			size = (size,20)
			sbn_size = ((size[0])//3, 20)

			sbn_img = generate_single_bond_img(n1, n2, size[0])

			surf = pg.surface.Surface(size, pg.locals.SRCALPHA)
			sbn_img = pg.transform.scale(sbn_img, sbn_size)

			surf.blit(sbn_img, (0, 0))
			surf.blit(sbn_img, (sbn_size[0], 0))
			surf.blit(sbn_img, (sbn_size[0]*2, 0))

			return surf

		anum = mol.unique_atom_numbers()

		at_imgs = {}
		for n in anum:
			im = generate_atom_img(n)
			at_imgs[n] = im

		single_bn_imgs = {}
		double_bn_imgs = {}
		triple_bn_imgs = {}
		for n1 in anum:
			for n2 in anum:
				single_bn_imgs[(n1,n2)] = generate_single_bond_img(n1, n2, 50)
				double_bn_imgs[(n1,n2)] = generate_double_bond_img(n1, n2, 50)
				triple_bn_imgs[(n1,n2)] = generate_triple_bond_img(n1, n2, 50)

		self.atom_imgs = at_imgs
		self.single_bond_imgs = single_bn_imgs
		self.double_bond_imgs = double_bn_imgs
		self.triple_bond_imgs = triple_bn_imgs


	def get_rotation_matrix(self, x=None, y=None, z=None):
		R = np.eye(3)

		if not x is None:
			c = cos(x)
			s = sin(x)
			R = R @ np.array(([ 1, 0, 0],
						      [ 0, c,-s],
						      [ 0, s, c]))

		if not y is None:
			c = cos(y)
			s = sin(y)
			R = R @ np.array(([ c, 0, s],
						      [ 0, 1, 0],
						      [-s, 0, c]))

		if not z is None:
			c = cos(z)
			s = sin(z)
			R = R @ np.array(([ c,-s, 0],
						   	  [ s, c, 0],
						   	  [ 0, 0, 1]))

		return R


	def project(self, points, settings={}):
		sg = settings.get
		mode = sg('projection_mode', 'perspective')

		if mode == 'orthographic':
			return (points[:,:2]*100 - r/2 + self.camera_position).astype(int)

		if mode == 'perspective':
			d = self.get_rotation_matrix(*self.camera_orientation) @ (points - np.asarray([*self.camera_position, self.camera_z])).T
			f = self.projection_plane @ d

			return np.vstack((f[0]/f[2], f[1]/f[2])).T.astype(int)


	def set_projection_plane(self):
		e = np.array([self.size[0]/2, self.size[1]/2, 600])
		self.projection_plane = np.array([[1, 0, e[0]/e[2]], [0, 1, e[1]/e[2]], [0, 0, 1/e[2]]])


	def draw_molecule(self, mol, settings={}):
		self.prepare_atom_bonds_imgs(mol)

		state = {}
		state['run'] = True
		state['main_mol'] = mol

		md = self.main_display
		ms = self.molecule_surf

		self.init_loop(state)

		while state['run']:
			md.fill(self.background_color)
			self.pre_update(state)
			self.update(state)
			self.post_update(state)

			md.blit(ms, (0,0))
			pg.display.update()


	def _prepare_molecule_surf(self, mol, smooth_bonds=False, atom_radius_factor=.5):
		def rotate_image(image, pos, originPos, angle):
			image_rect = image.get_rect(topleft = (pos[0]-originPos[0], pos[1]-originPos[1]))
			offset_center_to_pivot = pg.math.Vector2(list(pos)) - image_rect.center

			rotated_offset = offset_center_to_pivot.rotate(-angle)

			rotated_image_center = (pos[0] - rotated_offset.x, pos[1] - rotated_offset.y)

			rotated_image = pg.transform.rotate(image, angle)
			rotated_image_rect = rotated_image.get_rect(center = rotated_image_center)

			return rotated_image, rotated_image_rect

		if self.draw_mode.lower() == 'normal':
			draw_bond_rect = False
			draw_bond_center = False
			draw_bond_through_atom = False
			simple_atoms = False
			simple_bonds = False
		if self.draw_mode.lower() == 'simple':
			draw_bond_rect = False
			draw_bond_center = False
			draw_bond_through_atom = False
			simple_atoms = True
			simple_bonds = True
		if self.draw_mode.lower() == 'debug':
			draw_bond_rect = True
			draw_bond_center = True
			draw_bond_through_atom = True
			simple_atoms = True
			simple_bonds = True

		blits = []
		pos = mol.positions
		dist_to_cam = np.sqrt(np.sum((pos - (*self.camera_position, self.camera_z))**2, axis=1))
		idx = np.argsort(dist_to_cam)[::-1]
		atn = mol.atom_numbers
		atcolours = mol.colours

		# zero = self.project(np.array([0,0,0]))
		# r = (mol.radii/dist_to_cam * 800).astype(int)
		rad = np.zeros((len(mol.radii),3))
		rad[:,0] = 1
		rad = rad * mol.radii.reshape(-1,1) * atom_radius_factor
		offset_pos = rad + pos		
		
		mapped_positions = self.project(pos)
		mapped_radii = self.project(offset_pos)
		r = abs((mapped_radii - mapped_positions)[:,0])

		for i in idx:
			n = atn[i]
			a = mapped_positions[i]

			
			if simple_atoms: 
				pg.draw.circle(self.molecule_surf, atcolours[i], a, r[i])
			else:
				atom_img = self.atom_imgs[n]
				atom_img = pg.transform.scale(atom_img, (int(r[i]*2), int(r[i]*2)))
				blits.append((atom_img, a-r[i]))

		#moving on to bonds
		tuples = np.asarray(mol.bond_tuples)
		ra = mapped_positions[tuples[:,0]]
		rb = mapped_positions[tuples[:,1]]
		rabn = (rb-ra)
		rabn2 = np.clip(rabn**2, 0, 2**16)
		rabn2sum = np.sum(rabn2, axis=1)
		bond_dist = np.sqrt(rabn2sum).astype(int)
		B = mol.B

		atom_indices_in_blits = [i for i in idx]
		for i, bond in enumerate(mol.bond_tuples):
			try:
				a,b = bond

				n1 = atn[a]
				n2 = atn[b]

				radius = int((r[b]+r[a])/5)

				pab = (pos[b] - pos[a])
				npab = pab/np.linalg.norm(pab)
				X = pos[a] + npab * mol.radii[a] * atom_radius_factor
				mapped_on_sphere1 = self.project(X)[0]
				X = pos[b] - npab * mol.radii[b] * atom_radius_factor
				mapped_on_sphere2 = self.project(X)[0]



				# new_scale = (radius * B[a,b], np.clip(bond_dist[i], 0, self.size[0]))
				# new_scale = (radius * B[a,b], np.clip(int(np.linalg.norm(mapped_on_sphere2 - mapped_on_sphere1)), 0, self.size[0]))
				mapped_bond_dist = np.linalg.norm(mapped_on_sphere1 - mapped_on_sphere2)
				# Rab = (mapped_positions[a] - mapped_positions[b])
				# mapped_bond_dist = np.linalg.norm(Rab)
				bond_len = max(0, mapped_bond_dist)

				if bond_len > 0:
					bond_dir = (mapped_on_sphere1 - mapped_on_sphere2)/bond_len
					bond_center = mapped_positions[a] + (mapped_positions[b] - mapped_positions[a])/2
					new_scale = (radius * B[a,b], int(bond_len))
					bond_pos = mapped_positions[a] + bond_dir * r[a]

					if B[a,b] == 1:
						bond_img = self.single_bond_imgs[(n1,n2)].copy()
					if B[a,b] == 2:
						bond_img = self.double_bond_imgs[(n1,n2)].copy()
					if B[a,b] == 3:
						bond_img = self.triple_bond_imgs[(n1,n2)].copy()

					new_scale = (max(1,new_scale[0]), max(1,new_scale[1]))
					if smooth_bonds: 	
						bond_img = pg.transform.smoothscale(bond_img, new_scale)
					else: 				
						bond_img = pg.transform.scale(bond_img, new_scale)

					#insert the bond_img in the right place (after furthest atom)
					ai = atom_indices_in_blits.index(a)
					bi = atom_indices_in_blits.index(b)
					index = ai if ai < bi else bi
					atom_indices_in_blits.insert(index, None)

					angle = pg.math.Vector2([0,1]).angle_to(rabn[i])
					im, p = rotate_image(bond_img, mapped_on_sphere1, (bond_img.get_width()/2, 0), -angle)
					if draw_bond_rect:
						pg.draw.rect(self.molecule_surf, (0,255,0), p, width=2)

					if self.show_fig:

						plt.imshow(pg.PixelArray(im).transpose())
						plt.show()

					if not simple_bonds: 
						blits.insert(index+1, (im, p.topleft))
					else: 
						pg.draw.line(self.molecule_surf, (255,255,255), ra[i], rb[i], width=2)

				if draw_bond_through_atom: 
					pg.draw.line(self.molecule_surf, (255,0,255), ra[i], mapped_on_sphere1, width=5)
					pg.draw.line(self.molecule_surf, (255,0,255), rb[i], mapped_on_sphere2, width=5)
				if draw_bond_center:
					pg.draw.circle(self.molecule_surf, (255,255,255), bond_center, 10)

			except:
				raise

		self.molecule_surf.blits(blits)

		# for p, o in zip(mapped_positions, mapped_radii):
		# 	pg.draw.line(self.molecule_surf, (0,0,255), p, o, width=5)


	def init_loop(self, state):
		state['rot'] = np.array([0,0])
		state['rotation'] = np.array([0,0])
		state['zoom'] = 0
		state['fpss'] = []
		state['fps_num'] = 100
		state['dT'] = 0
		state['show_fps'] = True
		state['time'] = 0


	def pre_update(self, state):
		state['start_time'] = perf_counter()
		state['keys'] = pg.key.get_pressed()
		state['events'] = pg.event.get()
		self.molecule_surf.fill(self.background_color)
		self.handle_events(state)
		# state['main_mol'].shake(0.1)
		state['main_mol'].center()


	def update(self, state):
		if hasattr(state['main_mol'], 'frames'):
			i = state.get('mol_frame_i', 0)
			i = i % len(state['main_mol'].frames)
			state['main_mol'].positions = state['main_mol'].frames[i]
			state['mol_frame_i'] = i + 1

		state['main_mol'].rotate(x=state['rot'][0], y=state['rot'][1])
		self._prepare_molecule_surf(state['main_mol'])


	def post_update(self, state):
		state['rotation'] = state['rotation'] + state['rot']
		state['rot'] = state['rot'] * 0.8
		state['zoom'] = 0
		state['fpss'].append(1/(perf_counter()-state['start_time']))
		state['dT'] = perf_counter()-state['start_time']
		state['time'] += state['dT']

		if len(state['fpss']) > state['fps_num']:
			state['fpss'].pop(0)
		# if state['show_fps']: print(f"fps (avg. over {state['fps_num']}) = {sum(state['fpss'])/state['fps_num'] :.0f}")

		# self.camera_orientation = (0, state['time']/2, 0)
		# self.camera_position = [sin(state['time']*10), cos(state['time']*10)]


	def handle_events(self, state):
		if state['keys'][pg.K_ESCAPE]:
			state['run'] = False
		if state['keys'][pg.K_SPACE]:
			self.show_fig = True
		else:
			self.show_fig = False

		for e in state['events']:
			if e.type == pg.QUIT:
				state['run'] = False

			elif e.type == pg.MOUSEBUTTONDOWN:
				if e.button == 4:
						state['zoom'] = -state['dT'] * self.camera_z * 10
				elif e.button == 5:
						state['zoom'] = state['dT'] * self.camera_z * 10

		# print(state['zoom'], state['dT'])
		self.camera_z += state['zoom']

		move = pg.mouse.get_rel()
		if state['keys'][pg.K_LCTRL] or state['keys'][pg.K_RCTRL]:
			if pg.mouse.get_pressed()[2]:
				self.camera_position[0] += move[0]/100
				self.camera_position[1] += move[1]/100

			if pg.mouse.get_pressed()[0]:
				state['rot'] = np.array([move[1]/150, -move[0]/150])

	
	def draw_axes(self, surf):
		...


