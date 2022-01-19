import glfw
from OpenGL.GL import *
from OpenGL.GLU import *
from time import perf_counter
import molviewer.text as text


class Screen:
	def __init__(self, **kwargs):
		self.screen_width = kwargs.get('screen_width', 600)
		self.screen_height = kwargs.get('screen_height', 400)
		self.screen_title = kwargs.get('screen_title', 'screen')


	def mainloop(self, world):
		glfw.init()
		
		gluPerspective(60, 10, 0, 1)
		self.window = glfw.create_window(self.screen_width, self.screen_height, self.screen_title, None, None)
		glfw.make_context_current(self.window)

		# text.initialize()
		if not self.window:
			glfw.terminate()
			raise Exception("glfw window creation FAILED!")
		glfw.make_context_current(self.window)

		glClearColor(.2, .2, .2, 1)


		world.start()


		while not glfw.window_should_close(self.window):
			start = perf_counter()
			glfw.poll_events()

			glClear(GL_COLOR_BUFFER_BIT)

			ct = glfw.get_time()  # returns the elapsed time, since init was called

			for obj in world:
				obj.draw()


			glfw.swap_buffers(self.window)
			# print(f'FPS={1/(perf_counter()-start):.1f}')

		# terminate glfw, free up allocated resources
		glfw.terminate()



# class Camera():

#     def __init__(
#         self,
#         eye=None, target=None, up=None,
#         fov=None, near=0.1, far=100000
#     ):
#         self.eye = eye or glm.vec3(0, 0, 1)
#         self.target = target or glm.vec3(0, 0, 0)
#         self.up = up or glm.vec3(0, 1, 0)
#         self.original_up = glm.vec3(self.up)
#         self.fov = fov or glm.radians(45)
#         self.near = near
#         self.far = far

#     def update(self, aspect):
#         self.view = glm.lookAt(
#             self.eye, self.target, self.up
#         )
#         self.projection = glm.perspective(
#             self.fov, aspect, self.near, self.far
#         )

#     def rotate_target(self, delta):
#         right = glm.normalize(glm.cross(self.target - self.eye, self.up))
#         M = glm.mat4(1)
#         M = glm.translate(M, self.eye)
#         M = glm.rotate(M, delta.y, right)
#         M = glm.rotate(M, delta.x, self.up)
#         M = glm.translate(M, -self.eye)
#         self.target = glm.vec3(M * glm.vec4(self.target, 1.0))

#     def rotate_around_target(self, target, delta):
#         right = glm.normalize(glm.cross(self.target - self.eye, self.up))
#         amount = (right * delta.y + self.up * delta.x)
#         M = glm.mat4(1)
#         M = glm.rotate(M, amount.z, glm.vec3(0, 0, 1))
#         M = glm.rotate(M, amount.y, glm.vec3(0, 1, 0))
#         M = glm.rotate(M, amount.x, glm.vec3(1, 0, 0))
#         self.eye = glm.vec3(M * glm.vec4(self.eye, 1.0))
#         self.target = target
#         self.up = self.original_up

#     def rotate_around_origin(self, delta):
#         return self.rotate_around_target(glm.vec3(0), delta)


# class GlutController():

#     FPS = 0
#     ORBIT = 1

#     def __init__(self, camera, velocity=100, velocity_wheel=100):
#         self.velocity = velocity
#         self.velocity_wheel = velocity_wheel
#         self.camera = camera

#     def glut_mouse(self, button, state, x, y):
#         self.mouse_last_pos = glm.vec2(x, y)
#         self.mouse_down_pos = glm.vec2(x, y)

#         if button == GLUT_LEFT_BUTTON:
#             self.mode = self.FPS
#         elif button == GLUT_RIGHT_BUTTON:
#             self.mode = self.ORBIT

#     def glut_motion(self, x, y):
#         pos = glm.vec2(x, y)
#         move = self.mouse_last_pos - pos
#         self.mouse_last_pos = pos

#         if self.mode == self.FPS:
#             self.camera.rotate_target(move * 0.005)
#         elif self.mode == self.ORBIT:
#             self.camera.rotate_around_origin(move * 0.005)


# class Screen:

#     def __init__(self, **kwargs):
#         self.width = kwargs.get('width', 600)
# 		self.height = kwargs.get('height', 400)
# 		screen_title = kwargs.get('screen_title', 'screen')

#         glutInit()
#         glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
#         glutInitWindowSize(w, h)
#         glutCreateWindow(screen_title)

#         self.startup()

#         glutReshapeFunc(self.reshape)
#         glutDisplayFunc(self.display)
#         glutMouseFunc(self.controller.glut_mouse)
#         glutMotionFunc(self.controller.glut_motion)
#         glutIdleFunc(self.idle_func)

#     def startup(self):
#         glEnable(GL_DEPTH_TEST)

#         aspect = self.width / self.height
#         self.camera = Camera(
#             eye=glm.vec3(10, 10, 10),
#             target=glm.vec3(0, 0, 0),
#             up=glm.vec3(0, 1, 0)
#         )
#         self.model = glm.mat4(1)
#         self.controller = GlutController(self.camera)

#     def mainloop(self):
#         glutMainLoop()

#     def idle_func(self):
#         glutPostRedisplay()

#     def reshape(self, w, h):
#         glViewport(0, 0, w, h)
#         self.width = w
#         self.height = h

#     def display(self):
#         self.camera.update(self.width / self.height)

#         glClearColor(0.2, 0.3, 0.3, 1.0)
#         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

#         glMatrixMode(GL_PROJECTION)
#         glLoadIdentity()
#         gluPerspective(glm.degrees(self.camera.fov), self.width / self.height, self.camera.near, self.camera.far)
#         glMatrixMode(GL_MODELVIEW)
#         glLoadIdentity()
#         e = self.camera.eye
#         t = self.camera.target
#         u = self.camera.up
#         gluLookAt(e.x, e.y, e.z, t.x, t.y, t.z, u.x, u.y, u.z)
#         glColor3f(1, 1, 1)
#         glBegin(GL_LINES)
#         for i in range(-5, 6):
#             if i == 0:
#                 continue
#             glVertex3f(-5, 0, i)
#             glVertex3f(5, 0, i)
#             glVertex3f(i, 0, -5)
#             glVertex3f(i, 0, 5)
#         glEnd()

#         glBegin(GL_LINES)
#         glColor3f(1, 0, 0)
#         glVertex3f(-5, 0, 0)
#         glVertex3f(5, 0, 0)
#         glColor3f(0, 1, 0)
#         glVertex3f(0, -5, 0)
#         glVertex3f(0, 5, 0)
#         glColor3f(0, 0, 1)
#         glVertex3f(0, 0, -5)
#         glVertex3f(0, 0, 5)
#         glEnd()

#         glutSwapBuffers()