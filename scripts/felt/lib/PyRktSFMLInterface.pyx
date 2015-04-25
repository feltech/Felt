from libcpp cimport bool
from libcpp.string cimport string 
from libcpp.event cimport EventType
from libcpp.sfml cimport RenderWindow as cyWindow
from libcpp.sfml cimport Event as cyEvent
from pysfml.window cimport Window as pyWindow
from pysfml.window cimport Event as pyEvent


cdef extern from "gui/RktSFMLInterface.hpp" namespace "felt":
	cdef cppclass RktSFMLInterface:
		RktSFMLInterface(cyWindow* psfWindow)
		void rktInit(string sfContextName)
		void rktEvent(cyEvent* sfEvent)
		void rktInitDebugger()
	
	
cdef class PyRktSFMLInterface:

	cdef RktSFMLInterface* rktSFML

	def __cinit__(self, pyWindow sfWindow):
		self.rktSFML = new RktSFMLInterface(<cyWindow*>sfWindow.p_window)
		
	def __dealloc__(self):
		del self.rktSFML

	def rktInit(self, rktContext):
		self.rktSFML.rktInit(rktContext.name)
		
	def rktEvent(self, pyEvent pyevent):
		self.rktSFML.rktEvent(<cyEvent*>pyevent.p_this)
		
	def rktInitDebugger(self):
		self.rktSFML.rktInitDebugger()