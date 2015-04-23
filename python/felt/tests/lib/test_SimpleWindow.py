# import pydevd

import unittest2
from felt.lib.rocket import rocket
from felt.lib.sfml import sf

from felt.lib.PyRktSFMLInterface import PyRktSFMLInterface as RktSFML
#
class Test_rktSFML(unittest2.TestCase):
	
	def test_open_window_load_document(self):
		window = sf.RenderWindow(sf.VideoMode(640, 480), "pySFML Window")
		rktSFML = RktSFML(window)
		
		# Necessary to ensure context is unloaded before rocket-SFML interface.
		def doDefaultContext():
		
			context = rocket.CreateContext("default", rocket.Vector2i(640, 480))
			rktSFML.rktInit(context)
			
		# 	pydevd.settrace()
			rktSFML.rktInitDebugger()	
			document = context.LoadDocument('assets/blank.rml')
			
			self.assertIsNotNone(document)
			
			document.Show()
			
			self.assertTrue(window.is_open)
		
			context.Update()
			window.clear()
			context.Render()
			window.display()
			
			window.close()
			
			
			self.assertEqual(document.inner_rml.strip(), "This is blank")	
					
# 			self.assertTrue(window.is_open)
# 			while window.is_open:		
# 				for event in window.events:
# 	 				rktSFML.rktEvent(event)
# 					if type(event) is sf.CloseEvent:
# 						window.close()
# 	 			context.Update()
# 				window.clear()
# 	 			context.Render()
# 				window.display()				
			
					
		doDefaultContext()

	
if __name__ == "__main__":
	unittest2.main()