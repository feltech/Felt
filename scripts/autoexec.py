import unittest2
def main():
# 	import pydevd
# 	pydevd.settrace()
	unittest2.main(verbosity=2, module="felt.tests.lib.test_SimpleWindow")
	print("Done.")