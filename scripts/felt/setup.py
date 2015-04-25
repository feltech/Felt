from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

project_dir = os.environ.get("PROJECT_ROOT")

extensions = [
	Extension(
		"tests.lib.PyStubTemplate", 
		sources=[
			"tests/lib/PyStubTemplate.pyx"
		],
		include_dirs=[
			os.path.join(project_dir, "src/tests"), 
		],
		language="c++"
	),
	Extension(
		"lib.PyRktSFMLInterface", 
		sources=[
			"lib/PyRktSFMLInterface.pyx"
		],
		include_dirs=[
			os.path.join(project_dir, "include"),
			os.path.join(project_dir, "vendor/include")
		],
		library_dirs=[
			os.path.join(project_dir, "lib"),
			os.path.join(project_dir, "vendor/lib")
		],
		libraries=[
			"RktSFMLInterface", "RocketCore", "RocketControls", 
			"RocketDebugger", "sfml-window", 
			"sfml-system", "sfml-graphics"
		],
		language="c++",
		extra_compile_args=["-g"],
		extra_link_args=["-g"]
	),			
]


setup(
	name="felt",
	ext_modules=cythonize(
		extensions, gdb_debug=True, 
		output_dir=os.path.join(project_dir, "scripts/felt/lib")
	)
)