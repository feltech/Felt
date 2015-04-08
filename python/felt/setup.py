from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

project_dir = os.environ.get("PROJECT_DIR")

extensions = [
	Extension(
		"tests.lib.PyStubTemplate", 
		sources=[
			"tests/lib/PyStubTemplate.pyx"
		],
		include_dirs=[
			os.path.join(project_dir, "tests"), 
		],
		language="c++"
	),
]

setup(
	name="felt",
	ext_modules=cythonize(extensions)
)