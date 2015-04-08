cdef extern from "StubTemplate.cpp":
	cdef cppclass TestCppTemplate[T]:
		TestCppTemplate()
		TestCppTemplate(const T&)
		const unsigned getAbs() const
		void calcSize()
		void setVar(const T&)
		

def fn(val):
	cdef TestCppTemplate[float] cppObj
	cppObj.setVar(val)
	cppObj.calcSize()
	return cppObj.getAbs()