#include <cmath>

template <typename T>
class TestCppTemplate
{
protected:
	T m_val;
	unsigned m_size;
public:
	TestCppTemplate() : m_val(0), m_size(0) {}

	TestCppTemplate(const T& val) : m_val(val), m_size(0) {}

	const unsigned getAbs() const
	{
		return m_size;
	}

	void calcSize()
	{
		m_size = std::abs(m_val);
	}

	void setVar(const T& val)
	{
		m_val = val;
	}
};
