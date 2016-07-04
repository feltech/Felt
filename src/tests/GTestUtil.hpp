#ifndef GTESTUTIL_HPP
#define GTESTUTIL_HPP

// Borrowed from http://stackoverflow.com/a/29155677/535103

namespace testing
{
	namespace internal
	{
		enum GTestColor {
			COLOR_DEFAULT,
			COLOR_RED,
			COLOR_GREEN,
			COLOR_YELLOW
		};

		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
	}
}
#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__);  testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, "\n"); } while(0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
	~TestCout()
	{
		PRINTF("%s",str().c_str());
	}
};

#define TEST_COUT  TestCout()
#endif //GTESTUTIL_HPP
