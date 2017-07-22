#include "../catch.hpp"
#include <Felt/Impl/Common.hpp>
#include "Utils.hpp"
#include <Felt/Impl/Surface.hpp>

SCENARIO("Impl::Surface")
{

using namespace Felt;

GIVEN("a 2-layer 2D surface in a 7x7 isogrid with 3x3 spatial partitions")
{
	Surface<2, 2> surface(Vec2i(7, 7), Vec2i(3, 3));

	THEN("the isogrid is initialised correctly")
	{
		CHECK(surface.isogrid().size() == Vec2i(7, 7));
		CHECK(surface.isogrid().children().data().size() == 9);
		CHECK(surface.isogrid().children().get(Vec2i(0, 0)).size() == Vec2i(3, 3));
		CHECK(surface.isogrid().children().get(Vec2i(0, 0)).data().size() == 0);
		CHECK(surface.isogrid().size() == Vec2i(7, 7));
		// Grid is initialised to all points 'outside' the surface (since there is no surface yet).
		CHECK(surface.isogrid().get(Vec2i(0, 0)) == 3);
	}
}

}

