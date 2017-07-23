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

GIVEN("a 2-layer 2D surface in a 9x9 isogrid with a single 9x9 spatial partition")
{
	// 2D surface with 2x narrow band layers, respectively.
	using SurfaceType = Surface<2, 2>;

	// Construct the surface.
	SurfaceType surface(Vec2i(9, 9), Vec2i(9, 9));


	WHEN("a singularity seed is created at the centre")
	{
		// Create seed point in the centre
		surface.seed(Vec2i(0, 0));

		INFO(stringify_grid_slice(surface.isogrid()));

		// Trivially check centre of seed is indeed a zero-level point (i.e. point
		// on the surface).
		THEN("the value at the centre of the grid is 0")
		{
			const FLOAT val_centre = surface.isogrid().get(Vec2i(0, 0));
			CHECK(val_centre == 0);
		}

		// Grid to use for checking data.
		Impl::Grid::Snapshot<FLOAT, 2> isogrid_check(Vec2i(5, 5), Vec2i::Zero(), 0);

		THEN("the surface data matches a singularity seed point")
		{
			// A 2D 2-layer singularity (seed) point should look like the following.
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.vdata() -= surface.isogrid().snapshot()->vdata();

			const FLOAT diff = isogrid_check.vdata().sum();

			CHECK(diff == 0);

			// Check appropriate points have been added to narrow band layers.
			CHECK(surface.layer(0, -2).size() == 0);
			CHECK(surface.layer(0, -1).size() == 0);
			CHECK(surface.layer(0, 0).size() == 1);
			CHECK(surface.layer(0, 1).size() == 4);
			CHECK(surface.layer(0, 2).size() == 8);
		}
	}
}
}

