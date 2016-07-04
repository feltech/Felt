#include "catch.hpp"

#define _TESTING

#include <Felt/MultiTrackedGrid.hpp>

using namespace felt;

SCENARIO("MultiTrackedGrid")
{
	WHEN("initialisation")
	{
		// ==== Setup/Action ====
		MultiTrackedGrid<FLOAT, 3, 3> grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), 0);

		// ==== Confirm ====
		CHECK(grid.data().size() == 9*9*9);
		CHECK(grid.lookup().data().size() == 9*9*9);

		for (const FLOAT val : grid.data())
			CHECK(val == 0);

		for (const Vec3u& val : grid.lookup().data())
			CHECK(val == Vec3u::Constant(grid.lookup().NULL_IDX));

	}
}
