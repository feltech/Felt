#include "catch.hpp"

#define _TESTING

#include <Felt/MultiTrackedGrid.hpp>

using namespace felt;

SCENARIO("MultiTrackedGrid")
{
	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset and background value of 0")
	{
		MultiTrackedGrid<FLOAT, 3, 3> grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), 0);

		THEN("the grid size is as expected and is initialised to all zero")
		{
			CHECK(grid.data().size() == 9*9*9);
			for (const FLOAT val : grid.data())
				CHECK(val == 0);
		}

		THEN("the associated lookup grid's size is as expected and initialised to NULL indices")
		{
			CHECK(grid.lookup().data().size() == 9*9*9);
			for (const Vec3u& val : grid.lookup().data())
				CHECK(val == Vec3u::Constant(grid.lookup().NULL_IDX));
		}
	}
}
