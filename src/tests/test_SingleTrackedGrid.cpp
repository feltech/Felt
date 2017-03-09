#include "catch.hpp"

#define _TESTING

#include <Felt/TrackedGrid.hpp>

using namespace felt;

SCENARIO("LazyTrackedGrid")
{
	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 3")
	{
		LazyTrackedGrid<FLOAT, 3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3);

		const UINT NULL_IDX = LazyTrackedGrid<FLOAT, 3, 3>::Lookup::NULL_IDX;

		THEN("the data grid and associated lookup grid state is inactive")
		{
			CHECK(grid.is_active() == false);
			CHECK(grid.data().size() == 0);
			CHECK(grid.background() == 3);
			CHECK(grid.get(Vec3i(1,1,1)) == 3);
			CHECK(grid.lookup().is_active() == false);
			CHECK(grid.lookup().data().size() == 0);
			CHECK(grid.lookup().background() == NULL_IDX);
			CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
		}

		/// [LazySingleTrackedGrid activate]
		WHEN("the grid is activated")
		{
			grid.activate();

			THEN("the data grid and associated lookup grid state is active")
			{
				CHECK(grid.is_active() == true);
				CHECK(grid.data().size() == 3*3*3);
				CHECK(grid.get(Vec3i(1,1,1)) == 3);
				CHECK(grid.lookup().is_active() == true);
				CHECK(grid.lookup().data().size() == 3*3*3);
				CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
				/// [LazySingleTrackedGrid activate]
			}

			AND_WHEN("the grid is deactivated")
			{
				grid.deactivate();

				THEN("the data grid and associated lookup grid state is inactive")
				{
					CHECK(grid.is_active() == false);
					CHECK(grid.data().size() == 0);
					CHECK(grid.background() == 3);
					CHECK(grid.get(Vec3i(1,1,1)) == 3);
					CHECK(grid.lookup().is_active() == false);
					CHECK(grid.lookup().data().size() == 0);
					CHECK(grid.lookup().background() == NULL_IDX);
					CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
				}
			}
		}
	}
}
