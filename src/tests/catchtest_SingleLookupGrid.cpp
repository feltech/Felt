#include "catch.hpp"

#define _TESTING

#include <Felt/SingleLookupGrid.hpp>

using namespace felt;

SCENARIO("SingleLookupGrid")
{
	GIVEN("a grid and some locations")
	{
		using GridType = EagerSingleLookupGrid<3, 3>;
		GridType grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);
		const Vec3i pos7(7, 0, 0);

		WHEN("we add 4 locations to be tracked")
		{
			// Add the positions to the array and set index lookup values.
			grid.add(pos1, 0);
			grid.add(pos2, 1);
			grid.add(pos3, 1);
			grid.add(pos4, 2);


			THEN("the tracking lists contain the expected number of elements")
			{
				CHECK(grid.list(0).size() == 1);
				CHECK(grid.list(1).size() == 2);
				CHECK(grid.list(2).size() == 1);
			}

			THEN("the tracking list elements contain the position vectors")
			{
				CHECK(grid.list(0)[0] == pos1);
				CHECK(grid.list(1)[0] == pos2);
				CHECK(grid.list(1)[1] == pos3);
				CHECK(grid.list(2)[0] == pos4);
			}

			THEN("the grid contains the indices of the position vectors in the tracking list")
			{
				CHECK((UINT)grid(pos1) == 0);
				CHECK((UINT)grid(pos2) == 0);
				CHECK((UINT)grid(pos3) == 1);
				CHECK((UINT)grid(pos4) == 0);
			}
			AND_WHEN("we remove a position that is not tracked")
			{
				grid.remove(pos7, 1);

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 2);
					CHECK(grid.list(2).size() == 1);
				}

				THEN("the tracking list elements contain the position vectors")
				{
					CHECK(grid.list(0)[0] == pos1);
					CHECK(grid.list(1)[0] == pos2);
					CHECK(grid.list(1)[1] == pos3);
					CHECK(grid.list(2)[0] == pos4);
				}

				THEN("the grid contains the indices of the position vectors in the tracking list")
				{
					CHECK((UINT)grid(pos1) == 0);
					CHECK((UINT)grid(pos2) == 0);
					CHECK((UINT)grid(pos3) == 1);
					CHECK((UINT)grid(pos4) == 0);
				}
			}

			AND_WHEN("we remove a position vector from tracking in list 0")
			{
				grid.remove(pos2, 1);

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 1);
					CHECK(grid.list(2).size() == 1);
				}

				THEN(
				"the tracking list elements still contain the remaining position vectors,"
				" with the remaining position from list 1 having changed index"
				) {
					CHECK(grid.list(0)[0] == pos1);
					CHECK(grid.list(1)[0] == pos3);
					CHECK(grid.list(2)[0] == pos4);
				}

				THEN("the grid location corresponding to the removed point gives NULL index")
				{
					CHECK((UINT)grid(pos1) == 0);
					CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
					CHECK((UINT)grid(pos3) == 0);
					CHECK((UINT)grid(pos4) == 0);
				}

				AND_WHEN("we add two more points")
				{

					grid.add(pos5, 2);
					grid.add(pos6, 2);

					THEN("the tracking lists contain the expected number of elements")
					{
						CHECK(grid.list(0).size() == 1);
						CHECK(grid.list(1).size() == 1);
						CHECK(grid.list(2).size() == 3);
					}

					THEN("the tracking list elements contain the position vectors")
					{
						CHECK(grid.list(0)[0] == pos1);
						CHECK(grid.list(1)[0] == pos3);
						CHECK(grid.list(2)[0] == pos4);
						CHECK(grid.list(2)[1] == pos5);
						CHECK(grid.list(2)[2] == pos6);
					}

					THEN(
					"the grid contains the indices of the position vectors in the tracking"
					" list"
					) {
						CHECK((UINT)grid(pos1) == 0);
						CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
						CHECK((UINT)grid(pos3) == 0);
						CHECK((UINT)grid(pos4) == 0);
						CHECK((UINT)grid(pos5) == 1);
						CHECK((UINT)grid(pos6) == 2);
					}


					AND_WHEN("we remove a point by index and another point by location")
					{
						grid.remove(pos4, 2);
						grid.remove(0, 0);

						THEN("the tracking lists contain the expected number of elements")
						{
							CHECK(grid.list(0).size() == 0);
							CHECK(grid.list(1).size() == 1);
							CHECK(grid.list(2).size() == 2);
						}

						THEN("the tracking list elements contain the position vectors")
						{
							CHECK(grid.list(0)[0] == pos1);
							CHECK(grid.list(1)[0] == pos3);
							CHECK(grid.list(2)[1] == pos5);
							CHECK(grid.list(2)[0] == pos6);
						}

						THEN(
						"the grid contains the indices of the position vectors in the tracking"
						" list, with the two removed points having NULL index"
						) {
							CHECK((UINT)grid(pos1) == GridType::NULL_IDX);
							CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
							CHECK((UINT)grid(pos3) == 0);
							CHECK((UINT)grid(pos4) == GridType::NULL_IDX);
							CHECK((UINT)grid(pos5) == 1);
							CHECK((UINT)grid(pos6) == 0);
						}

					}
				}
			}

			AND_WHEN("we reset list 1")
			{
				grid.reset(1);

				THEN("list 1 is empty but the other lists are unaffected")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 0);
					CHECK(grid.list(2).size() == 1);
				}

				THEN("the locations in the grid that were in list 2 are now NULL index")
				{
					CHECK(grid(pos1) == 0);
					CHECK(grid(pos2) == GridType::NULL_IDX);
					CHECK(grid(pos3) == GridType::NULL_IDX);
					CHECK(grid(pos4) == 0);
				}
			}
		}

	}
}


SCENARIO("LazySingleLookupGrid")
{
	WHEN("initialisation")
	{
		/// [LazySingleLookupGrid initialisation]
		// ==== Setup ====
		LazySingleLookupGrid<3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
		const UINT NULL_IDX_DATA = LazySingleLookupGrid<3, 3>::NULL_IDX;

		// ==== Confirm ====
		CHECK(grid.is_active() == false);
		CHECK(grid.data().size() == 0);
		CHECK(grid.background() == NULL_IDX_DATA);
		CHECK(grid.get(Vec3i(1,1,1)) == NULL_IDX_DATA);
		/// [LazySingleLookupGrid initialisation]
	}
}
