#include "catch.hpp"

#define _TESTING

#include <Felt/MultiLookupGrid.hpp>

using namespace felt;

SCENARIO("MultiLookupGrid")
{
	GIVEN("a 10x10x10 EagerMultiLookupGrid with 3 tracking lists, and some locations")
	{
		using GridType = EagerMultiLookupGrid<3, 3>;
		GridType grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		THEN("the grid is initialised with NULL indices and the tracking lists are empty")
		{
			REQUIRE(grid.list(0).size() == 0);
			REQUIRE(grid.list(1).size() == 0);
			REQUIRE(grid.list(2).size() == 0);
			CHECK(grid(pos1)(0) == GridType::NULL_IDX);
			CHECK(grid(pos2)(0) == GridType::NULL_IDX);
			CHECK(grid(pos3)(0) == GridType::NULL_IDX);
			CHECK(grid(pos4)(0) == GridType::NULL_IDX);
			CHECK(grid(pos5)(0) == GridType::NULL_IDX);
			CHECK(grid(pos6)(0) == GridType::NULL_IDX);
		}

		WHEN("we append 4 locations to be tracked in list 0")
		{
			grid.add(pos1, 0);
			grid.add(pos2, 0);
			grid.add(pos3, 0);
			grid.add(pos4, 0);

			THEN(
				"the tracking list is populated and the grid locations hold their indices in the"
				" list"
			) {
				REQUIRE(grid.list(0).size() == 4);
				CHECK(grid.list(0)[0] == pos1);
				CHECK(grid.list(0)[1] == pos2);
				CHECK(grid.list(0)[2] == pos3);
				CHECK(grid.list(0)[3] == pos4);
				CHECK(grid(pos1)(0) == 0);
				CHECK(grid(pos2)(0) == 1);
				CHECK(grid(pos3)(0) == 2);
				CHECK(grid(pos4)(0) == 3);
			}

			AND_WHEN("we add a location that is already tracked")
			{
				grid.add(pos2, 0);

				THEN("the grid state is just as if the final point was not added")
				{
					REQUIRE(grid.list(0).size() == 4);
					CHECK(grid.list(0)[0] == pos1);
					CHECK(grid.list(0)[1] == pos2);
					CHECK(grid.list(0)[2] == pos3);
					CHECK(grid.list(0)[3] == pos4);
					CHECK(grid(pos1)(0) == 0);
					CHECK(grid(pos2)(0) == 1);
					CHECK(grid(pos3)(0) == 2);
					CHECK(grid(pos4)(0) == 3);
				}
			}

			AND_WHEN("we remove a location by it's index in the tracking list")
			{
				grid.remove(1, 0);

				THEN("the tracking list is updated and the grid location now holds NULL")
				{
					REQUIRE(grid.list(0).size() == 3);
					CHECK(grid.list(0)[0] == pos1);
					CHECK(grid.list(0)[1] == pos4);
					CHECK(grid.list(0)[2] == pos3);
					CHECK(grid(pos1)(0) == 0);
					CHECK(grid(pos2)(0) == GridType::NULL_IDX);
					CHECK(grid(pos3)(0) == 2);
					CHECK(grid(pos4)(0) == 1);
				}

				AND_WHEN("we remove a location by position in the grid")
				{
					grid.remove(pos1, 0);

					THEN("the tracking list is updated and the grid location now holds NULL")
					{
						REQUIRE(grid.list(0).size() == 2);
						CHECK(grid.list(0)[0] == pos3);
						CHECK(grid.list(0)[1] == pos4);
						CHECK(grid(pos1)(0) == GridType::NULL_IDX);
						CHECK(grid(pos2)(0) == GridType::NULL_IDX);
						CHECK(grid(pos3)(0) == 0);
						CHECK(grid(pos4)(0) == 1);
					}
				}
			}

			AND_WHEN("we reset the grid for tracking list 0")
			{
				grid.reset(0);

				THEN(
					"all added location are NULL in the grid and the tracking list contains no"
					" points"
				) {
					CHECK(grid.list(0).size() == 0);
					CHECK(grid(pos1)(0) == GridType::NULL_IDX);
					CHECK(grid(pos2)(0) == GridType::NULL_IDX);
					CHECK(grid(pos3)(0) == GridType::NULL_IDX);
					CHECK(grid(pos4)(0) == GridType::NULL_IDX);
				}
			}
		}

		WHEN("we append 4 locations to be tracked spread across all 3 lists")
		{
			grid.add(pos1, 0);
			grid.add(pos2, 1);
			grid.add(pos3, 1);
			grid.add(pos4, 2);

			THEN("the tracking lists and index tuples within the grid are updated")
			{
				CHECK(grid.list(0).size() == 1);
				CHECK(grid.list(1).size() == 2);
				CHECK(grid.list(2).size() == 1);
				CHECK(grid.list(0)[0] == pos1);
				CHECK(grid.list(1)[0] == pos2);
				CHECK(grid.list(1)[1] == pos3);
				CHECK(grid.list(2)[0] == pos4);
				CHECK(grid(pos1)(0) == 0);
				CHECK(grid(pos2)(1) == 0);
				CHECK(grid(pos3)(1) == 1);
				CHECK(grid(pos4)(2) == 0);
			}

			AND_WHEN("we remove a location from tracking list 1")
			{
				grid.remove(pos2, 1);

				THEN("the grid and tracking list is updated")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 1);
					CHECK(grid.list(2).size() == 1);
					CHECK(grid.list(0)[0] == pos1);
					CHECK(grid.list(1)[0] == pos3);
					CHECK(grid.list(2)[0] == pos4);
					CHECK(grid(pos1)(0) == 0);
					CHECK(grid(pos2)(1) == GridType::NULL_IDX);
					CHECK(grid(pos3)(1) == 0);
					CHECK(grid(pos4)(2) == 0);
				}

				AND_WHEN("we add 2 more points to tracking list 2")
				{
					grid.add(pos5, 2);
					grid.add(pos6, 2);

					THEN("the grid and tracking list is updated")
					{
						CHECK(grid.list(0).size() == 1);
						CHECK(grid.list(1).size() == 1);
						CHECK(grid.list(2).size() == 3);
						CHECK(grid.list(0)[0] == pos1);
						CHECK(grid.list(1)[0] == pos3);
						CHECK(grid.list(2)[0] == pos4);
						CHECK(grid.list(2)[1] == pos5);
						CHECK(grid.list(2)[2] == pos6);
						CHECK(grid(pos1)(0) == 0);
						CHECK(grid(pos2)(1) == GridType::NULL_IDX);
						CHECK(grid(pos3)(1) == 0);
						CHECK(grid(pos4)(2) == 0);
						CHECK(grid(pos5)(2) == 1);
						CHECK(grid(pos6)(2) == 2);
					}

					AND_WHEN("we remove 2 points from different tracking lists")
					{
						grid.remove(pos4, 2);
						grid.remove(0, 0);

						THEN("the grid and tracking lists are updated")
						{
							CHECK(grid.list(0).size() == 0);
							CHECK(grid.list(1).size() == 1);
							CHECK(grid.list(2).size() == 2);
							CHECK(grid.list(0)[0] == pos1);
							CHECK(grid.list(1)[0] == pos3);
							CHECK(grid.list(2)[1] == pos5);
							CHECK(grid.list(2)[0] == pos6);
							CHECK(grid(pos1)(0) == GridType::NULL_IDX);
							CHECK(grid(pos2)(1) == GridType::NULL_IDX);
							CHECK(grid(pos3)(1) == 0);
							CHECK(grid(pos4)(2) == GridType::NULL_IDX);
							CHECK(grid(pos5)(2) == 1);
							CHECK(grid(pos6)(2) == 0);
						}

						AND_WHEN("we reset tracking list 2")
						{
							grid.reset(2);

							THEN("only tracking list 2 is affected")
							{
								CHECK(grid.list(0).size() == 0);
								CHECK(grid.list(1).size() == 1);
								CHECK(grid.list(2).size() == 0);
								CHECK(grid.list(1)[0] == pos3);
								CHECK(grid(pos1)(0) == GridType::NULL_IDX);
								CHECK(grid(pos2)(1) == GridType::NULL_IDX);
								CHECK(grid(pos3)(1) == 0);
								CHECK(grid(pos4)(2) == GridType::NULL_IDX);
								CHECK(grid(pos5)(2) == GridType::NULL_IDX);
								CHECK(grid(pos6)(2) == GridType::NULL_IDX);
							}
						}
					}

				}
			}
		}
	}
}
