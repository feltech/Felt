#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Tracked.hpp>

#include "../catch.hpp"

using namespace Felt;

/**
 * Test the Grid class.
 */
SCENARIO("Grid::Simple")
{
	using GridType = Impl::Grid::Simple<FLOAT, 3>;
	/// [Grid - basics: GIVEN 3x7x11]
	GIVEN("a 3x7x11 grid with no offset and background value of 0")
	{
		GridType grid(Vec3i(3, 7, 11), Vec3i::Zero(), 0);
		/// [Grid - basics: GIVEN 3x7x11]

		/// [Grid - basics: THEN memory is allocated]
		THEN("memory is allocated and the size is reported correctly")
		{
			CHECK(grid.size()(0) == 3);
			CHECK(grid.size()(1) == 7);
			CHECK(grid.size()(2) == 11);
			CHECK(grid.data().size() == (3 * 7 * 11));
		}
		/// [Grid - basics: THEN memory is allocated]

		THEN("we can test if locations lie within the grid")
		{
			CHECK(grid.inside(Vec3i(-1,0,0)) == false);
			CHECK(grid.inside(Vec3i(0,0,0)) == true);
			CHECK(grid.inside(Vec3i(1,2,3)) == true);
			CHECK(grid.inside(Vec3i(3,7,11)) == false);
			CHECK(grid.inside(Vec3f(0,-0.00001,0)) == false);
			CHECK(grid.inside(Vec3f(0,0,9.99999)) == true);
		}

		WHEN("some positions values are set")
		{
			grid.set(Vec3i(0, 0, 0), 13.0f);
			grid.set(Vec3i(1, 2, 3), 17.0f);
			grid.set(Vec3i(2, 6, 10), 19.0f);

			//! [Get and set]
			THEN("querying those positions returns the same values")
			{
				CHECK(grid.get(Vec3i(1, 2, 3)) == 17.0f);
			}

			THEN("expected elements of the underlying array contain those values")
			{
				CHECK(grid.data()[0] == 13.0f);
				CHECK(grid.data()[grid.data().size() - 1] == 19.0f);
			}
			//! [Get and set]
		}
	}

	GIVEN("a 7x11x13 grid with (-3,-3,-3) offset and background value of 0")
	{
		Vec3i size(7, 11, 13);
		Vec3i offset(-3, -3, -3);
		GridType grid(size, offset, 0);

		//! [Position index]
		THEN("the index of a point in the data array is reported correctly")
		{
			CHECK(GridType::index(Vec3i(1, 0, -1), size, offset) == 613);
			CHECK(grid.index(Vec3i(1, 0, -1)) == 613);
		}

		THEN("the point represented by an index in the data array is reported correctly")
		{
			CHECK(grid.index(613) == Vec3i(1, 0, -1));
			CHECK(GridType::index(613, size, offset) == Vec3i(1, 0, -1));
		}
		//! [Position index]

		THEN("we can test if locations lie within the offset grid")
		{
			CHECK(grid.inside(Vec3i(-2,0,0)) == true);
			CHECK(grid.inside(Vec3i(-4,0,0)) == false);

		}

		WHEN("editing points in the offset grid")
		{
			grid.set(Vec3i(-3,-3,-3), 21.0f);
			grid.set(Vec3i(-1,0,-1), 23.0f);

			THEN("we can retrieve the values from the offset positions")
			{
				CHECK(grid.data()[0] == 21.0f);
				CHECK(grid.get(Vec3i(-1,0,-1)) == 23.0f);
			}
		}
	}
}


SCENARIO("Lookup::Simple")
{
	using GridType = Impl::Lookup::Simple<3>;
	GIVEN("a grid and some locations")
	{
		GridType grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

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
			grid.add(pos1);
			grid.add(pos2);
			grid.add(pos3);
			grid.add(pos4);


			THEN("the tracking lists contain the expected number of elements")
			{
				CHECK(grid.list().size() == 4);
			}

			THEN("the grid reports the active state of positions correctly")
			{
				CHECK(grid.is_active(pos1) == true);
				CHECK(grid.is_active(pos2) == true);
				CHECK(grid.is_active(pos3) == true);
				CHECK(grid.is_active(pos4) == true);
				CHECK(grid.is_active(pos5) == false);
			}

			THEN("the tracking list elements contain the positions")
			{
				CHECK(grid.list()[0] == pos1);
				CHECK(grid.list()[1] == pos2);
				CHECK(grid.list()[2] == pos3);
				CHECK(grid.list()[3] == pos4);
			}

			THEN("the grid contains the indices of the positions in the tracking list")
			{
				CHECK(grid.get(pos1) == 0);
				CHECK(grid.get(pos2) == 1);
				CHECK(grid.get(pos3) == 2);
				CHECK(grid.get(pos4) == 3);
			}
			AND_WHEN("we remove a position that is not tracked")
			{
				grid.remove(pos7);

				THEN("the tracking lists contain the same number of elements")
				{
					CHECK(grid.list().size() == 4);
				}

				THEN("the tracking list elements contain the same position vectors")
				{
					CHECK(grid.list()[0] == pos1);
					CHECK(grid.list()[1] == pos2);
					CHECK(grid.list()[2] == pos3);
					CHECK(grid.list()[3] == pos4);
				}

				THEN("the grid contains the same indices of the positions in the tracking list")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == 1);
					CHECK(grid.get(pos3) == 2);
					CHECK(grid.get(pos4) == 3);
				}
			}

			AND_WHEN("we remove a position from tracking")
			{
				grid.remove(pos2);

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list().size() == 3);
				}

				THEN(
				"the removed position is gone from the list to be replaced by the final position"
				) {
					CHECK(grid.list()[0] == pos1);
					CHECK(grid.list()[1] == pos4);
					CHECK(grid.list()[2] == pos3);
				}

				THEN("the grid location corresponding to the removed point gives NULL index")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == NULL_IDX);
					CHECK(grid.get(pos3) == 2);
					CHECK(grid.get(pos4) == 1);
				}

				AND_WHEN("we add two more points")
				{

					grid.add(pos5);
					grid.add(pos6);

					THEN("the tracking lists contain an extra of element")
					{
						CHECK(grid.list().size() == 5);
					}

					THEN("the tracking list elements contain the position vectors")
					{
						CHECK(grid.list()[0] == pos1);
						CHECK(grid.list()[1] == pos4);
						CHECK(grid.list()[2] == pos3);
						CHECK(grid.list()[3] == pos5);
						CHECK(grid.list()[4] == pos6);
					}

					THEN(
						"the grid contains the indices of the position vectors in the tracking"
						" list"
					) {
						CHECK(grid.get(pos1) == 0);
						CHECK(grid.get(pos2) == NULL_IDX);
						CHECK(grid.get(pos3) == 2);
						CHECK(grid.get(pos4) == 1);
						CHECK(grid.get(pos5) == 3);
						CHECK(grid.get(pos6) == 4);
					}
				}
			}

			AND_WHEN("we reset the list")
			{
				grid.reset();

				THEN("list is empty")
				{
					CHECK(grid.list().size() == 0);
				}

				THEN("the locations in the grid are now NULL index")
				{
					CHECK(grid.get(pos1) == NULL_IDX);
					CHECK(grid.get(pos2) == NULL_IDX);
					CHECK(grid.get(pos3) == NULL_IDX);
					CHECK(grid.get(pos4) == NULL_IDX);
				}

				THEN("the grid reports the active state of positions correctly")
				{
					CHECK(grid.is_active(pos1) == false);
					CHECK(grid.is_active(pos2) == false);
					CHECK(grid.is_active(pos3) == false);
					CHECK(grid.is_active(pos4) == false);
					CHECK(grid.is_active(pos5) == false);
				}
			}
		}

	}
}



SCENARIO("Lookup::Single")
{
	using GridType = Impl::Lookup::Single<3, 3>;

	GIVEN("a grid and some locations")
	{
		GridType grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

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

			THEN("the tracking list elements contain the positions")
			{
				CHECK(grid.list(0)[0] == pos1);
				CHECK(grid.list(1)[0] == pos2);
				CHECK(grid.list(1)[1] == pos3);
				CHECK(grid.list(2)[0] == pos4);
			}

			THEN("the grid contains the indices of the positions in the tracking list")
			{
				CHECK(grid.get(pos1) == 0);
				CHECK(grid.get(pos2) == 0);
				CHECK(grid.get(pos3) == 1);
				CHECK(grid.get(pos4) == 0);
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
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == 0);
					CHECK(grid.get(pos3) == 1);
					CHECK(grid.get(pos4) == 0);
				}
			}

			AND_WHEN("we remove a position from tracking in list 0")
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
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == Felt::NULL_IDX);
					CHECK(grid.get(pos3) == 0);
					CHECK(grid.get(pos4) == 0);
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
						"the grid contains the indices of the position vectors in the tracking list"
					) {
						CHECK(grid.get(pos1) == 0);
						CHECK(grid.get(pos2) == Felt::NULL_IDX);
						CHECK(grid.get(pos3) == 0);
						CHECK(grid.get(pos4) == 0);
						CHECK(grid.get(pos5) == 1);
						CHECK(grid.get(pos6) == 2);
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
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == Felt::NULL_IDX);
					CHECK(grid.get(pos3) == Felt::NULL_IDX);
					CHECK(grid.get(pos4) == 0);
				}
			}
		}

	}
}



SCENARIO("Lookup::Multi")
{
	using GridType = Impl::Lookup::Multi<3, 3>;

	GIVEN("a 10x10x10 EagerMultiLookupGrid with 3 tracking lists, and some locations")
	{
		GridType grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

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
			CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
			CHECK(grid.get(pos2)(0) == Felt::NULL_IDX);
			CHECK(grid.get(pos3)(0) == Felt::NULL_IDX);
			CHECK(grid.get(pos4)(0) == Felt::NULL_IDX);
			CHECK(grid.get(pos5)(0) == Felt::NULL_IDX);
			CHECK(grid.get(pos6)(0) == Felt::NULL_IDX);
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
				CHECK(grid.get(pos1)(0) == 0);
				CHECK(grid.get(pos2)(0) == 1);
				CHECK(grid.get(pos3)(0) == 2);
				CHECK(grid.get(pos4)(0) == 3);
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
					CHECK(grid.get(pos1)(0) == 0);
					CHECK(grid.get(pos2)(0) == 1);
					CHECK(grid.get(pos3)(0) == 2);
					CHECK(grid.get(pos4)(0) == 3);
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
					CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
					CHECK(grid.get(pos2)(0) == Felt::NULL_IDX);
					CHECK(grid.get(pos3)(0) == Felt::NULL_IDX);
					CHECK(grid.get(pos4)(0) == Felt::NULL_IDX);
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
				CHECK(grid.get(pos1)(0) == 0);
				CHECK(grid.get(pos2)(1) == 0);
				CHECK(grid.get(pos3)(1) == 1);
				CHECK(grid.get(pos4)(2) == 0);
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
					CHECK(grid.get(pos1)(0) == 0);
					CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
					CHECK(grid.get(pos3)(1) == 0);
					CHECK(grid.get(pos4)(2) == 0);
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
						CHECK(grid.get(pos1)(0) == 0);
						CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
						CHECK(grid.get(pos3)(1) == 0);
						CHECK(grid.get(pos4)(2) == 0);
						CHECK(grid.get(pos5)(2) == 1);
						CHECK(grid.get(pos6)(2) == 2);
					}

					AND_WHEN("we remove 2 points from different tracking lists")
					{
						grid.remove(pos4, 2);
						grid.remove(pos1, 0);

						THEN("the grid and tracking lists are updated")
						{
							CHECK(grid.list(0).size() == 0);
							CHECK(grid.list(1).size() == 1);
							CHECK(grid.list(2).size() == 2);
							CHECK(grid.list(0)[0] == pos1);
							CHECK(grid.list(1)[0] == pos3);
							CHECK(grid.list(2)[1] == pos5);
							CHECK(grid.list(2)[0] == pos6);
							CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
							CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
							CHECK(grid.get(pos3)(1) == 0);
							CHECK(grid.get(pos4)(2) == Felt::NULL_IDX);
							CHECK(grid.get(pos5)(2) == 1);
							CHECK(grid.get(pos6)(2) == 0);
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
								CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
								CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
								CHECK(grid.get(pos3)(1) == 0);
								CHECK(grid.get(pos4)(2) == Felt::NULL_IDX);
								CHECK(grid.get(pos5)(2) == Felt::NULL_IDX);
								CHECK(grid.get(pos6)(2) == Felt::NULL_IDX);
							}
						}
					}

				}
			}
		}
	}
}


SCENARIO("Lookup::LazySingle")
{
	using GridType = Impl::Lookup::LazySingle<3, 3>;

	GIVEN("a 3x3x3 lazy single-index lookup grid with 3 tracking lists")
	{
		/// [LazySingleLookupGrid initialisation]

		GridType grid(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

		THEN("the grid is initially inactive")
		{
			CHECK(grid.data().size() == 0);
			CHECK(grid.list(0).size() == 0);
			CHECK(grid.list(1).size() == 0);
			CHECK(grid.list(2).size() == 0);
		}
		/// [LazySingleLookupGrid initialisation]

		THEN("queries return the NULL background value")
		{
			CHECK(grid.get(Vec3i(1,1,1)) == Felt::NULL_IDX);
			CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::NULL_IDX);
		}

		/// [LazySingleLookupGrid activation]
		WHEN("the grid is activated")
		{
			grid.activate();

			THEN("memory is allocated and the grid filled with background value")
			{
				CHECK(grid.data().size() == 3 * 3 * 3);
			}

			THEN("queries still return the NULL background value")
			{
				CHECK(grid.get(Vec3i(1, 1, 1)) == Felt::NULL_IDX);
				CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::NULL_IDX);
			}

			AND_WHEN("we add a position to be tracked to list 1")
			{
				grid.add(Vec3i(1, 1, 1), 1);

				THEN("that position's value is updated and is added to the tracking list")
				{
					CHECK(grid.get(Vec3i(1, 1, 1)) == 0);
					CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::NULL_IDX);
					CHECK(grid.list(1)[0] == Vec3i(1, 1, 1));
				}

				AND_WHEN("the grid is deactivated")
				{
					grid.deactivate();

					THEN("the grid is once again inactive")
					{
						CHECK(grid.data().size() == 0);
						CHECK(grid.list(0).size() == 0);
						CHECK(grid.list(1).size() == 0);
						CHECK(grid.list(2).size() == 0);
					}

					THEN("queries once again return the NULL background value")
					{
						CHECK(grid.get(Vec3i(1,1,1)) == Felt::NULL_IDX);
						CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::NULL_IDX);
					}
				}
			}
			/// [LazySingleLookupGrid activation]
		}
	}
}


SCENARIO("Tracked::LazySingle")
{
	using GridType = Tracked::LazySingle<FLOAT, 3, 3>;
	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 3")
	{
		GridType grid(Vec3i(3, 3, 3), Vec3i(-1,-1,-1), 3.14159);

		const UINT NULL_IDX = Felt::NULL_IDX;

		THEN("the data grid and associated lookup grid state is inactive")
		{
			CHECK(grid.data().size() == 0);
			CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
			CHECK(grid.lookup().data().size() == 0);
			CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
		}

		/// [LazySingleTrackedGrid activate]
		WHEN("the grid is activated")
		{
			grid.activate();

			THEN("the data grid and associated lookup grid state is active")
			{
				CHECK(grid.data().size() == 3*3*3);
				CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
				CHECK(grid.lookup().data().size() == 3*3*3);
				CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
				/// [LazySingleTrackedGrid activate]
			}

			AND_WHEN("the grid is deactivated")
			{
				grid.deactivate();

				THEN("the data grid and associated lookup grid state is inactive")
				{
					CHECK(grid.data().size() == 0);
					CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
					CHECK(grid.lookup().data().size() == 0);
					CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
				}
			}
		}
	}
}

