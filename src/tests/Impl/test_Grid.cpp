#include <type_traits>
#include <experimental/type_traits>

#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Partitioned.hpp>
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
				CHECK(grid.get(Vec3i(0, 0, 0)) == 13.0f);
				CHECK(grid.get(0) == 13.0f);
				CHECK(grid.get(3 * 7 * 11 - 1) == 19.0f);
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
		const UINT pos1_idx = grid.index(pos1);
		const UINT pos2_idx = grid.index(pos2);
		const UINT pos3_idx = grid.index(pos3);
		const UINT pos4_idx = grid.index(pos4);
		const UINT pos5_idx = grid.index(pos5);
		const UINT pos6_idx = grid.index(pos6);
		const UINT pos7_idx = grid.index(pos7);

		WHEN("we track 4 locations to be tracked")
		{
			// Add the positions to the array and set index lookup values.
			grid.track(pos1_idx);
			grid.track(pos2_idx);
			grid.track(pos3_idx);
			grid.track(pos4_idx);


			THEN("the tracking lists contain the expected number of elements")
			{
				CHECK(grid.list().size() == 4);
			}

			THEN("the grid reports the active state of positions correctly")
			{
				CHECK(grid.is_tracked(pos1_idx) == true);
				CHECK(grid.is_tracked(pos2_idx) == true);
				CHECK(grid.is_tracked(pos3_idx) == true);
				CHECK(grid.is_tracked(pos4_idx) == true);
				CHECK(grid.is_tracked(pos5_idx) == false);
			}

			THEN("the tracking list elements contain the positions")
			{
				CHECK(grid.list()[0] == pos1_idx);
				CHECK(grid.list()[1] == pos2_idx);
				CHECK(grid.list()[2] == pos3_idx);
				CHECK(grid.list()[3] == pos4_idx);
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
				grid.remove(grid.index(pos7));

				THEN("the tracking lists contain the same number of elements")
				{
					CHECK(grid.list().size() == 4);
				}

				THEN("the tracking list elements contain the same position vectors")
				{
					CHECK(grid.list()[0] == pos1_idx);
					CHECK(grid.list()[1] == pos2_idx);
					CHECK(grid.list()[2] == pos3_idx);
					CHECK(grid.list()[3] == pos4_idx);
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
				grid.remove(grid.index(pos2));

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list().size() == 3);
				}

				THEN(
				"the removed position is gone from the list to be replaced by the final position"
				) {
					CHECK(grid.list()[0] == pos1_idx);
					CHECK(grid.list()[1] == pos4_idx);
					CHECK(grid.list()[2] == pos3_idx);
				}

				THEN("the grid location corresponding to the removed point gives NULL index")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == NULL_IDX);
					CHECK(grid.get(pos3) == 2);
					CHECK(grid.get(pos4) == 1);
				}

				AND_WHEN("we track two more points")
				{

					grid.track(pos5_idx);
					grid.track(pos6_idx);

					THEN("the tracking lists contain an extra of element")
					{
						CHECK(grid.list().size() == 5);
					}

					THEN("the tracking list elements contain the position vectors")
					{
						CHECK(grid.list()[0] == pos1_idx);
						CHECK(grid.list()[1] == pos4_idx);
						CHECK(grid.list()[2] == pos3_idx);
						CHECK(grid.list()[3] == pos5_idx);
						CHECK(grid.list()[4] == pos6_idx);
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
					CHECK(grid.is_tracked(pos1_idx) == false);
					CHECK(grid.is_tracked(pos2_idx) == false);
					CHECK(grid.is_tracked(pos3_idx) == false);
					CHECK(grid.is_tracked(pos4_idx) == false);
					CHECK(grid.is_tracked(pos5_idx) == false);
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
		const UINT pos1_idx = grid.index(pos1);
		const UINT pos2_idx = grid.index(pos2);
		const UINT pos3_idx = grid.index(pos3);
		const UINT pos4_idx = grid.index(pos4);
		const UINT pos5_idx = grid.index(pos5);
		const UINT pos6_idx = grid.index(pos6);
		const UINT pos7_idx = grid.index(pos7);


		WHEN("we track 4 locations to be tracked")
		{
			// Add the positions to the array and set index lookup values.
			grid.track(pos1_idx, 0);
			grid.track(pos2_idx, 1);
			grid.track(pos3_idx, 1);
			grid.track(pos4_idx, 2);


			THEN("the tracking lists contain the expected number of elements")
			{
				CHECK(grid.list(0).size() == 1);
				CHECK(grid.list(1).size() == 2);
				CHECK(grid.list(2).size() == 1);
			}

			THEN("the tracking list elements contain the positions")
			{
				CHECK(grid.list(0)[0] == pos1_idx);
				CHECK(grid.list(1)[0] == pos2_idx);
				CHECK(grid.list(1)[1] == pos3_idx);
				CHECK(grid.list(2)[0] == pos4_idx);
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
				grid.remove(pos7_idx, 1);

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 2);
					CHECK(grid.list(2).size() == 1);
				}

				THEN("the tracking list elements contain the position vectors")
				{
					CHECK(grid.list(0)[0] == pos1_idx);
					CHECK(grid.list(1)[0] == pos2_idx);
					CHECK(grid.list(1)[1] == pos3_idx);
					CHECK(grid.list(2)[0] == pos4_idx);
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
				grid.remove(pos2_idx, 1);

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
					CHECK(grid.list(0)[0] == pos1_idx);
					CHECK(grid.list(1)[0] == pos3_idx);
					CHECK(grid.list(2)[0] == pos4_idx);
				}

				THEN("the grid location corresponding to the removed point gives NULL index")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == Felt::NULL_IDX);
					CHECK(grid.get(pos3) == 0);
					CHECK(grid.get(pos4) == 0);
				}

				AND_WHEN("we track two more points")
				{

					grid.track(pos5_idx, 2);
					grid.track(pos6_idx, 2);

					THEN("the tracking lists contain the expected number of elements")
					{
						CHECK(grid.list(0).size() == 1);
						CHECK(grid.list(1).size() == 1);
						CHECK(grid.list(2).size() == 3);
					}

					THEN("the tracking list elements contain the position vectors")
					{
						CHECK(grid.list(0)[0] == pos1_idx);
						CHECK(grid.list(1)[0] == pos3_idx);
						CHECK(grid.list(2)[0] == pos4_idx);
						CHECK(grid.list(2)[1] == pos5_idx);
						CHECK(grid.list(2)[2] == pos6_idx);
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

			AND_WHEN("we reset the grid")
			{
				grid.reset();

				THEN("lists are empty ")
				{
					CHECK(grid.list(0).size() == 0);
					CHECK(grid.list(1).size() == 0);
					CHECK(grid.list(2).size() == 0);
				}

				THEN("the locations in the grid are now NULL index")
				{
					CHECK(grid.get(pos1) == Felt::NULL_IDX);
					CHECK(grid.get(pos2) == Felt::NULL_IDX);
					CHECK(grid.get(pos3) == Felt::NULL_IDX);
					CHECK(grid.get(pos4) == Felt::NULL_IDX);
				}
			}
		}

	}
}



SCENARIO("Lookup::Multi")
{
	using GridType = Impl::Lookup::Multi<3, 3>;

	GIVEN("a 10x10x10 grid with 3 tracking lists, and some locations")
	{
		GridType grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);
		const UINT pos1_idx = grid.index(pos1);
		const UINT pos2_idx = grid.index(pos2);
		const UINT pos3_idx = grid.index(pos3);
		const UINT pos4_idx = grid.index(pos4);
		const UINT pos5_idx = grid.index(pos5);
		const UINT pos6_idx = grid.index(pos6);

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
			grid.track(pos1_idx, 0);
			grid.track(pos2_idx, 0);
			grid.track(pos3_idx, 0);
			grid.track(pos4_idx, 0);

			THEN(
				"the tracking list is populated and the grid locations hold their indices in the"
				" list"
			) {
				REQUIRE(grid.list(0).size() == 4);
				CHECK(grid.list(0)[0] == pos1_idx);
				CHECK(grid.list(0)[1] == pos2_idx);
				CHECK(grid.list(0)[2] == pos3_idx);
				CHECK(grid.list(0)[3] == pos4_idx);
				CHECK(grid.get(pos1)(0) == 0);
				CHECK(grid.get(pos2)(0) == 1);
				CHECK(grid.get(pos3)(0) == 2);
				CHECK(grid.get(pos4)(0) == 3);
			}

			AND_WHEN("we track a location that is already tracked")
			{
				grid.track(pos2_idx, 0);

				THEN("the grid state is just as if the final point was not added")
				{
					REQUIRE(grid.list(0).size() == 4);
					CHECK(grid.list(0)[0] == pos1_idx);
					CHECK(grid.list(0)[1] == pos2_idx);
					CHECK(grid.list(0)[2] == pos3_idx);
					CHECK(grid.list(0)[3] == pos4_idx);
					CHECK(grid.get(pos1)(0) == 0);
					CHECK(grid.get(pos2)(0) == 1);
					CHECK(grid.get(pos3)(0) == 2);
					CHECK(grid.get(pos4)(0) == 3);
				}
			}

			AND_WHEN("we reset the grid")
			{
				grid.reset();

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
			grid.track(pos1_idx, 0);
			grid.track(pos2_idx, 1);
			grid.track(pos3_idx, 1);
			grid.track(pos4_idx, 2);
			grid.track(pos3_idx, 2);

			THEN("the tracking lists and index tuples within the grid are updated")
			{
				CHECK(grid.list(0).size() == 1);
				CHECK(grid.list(1).size() == 2);
				CHECK(grid.list(2).size() == 2);
				CHECK(grid.list(0)[0] == pos1_idx);
				CHECK(grid.list(1)[0] == pos2_idx);
				CHECK(grid.list(1)[1] == pos3_idx);
				CHECK(grid.list(2)[0] == pos4_idx);
				CHECK(grid.list(2)[1] == pos3_idx);
				CHECK(grid.get(pos1)(0) == 0);
				CHECK(grid.get(pos2)(1) == 0);
				CHECK(grid.get(pos3)(1) == 1);
				CHECK(grid.get(pos3)(2) == 1);
				CHECK(grid.get(pos4)(2) == 0);
			}

			AND_WHEN("we remove a location from tracking list 1")
			{
				grid.remove(pos2_idx, 1);

				THEN("the grid and tracking list is updated")
				{
					CHECK(grid.list(0).size() == 1);
					CHECK(grid.list(1).size() == 1);
					CHECK(grid.list(2).size() == 2);
					CHECK(grid.list(0)[0] == pos1_idx);
					CHECK(grid.list(1)[0] == pos3_idx);
					CHECK(grid.list(2)[0] == pos4_idx);
					CHECK(grid.list(2)[1] == pos3_idx);
					CHECK(grid.get(pos1)(0) == 0);
					CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
					CHECK(grid.get(pos3)(1) == 0);
					CHECK(grid.get(pos4)(2) == 0);
					CHECK(grid.get(pos3)(2) == 1);
				}

				AND_WHEN("we track 2 more points to tracking list 2")
				{
					grid.track(pos5_idx, 2);
					grid.track(pos6_idx, 2);

					THEN("the grid and tracking list is updated")
					{
						CHECK(grid.list(0).size() == 1);
						CHECK(grid.list(1).size() == 1);
						CHECK(grid.list(2).size() == 4);
						CHECK(grid.list(0)[0] == pos1_idx);
						CHECK(grid.list(1)[0] == pos3_idx);
						CHECK(grid.list(2)[0] == pos4_idx);
						CHECK(grid.list(2)[1] == pos3_idx);
						CHECK(grid.list(2)[2] == pos5_idx);
						CHECK(grid.list(2)[3] == pos6_idx);
						CHECK(grid.get(pos1)(0) == 0);
						CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
						CHECK(grid.get(pos3)(1) == 0);
						CHECK(grid.get(pos4)(2) == 0);
						CHECK(grid.get(pos3)(2) == 1);
						CHECK(grid.get(pos5)(2) == 2);
						CHECK(grid.get(pos6)(2) == 3);
					}

					AND_WHEN("we remove 2 points from different tracking lists")
					{
						grid.remove(pos4_idx, 2);
						grid.remove(pos1_idx, 0);

						THEN("the grid and tracking lists are updated")
						{
							CHECK(grid.list(0).size() == 0);
							CHECK(grid.list(1).size() == 1);
							CHECK(grid.list(2).size() == 3);
							CHECK(grid.list(0)[0] == pos1_idx);
							CHECK(grid.list(1)[0] == pos3_idx);
							CHECK(grid.list(2)[0] == pos6_idx);
							CHECK(grid.list(2)[1] == pos3_idx);
							CHECK(grid.list(2)[2] == pos5_idx);
							CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
							CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
							CHECK(grid.get(pos3)(1) == 0);
							CHECK(grid.get(pos3)(2) == 1);
							CHECK(grid.get(pos4)(2) == Felt::NULL_IDX);
							CHECK(grid.get(pos5)(2) == 2);
							CHECK(grid.get(pos6)(2) == 0);
						}

						AND_WHEN("the grid is reset")
						{
							grid.reset();

							THEN("the lists are reset and the grid contains null indices")
							{
								CHECK(grid.list(0).size() == 0);
								CHECK(grid.list(1).size() == 0);
								CHECK(grid.list(2).size() == 0);
								CHECK(grid.get(pos1)(0) == Felt::NULL_IDX);
								CHECK(grid.get(pos2)(1) == Felt::NULL_IDX);
								CHECK(grid.get(pos3)(1) == Felt::NULL_IDX);
								CHECK(grid.get(pos3)(2) == Felt::NULL_IDX);
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

		GridType grid = GridType();

		THEN("the grid is initially inactive")
		{
			CHECK(grid.data().size() == 0);
			CHECK(grid.list(0).size() == 0);
			CHECK(grid.list(1).size() == 0);
			CHECK(grid.list(2).size() == 0);
		}
		/// [LazySingleLookupGrid initialisation]

		WHEN("the grid size is set")
		{
			grid.resize(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

			THEN("the grid is still inactive but now reports new size")
			{
				CHECK(grid.data().size() == 0);
				CHECK(grid.list(0).size() == 0);
				CHECK(grid.list(1).size() == 0);
				CHECK(grid.list(2).size() == 0);
				CHECK(grid.size() == Vec3i(3, 3, 3));
				CHECK(grid.offset() == Vec3i(-1,-1,-1));
			}

			THEN("queries return the NULL background value")
			{
				CHECK(grid.get(Vec3i(1, 1, 1)) == Felt::NULL_IDX);
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

				AND_WHEN("we track a position to be tracked to list 1")
				{
					grid.track(grid.index(Vec3i(1, 1, 1)), 1);

					THEN("that position's value is updated and is added to the tracking list")
					{
						CHECK(grid.get(Vec3i(1, 1, 1)) == 0);
						CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::NULL_IDX);
						CHECK(grid.list(1)[0] == grid.index(Vec3i(1, 1, 1)));
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
}


/// Utility for static_assert of presence of reset method.
template <typename T> using has_reset_t = decltype(&T::reset);


SCENARIO("Tracked::LazySingle")
{
	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 3.14159")
	{
		using GridType = Impl::Tracked::LazySingleByValue<FLOAT, 3, 3>;

		static_assert(
			std::experimental::is_detected<has_reset_t, GridType>::value,
			"Tracked grids with a single lookup index per grid node should have a reset method."
		);

		GridType grid = GridType(3.14159);

		const PosIdx NULL_IDX = Felt::NULL_IDX;

		THEN("the data grid and associated lookup grid state is zero size and inactive")
		{
			CHECK(grid.data().size() == 0);
			CHECK(grid.lookup().data().size() == 0);
		}

		WHEN("the grid is resized")
		{
			grid.resize(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

			THEN("the data grid and associated lookup grid report new size and remains inactive")
			{
				CHECK(grid.data().size() == 0);
				CHECK(grid.lookup().data().size() == 0);

				CHECK(grid.size() == Vec3i(3,3,3));
				CHECK(grid.offset() == Vec3i(-1,-1,-1));
				CHECK(grid.lookup().size() == Vec3i(3,3,3));
				CHECK(grid.lookup().offset() == Vec3i(-1,-1,-1));

				CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
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

				AND_WHEN("a location is updated and tracked in list 1")
				{
					grid.track(42.0f, grid.index(Vec3i(1,1,1)), 1);

					THEN("the data grid is updated and the lookup grid tracks the point")
					{
						CHECK(grid.get(Vec3i(1,1,1)) == 42);
						CHECK(grid.lookup().get(Vec3i(1,1,1)) == 0);
						CHECK(grid.lookup().list(0).size() == 0);
						CHECK(grid.lookup().list(1).size() == 1);
						CHECK(grid.lookup().list(1)[0] == grid.index(Vec3i(1,1,1)));
						CHECK(grid.lookup().list(2).size() == 0);
					}

					AND_WHEN("the grid is reset")
					{
						grid.reset();

						THEN(
							"the value in the data grid is reset to background value and locations"
							" are no longer tracked"
						) {
							CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
							CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
							CHECK(grid.lookup().list(0).size() == 0);
							CHECK(grid.lookup().list(1).size() == 0);
							CHECK(grid.lookup().list(2).size() == 0);
						}
					}
				}
			}
		}
	}

	GIVEN("a 9x9x9 grid of std::vectors with (-4,-4,-4) offset")
	{
		using LeafType = std::vector<int>;
		using GridType = Impl::Tracked::LazySingleByValue<LeafType, 3, 3>;

		GridType grid(LeafType{1,2,3});
		grid.resize(Vec3i(9,9,9), Vec3i(-4,-4,-4));
		grid.activate();

		WHEN("a value is set with an lvalue reference")
		{
			LeafType move_me{5,6,7};
			int* pdata = &move_me[0];
			grid.set(Vec3i(2,2,2), move_me);

			THEN("the data from the input has been copied into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == move_me);
				CHECK(&grid.get(Vec3i(2,2,2))[0] != pdata);
			}
		}

		WHEN("a value is set with an rvalue reference")
		{
			LeafType move_me{5,6,7};
			LeafType copied = move_me;
			int* pdata = &move_me[0];
			grid.set(Vec3i(2,2,2), std::move(move_me));

			THEN("the data from the input has been copied into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == copied);
				CHECK(&grid.get(Vec3i(2,2,2))[0] != pdata);
			}
		}
	}
}


SCENARIO("Tracked::MultiByRef")
{
	GIVEN("a 9x9x9 grid of floats with (-4,-4,-4) offset and background value of 0")
	{
		using GridType = Impl::Tracked::MultiByRef<FLOAT, 3, 3>;

		static_assert(
			!std::experimental::is_detected<has_reset_t, GridType>::value,
			"Tracked grids with multiple lookup indices per grid node should not have a reset"
			" method."
		);

		GridType grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), 0);

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
				CHECK(val == Vec3u::Constant(Felt::NULL_IDX));
		}

		WHEN("a simple value is added to the grid to be tracked by list 1 and 2")
		{
			grid.track(42.0f, grid.index(Vec3i(2,2,2)), 1);
			grid.lookup().track(grid.index(Vec3i(2,2,2)), 2);

			THEN("the value stored in the grid is correct")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == 42.0f);
			}

			THEN("the lookup grid is tracking the location just added")
			{
				CHECK(grid.lookup().get(Vec3i(2,2,2)) == Vec3u(Felt::NULL_IDX, 0, 0));
				CHECK(grid.lookup().list(0).size() == 0);
				CHECK(grid.lookup().list(1).size() == 1);
				CHECK(grid.lookup().list(2).size() == 1);
				CHECK(grid.lookup().list(1)[0] == grid.index(Vec3i(2,2,2)));
				CHECK(grid.lookup().list(2)[0] == grid.index(Vec3i(2,2,2)));
			}

			AND_WHEN("the value is modified by reference")
			{
				grid.get(Vec3i(2,2,2)) = 3;

				THEN("the value in the grid is updated")
				{
					CHECK(grid.get(Vec3i(2,2,2)) == 3.0f);
				}
			}
		}
	}

	GIVEN("a 9x9x9 grid of std::vectors with (-4,-4,-4) offset")
	{
		using LeafType = std::vector<int>;
		using GridType = Impl::Tracked::MultiByRef<LeafType, 3, 3>;
		GridType grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), LeafType{1,2,3});

		WHEN("a value is set with an lvalue reference")
		{
			LeafType move_me{5,6,7};
			int* pdata = &move_me[0];
			grid.get(Vec3i(2,2,2)) = move_me;

			THEN("the data from the input has been copied into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == move_me);
				CHECK(&grid.get(Vec3i(2,2,2))[0] != pdata);
			}
		}

		WHEN("a value is set with an rvalue reference")
		{
			LeafType move_me{5,6,7};
			LeafType copied = move_me;
			int* pdata = &move_me[0];
			grid.get(Vec3i(2,2,2)) = std::move(move_me);

			THEN("the data from the input has been moved into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == copied);
				CHECK(&grid.get(Vec3i(2,2,2))[0] == pdata);
			}
		}
	}
}


SCENARIO("Paritioned::Lookup")
{
	using GridType = Impl::Partitioned::Lookup<3, 3>;
	using ChildrenGrid = GridType::ChildrenGrid;

	static_assert(
		std::is_same<
			ChildrenGrid, Impl::Tracked::MultiByRef<Impl::Lookup::LazySingle<3, 3>, 3, 3>
		>::value,
		"Children grid of partitioned lookup must be a multi-list multi-index tracked grid with"
		" lazily activated lookup sub-grids as the leaf type."
	);

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions")
	{
		GridType grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3));

		THEN("the children tracking grid has been initialised with lazy subgrids")
		{
			CHECK(grid.children().data().size() == 3*3*3);
			CHECK(grid.children().lookup().data().size() == 3*3*3);

			CHECK(grid.children().get(Vec3i(-1,-1,-1)).data().size() == 0);
			CHECK(grid.children().get(Vec3i(-1,-1,-1)).size() == Vec3i(3,3,3));
			CHECK(grid.children().get(Vec3i(-1,-1,-1)).offset() == Vec3i(-4,-4,-4));
			CHECK(grid.children().get(Vec3i( 1, 1, 1)).size() == Vec3i(3,3,3));
			CHECK(grid.children().get(Vec3i( 1, 1, 1)).offset() == Vec3i( 2, 2, 2));
		}

		WHEN("some locations are tracked")
		{
			const Vec3i pos1(1, -4, -1);
			const Vec3i pos2(2, -3, -2);
			const Vec3i pos3(3, -2, -3);
			const Vec3i pos4(4, -1, -4);
			const Vec3i part1(0, -1,  0);
			const Vec3i part2_3(1, -1, -1);
			const Vec3i part4(1,  0, -1);
			const PosIdx part1_idx = grid.children().index(part1);
			const PosIdx part2_3_idx = grid.children().index(part2_3);
			const PosIdx part4_idx = grid.children().index(part4);

			grid.track(pos1, 0);
			grid.track(pos2, 0);
			// Will not activate partition, but already activated by previous line.
			grid.track(part2_3_idx, grid.children().get(part2_3_idx).index(pos3), 0);
			grid.track(pos4, 2);

			THEN("the underlying data is as expected")
			{
				CHECK(grid.children().lookup().list(0).size() == 2);
				CHECK(grid.children().lookup().list(2).size() == 1);
				CHECK(grid.children().lookup().list(0)[0] == part1_idx);
				CHECK(grid.children().lookup().list(0)[1] == part2_3_idx);
				CHECK(grid.children().lookup().list(2)[0] == part4_idx);
				CHECK(grid.children().lookup().get(part1)(0) == 0);
				CHECK(grid.children().lookup().get(part2_3)(0) == 1);
				CHECK(grid.children().lookup().get(part4)(2) == 0);

				CHECK(grid.children().get(part1).list(0).size() == 1);
				CHECK(grid.children().get(part2_3).list(0).size() == 2);
				CHECK(grid.children().get(part4).list(2).size() == 1);
				CHECK(grid.children().get(part1).get(pos1) == 0);
				CHECK(grid.children().get(part2_3).get(pos2) == 0);
				CHECK(grid.children().get(part2_3).get(pos3) == 1);
				CHECK(grid.children().get(part4).get(pos4) == 0);
			}
		}

		WHEN("some points are tracked which overlap points tracked in a masking grid")
		{
			GridType grid_master(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3));

			const Vec3i pos_list_0(0, 0, 0);
			const Vec3i pos_active_because_master(-4, 0, 4);
			const Vec3i pos_list_1(4, 0, 0);
			const Vec3i pos_child_list_0(0, 0, 0);
			const Vec3i pos_child_active_because_master(-1, 0, 1);
			const Vec3i pos_child_list_1(1, 0, 0);

			grid_master.track(pos_active_because_master, 0);
			grid.track(pos_active_because_master, 0);
			grid.track(pos_list_0, 0);
			grid.track(pos_list_1, 1);

			AND_WHEN("resetting the grid")
			{
				grid.reset(grid_master);

				THEN("all children are reset")
				{
					CHECK(
						grid.children().get(
							pos_child_active_because_master
						).get(pos_active_because_master) == NULL_IDX
					);
					CHECK(grid.children().get(pos_child_list_0).get(pos_list_0) == NULL_IDX);
					CHECK(grid.children().get(pos_child_list_1).get(pos_list_1) == NULL_IDX);
					CHECK(
						grid.children().get(pos_child_active_because_master).list(0).size() == 0
					);
					CHECK(grid.children().get(pos_child_list_0).list(0).size() == 0);
					CHECK(grid.children().get(pos_child_list_1).list(1).size() == 0);
				}

				THEN("all children not tracked by mask grid are deactivated")
				{
					CHECK(grid.children().get(Vec3i(pos_child_list_0)).is_active() == false);
					CHECK(grid.children().get(Vec3i(pos_child_list_0)).data().size() == 0);
					CHECK(grid.children().get(pos_child_list_1).is_active() == false);
					CHECK(grid.children().get(pos_child_list_1).data().size() == 0);
					CHECK(grid.children().lookup().list(1).size() == 0);
				}

				THEN("children that are tracked by mask grid remain active")
				{
					CHECK(grid.children().get(pos_child_active_because_master).is_active() == true);
					CHECK(
						grid.children().get(pos_child_active_because_master).data().size() == 3*3*3
					);
					CHECK(grid.children().lookup().list(0).size() == 0);
				}
			}
		}
	}
}


SCENARIO("Paritioned::Tracked::Simple")
{
	using GridType = Impl::Partitioned::Tracked::Simple<INT, 3, 3>;

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions with background value -42")
	{
		GridType grid(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), -42);

		const Vec3i pos12_child(-1, -1, -1);
		const Vec3i pos3_child(0, 0, 0);
		const Vec3i pos4_child(1, 1, 1);
		const PosIdx pos12_child_idx = grid.children().index(pos12_child);
		const PosIdx pos3_child_idx = grid.children().index(pos3_child);
		const PosIdx pos4_child_idx = grid.children().index(pos4_child);

		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos3(0, 0, 0);
		const Vec3i pos4(4, 4, 4);
		const PosIdx pos1_idx = grid.children().get(pos12_child_idx).index(pos1);
		const PosIdx pos2_idx = grid.children().get(pos12_child_idx).index(pos2);
		const PosIdx pos3_idx = grid.children().get(pos3_child_idx).index(pos3);
		const PosIdx pos4_idx = grid.children().get(pos4_child_idx).index(pos4);

		THEN("grid is initialised as inactive with no tracking and queries give background value")
		{
			CHECK(grid.children().get(pos12_child_idx).is_active() == false);
			CHECK(grid.children().get(pos12_child_idx).data().size() == 0);
			CHECK(grid.children().get(pos12_child_idx).get(pos1_idx) == -42);
			CHECK(grid.children().get(pos12_child_idx).lookup().data().size() == 0);
			CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == NULL_IDX);
			CHECK(grid.children().get(pos12_child_idx).get(pos2_idx) == -42);
			CHECK(grid.children().get(pos12_child_idx).lookup().get(pos2_idx) == NULL_IDX);

			CHECK(grid.children().get(pos3_child_idx).is_active() == false);
			CHECK(grid.children().get(pos3_child_idx).data().size() == 0);
			CHECK(grid.children().get(pos3_child_idx).get(pos3_idx) == -42);
			CHECK(grid.children().get(pos3_child_idx).lookup().data().size() == 0);
			CHECK(grid.children().get(pos3_child_idx).lookup().get(pos3_idx) == NULL_IDX);

			CHECK(grid.children().get(pos4_child_idx).is_active() == false);
			CHECK(grid.children().get(pos4_child_idx).data().size() == 0);
			CHECK(grid.children().get(pos4_child_idx).get(pos4_idx) == -42);
			CHECK(grid.children().get(pos4_child_idx).lookup().data().size() == 0);
			CHECK(grid.children().get(pos4_child_idx).lookup().get(pos4_idx) == NULL_IDX);
		}

		AND_WHEN("a mask grid is tracking some partitions")
		{
			GridType grid_master(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), 3.14159);
			grid_master.track(1234, pos1, 0);
			grid_master.track(1234, pos3, 0);

			AND_WHEN("children are added based on the mask grid")
			{
				grid.track_children(grid_master);

				THEN("those children are now active and initialised to background value")
				{
					CHECK(grid.children().get(pos12_child_idx).is_active() == true);
					CHECK(grid.children().get(pos12_child_idx).data().size() == 3*3*3);
					CHECK(grid.children().get(pos12_child_idx).get(pos1_idx) == -42);
					CHECK(grid.children().get(pos12_child_idx).lookup().data().size() == 3*3*3);
					CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == NULL_IDX);
					CHECK(grid.children().get(pos12_child_idx).get(pos2_idx) == -42);
					CHECK(grid.children().get(pos12_child_idx).lookup().get(pos2_idx) == NULL_IDX);

					CHECK(grid.children().get(pos3_child_idx).is_active() == true);
					CHECK(grid.children().get(pos3_child_idx).data().size() == 3*3*3);
					CHECK(grid.children().get(pos3_child_idx).get(pos3_idx) == -42);
					CHECK(grid.children().get(pos3_child_idx).lookup().data().size() == 3*3*3);
					CHECK(grid.children().get(pos3_child_idx).lookup().get(pos3_idx) == NULL_IDX);

					CHECK(grid.children().lookup().list(0).size() == 2);
					CHECK(grid.children().lookup().list(1).size() == 0);
					CHECK(grid.children().lookup().list(2).size() == 0);
				}

				THEN("other children are remain inactive")
				{
					CHECK(grid.children().get(pos4_child_idx).is_active() == false);
					CHECK(grid.children().get(pos4_child_idx).data().size() == 0);
					CHECK(grid.children().get(pos4_child_idx).get(pos4_idx) == -42);
					CHECK(grid.children().get(pos4_child_idx).lookup().data().size() == 0);
					CHECK(grid.children().get(pos4_child_idx).lookup().get(pos4_idx) == NULL_IDX);
				}

				AND_WHEN("a position in an active partition is tracked by index for list 1")
				{
					grid.track(345, pos12_child_idx, pos1_idx, 1);

					THEN("the grid value is updated")
					{
						CHECK(grid.children().get(pos12_child_idx).get(pos1_idx) == 345);
					}

					THEN("list 1 is *not* tracked by the children grid")
					{
						CHECK(grid.children().lookup().list(1).size() == 0);
					}

					THEN("list 1 tracks the leaf position in the child's lookup grid")
					{
						CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == 0);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(1).size() == 1);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(1)[0] == pos1_idx);

						CHECK(grid.children().get(pos12_child_idx).lookup().list(0).size() == 0);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(2).size() == 0);

					}
				}

				AND_WHEN("a position in an active partition is tracked by location for list 1")
				{
					grid.track(345, pos1, 1);

					THEN("the grid value is updated")
					{
						CHECK(grid.children().get(pos12_child_idx).get(pos1_idx) == 345);
					}

					THEN("list 1 is tracked by the children grid")
					{
						CHECK(grid.children().lookup().list(1).size() == 1);
						CHECK(grid.children().lookup().get(pos12_child_idx)[1] == 0);
					}

					THEN("list 1 tracks the leaf position in the child's lookup grid")
					{
						CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == 0);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(1).size() == 1);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(1)[0] == pos1_idx);
					}
				}

				AND_WHEN("a position in an inactive partition is tracked by location")
				{
					grid.track(345, pos4, 0);

					THEN("the partition is activated")
					{
						CHECK(grid.children().get(pos4_child_idx).is_active() == true);
						CHECK(grid.children().get(pos4_child_idx).data().size() == 3*3*3);
					}

					THEN("the grid value is updated")
					{
						CHECK(grid.children().get(pos4_child_idx).get(pos4_idx) == 345);
					}

					THEN("child is tracked by the children grid")
					{
						CHECK(grid.children().lookup().list(0).size() == 3);
						CHECK(grid.children().lookup().list(0)[2] == pos4_child_idx);
						CHECK(grid.children().lookup().get(pos4_child_idx)[0] == 2);
					}

					THEN("child tracks the leaf position in the child's lookup grid")
					{
						CHECK(grid.children().get(pos4_child_idx).lookup().get(pos4_idx) == 0);
						CHECK(grid.children().get(pos4_child_idx).lookup().list(0).size() == 1);
						CHECK(grid.children().get(pos4_child_idx).lookup().list(0)[0] == pos4_idx);
					}

					AND_WHEN("the grid is reset")
					{
						grid.reset(grid_master);

						THEN("partitions tracked in master remain active")
						{
							CHECK(grid.children().get(pos12_child_idx).is_active() == true);
							CHECK(grid.children().get(pos12_child_idx).data().size() == 3*3*3);
							CHECK(grid.children().get(pos3_child_idx).is_active() == true);
							CHECK(grid.children().get(pos3_child_idx).data().size() == 3*3*3);
						}

						THEN("partitions not tracked in master are deactivated")
						{
							CHECK(grid.children().get(pos4_child_idx).is_active() == false);
							CHECK(grid.children().get(pos4_child_idx).data().size() == 0);
						}

						THEN("children grid lookup is reset")
						{
							auto NULL_IDX_TUPLE = grid.children().lookup().NULL_IDX_TUPLE;
							CHECK(grid.children().lookup().list(0).size() == 0);
							CHECK(grid.children().lookup().list(1).size() == 0);
							CHECK(grid.children().lookup().list(2).size() == 0);
							CHECK(grid.children().lookup().get(pos12_child_idx) == NULL_IDX_TUPLE);
							CHECK(grid.children().lookup().get(pos3_child_idx) == NULL_IDX_TUPLE);
							CHECK(grid.children().lookup().get(pos4_child_idx) == NULL_IDX_TUPLE);
						}

						THEN("child grids' lookups are reset")
						{
							CHECK(grid.children().get(pos12_child_idx).lookup().list(0).size() == 0);
							CHECK(grid.children().get(pos12_child_idx).lookup().list(1).size() == 0);
							CHECK(grid.children().get(pos12_child_idx).lookup().list(2).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(0).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(1).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(2).size() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(0).size() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(1).size() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(2).size() == 0);
						}
					}
				}
			}
		}
	}
}


SCENARIO("Paritioned::Tracked::Numeric")
{
	using GridType = Impl::Partitioned::Tracked::Numeric<INT, 3, 3>;

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions with background value -42")
	{
		GridType grid(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), -42);

		const Vec3i pos123_child(-1, -1, -1);
		const Vec3i pos3_child(-1, -1, -1);
		const Vec3i pos4_child(1, 1, 1);
		const PosIdx pos123_child_idx = grid.children().index(pos123_child);
		const PosIdx pos3_child_idx = grid.children().index(pos3_child);
		const PosIdx pos4_child_idx = grid.children().index(pos4_child);

		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos3(-4, -3, -4);
		const Vec3i pos4(4, 4, 4);
		const PosIdx pos1_idx = grid.children().get(pos123_child_idx).index(pos1);
		const PosIdx pos2_idx = grid.children().get(pos123_child_idx).index(pos2);
		const PosIdx pos3_idx = grid.children().get(pos3_child_idx).index(pos3);
		const PosIdx pos4_idx = grid.children().get(pos4_child_idx).index(pos4);

		const GridType::ChildType& child = grid.children().get(pos123_child_idx);

		WHEN("two positions in the same partition are tracked")
		{
			grid.track(345, pos1, 0);
			grid.track(789, pos2, 0);
			grid.track(123, pos3, 1);

			THEN("the grid values are updated")
			{
				CHECK(child.get(pos1_idx) == 345);
				CHECK(child.get(pos2_idx) == 789);
				CHECK(child.get(pos3_idx) == 123);
			}

			THEN("partition is tracked by the children grid")
			{
				CHECK(grid.children().lookup().list(0).size() == 1);
				CHECK(grid.children().lookup().list(1).size() == 1);
				CHECK(grid.children().lookup().list(2).size() == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[0] == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[2] == NULL_IDX);
			}

			THEN("child tracks the leaf positions")
			{
				CHECK(child.lookup().list(0).size() == 2);
				CHECK(child.lookup().list(0)[0] == pos1_idx);
				CHECK(child.lookup().list(0)[1] == pos2_idx);
				CHECK(child.lookup().list(1).size() == 1);
				CHECK(child.lookup().list(1)[0] == pos3_idx);
				CHECK(child.lookup().get(pos1_idx) == 0);
				CHECK(child.lookup().get(pos2_idx) == 1);
				CHECK(child.lookup().get(pos3_idx) == 0);
			}

			AND_WHEN("a position is removed")
			{
				grid.remove(-100, pos123_child_idx, pos1_idx, 0);

				THEN("the partition is still tracked by the children grid")
				{
					CHECK(grid.children().lookup().list(0).size() == 1);
					CHECK(grid.children().lookup().list(1).size() == 1);
					CHECK(grid.children().lookup().list(2).size() == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[0] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[2] == NULL_IDX);

				}

				THEN("the child now tracks just the other positions")
				{
					CHECK(child.lookup().list(0).size() == 1);
					CHECK(child.lookup().list(1).size() == 1);
					CHECK(child.lookup().list(2).size() == 0);
					CHECK(child.lookup().list(0)[0] == pos2_idx);
					CHECK(child.lookup().list(1)[0] == pos3_idx);
					CHECK(child.lookup().get(pos1_idx) == NULL_IDX);
					CHECK(child.lookup().get(pos2_idx) == 0);
					CHECK(child.lookup().get(pos3_idx) == 0);
				}

				THEN("the grid value of the removed position is set to the passed background value")
				{
					CHECK(child.get(pos1_idx) == -100);
					CHECK(child.get(pos2_idx) == 789);
					CHECK(child.get(pos3_idx) == 123);
				}

				AND_WHEN("another postion is removed from the same tracking list")
				{
					grid.remove(-102, pos123_child_idx, pos2_idx, 0);

					THEN("the partition is still tracked by the children grid")
					{
						CHECK(grid.children().lookup().list(0).size() == 0);
						CHECK(grid.children().lookup().list(1).size() == 1);
						CHECK(grid.children().lookup().list(2).size() == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[0] == NULL_IDX);
						CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[2] == NULL_IDX);
					}

					THEN("the child now tracks just the final position")
					{
						CHECK(child.lookup().list(0).size() == 0);
						CHECK(child.lookup().list(1).size() == 1);
						CHECK(child.lookup().list(2).size() == 0);
						CHECK(child.lookup().list(1)[0] == pos3_idx);
						CHECK(child.lookup().get(pos1_idx) == NULL_IDX);
						CHECK(child.lookup().get(pos2_idx) == NULL_IDX);
						CHECK(child.lookup().get(pos3_idx) == 0);
					}

					AND_WHEN("the final postion is removed from tracking")
					{
						grid.remove(-999, pos123_child_idx, pos3_idx, 1);

						THEN("the partition is no longer tracked by the children grid")
						{
							CHECK(grid.children().lookup().list(0).size() == 0);
							CHECK(grid.children().lookup().list(1).size() == 0);
							CHECK(grid.children().lookup().list(2).size() == 0);
							CHECK(grid.children().lookup().get(pos123_child_idx)[0] == NULL_IDX);
							CHECK(grid.children().lookup().get(pos123_child_idx)[1] == NULL_IDX);
							CHECK(grid.children().lookup().get(pos123_child_idx)[2] == NULL_IDX);
						}

						THEN("the child no longer tracks any final positions")
						{
							CHECK(child.lookup().list(0).size() == 0);
							CHECK(child.lookup().list(1).size() == 0);
							CHECK(child.lookup().list(2).size() == 0);
							CHECK(child.lookup().get(pos1_idx) == NULL_IDX);
							CHECK(child.lookup().get(pos2_idx) == NULL_IDX);
							CHECK(child.lookup().get(pos3_idx) == NULL_IDX);
						}

						THEN("the child has been deactivated")
						{
							CHECK(child.is_active() == false);
							CHECK(child.data().size() == 0);
						}

						THEN("all positions report the final background value")
						{
							CHECK(child.get(pos1_idx) == -999);
							CHECK(child.get(pos2_idx) == -999);
							CHECK(child.get(pos3_idx) == -999);
						}
					}
				}
			}
		}
	}
}

