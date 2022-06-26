#include <type_traits>
#include <experimental/type_traits>
#include <fstream>

#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Partitioned.hpp>
#include <Felt/Impl/Tracked.hpp>

#include "catch.hpp"
#include "Utils.hpp"

using namespace Felt;


/**
 * Test the Grid class.
 */
SCENARIO("Grid::Simple")
{
	using Grid = Impl::Grid::Simple<Distance, 3>;
	/// [Grid - basics: GIVEN 3x7x11]
	GIVEN("a 3x7x11 grid with no offset and background value of 0")
	{
		Grid grid(Vec3i(3, 7, 11), Vec3i::Zero(), 0);
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
			CHECK(grid.inside(Vec3f(0,-0.00001f,0)) == false);
			CHECK(grid.inside(Vec3f(0,0,9.99999f)) == true);
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
		Grid grid(size, offset, 0);

		//! [Position index]
		THEN("the index of a point in the data array is reported correctly")
		{
			CHECK(Felt::index<3>(Vec3i(1, 0, -1), size, offset) == 613);
			CHECK(grid.index(Vec3i(1, 0, -1)) == 613);
		}

		THEN("the point represented by an index in the data array is reported correctly")
		{
			CHECK(grid.index(613) == Vec3i(1, 0, -1));
			CHECK(Felt::index<3>(613, size, offset) == Vec3i(1, 0, -1));
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


SCENARIO("Grid::Snapshot")
{
	using Grid = Impl::Grid::Snapshot<Distance, 3>;
	GIVEN("a 7x11x13 grid with (-3,-3,-3) offset and background value of 2")
	{
		Grid grid(Vec3i(7, 11, 13), Vec3i(-3, -3, -3), 2);

		static_assert(
			std::is_same<
				Grid::VArrayData, Eigen::Map< Eigen::Array<Distance, 1, Eigen::Dynamic> >
			>::value,
			"Vector form of underlying data should be an Eigen::Map to an Eigen::Array type"
		);

		WHEN("the underlying data is wrapped in an Eigen vector")
		{
			Grid::VArrayData array = grid.array();

			THEN("the data contains the expected (background) values")
			{
				CHECK(array(123) == 2);
			}

			AND_WHEN("values in the wrapped data are modified")
			{
				array(123) = 7;

				THEN("the underlying data in the grid is also modified")
				{
					CHECK(grid.get(123) == 7);
				}
			}
		}
	}
}


SCENARIO("Lookup::SingleListSingleIdx")
{
	using Grid = Impl::Lookup::SingleListSingleIdx<3>;
	GIVEN("a 10x10x10 grid with (0,-5,-5) offset")
	{
		Grid grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);
		const Vec3i pos7(7, 0, 0);
		const PosIdx pos1_idx = grid.index(pos1);
		const PosIdx pos2_idx = grid.index(pos2);
		const PosIdx pos3_idx = grid.index(pos3);
		const PosIdx pos4_idx = grid.index(pos4);
		const PosIdx pos5_idx = grid.index(pos5);
		const PosIdx pos6_idx = grid.index(pos6);
		const PosIdx pos7_idx = grid.index(pos7);

		WHEN("we track 4 locations")
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

			THEN("the grid reports the tracked state of positions correctly")
			{
				CHECK(grid.is_tracked(pos1_idx) == true);
				CHECK(grid.is_tracked(pos2_idx) == true);
				CHECK(grid.is_tracked(pos3_idx) == true);
				CHECK(grid.is_tracked(pos4_idx) == true);
				CHECK(grid.is_tracked(pos5_idx) == false);
			}

			THEN("the tracking list elements contain the position indices")
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
			AND_WHEN("we untrack a position that is not tracked")
			{
				grid.untrack(grid.index(pos7));

				THEN("the tracking lists contain the same number of elements")
				{
					CHECK(grid.list().size() == 4);
				}

				THEN("the tracking list elements contain the same position indices")
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

			AND_WHEN("we untrack a tracked position")
			{
				grid.untrack(grid.index(pos2));

				THEN("the tracking lists contain the expected number of elements")
				{
					CHECK(grid.list().size() == 3);
				}

				THEN(
					"the untracked position is gone from the list to be replaced by the final"
					" position"
				) {
					CHECK(grid.list()[0] == pos1_idx);
					CHECK(grid.list()[1] == pos4_idx);
					CHECK(grid.list()[2] == pos3_idx);
				}

				THEN("the grid location corresponding to the untracked point gives NULL index")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == null_idx);
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
						CHECK(grid.get(pos2) == null_idx);
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
					CHECK(grid.get(pos1) == null_idx);
					CHECK(grid.get(pos2) == null_idx);
					CHECK(grid.get(pos3) == null_idx);
					CHECK(grid.get(pos4) == null_idx);
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



SCENARIO("Lookup::MultiListSingleIdx")
{
	using Grid = Impl::Lookup::MultiListSingleIdx<3, 3>;

	GIVEN("a grid and some locations")
	{
		Grid grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);
		const Vec3i pos7(7, 0, 0);
		const PosIdx pos1_idx = grid.index(pos1);
		const PosIdx pos2_idx = grid.index(pos2);
		const PosIdx pos3_idx = grid.index(pos3);
		const PosIdx pos4_idx = grid.index(pos4);
		const PosIdx pos5_idx = grid.index(pos5);
		const PosIdx pos6_idx = grid.index(pos6);
		const PosIdx pos7_idx = grid.index(pos7);


		WHEN("we add 4 locations to be tracked")
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
			AND_WHEN("we untrack a position that is not tracked")
			{
				grid.untrack(pos7_idx, 1);

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

			AND_WHEN("we untrack a position from tracking in list 0")
			{
				grid.untrack(pos2_idx, 1);

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

				THEN("the grid location corresponding to the untrackd point gives NULL index")
				{
					CHECK(grid.get(pos1) == 0);
					CHECK(grid.get(pos2) == Felt::null_idx);
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
						CHECK(grid.get(pos2) == Felt::null_idx);
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
					CHECK(grid.get(pos1) == Felt::null_idx);
					CHECK(grid.get(pos2) == Felt::null_idx);
					CHECK(grid.get(pos3) == Felt::null_idx);
					CHECK(grid.get(pos4) == Felt::null_idx);
				}
			}
		}

	}
}



SCENARIO("Lookup::Multi")
{

	GIVEN("a 2D grid with 5 tracking lists")
	{
		using Grid = Impl::Lookup::MultiListMultiIdx<2, 5>;
		Grid grid(Vec2i(3,3), Vec2i(-1,-1));

		THEN("the grid has the correct dimension")
		{
			const Vec2i pos1(-1, -1);
			const PosIdx pos1_idx = grid.index(pos1);
			CHECK(pos1_idx == 0);

			AND_THEN("the grid holds index tuples of the correct size")
			{
				CHECK(grid.get(pos1_idx).size() == 5);
			}
		}
	}

	using Grid = Impl::Lookup::MultiListMultiIdx<3, 3>;

	GIVEN("a 10x10x10 grid with 3 tracking lists, and some locations")
	{
		Grid grid(Vec3i(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);
		const PosIdx pos1_idx = grid.index(pos1);
		const PosIdx pos2_idx = grid.index(pos2);
		const PosIdx pos3_idx = grid.index(pos3);
		const PosIdx pos4_idx = grid.index(pos4);
		const PosIdx pos5_idx = grid.index(pos5);
		const PosIdx pos6_idx = grid.index(pos6);

		THEN("the grid is initialised with NULL indices and the tracking lists are empty")
		{
			REQUIRE(grid.list(0).capacity() == 0);
			REQUIRE(grid.list(1).capacity() == 0);
			REQUIRE(grid.list(2).capacity() == 0);
			CHECK(grid.get(pos1)(0) == Felt::null_idx);
			CHECK(grid.get(pos2)(0) == Felt::null_idx);
			CHECK(grid.get(pos3)(0) == Felt::null_idx);
			CHECK(grid.get(pos4)(0) == Felt::null_idx);
			CHECK(grid.get(pos5)(0) == Felt::null_idx);
			CHECK(grid.get(pos6)(0) == Felt::null_idx);
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
					CHECK(grid.get(pos1)(0) == Felt::null_idx);
					CHECK(grid.get(pos2)(0) == Felt::null_idx);
					CHECK(grid.get(pos3)(0) == Felt::null_idx);
					CHECK(grid.get(pos4)(0) == Felt::null_idx);
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

			AND_WHEN("we untrack a location from tracking list 1")
			{
				grid.untrack(pos2_idx, 1);

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
					CHECK(grid.get(pos2)(1) == Felt::null_idx);
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
						CHECK(grid.get(pos2)(1) == Felt::null_idx);
						CHECK(grid.get(pos3)(1) == 0);
						CHECK(grid.get(pos4)(2) == 0);
						CHECK(grid.get(pos3)(2) == 1);
						CHECK(grid.get(pos5)(2) == 2);
						CHECK(grid.get(pos6)(2) == 3);
					}

					AND_WHEN("we untrack 2 points from different tracking lists")
					{
						grid.untrack(pos4_idx, 2);
						grid.untrack(pos1_idx, 0);

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
							CHECK(grid.get(pos1)(0) == Felt::null_idx);
							CHECK(grid.get(pos2)(1) == Felt::null_idx);
							CHECK(grid.get(pos3)(1) == 0);
							CHECK(grid.get(pos3)(2) == 1);
							CHECK(grid.get(pos4)(2) == Felt::null_idx);
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
								CHECK(grid.get(pos1)(0) == Felt::null_idx);
								CHECK(grid.get(pos2)(1) == Felt::null_idx);
								CHECK(grid.get(pos3)(1) == Felt::null_idx);
								CHECK(grid.get(pos3)(2) == Felt::null_idx);
								CHECK(grid.get(pos4)(2) == Felt::null_idx);
								CHECK(grid.get(pos5)(2) == Felt::null_idx);
								CHECK(grid.get(pos6)(2) == Felt::null_idx);
							}
						}
					}

				}
			}
		}
	}
}


SCENARIO("Lookup::LazySimple")
{
	using Grid = Impl::Lookup::LazySingleListSingleIdx<3>;

	GIVEN("a 3x3x3 lazy single-index lookup grid with 3 tracking lists")
	{
		/// [LazySingleLookupGrid initialisation]

		Grid grid = Grid();

		THEN("the grid is initially inactive")
		{
			CHECK(grid.data().capacity() == 0);
			CHECK(grid.list().capacity() == 0);
		}
		/// [LazySingleLookupGrid initialisation]

		WHEN("the grid size is set")
		{
			grid.resize(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

			THEN("the grid is still inactive but now reports new size")
			{
				CHECK(grid.data().capacity() == 0);
				CHECK(grid.list().capacity() == 0);
				CHECK(grid.size() == Vec3i(3, 3, 3));
				CHECK(grid.offset() == Vec3i(-1,-1,-1));
			}

			/// [LazySingleLookupGrid activation]
			WHEN("the grid is activated")
			{
				grid.activate();

				THEN("memory is allocated and the grid filled with background value")
				{
					CHECK(grid.data().size() == 3 * 3 * 3);
				}

				THEN("queries return the NULL background value")
				{
					CHECK(grid.get(Vec3i(1, 1, 1)) == Felt::null_idx);
					CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
				}

				AND_WHEN("we track a position to be tracked to list 1")
				{
					grid.track(grid.index(Vec3i(1, 1, 1)));

					THEN("that position's value is updated and is added to the tracking list")
					{
						CHECK(grid.get(Vec3i(1, 1, 1)) == 0);
						CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
						CHECK(grid.list()[0] == grid.index(Vec3i(1, 1, 1)));
					}

					AND_WHEN("the grid is deactivated")
					{
						grid.deactivate();

						THEN("the grid is once again inactive")
						{
							CHECK(grid.data().capacity() == 0);
							CHECK(grid.list().capacity() == 0);
							CHECK(grid.data().capacity() == 0);
							CHECK(grid.list().capacity() == 0);
						}
					}
				}
				/// [LazySingleLookupGrid activation]
			}
		}
	}
}


SCENARIO("Lookup::LazySingle")
{
	using Grid = Impl::Lookup::LazyMultiListSingleIdx<3, 3>;

	GIVEN("a 3x3x3 lazy single-index lookup grid with 3 tracking lists")
	{
		/// [LazySingleLookupGrid initialisation]

		Grid grid = Grid();

		THEN("the grid is initially inactive")
		{
			CHECK(grid.data().capacity() == 0);
			CHECK(grid.list(0).capacity() == 0);
			CHECK(grid.list(1).capacity() == 0);
			CHECK(grid.list(2).capacity() == 0);
		}
		/// [LazySingleLookupGrid initialisation]

		WHEN("the grid size is set")
		{
			grid.resize(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

			THEN("the grid is still inactive but now reports new size")
			{
				CHECK(grid.data().capacity() == 0);
				CHECK(grid.list(0).capacity() == 0);
				CHECK(grid.list(1).capacity() == 0);
				CHECK(grid.list(2).capacity() == 0);
				CHECK(grid.size() == Vec3i(3, 3, 3));
				CHECK(grid.offset() == Vec3i(-1,-1,-1));
			}

			THEN("queries return the NULL background value")
			{
				CHECK(grid.get(Vec3i(1, 1, 1)) == Felt::null_idx);
				CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
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
					CHECK(grid.get(Vec3i(1, 1, 1)) == Felt::null_idx);
					CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
				}

				AND_WHEN("we track a position to be tracked to list 1")
				{
					grid.track(grid.index(Vec3i(1, 1, 1)), 1);

					THEN("that position's value is updated and is added to the tracking list")
					{
						CHECK(grid.get(Vec3i(1, 1, 1)) == 0);
						CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
						CHECK(grid.list(1)[0] == grid.index(Vec3i(1, 1, 1)));
					}

					AND_WHEN("the grid is deactivated")
					{
						grid.deactivate();

						THEN("the grid is once again inactive")
						{
							CHECK(grid.data().capacity() == 0);
							CHECK(grid.list(0).capacity() == 0);
							CHECK(grid.list(1).capacity() == 0);
							CHECK(grid.list(2).capacity() == 0);
						}

						THEN("queries once again return the NULL background value")
						{
							CHECK(grid.get(Vec3i(1,1,1)) == Felt::null_idx);
							CHECK(grid.get(Vec3i(0, 1, 1)) == Felt::null_idx);
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


SCENARIO("Tracked::LazySingleByValue")
{
	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 3.14159")
	{
		using Grid = Impl::Tracked::LazyMultiListSingleIdxByValue<Distance, 3, 3>;

//		static_assert(
//			std::experimental::is_detected<has_reset_t, GridType>::value,
//			"Tracked grids with a single lookup index per grid node should have a reset method."
//		);

		Grid grid = Grid(3.14159f);

		const PosIdx NULL_IDX = Felt::null_idx;

		THEN("the data grid and associated lookup grid state is zero size and inactive")
		{
			CHECK(grid.data().size() == 0);
			CHECK(grid.lookup().data().size() == 0);
		}

		THEN("there is an alias method to the lookup grid's lists")
		{
			CHECK(&grid.lookup().list(0) == &grid.list(0));
			CHECK(&grid.lookup().list(1) == &grid.list(1));
			CHECK(&grid.lookup().list(2) == &grid.list(2));
		}

		WHEN("the grid is resized")
		{
			grid.resize(Vec3i(3, 3, 3), Vec3i(-1,-1,-1));

			THEN("the data grid and associated lookup grid report new size and remains inactive")
			{
				CHECK(grid.data().capacity() == 0);
				CHECK(grid.lookup().data().capacity() == 0);

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
						CHECK(grid.data().capacity() == 0);
						CHECK(grid.get(Vec3i(1,1,1)) == 3.14159f);
						CHECK(grid.lookup().data().capacity() == 0);
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
						CHECK(grid.lookup().list(0).capacity() == 0);
						CHECK(grid.lookup().list(1).size() == 1);
						CHECK(grid.lookup().list(1)[0] == grid.index(Vec3i(1,1,1)));
						CHECK(grid.lookup().list(2).capacity() == 0);
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
		using Leaf = std::vector<int>;
		using Grid = Impl::Tracked::LazyMultiListSingleIdxByValue<Leaf, 3, 3>;

		Grid grid(Leaf{1,2,3});
		grid.resize(Vec3i(9,9,9), Vec3i(-4,-4,-4));
		grid.activate();

		WHEN("a value is set with an lvalue reference")
		{
			Leaf move_me{5,6,7};
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
			Leaf move_me{5,6,7};
			Leaf copied = move_me;
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
		using Grid = Impl::Tracked::MultiListMultiIdxByRef<Distance, 3, 3>;
		using IndexTuple = Tuple<ListIdx, 3>;

//		static_assert(
//			!std::experimental::is_detected<has_reset_t, GridType>::value,
//			"Tracked grids with multiple lookup indices per grid node should not have a reset"
//			" method."
//		);

		Grid grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), 0);

		THEN("the grid size is as expected and is initialised to all zero")
		{
			CHECK(grid.data().size() == 9*9*9);
			for (const Distance val : grid.data())
				CHECK(val == 0);
		}

		THEN("the associated lookup grid's size is as expected and initialised to NULL indices")
		{
			CHECK(grid.lookup().data().size() == 9*9*9);
			for (const IndexTuple& val : grid.lookup().data())
				CHECK(val == IndexTuple::Constant(Felt::null_idx));
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
				CHECK(grid.lookup().get(Vec3i(2,2,2)) == IndexTuple(Felt::null_idx, 0, 0));
				CHECK(grid.lookup().list(0).capacity() == 0);
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
		using Leaf = std::vector<int>;
		using Grid = Impl::Tracked::MultiListMultiIdxByRef<Leaf, 3, 3>;
		Grid grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), Leaf{1,2,3});

		WHEN("a value is set with an lvalue reference")
		{
			Leaf move_me{5,6,7};
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
			Leaf move_me{5,6,7};
			Leaf copied = move_me;
			int* pdata = &move_me[0];
			grid.get(Vec3i(2,2,2)) = std::move(move_me);

			THEN("the data from the input has been moved into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == copied);
				CHECK(&grid.get(Vec3i(2,2,2))[0] == pdata);
			}
		}

		WHEN("a value is tracked with an lvalue reference")
		{
			Leaf move_me{5,6,7};
			int* pdata = &move_me[0];
			grid.track(move_me, grid.index(Vec3i(2,2,2)), 1);

			THEN("the data from the input has been copied into the grid")
			{
				CHECK(grid.get(Vec3i(2,2,2)) == move_me);
				CHECK(&grid.get(Vec3i(2,2,2))[0] != pdata);
			}
		}

		WHEN("a value is tracked with an rvalue reference")
		{
			Leaf move_me{5,6,7};
			Leaf copied = move_me;
			int* pdata = &move_me[0];
			grid.track(std::move(move_me), grid.index(Vec3i(2,2,2)), 1);

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
	using Grid = Impl::Partitioned::Lookup<3, 3>;
	using ChildrenGrid = Grid::ChildrenGrid;

	static_assert(
		std::is_same<
			ChildrenGrid, Impl::Tracked::MultiListMultiIdxByRef<Impl::Lookup::LazyMultiListSingleIdx<3, 3>, 3, 3>
		>::value,
		"Children grid of partitioned lookup must be a multi-list multi-index tracked grid with"
		" lazily activated lookup sub-grids as the leaf type."
	);

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions")
	{
		Grid grid(Vec3i(9,9,9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3));

		THEN("the children tracking grid has been initialised with lazy subgrids")
		{
			CHECK(grid.children().data().size() == 3*3*3);
			CHECK(grid.children().lookup().data().size() == 3*3*3);

			CHECK(grid.children().get(Vec3i(-1,-1,-1)).data().capacity() == 0);
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
			Grid grid_master(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3));

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
						).get(pos_active_because_master) == null_idx
					);
					CHECK(grid.children().get(pos_child_list_0).get(pos_list_0) == null_idx);
					CHECK(grid.children().get(pos_child_list_1).get(pos_list_1) == null_idx);
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
	using Grid = Impl::Partitioned::Tracked::Simple<int, 3, 3>;

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions with background value -42")
	{
		Grid grid(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), -42);

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
			CHECK(grid.children().get(pos12_child_idx).data().capacity() == 0);
			CHECK(grid.children().get(pos12_child_idx).get(pos1_idx) == -42);
			CHECK(grid.children().get(pos12_child_idx).lookup().data().capacity() == 0);
			CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == null_idx);
			CHECK(grid.children().get(pos12_child_idx).get(pos2_idx) == -42);
			CHECK(grid.children().get(pos12_child_idx).lookup().get(pos2_idx) == null_idx);

			CHECK(grid.children().get(pos3_child_idx).is_active() == false);
			CHECK(grid.children().get(pos3_child_idx).data().capacity() == 0);
			CHECK(grid.children().get(pos3_child_idx).get(pos3_idx) == -42);
			CHECK(grid.children().get(pos3_child_idx).lookup().data().capacity() == 0);
			CHECK(grid.children().get(pos3_child_idx).lookup().get(pos3_idx) == null_idx);

			CHECK(grid.children().get(pos4_child_idx).is_active() == false);
			CHECK(grid.children().get(pos4_child_idx).data().capacity() == 0);
			CHECK(grid.children().get(pos4_child_idx).get(pos4_idx) == -42);
			CHECK(grid.children().get(pos4_child_idx).lookup().data().capacity() == 0);
			CHECK(grid.children().get(pos4_child_idx).lookup().get(pos4_idx) == null_idx);
		}

		AND_WHEN("a mask grid is tracking some partitions")
		{
			Grid grid_master(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), 0);
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
					CHECK(grid.children().get(pos12_child_idx).lookup().get(pos1_idx) == null_idx);
					CHECK(grid.children().get(pos12_child_idx).get(pos2_idx) == -42);
					CHECK(grid.children().get(pos12_child_idx).lookup().get(pos2_idx) == null_idx);

					CHECK(grid.children().get(pos3_child_idx).is_active() == true);
					CHECK(grid.children().get(pos3_child_idx).data().size() == 3*3*3);
					CHECK(grid.children().get(pos3_child_idx).get(pos3_idx) == -42);
					CHECK(grid.children().get(pos3_child_idx).lookup().data().size() == 3*3*3);
					CHECK(grid.children().get(pos3_child_idx).lookup().get(pos3_idx) == null_idx);

					CHECK(grid.children().lookup().list(0).size() == 2);
					CHECK(grid.children().lookup().list(1).capacity() == 0);
					CHECK(grid.children().lookup().list(2).capacity() == 0);
				}

				THEN("other children are remain inactive")
				{
					CHECK(grid.children().get(pos4_child_idx).is_active() == false);
					CHECK(grid.children().get(pos4_child_idx).data().capacity() == 0);
					CHECK(grid.children().get(pos4_child_idx).get(pos4_idx) == -42);
					CHECK(grid.children().get(pos4_child_idx).lookup().data().capacity() == 0);
					CHECK(grid.children().get(pos4_child_idx).lookup().get(pos4_idx) == null_idx);
				}

				AND_WHEN("a position in an active partition is tracked by index for list 1")
				{
					grid.track(345, pos12_child_idx, pos1_idx, 1);

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

						CHECK(grid.children().get(pos12_child_idx).lookup().list(0).capacity() == 0);
						CHECK(grid.children().get(pos12_child_idx).lookup().list(2).capacity() == 0);

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
							CHECK(grid.children().get(pos4_child_idx).data().capacity() == 0);
						}

						THEN("children grid lookup is reset")
						{
							auto s_null_idxs = grid.children().lookup().s_null_idxs;
							CHECK(grid.children().lookup().list(0).size() == 0);
							CHECK(grid.children().lookup().list(1).size() == 0);
							CHECK(grid.children().lookup().list(2).size() == 0);
							CHECK(grid.children().lookup().get(pos12_child_idx) == s_null_idxs);
							CHECK(grid.children().lookup().get(pos3_child_idx) == s_null_idxs);
							CHECK(grid.children().lookup().get(pos4_child_idx) == s_null_idxs);
						}

						THEN("child grids' lookups are reset")
						{
							CHECK(grid.children().get(pos12_child_idx).lookup().list(0).size() == 0);
							CHECK(grid.children().get(pos12_child_idx).lookup().list(1).size() == 0);
							CHECK(grid.children().get(pos12_child_idx).lookup().list(2).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(0).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(1).size() == 0);
							CHECK(grid.children().get(pos3_child_idx).lookup().list(2).size() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(0).capacity() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(1).capacity() == 0);
							CHECK(grid.children().get(pos4_child_idx).lookup().list(2).capacity() == 0);
						}
					}
				}
			}
		}
	}
}


SCENARIO("Paritioned::Tracked::Numeric")
{
	using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 3, 3>;

	GIVEN("a 9x9x9 grid with (-4,-4,-4) offset in 3x3x3 partitions with background value -42")
	{
		Grid grid(Vec3i(9, 9, 9), Vec3i(-4,-4,-4), Vec3i(3, 3, 3), -42);

		const Vec3i pos123_child(-1, -1, -1);
		const Vec3i pos4_child(1, 1, 1);
		const Vec3i pos_inactive_child(1, 0, -1);
		const PosIdx pos123_child_idx = grid.children().index(pos123_child);
		const PosIdx pos4_child_idx = grid.children().index(pos4_child);

		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos3(-4, -3, -4);
		const Vec3i pos4(4, 4, 4);
		const PosIdx pos1_idx = grid.children().get(pos123_child_idx).index(pos1);
		const PosIdx pos2_idx = grid.children().get(pos123_child_idx).index(pos2);
		const PosIdx pos3_idx = grid.children().get(pos123_child_idx).index(pos3);
		const PosIdx pos4_idx = grid.children().get(pos4_child_idx).index(pos4);

		Grid::Child& child = grid.children().get(pos123_child_idx);

		THEN("querying outside of the grid returns background value")
		{
			CHECK(grid.get(Vec3i(20,20,20)) == -42);
			CHECK(grid.get(Vec3i(-20,-20,-20)) == -42);
		}

		WHEN("positions in the same partition are tracked")
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
				CHECK(grid.children().lookup().list(2).capacity() == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[0] == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
				CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
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

			AND_WHEN("all grid points in list 0 are iterated over and their values stored")
			{
				std::vector<Distance> aval;

				grid.leafs(0, [&aval, &grid](const auto& pos_) {
					aval.push_back(grid.get(pos_));
				});

				THEN("array of recorded values is correct and in expected order")
				{
					std::vector<Distance>aval_check {345, 789};
					CHECK(aval_check == aval);
				}
			}


			AND_WHEN("a position is untracked")
			{
				grid.untrack(-100, pos123_child_idx, pos1_idx, 0);

				THEN("the partition is still tracked by the children grid")
				{
					CHECK(grid.children().lookup().list(0).size() == 1);
					CHECK(grid.children().lookup().list(1).size() == 1);
					CHECK(grid.children().lookup().list(2).capacity() == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[0] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
				}

				THEN("the child now tracks just the other positions")
				{
					CHECK(child.lookup().list(0).size() == 1);
					CHECK(child.lookup().list(1).size() == 1);
					CHECK(child.lookup().list(2).capacity() == 0);
					CHECK(child.lookup().list(0)[0] == pos2_idx);
					CHECK(child.lookup().list(1)[0] == pos3_idx);
					CHECK(child.lookup().get(pos1_idx) == null_idx);
					CHECK(child.lookup().get(pos2_idx) == 0);
					CHECK(child.lookup().get(pos3_idx) == 0);
				}

				THEN("value at the untracked position is set to the passed background value")
				{
					CHECK(child.get(pos1_idx) == -100);
					CHECK(child.get(pos2_idx) == 789);
					CHECK(child.get(pos3_idx) == 123);
				}

				AND_WHEN("another postion is untracked from the same tracking list")
				{
					grid.untrack(-102, pos123_child_idx, pos2_idx, 0);

					THEN("the partition is still tracked by the children grid")
					{
						CHECK(grid.children().lookup().list(0).size() == 0);
						CHECK(grid.children().lookup().list(1).size() == 1);
						CHECK(grid.children().lookup().list(2).capacity() == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[0] == null_idx);
						CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
					}

					THEN("the child now tracks just the final position")
					{
						CHECK(child.lookup().list(0).size() == 0);
						CHECK(child.lookup().list(1).size() == 1);
						CHECK(child.lookup().list(2).capacity() == 0);
						CHECK(child.lookup().list(1)[0] == pos3_idx);
						CHECK(child.lookup().get(pos1_idx) == null_idx);
						CHECK(child.lookup().get(pos2_idx) == null_idx);
						CHECK(child.lookup().get(pos3_idx) == 0);
					}

					AND_WHEN("the final postion is untracked from tracking")
					{
						grid.untrack(-999, pos123_child_idx, pos3_idx, 1);

						THEN("the partition is no longer tracked by the children grid")
						{
							CHECK(grid.children().lookup().list(0).size() == 0);
							CHECK(grid.children().lookup().list(1).size() == 0);
							CHECK(grid.children().lookup().list(2).capacity() == 0);
							CHECK(grid.children().lookup().get(pos123_child_idx)[0] == null_idx);
							CHECK(grid.children().lookup().get(pos123_child_idx)[1] == null_idx);
							CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
						}

						THEN("the child no longer tracks any positions")
						{
							CHECK(child.lookup().list(0).capacity() == 0);
							CHECK(child.lookup().list(1).capacity() == 0);
							CHECK(child.lookup().list(2).capacity() == 0);
							CHECK(child.lookup().get(pos1_idx) == null_idx);
							CHECK(child.lookup().get(pos2_idx) == null_idx);
							CHECK(child.lookup().get(pos3_idx) == null_idx);
						}

						THEN("the child has been deactivated")
						{
							CHECK(child.is_active() == false);
							CHECK(child.data().capacity() == 0);
						}

						THEN("all positions report the final background value")
						{
							CHECK(child.get(pos1_idx) == -999);
							CHECK(child.get(pos2_idx) == -999);
							CHECK(child.get(pos3_idx) == -999);
						}
					}
				}
			} // End WHEN "a position is untracked"

			AND_WHEN("a position is moved from one tracking list to another")
			{
				grid.retrack(pos123_child_idx, pos1_idx, 0, 1);

				THEN("the value at the retracked position is unchanged")
				{
					CHECK(child.get(pos1_idx) == 345);
				}

				THEN("the children grid still tracks the list")
				{
					CHECK(grid.children().lookup().list(0).size() == 1);
					CHECK(grid.children().lookup().list(1).size() == 1);
					CHECK(grid.children().lookup().list(2).capacity() == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[0] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
					CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
				}

				THEN("the position has been appended to the other list")
				{
					CHECK(child.lookup().list(0).size() == 1);
					CHECK(child.lookup().list(1).size() == 2);
					CHECK(child.lookup().list(2).capacity() == 0);
					CHECK(child.lookup().list(0)[0] == pos2_idx);
					CHECK(child.lookup().list(1)[0] == pos3_idx);
					CHECK(child.lookup().list(1)[1] == pos1_idx);
					CHECK(child.lookup().get(pos1_idx) == 1);
					CHECK(child.lookup().get(pos2_idx) == 0);
					CHECK(child.lookup().get(pos3_idx) == 0);
				}

				AND_WHEN("the final position in a list is moved to another list")
				{
					grid.retrack(pos123_child_idx, pos2_idx, 0, 1);

					THEN("the children grid no longer tracks the empty list")
					{
						CHECK(grid.children().lookup().list(0).size() == 0);
						CHECK(grid.children().lookup().list(1).size() == 1);
						CHECK(grid.children().lookup().list(2).capacity() == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[0] == null_idx);
						CHECK(grid.children().lookup().get(pos123_child_idx)[1] == 0);
						CHECK(grid.children().lookup().get(pos123_child_idx)[2] == null_idx);
					}
				}
			}
		} // End WHEN positions in the same partition are tracked

		WHEN("values are changed in a child sub-grid")
		{
			child.activate();
			child.set(pos1_idx, 3.14f);

			THEN("the value can be retrieved by position vector through the parent grid")
			{
				CHECK(grid.get(pos1) == 3.14f);
			}

			AND_WHEN("a snapshot is taken of the grid")
			{
				static_assert(
					std::is_same<
						Grid::SnapshotPtr, std::unique_ptr< Impl::Grid::Snapshot<Distance, 3> >
					>::value,
					"Snapshot grid must be smart pointer to a simple Grid::Snapshot."
				);
				Grid::SnapshotPtr psnapshot = grid.snapshot();

				THEN("the snapshot grid is of the correct size")
				{
					CHECK(psnapshot->size() == Vec3i(9,9,9));
					CHECK(psnapshot->offset() == Vec3i(-4,-4,-4));
				}

				THEN("the snapshot reports the same values as the main grid")
				{
					CHECK(psnapshot->get(pos1_idx) == 3.14f);
					CHECK(psnapshot->get(pos2_idx) == -42);
					CHECK(psnapshot->get(pos3_idx) == -42);
					CHECK(psnapshot->get(pos4_idx) == -42);
				}

				AND_WHEN("a value is changed in the snapshot and the snapshot is flushed")
				{
					psnapshot->set(pos4, 567);
					grid.snapshot(psnapshot);

					THEN("the main grid value is updated")
					{
						CHECK(grid.get(pos4) == 567);
					}

					THEN("only required children are activated")
					{
						CHECK(child.is_active() == true);
						CHECK(grid.children().get(pos4_child_idx).is_active() == true);
						CHECK(grid.children().get(pos_inactive_child).is_active() == false);
					}

				}
			}
		} // WHEN("values are changed in a child sub-grid")

		WHEN("grid is modified and serialised to disk")
		{
			grid.track(345, pos1, 0);
			grid.track(789, pos2, 1);
			grid.track(123, pos3, 2);
			grid.track(234, pos4, 1);

			std::ofstream writer{"/tmp/grid.bin", std::ios::binary};
			grid.write(writer);
			writer.flush();

			AND_WHEN("grid is loaded")
			{
				std::ifstream reader{"/tmp/grid.bin", std::ios::binary};
				Grid grid_loaded{Grid::read(reader)};

				THEN("loaded grid size matches")
				{
					CHECK(grid_loaded.size() == grid.size());
					CHECK(grid_loaded.offset() == grid.offset());
				}

				THEN("loaded child subgrid background value matches")
				{
					for (
						PosIdx pos_idx_child = 0;
						pos_idx_child < grid.children().data().size(); pos_idx_child++
					) {
						CHECK(
							grid.children().get(pos_idx_child).background() ==
								grid_loaded.children().get(pos_idx_child).background()
						);
					}
				}

				THEN("loaded child data matches")
				{
					for (
						PosIdx pos_idx_child = 0;
						pos_idx_child < grid.children().data().size(); pos_idx_child++
					) {
						CHECK(
							grid.children().get(pos_idx_child).data() ==
								grid_loaded.children().get(pos_idx_child).data()
						);
					}
				}

				THEN("loaded child lists match")
				{
					for (TupleIdx list_idx = 0; list_idx < Grid::num_lists; list_idx++)
					{
						CHECK(
							grid.children().lookup().list(list_idx) ==
								grid_loaded.children().lookup().list(list_idx)
						);
						for (
							PosIdx pos_idx_child = 0;
							pos_idx_child < grid.children().data().size(); pos_idx_child++
						) {
							CHECK(
								grid.children().get(pos_idx_child).list(list_idx) ==
									grid_loaded.children().get(pos_idx_child).list(list_idx)
							);
						}
					}
				}
			}
		}
	} // End GIVEN a 9x9x9 grid.

	GIVEN("A 1D grid type and input vector of values")
	{
		using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 1, 3>;
		std::vector<Distance> input = { 1.0f, 0 };

		WHEN("we interpolate a distance of 0.3 between the vector values")
		{
			using Vec1f = Eigen::Matrix<Distance, 1, 1>;
			const Vec1f pos(0.3);

			Grid::interp(input, pos);

			THEN("the input vector now contains the single interpolated value")
			{
				CHECK(input.size() == 1);
				CHECK(input[0] == 0.7f);
			}
		}
	}

	GIVEN("A 2D grid type and input vector of values")
	{
		using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 2, 3>;

		std::vector<Distance> input = std::vector<Distance>(4);
		input[0 /*00*/] = 2.0f;
		input[1 /*01*/] = 0;
		input[2 /*10*/] = 0.0f;
		input[3 /*11*/] = 1.0;

		WHEN("we bilinearly interpolate a distance of (0.8, 0.5) between the vector values")
		{
			const Vec2f pos(0.8f, 0.5f);

			Grid::interp(input, pos);

			THEN("the input vector now contains the correct interpolated values")
			{
				CHECK(input.size() == 2);
				CHECK(input[0] == Approx(0.4f));
				CHECK(input[1] == Approx(0.8f));

				AND_WHEN("we interpolate along this line")
				{
					Grid::interp(input, pos);

					THEN("the input vector now contains the final interpolated value")
					{
						CHECK(input.size() == 1);
						CHECK(input[0] == Approx(0.6f));
					}
				}
			}
		}
	}

	GIVEN("A 3D grid type and input vector of values")
	{
		/**
				  011----111
				 /|		  /|
				010----011 |
				| 100----|101
				|/		 |/
				000----001
		*/
		using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 3, 3>;

		std::vector<Distance> input = std::vector<Distance>(8);
		input[0 /**000*/] = 0.0f;
		input[1 /**001*/] = 0.8f;
		input[2 /**010*/] = 1.0f;
		input[3 /**011*/] = 1.0f;
		input[4 /**100*/] = 0.0f;
		input[5 /**101*/] = 0.0f;
		input[6 /**110*/] = 1.0f;
		input[7 /**111*/] = 1.0f;

		WHEN("we trilinearly interpolate a distance of (0.5, 0.75, 0.5) between the vector values")
		{
			Vec3f pos(0.5f, 0.75f, 0.5f);

			Grid::interp(input, pos);

			THEN("the input vector now contains the correct interpolated values")
			{
				CHECK(input.size() == 4);
				CHECK(input[0 /**00x*/] == 0.4f);
				CHECK(input[1 /**01x*/] == 1.0f);
				CHECK(input[2 /**10x*/] == 0.0f);
				CHECK(input[3 /**11x*/] == 1.0f);

				AND_WHEN(
					"we bilinearly interpolate a distance of (0.8, 0.5) between these vector values"
				) {
					Grid::interp(input, pos);

					THEN("the input vector now contains the correct interpolated values")
					{
						CHECK(input.size() == 2);
						CHECK(input[0 /**0yx*/] == Approx(0.85f));
						CHECK(input[1 /**1yx*/] == Approx(0.75f));

						AND_WHEN("we interpolate along this line")
						{
							Grid::interp(input, pos);

							THEN("the input vector now contains the final interpolated value")
							{
								CHECK(input.size() == 1);
								CHECK(input[0 /**zyx*/] == Approx(0.8f));
							}
						}
					}
				}
			}
		}
	} // End GIVEN A 3D grid type and input vector of values.


	GIVEN("a 3x3 grid with (-1,-1) offset and background value of 0")
	{
		using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 2, 3>;
		Grid grid(Vec2i(3, 3), Vec2i(-1,-1), Vec2i(3,3), 0);
		// Only a single partition in this case.
		grid.children().get(0).activate();

		THEN("the spatial resolution has a default of 1")
		{
			CHECK(grid.dx() == 1.0f);

			AND_WHEN("we change the spatial resolution to 2")
			{
				grid.dx(2.0f);

				THEN("the resolution is reported as 2")
				{
					CHECK(grid.dx() == 2.0f);
				}
			}
		}


		WHEN("we modify some values near the centre")
		{
			grid.set(Vec2i(-1,-1), 1.0f);
			grid.set(Vec2i(-1,0), 1.0f);
			grid.set(Vec2i(0,1), 2.0f);
			grid.set(Vec2i(1,1), 2.0f);

			AND_WHEN("we interpolate at some real locations using explicit function calls")
			{
				const Distance val1 = grid.interp(Vec2f(0.0f, 0.0f));
				const Distance val2 = grid.interp(Vec2f(-0.5f, -0.5f));
				const Distance val3 = grid.interp(Vec2f(0.5f, 0.5f));

				THEN("the interpolated values are correct")
				{
					CHECK(val1 == Approx(0.0f));
					CHECK(val2 == Approx(0.5f));
					CHECK(val3 == Approx(1.0f));
				}
			}

			AND_WHEN("we implicitly interpolate at a real location using value getter")
			{
				const Distance val = grid.get(Vec2f(0.5f, 0.5f));

				THEN("the interpolated value is correct")
				{
					CHECK(val == Approx(1.0f));
				}
			}
		}


		WHEN("we calculate the forward difference gradient at the centre")
		{
			const Vec2f grad = grid.gradF(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we calculate the backward difference gradient at the centre")
		{
			const Vec2f grad = grid.gradB(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we calculate the central difference gradient at the centre")
		{
			const Vec2f grad = grid.gradC(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we set the central grid location to 1")
		{
			grid.set(Vec2i(0,0), 1.0f);

			AND_WHEN("we calculate the forward difference gradient at the centre")
			{
				const Vec2f grad = grid.gradF(Vec2i(0,0));

				THEN("the gradient is negative")
				{
					CHECK(grad == Vec2f(-1, -1));
				}
			}

			AND_WHEN("we calculate the backward difference gradient at the centre")
			{
				const Vec2f grad = grid.gradB(Vec2i(0,0));

				THEN("the gradient is positive")
				{
					CHECK(grad == Vec2f(1, 1));
				}
			}

			AND_WHEN("we calculate the central difference gradient at the centre")
			{
				const Vec2f grad = grid.gradC(Vec2i(0,0));

				THEN("the gradient is zero")
				{
					CHECK(grad == Vec2f(0, 0));
				}
			}
		}


		WHEN("grid data is set using an initialiser list")
		{
			grid = {
			// <  -	y  + >
				1,	2,	3,	// -
				4,	5,	6,	// x
				7,	8,	9	// +
			};

			THEN("the underlying data has been updated")
			{
				CHECK(grid.get(Vec2i(-1,-1)) == 1);
				CHECK(grid.get(Vec2i(-1, 0)) == 2);
				CHECK(grid.get(Vec2i(-1, 1)) == 3);
				CHECK(grid.get(Vec2i( 0,-1)) == 4);
				CHECK(grid.get(Vec2i( 0, 0)) == 5);
				CHECK(grid.get(Vec2i( 0, 1)) == 6);
				CHECK(grid.get(Vec2i( 1,-1)) == 7);
				CHECK(grid.get(Vec2i( 1, 0)) == 8);
				CHECK(grid.get(Vec2i( 1, 1)) == 9);
			}
		}

		WHEN("we have a (positively directed) entropic flow")
		{
			grid = {
			//<	  -	y  + >
				0,	0,	0,	// -
				0,	1,	3,	// x
				0,	3,	0	// +
			};

			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(1, 1));
			}
		}

		WHEN("we have a (negatively directed) entropic flow")
		{
			grid = {
			//<	  -	y  + >
				0,	3,	0,	// -
				3,	1,	0,	// x
				0,	0,	0	// +
			};

			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(-1, -1));
			}
		}

		WHEN("we have a positively divergent flow")
		{
			grid = {
			//<	  -	y  + >
				0,	2,	0,	// -
				2,	1,	3,	// x
				0,	3,	0	// +
			};
			THEN("the entropy satisfying gradient is zero")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(0, 0));
			}
		}

		WHEN("we have a negatively divergent flow")
		{
			grid = {
			//<	  -	y  + >
				0,	6,	0,	// -
				6,	9,	1,	// x
				0,	1,	0	// +
			};
			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(-5, -5));
			}
		}

		WHEN("we set up a positive divergence at the centre")
		{
			grid = {
				1,	1,	1,
				1,	0,	1,
				1,	1,	1
			};

			THEN("the curvature is at its maximum")
			{
				CHECK(grid.curv(Vec2i(0,0)) == 2);
			}
		}

		WHEN("we set up a 'corner' at the centre")
		{
			grid = {
				 1,	 1,	 1,
				 0,	 0,  1,
				-1,	 0,	 1
			};

			THEN("the curvature is 1")
			{
				CHECK(grid.curv(Vec2i(0,0)) == 1);
			}
		}

	} // End GIVEN("a 3x3 grid with (-1,-1) offset and background value of 0")


	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 0")
	{
		using Grid = Impl::Partitioned::Tracked::Numeric<Distance, 3, 3>;

		Grid grid(Vec3i(3, 3, 3), Vec3i(-1,-1,-1), Vec3i(3,3,3), 0);
		grid.children().get(0).activate();

		WHEN("we set up a gradient about the centre")
		{
			grid = {
				0,	0,	0,
				0,	2,	0,
				0,	0,	0,

				0,  0,	0,
				0,	1,	2,
				0,	0,	0,

				0,	0,	0,
				0,	0,	0,
				0,	0,	0
			};

			THEN("the forward difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradF(Vec3i(0,0,0));
				CHECK(grad == Vec3f(-1, -1, 1));
			}

			THEN("the backward difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradB(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 1, 1));
			}

			THEN("the central difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradC(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 0, 1));
			}

			THEN("the safe gradient at the centre uses the central difference")
			{
				const Vec3f grad = grid.grad(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 0, 1));
			}

			THEN("the safe gradient at the bottom face uses the forward difference")
			{
				const Vec3f grad = grid.grad(Vec3i(0,-1,0));

				CHECK(grad == Vec3f(0, 1, 0));
			}

			THEN("the safe gradient at the right-forward edge uses the backward difference")
			{
				const Vec3f grad = grid.grad(Vec3i(1,0,1));

				CHECK(grad == Vec3f(-2, 0, 0));
			}

			AND_WHEN("we decrease the spatial resolution")
			{
				grid.dx(2);

				THEN("the forward difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradF(Vec3i(0,0,0));
					CHECK(grad == Vec3f(-0.5, -0.5, 0.5));
				}

				THEN("the backward difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradB(Vec3i(0,0,0));

					CHECK(grad == Vec3f(-0.5, 0.5, 0.5));
				}

				THEN("the central difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradC(Vec3i(0,0,0));

					CHECK(grad == Vec3f(-0.5, 0, 0.5));
				}
			}
		}

		WHEN("we set up a positive divergence at the centre")
		{
			grid = {
				1,	1,	1,
				1,	1,	1,
				1,	1,	1,

				1,	1,	1,
				1,	0,	1,
				1,	1,	1,

				1,	1,	1,
				1,	1,	1,
				1,	1,	1
			};

			THEN("the divergence is calculated correctly")
			{
				CHECK(grid.divergence(Vec3i(0,0,0)) == Approx(6));
			}

			AND_WHEN("we decrease the spatial resolution")
			{
				grid.dx(2);

				THEN("the divergence is calculated correctly")
				{
					CHECK(grid.divergence(Vec3i(0,0,0)) == Approx(0.75));
				}
			}

			THEN("the mean curvature is at the maximum")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 3);
			}
		}

		WHEN("we set up a 'corner' along two dimensions")
		{
			grid = {
				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1
			};

			THEN("the mean curvature is 1")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 1);
			}
		}

		WHEN("we set up a sharp 'corner' along all three dimensions")
		{
			grid = {
				 1,	 1,	 1,
				 1,	 1,	 1,
				 1,	 1,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				 0,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1
			};

			THEN("the mean curvature is 1.5")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 1.5);
			}
		}

		WHEN("neighbours of an arbitrary position vector are cycled")
		{
			std::vector<Vec3i> apos;
			Grid::neighs(Vec3i(-67, 54, 3), [&apos](const Vec3i& pos) {
				apos.push_back(pos);
			});

			THEN("the positions are cycled in the expected order")
			{
				CHECK(
					apos == (std::vector<Vec3i>{
						Vec3i(-68, 54, 3), Vec3i(-66, 54, 3),
						Vec3i(-67, 53, 3), Vec3i(-67, 55, 3),
						Vec3i(-67, 54, 2), Vec3i(-67, 54, 4)
					})
				);
			}
		}
	}

}

