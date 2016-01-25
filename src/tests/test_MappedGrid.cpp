#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "Felt/LookupGrid.hpp"
#include "Felt/SharedLookupGrid.hpp"
#include "Felt/TrackedGrid.hpp"


using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_LookupGrid)

	BOOST_AUTO_TEST_CASE(initialise_and_populate_single_tracking_list)
	{
		typedef LookupGrid<3> Grid_t;
		Grid_t grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		Vec3i pos1(1, 0, -1);
		Vec3i pos2(2, 1, 0);
		Vec3i pos3(3, -1, 0);
		Vec3i pos4(4, -1, 2);

		// Check initialised to zero length with null index references.
		BOOST_REQUIRE_EQUAL(grid.list().size(), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), Grid_t::NULL_IDX);

		// Add the positions to the array and set index lookup values.
		grid.add(pos1);
		grid.add(pos2);
		grid.add(pos3);
		grid.add(pos4);


		// Check the positions were added to the array and their respective
		// index lookups are as expected.
		BOOST_REQUIRE_EQUAL(grid.list().size(), 4);
		BOOST_CHECK_EQUAL(grid.list()[0], pos1);
		BOOST_CHECK_EQUAL(grid.list()[1], pos2);
		BOOST_CHECK_EQUAL(grid.list()[2], pos3);
		BOOST_CHECK_EQUAL(grid.list()[3], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), 2);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), 3);

		// Attempt to add the same position to the array again (i.e. duplicate).
		grid.add(pos2);

		// Ensure nothing changed.
		BOOST_REQUIRE_EQUAL(grid.list().size(), 4);
		BOOST_CHECK_EQUAL(grid.list()[0], pos1);
		BOOST_CHECK_EQUAL(grid.list()[1], pos2);
		BOOST_CHECK_EQUAL(grid.list()[2], pos3);
		BOOST_CHECK_EQUAL(grid.list()[3], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), 2);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), 3);

		// Remove a position by index.
		grid.remove(1);

		// Ensure position is removed from the array, the associated index
		// lookup is set to null, and the array restructured as expected.
		BOOST_REQUIRE_EQUAL(grid.list().size(), 3);
		BOOST_CHECK_EQUAL(grid.list()[0], pos1);
		BOOST_CHECK_EQUAL(grid.list()[1], pos4);
		BOOST_CHECK_EQUAL(grid.list()[2], pos3);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), 2);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), 1);

		// Remove a position by position (using index lookup).
		grid.remove(pos1);

		// Ensure as above that position is removed and lookup nulled.
		BOOST_REQUIRE_EQUAL(grid.list().size(), 2);
		BOOST_CHECK_EQUAL(grid.list()[0], pos3);
		BOOST_CHECK_EQUAL(grid.list()[1], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), 1);

		// Reset the array.
		grid.reset();

		// Ensure array is zero size and remaining associated index lookups
		// have null value.
		BOOST_CHECK_EQUAL(grid.list().size(), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(0), Grid_t::NULL_IDX);
	}

	BOOST_AUTO_TEST_CASE(initialise_and_populatate_multiple_tracking_lists)
	{
		typedef LookupGrid<3, 3> Grid_t;
		Grid_t grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		// Add the positions to the array and set index lookup values.
		grid.add(pos1, 0);
		grid.add(pos1, 0); // Shouldn't do anything.
		grid.add(pos2, 1);
		grid.add(pos3, 1);
		grid.add(pos4, 2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 2);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos2);
		BOOST_CHECK_EQUAL(grid.list(1)[1], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(1), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), 0);

		grid.remove(pos2, 1);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), 0);

		grid.add(pos5, 2);
		grid.add(pos6, 2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 3);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL(grid.list(2)[1], pos5);
		BOOST_CHECK_EQUAL(grid.list(2)[2], pos6);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos5)(2), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos6)(2), 2);

		grid.remove(pos4, 2);
		grid.remove(0, 0);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 2);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[1], pos5);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos6);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5)(2), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos6)(2), 0);

		grid.reset(2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5)(2), Grid_t::NULL_IDX);
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_SharedLookupGrid)
	BOOST_AUTO_TEST_CASE(intialise_and_populate)
	{
		using GridType = SharedLookupGrid<3, 3>;
		GridType grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		// Add the positions to the array and set index lookup values.
		grid.add(pos1, 0);
		grid.add(pos2, 1);
		grid.add(pos3, 1);
		grid.add(pos4, 2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 2);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos2);
		BOOST_CHECK_EQUAL(grid.list(1)[1], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), 0);

		grid.remove(pos2, 1);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), 0);

		grid.add(pos5, 2);
		grid.add(pos6, 2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 3);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL(grid.list(2)[1], pos5);
		BOOST_CHECK_EQUAL(grid.list(2)[2], pos6);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos5), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos6), 2);

		grid.remove(pos4, 2);
		grid.remove(0, 0);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 2);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[1], pos5);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos6);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos6), 0);

		grid.reset(2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), GridType::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5), GridType::NULL_IDX);

	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_TrackedGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		TrackedGrid<FLOAT, 3, 3> grid(Vec3u(9,9,9), Vec3i(-4,-4,-4));
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_LazyLookupGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazyLookupGrid initialisation]
		// ==== Setup ====
		LazyLookupGrid<3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
		const Vec3u NULL_IDX_DATA = LazyLookupGrid<3, 3>::Traits::NULL_IDX_DATA;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), NULL_IDX_DATA);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), NULL_IDX_DATA);
		/// [LazyLookupGrid initialisation]
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_LazySharedLookupGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazySharedLookupGrid initialisation]
		// ==== Setup ====
		LazySharedLookupGrid<3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
		const UINT NULL_IDX_DATA = LazySharedLookupGrid<3, 3>::NULL_IDX;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), NULL_IDX_DATA);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), NULL_IDX_DATA);
		/// [LazySharedLookupGrid initialisation]
	}
BOOST_AUTO_TEST_SUITE_END()

