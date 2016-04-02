#include <boost/test/unit_test.hpp>

#define _TESTING

#include "Felt/SharedLookupGrid.hpp"

using namespace felt;

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
