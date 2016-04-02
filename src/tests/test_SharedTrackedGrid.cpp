#include <boost/test/unit_test.hpp>

#define _TESTING

#include "Felt/SharedTrackedGrid.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_LazySharedTrackedGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazySharedTrackedGrid initialisation]
		// ==== Setup ====
		LazySharedTrackedGrid<FLOAT, 3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3);
		const UINT NULL_IDX = LazySharedTrackedGrid<FLOAT, 3, 3>::Lookup::NULL_IDX;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), 3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3);
		BOOST_CHECK_EQUAL(grid.lookup().is_active(), false);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 0);
		BOOST_CHECK_EQUAL(grid.lookup().background(), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.lookup().get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySharedTrackedGrid initialisation]
	}

	struct Fixture {
		const UINT NULL_IDX = LazySharedTrackedGrid<FLOAT, 3, 3>::Lookup::NULL_IDX;
		LazySharedTrackedGrid<FLOAT, 3, 3> grid;
		Fixture()
			: grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3)
		{}
	};

	BOOST_FIXTURE_TEST_CASE(activate_should_activate_lookup, Fixture)
	{
		/// [LazySharedTrackedGrid activate]
		// ==== Action ====
		grid.activate();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), true);
		BOOST_CHECK_EQUAL(grid.data().size(), 3*3*3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3);
		BOOST_CHECK_EQUAL(grid.lookup().is_active(), true);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 3*3*3);
		BOOST_CHECK_EQUAL(grid.lookup().get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySharedTrackedGrid activate]
	}
BOOST_AUTO_TEST_SUITE_END()

