#include <boost/test/unit_test.hpp>

#define _TESTING

#include "Felt/SingleTrackedGrid.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_LazySingleTrackedGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		// ==== Setup ====
		LazySingleTrackedGrid<FLOAT, 3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3);
		const UINT NULL_IDX = LazySingleTrackedGrid<FLOAT, 3, 3>::MultiLookup::NULL_IDX;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), 3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3);
		BOOST_CHECK_EQUAL(grid.lookup().is_active(), false);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 0);
		BOOST_CHECK_EQUAL(grid.lookup().background(), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.lookup().get(Vec3i(1,1,1)), NULL_IDX);
	}

	struct Fixture {
		const UINT NULL_IDX = LazySingleTrackedGrid<FLOAT, 3, 3>::MultiLookup::NULL_IDX;
		LazySingleTrackedGrid<FLOAT, 3, 3> grid;
		Fixture()
			: grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3)
		{}
	};

	BOOST_FIXTURE_TEST_CASE(activate_should_activate_lookup, Fixture)
	{
		/// [LazySingleTrackedGrid activate]
		// ==== Action ====
		grid.activate();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), true);
		BOOST_CHECK_EQUAL(grid.data().size(), 3*3*3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3);
		BOOST_CHECK_EQUAL(grid.lookup().is_active(), true);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 3*3*3);
		BOOST_CHECK_EQUAL(grid.lookup().get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySingleTrackedGrid activate]
	}

	BOOST_FIXTURE_TEST_CASE(deactivate_should_deactivate_lookup, Fixture)
	{
		/// [LazySingleTrackedGrid deactivate]
		// ==== Action ====
		grid.activate();
		grid.deactivate();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), 3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3);
		BOOST_CHECK_EQUAL(grid.lookup().is_active(), false);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 0);
		BOOST_CHECK_EQUAL(grid.lookup().background(), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.lookup().get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySingleTrackedGrid deactivate]
	}
BOOST_AUTO_TEST_SUITE_END()

