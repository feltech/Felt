#include <boost/test/unit_test.hpp>

#define _TESTING

#include "Felt/TrackedGrid.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_TrackedGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		// ==== Setup/Action ====
		TrackedGrid<FLOAT, 3, 3> grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), 0);

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.data().size(), 9*9*9);
		BOOST_CHECK_EQUAL(grid.lookup().data().size(), 9*9*9);

		for (const FLOAT val : grid.data())
			BOOST_CHECK_EQUAL(val, 0);

		for (const Vec3u& val : grid.lookup().data())
			BOOST_CHECK_EQUAL(val, Vec3u::Constant(grid.lookup().NULL_IDX));

	}
BOOST_AUTO_TEST_SUITE_END()
