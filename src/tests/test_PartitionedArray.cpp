#include <boost/test/unit_test.hpp>
#include "Felt/PartitionedArray.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_PartitionedArray)
	BOOST_AUTO_TEST_CASE(initialise_and_populate)
	{
		typedef PartitionedArray<FLOAT, 3> ArrayGrid;

		ArrayGrid grid(Vec3u(9,9,9), Vec3i(-4, -4, -4), Vec3u(3, 3, 3));

		const Vec3i pos1(1, -4, -1);
		const Vec3i pos2(2, -3, -2);
		const Vec3i pos3(3, -2, -3);
		const Vec3i pos4(4, -1, -4);
		const Vec3i part1(0, -1,  0);
		const Vec3i part2_3(1, -1, -1);
		const Vec3i part4(1,  0, -1);

		grid.add(pos1, 1.0f);
		grid.add(pos2, 2.0f);
		grid.add(pos3, 3.0f);
		grid.add(pos4, 4.0f);

		BOOST_CHECK_EQUAL(grid.children().list().size(), 3);
		BOOST_CHECK_EQUAL(grid.children().get(part1).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(part2_3).size(), 2);
		BOOST_CHECK_EQUAL(grid.children().get(part4).size(), 1);

		BOOST_CHECK_EQUAL(grid.children().get(part1)[0], 1.0f);
		BOOST_CHECK_EQUAL(grid.children().get(part2_3)[0], 2.0f);
		BOOST_CHECK_EQUAL(grid.children().get(part2_3)[1], 3.0f);
		BOOST_CHECK_EQUAL(grid.children().get(part4)[0], 4.0f);

		grid.reset();

		BOOST_CHECK_EQUAL(grid.children().list().size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(part1).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(part2_3).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(part4).size(), 0);
	}
BOOST_AUTO_TEST_SUITE_END()
