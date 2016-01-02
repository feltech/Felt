#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include <Felt/LazyPartitionedGrid.hpp>

using namespace felt;

/**
 * @addtogroup Tests
 * @addtogroup PartitionedGridTests Spatially Partitioned Grid Tests
 *
 * @{
 */
BOOST_AUTO_TEST_SUITE(test_LazyPartitionedGrid)
/**
 * @name LazyPartitionedGrid
 * @ref felt::LazyPartitionedGrid
 */

BOOST_AUTO_TEST_CASE(basic_usage)
{
	//! [Lazy grid basic usage]
	// ==== Setup ====
	using LazyGrid = LazyPartitionedGrid<FLOAT, 3>;
	LazyGrid grid(Vec3u(9,9,9), Vec3i(-4, -4, -4), -3.0f, Vec3u(3, 3, 3));

	const Vec3i pos1(1, -4, -1);
	const Vec3i pos2(2, -3, -2);
	const Vec3i pos3(4, -1, -4);
	const Vec3i part1(0, -1,  0);
	const Vec3i part2(1, -1, -1);
	const Vec3i part3(1,  0, -1);

	// ==== Action ====
	// Initialise child grid
	grid.add_child(part2);
	// Set value inside child grid
	grid.get(pos2) = 2.0f;

	// Add then remove child.
	grid.add_child(part3);
	grid.get(pos3) = 5.0f;
	grid.remove_child(part3);

	// ==== Confirm ====
	BOOST_CHECK_EQUAL(grid.get(pos1), -3.0f);
	BOOST_CHECK_EQUAL(grid.child(part1).data().size(), 0);
	BOOST_CHECK_EQUAL(grid.get(pos2), 2.0f);
	BOOST_CHECK_EQUAL(grid.child(part2).data().size(), 3*3*3);
	BOOST_CHECK_EQUAL(grid.get(pos3), -3.0f);
	BOOST_CHECK_EQUAL(grid.child(part3).data().size(), 0);
	//! [Lazy grid basic usage]
}

BOOST_AUTO_TEST_SUITE_END()

/** @} */ // End group Tests.

/**
 *  @class felt::LazyPartitionedGrid
 *  @test see @ref Tests
 */
