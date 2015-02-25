#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "ArrayMappedGrid.hpp"

using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_MappedGrid)
	/*
	 * ArrayMappedGrid.
	 */
	BOOST_AUTO_TEST_CASE(main_ArrayMappedGrid)
	{
		typedef ArrayMappedGrid<FLOAT, 3> Grid_t;

		Grid_t grid(Vec3u(5,5,5), Vec3i(-2,-2,-2));

		Vec3i pos1(-1, 0, 1);
		Vec3i pos2(-1, 1, 0);
		grid.add(pos1, 3.0f);
		grid.add(pos2);
		grid(pos2) = 5.0f;

		BOOST_CHECK_EQUAL(grid(pos1), 3.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 5.0f);
		BOOST_CHECK_EQUAL(grid.list().size(), 2);

		for (Vec3i pos : grid.list())
			grid(pos) = 4.0f;

		BOOST_CHECK_EQUAL(grid(pos1), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 4.0f);

		grid.reset(-1.0f);
		BOOST_CHECK_EQUAL(grid(pos1), -1.0f);
		BOOST_CHECK_EQUAL(grid(pos2), -1.0f);
		BOOST_CHECK_EQUAL(grid.list().size(), 0);
	}


	BOOST_AUTO_TEST_CASE(main_GridMappedArray)
	{
		typedef GridMappedArray<FLOAT, 3> Array_t;
		Array_t arr(Vec3u(10,10,10), Vec3i(-5, -5, -5));

		Vec3i pos1(-1, 0, 1);
		Vec3i pos2(-1, 1, 0);

		BOOST_CHECK_EQUAL(arr.idx(pos1), Array_t::NULL_IDX);
		BOOST_CHECK_EQUAL(arr.idx(pos2), Array_t::NULL_IDX);

		arr.add(pos1, 3.0f);
		arr.add(pos2, 5.0f);

		BOOST_CHECK_EQUAL(arr.size(), 2);
		BOOST_CHECK_EQUAL(arr.idx(pos1), 0);
		BOOST_CHECK_EQUAL(arr.idx(pos2), 1);
	}
BOOST_AUTO_TEST_SUITE_END()
