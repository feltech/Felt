#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "MappedGrid.hpp"

using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_MappedGrid)
	/*
	 * ArrayMappedGrid.
	 */
	BOOST_AUTO_TEST_CASE(test_ArrayMappedGrid)
	{
		typedef ArrayMappedGrid<FLOAT, 3> Grid_t;

		Grid_t grid(Vec3u(5,5,5), Vec3i(-2,-2,-2));

		Vec3i pos1(0, 0, 1);
		Vec3i pos2(1, 1, 0);
		Vec3i pos3(2, 0, -1);
		UINT idx1 = grid.add(pos1, 3.0f);
		UINT idx2 = grid.add(pos2, -1.0f);
		grid(pos2) = 5.0f;
		UINT idx3 = grid.add(pos3, 7.0f);


		BOOST_CHECK_EQUAL(grid.list().size(), 3);
		BOOST_CHECK_EQUAL(idx1, 0);
		BOOST_CHECK_EQUAL(idx2, 1);
		BOOST_CHECK_EQUAL(idx3, 2);
		BOOST_CHECK_EQUAL(grid(pos1), 3.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 5.0f);
		BOOST_CHECK_EQUAL(grid(pos3), 7.0f);
		BOOST_CHECK_EQUAL(grid.list()[0], pos1);
		BOOST_CHECK_EQUAL(grid.list()[1], pos2);
		BOOST_CHECK_EQUAL(grid.list()[2], pos3);

		for (Vec3i pos : grid.list())
			grid(pos) = 4.0f;

		BOOST_CHECK_EQUAL(grid(pos1), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos3), 4.0f);

		grid.remove(1);
		BOOST_CHECK_EQUAL(grid.list().size(), 2);
		BOOST_CHECK_EQUAL(grid.list()[0], pos1);
		BOOST_CHECK_EQUAL(grid.list()[1], pos3);

		grid.reset(-1.0f);
		BOOST_CHECK_EQUAL(grid.list().size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), -1.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos3), -1.0f);

		grid.add(pos1, 3.0f);
		grid.add(pos2, 5.0f);
		grid.list().clear();

		BOOST_CHECK_EQUAL(grid.list().size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), 3.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 5.0f);
	}


	BOOST_AUTO_TEST_CASE(test_PosArrayMappedGrid)
	{
		typedef PosArrayMappedGrid<3> Grid_t;
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
		const UINT& idx1 = grid.add(pos1);
		const UINT& idx2 = grid.add(pos2);
		const UINT& idx3 = grid.add(pos3);
		const UINT& idx4 = grid.add(pos4);

		BOOST_CHECK_EQUAL(idx1, 0);
		BOOST_CHECK_EQUAL(idx2, 1);
		BOOST_CHECK_EQUAL(idx3, 2);
		BOOST_CHECK_EQUAL(idx4, 3);


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
		UINT idx5 = grid.add(pos2);

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
		BOOST_CHECK_EQUAL(idx5, 1);

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


	BOOST_AUTO_TEST_CASE(test_multi_PosArrayMappedGrid)
	{
		typedef PosArrayMappedGrid<3, 3> Grid_t;
		Grid_t grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		BOOST_CHECK_EQUAL(Grid_t::num_lists(), 3);


		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		// Add the positions to the array and set index lookup values.
		const UINT& idx1 = grid.add(pos1, 0);
		const UINT& idx2 = grid.add(pos2, 1);
		const UINT& idx3 = grid.add(pos3, 1);
		const UINT& idx4 = grid.add(pos4, 2);

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
		BOOST_CHECK_EQUAL(idx1, 0);
		BOOST_CHECK_EQUAL(idx2, 0);
		BOOST_CHECK_EQUAL(idx3, 1);
		BOOST_CHECK_EQUAL(idx4, 0);

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


	BOOST_AUTO_TEST_CASE(test_multi_shared_PosArrayMappedGrid)
	{
		typedef PosArrayMappedSharedGrid<3, 3> Grid_t;
		Grid_t grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		// Add the positions to the array and set index lookup values.
		const UINT& idx1 = grid.add(pos1, 0);
		const UINT& idx2 = grid.add(pos2, 1);
		const UINT& idx3 = grid.add(pos3, 1);
		const UINT& idx4 = grid.add(pos4, 2);

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
		BOOST_CHECK_EQUAL(idx1, 0);
		BOOST_CHECK_EQUAL(idx2, 0);
		BOOST_CHECK_EQUAL(idx3, 1);
		BOOST_CHECK_EQUAL(idx4, 0);

		grid.remove(pos2, 1);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(0)[0], pos1);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL(grid.list(2)[0], pos4);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), Grid_t::NULL_IDX);
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
		BOOST_CHECK_EQUAL((UINT)grid(pos2), Grid_t::NULL_IDX);
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
		BOOST_CHECK_EQUAL((UINT)grid(pos1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos6), 0);

		grid.reset(2);

		BOOST_CHECK_EQUAL(grid.list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid.list(1)[0], pos3);
		BOOST_CHECK_EQUAL((UINT)grid(pos1), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), Grid_t::NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)grid(pos5), Grid_t::NULL_IDX);

	}
BOOST_AUTO_TEST_SUITE_END()
