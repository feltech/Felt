#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "PartitionedGrid.hpp"

using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_PartitionedGrid)

	BOOST_AUTO_TEST_CASE(init)
	{
		{
			typedef PartitionedGrid<FLOAT, 3, 3> Grid_t;
			// Default constructor.
			Grid_t grid;
			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(0,0,0));
		}

		{
			typedef PartitionedGrid<FLOAT, 3, 3> Grid_t;
			// Invalid dimensions (not divisible by partition size).
			BOOST_CHECK_THROW(
				Grid_t grid(Vec3u(9,8,9), Vec3i(-5, -5,-5)),
				std::invalid_argument
			);
		}

		{
			typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;
			typedef Grid_t::PartGrid_t PartGrid_t;

			Grid_t grid(Vec3u(4,4,4), Vec3i(-2, -2,-2));
			PartGrid_t& parts = grid.parts();

			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i(-1, -1, -1)), (void*)&parts.data()[0]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i(-1, -1,  0)), (void*)&parts.data()[1]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i(-1,  0, -1)), (void*)&parts.data()[2]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i(-1,  0,  0)), (void*)&parts.data()[3]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i( 0, -1, -1)), (void*)&parts.data()[4]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i( 0, -1,  0)), (void*)&parts.data()[5]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i( 0,  0, -1)), (void*)&parts.data()[6]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parts(Vec3i( 0,  0,  0)), (void*)&parts.data()[7]
			);

			BOOST_CHECK_EQUAL(
				parts(Vec3i(-1, -1, -1)).offset(), Vec3i(-2, -2, -2)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i(-1, -1,  0)).offset(), Vec3i(-2, -2,  0)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i(-1,  0, -1)).offset(), Vec3i(-2,  0, -2)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i(-1,  0,  0)).offset(), Vec3i(-2,  0,  0)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 0, -1, -1)).offset(), Vec3i( 0, -2, -2)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 0, -1,  0)).offset(), Vec3i( 0, -2,  0)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 0,  0, -1)).offset(), Vec3i( 0,  0, -2)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 0,  0,  0)).offset(), Vec3i( 0,  0,  0)
			);
		}

		{
			typedef PartitionedGrid<FLOAT, 3, 3> Grid_t;
			typedef Grid_t::PartGrid_t PartGrid_t;

			Grid_t grid(Vec3u(9,9,9), Vec3i(-4, -4, -4));
			PartGrid_t& parts = grid.parts();

			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(9,9,9));

			BOOST_CHECK_EQUAL(grid.data().size(), 0);
			BOOST_CHECK_EQUAL(parts.data().size(), 27);

			const Vec3u part_dims(3,3,3);
			BOOST_CHECK_EQUAL(parts(Vec3i(-1, -1, -1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parts(Vec3i( 1,  1,  1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parts(Vec3i( 0,  0,  0)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parts(Vec3i(-1,  0,  1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(
				parts(Vec3i(-1, -1, -1)).offset(), Vec3i(-4, -4, -4)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 0,  0,  0)).offset(), Vec3i(-1, -1, -1)
			);
			BOOST_CHECK_EQUAL(
				parts(Vec3i( 1,  1,  1)).offset(), Vec3i( 2,  2,  2)
			);
		}

		{
			typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;

			Grid_t grid(Vec3u(8,8,8), Vec3i(-3,-3,-3));

			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(8,8,8));

			BOOST_CHECK_EQUAL(grid.parts().data().size(), 64);

			const Vec3u part_dims(2,2,2);
			BOOST_CHECK_EQUAL(grid.parts()(Vec3i(-1, -1, -1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(grid.parts()(Vec3i( 0,  0,  0)).dims(), part_dims);
			BOOST_CHECK_EQUAL(grid.parts()(Vec3i( 1,  1,  1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(grid.parts()(Vec3i( 2,  2,  2)).dims(), part_dims);

			BOOST_CHECK_EQUAL(
				grid.parts()(Vec3i(-1, -1, -1)).offset(), Vec3i(-3, -3, -3)
			);
			BOOST_CHECK_EQUAL(
				grid.parts()(Vec3i( 0,  0,  0)).offset(), Vec3i(-1, -1, -1)
			);
			BOOST_CHECK_EQUAL(
				grid.parts()(Vec3i( 1,  1,  1)).offset(), Vec3i( 1,  1,  1)
			);
			BOOST_CHECK_EQUAL(
				grid.parts()(Vec3i( 2,  2,  2)).offset(), Vec3i( 3,  3,  3)
			);
		}
	}


	BOOST_AUTO_TEST_CASE(get_and_set)
	{
		typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;
		typedef Grid_t::PartGrid_t PartGrid_t;

		Grid_t grid(Vec3u(4,4,4), Vec3i(-2, -2,-2));
		PartGrid_t& parts = grid.parts();

		grid.fill(-1.0f);

		BOOST_CHECK_EQUAL(grid(Vec3i(-2,-2,1)), -1.0f);

		for (INT x = -2; x <= 1; x++)
			for (INT y = -2; y <= 1; y++)
				for (INT z = -2; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), -1.0f);
				}
	}

BOOST_AUTO_TEST_SUITE_END()
