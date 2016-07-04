#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include <Felt/PartitionedGrid.hpp>

using namespace felt;

SCENARIO("PartitionedGrid")
{
	/**
	 * Basic initialisation.
	 */
	WHEN("init_simple")
	{
		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid;

			// ==== Confirm ====
			CHECK(grid.size() == Vec3u(0,0,0));
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid(Vec3u(4,4,4), Vec3i(-2, -2,-2), 0, Vec3u(2, 2, 2));
			PartitionedGrid<FLOAT, 3>::ChildrenGrid& parent = grid.children();

			// ==== Confirm ====
			CHECK(
				(void*)&parent(Vec3i(-1, -1, -1)) == (void*)&parent.data()[0]
			);
			CHECK(
				(void*)&parent(Vec3i(-1, -1,  0)) == (void*)&parent.data()[1]
			);
			CHECK(
				(void*)&parent(Vec3i(-1,  0, -1)) == (void*)&parent.data()[2]
			);
			CHECK(
				(void*)&parent(Vec3i(-1,  0,  0)) == (void*)&parent.data()[3]
			);
			CHECK(
				(void*)&parent(Vec3i( 0, -1, -1)) == (void*)&parent.data()[4]
			);
			CHECK(
				(void*)&parent(Vec3i( 0, -1,  0)) == (void*)&parent.data()[5]
			);
			CHECK(
				(void*)&parent(Vec3i( 0,  0, -1)) == (void*)&parent.data()[6]
			);
			CHECK(
				(void*)&parent(Vec3i( 0,  0,  0)) == (void*)&parent.data()[7]
			);

			CHECK(
				parent(Vec3i(-1, -1, -1)).offset() == Vec3i(-2, -2, -2)
			);
			CHECK(
				parent(Vec3i(-1, -1,  0)).offset() == Vec3i(-2, -2,  0)
			);
			CHECK(
				parent(Vec3i(-1,  0, -1)).offset() == Vec3i(-2,  0, -2)
			);
			CHECK(
				parent(Vec3i(-1,  0,  0)).offset() == Vec3i(-2,  0,  0)
			);
			CHECK(
				parent(Vec3i( 0, -1, -1)).offset() == Vec3i( 0, -2, -2)
			);
			CHECK(
				parent(Vec3i( 0, -1,  0)).offset() == Vec3i( 0, -2,  0)
			);
			CHECK(
				parent(Vec3i( 0,  0, -1)).offset() == Vec3i( 0,  0, -2)
			);
			CHECK(
				parent(Vec3i( 0,  0,  0)).offset() == Vec3i( 0,  0,  0)
			);
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid(Vec3u(9,9,9), Vec3i(-4, -4, -4), 0, Vec3u(3, 3, 3));
			PartitionedGrid<FLOAT, 3>::ChildrenGrid& parent = grid.children();

			// ==== Confirm ====
			CHECK(grid.size() == Vec3u(9,9,9));

			CHECK((UINT)parent.data().size() == 27u);

			const Vec3u part_size(3,3,3);
			CHECK(parent(Vec3i(-1, -1, -1)).size() == part_size);
			CHECK(parent(Vec3i( 1,  1,  1)).size() == part_size);
			CHECK(parent(Vec3i( 0,  0,  0)).size() == part_size);
			CHECK(parent(Vec3i(-1,  0,  1)).size() == part_size);
			CHECK(
				parent(Vec3i(-1, -1, -1)).offset() == Vec3i(-4, -4, -4)
			);
			CHECK(
				parent(Vec3i( 0,  0,  0)).offset() == Vec3i(-1, -1, -1)
			);
			CHECK(
				parent(Vec3i( 1,  1,  1)).offset() == Vec3i( 2,  2,  2)
			);
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3>  grid(Vec3u(8,8,8), Vec3i(-3,-3,-3), 0, Vec3u(2, 2, 2));
			const PartitionedGrid<FLOAT, 3> ::ChildrenGrid& parent = grid.children();

			// ==== Confirm ====
			CHECK(grid.size() == Vec3u(8,8,8));

			CHECK((UINT)parent.data().size() == 64);

			const Vec3u part_size(2,2,2);
			CHECK(parent(Vec3i(-1, -1, -1)).size() == part_size);
			CHECK(parent(Vec3i( 0,  0,  0)).size() == part_size);
			CHECK(parent(Vec3i( 1,  1,  1)).size() == part_size);
			CHECK(parent(Vec3i( 2,  2,  2)).size() == part_size);

			CHECK(
				parent(Vec3i(-1, -1, -1)).offset() == Vec3i(-3, -3, -3)
			);
			CHECK(
				parent(Vec3i( 0,  0,  0)).offset() == Vec3i(-1, -1, -1)
			);
			CHECK(
				parent(Vec3i( 1,  1,  1)).offset() == Vec3i( 1,  1,  1)
			);
			CHECK(
				parent(Vec3i( 2,  2,  2)).offset() == Vec3i( 3,  3,  3)
			);
		}
	}

	/**
	 * Simple get and set values.
	 */
	WHEN("get_and_set_simple")
	{
		// ==== Setup ====
		PartitionedGrid<FLOAT, 3> grid(Vec3u(4,4,4), Vec3i(-2, -2,-2), 0, Vec3u(4, 4, 4));
		PartitionedGrid<FLOAT, 3>::ChildrenGrid& parent = grid.children();

		// ==== Action ====
		grid.fill(-1.0f);

		// ==== Confirm ====
		for (INT x = -2; x <= 1; x++)
			for (INT y = -2; y <= 1; y++)
				for (INT z = -2; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					CHECK(grid(pos) == -1.0f);
				}

		// ==== Setup ====
		const Vec3i pos1(-2, -2, -2);
		const Vec3i pos2(-1, -1, -1);
		const Vec3i pos3( 0,  0,  0);
		const Vec3i pos4( 1,  1,  1);
		const Vec3i pos5(-2, -1,  1);
		const Vec3i pos6( 0,  1,  0);

		// ==== Action ====
		grid(pos1) = 1.0f;
		grid(pos2) = 2.0f;
		grid(pos3) = 3.0f;
		grid(pos4) = 4.0f;
		grid(pos5) = 5.0f;
		grid(pos6) = 6.0f;

		const Vec3f grad = grid.grad(pos3);

		// ==== Confirm ====
		CHECK(grid(pos1) == 1.0f);
		CHECK(grid(pos2) == 2.0f);
		CHECK(grid(pos3) == 3.0f);
		CHECK(grid(pos4) == 4.0f);
		CHECK(grid(pos5) == 5.0f);
		CHECK((FLOAT)grad(0) == 0.0f);
		CHECK((FLOAT)grad(1) == 3.5f);
		CHECK((FLOAT)grad(2) == 0.0f);
	}
}
