#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include <Felt/PartitionedGrid.hpp>

using namespace felt;

/**
 * @addtogroup Tests
 * @addtogroup PartitionedGridTests Spatially Partitioned Grid Tests
 *
 * Test the various spatially partitioned wrappers around felt::Grid based classes.
 *
 * @{
 */
BOOST_AUTO_TEST_SUITE(test_PartitionedGrid)
	/**
	 * @name PartitionedGrid
	 * @ref felt::PartitionedGrid
	 */

	/**
	 * Basic initialisation.
	 */
	BOOST_AUTO_TEST_CASE(init_simple)
	{
		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid;

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(grid.size(), Vec3u(0,0,0));
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid(Vec3u(4,4,4), Vec3i(-2, -2,-2), Vec3u(2, 2, 2));
			PartitionedGrid<FLOAT, 3>::BranchGrid& parent = grid.branch();

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i(-1, -1, -1)), (void*)&parent.data()[0]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i(-1, -1,  0)), (void*)&parent.data()[1]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i(-1,  0, -1)), (void*)&parent.data()[2]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i(-1,  0,  0)), (void*)&parent.data()[3]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i( 0, -1, -1)), (void*)&parent.data()[4]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i( 0, -1,  0)), (void*)&parent.data()[5]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i( 0,  0, -1)), (void*)&parent.data()[6]
			);
			BOOST_CHECK_EQUAL(
				(void*)&parent(Vec3i( 0,  0,  0)), (void*)&parent.data()[7]
			);

			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1, -1, -1)).offset(), Vec3i(-2, -2, -2)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1, -1,  0)).offset(), Vec3i(-2, -2,  0)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1,  0, -1)).offset(), Vec3i(-2,  0, -2)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1,  0,  0)).offset(), Vec3i(-2,  0,  0)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0, -1, -1)).offset(), Vec3i( 0, -2, -2)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0, -1,  0)).offset(), Vec3i( 0, -2,  0)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0,  0, -1)).offset(), Vec3i( 0,  0, -2)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0,  0,  0)).offset(), Vec3i( 0,  0,  0)
			);
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3> grid(Vec3u(9,9,9), Vec3i(-4, -4, -4), Vec3u(3, 3, 3));
			PartitionedGrid<FLOAT, 3>::BranchGrid& parent = grid.branch();

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(grid.size(), Vec3u(9,9,9));

			BOOST_CHECK_EQUAL((UINT)parent.data().size(), 27u);

			const Vec3u part_size(3,3,3);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1, -1, -1)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i( 1,  1,  1)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i( 0,  0,  0)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1,  0,  1)).size(), part_size);
			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1, -1, -1)).offset(), Vec3i(-4, -4, -4)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0,  0,  0)).offset(), Vec3i(-1, -1, -1)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 1,  1,  1)).offset(), Vec3i( 2,  2,  2)
			);
		}

		{
			// ==== Setup/action ====
			PartitionedGrid<FLOAT, 3>  grid(Vec3u(8,8,8), Vec3i(-3,-3,-3), Vec3u(2, 2, 2));
			const PartitionedGrid<FLOAT, 3> ::BranchGrid& parent = grid.branch();

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(grid.size(), Vec3u(8,8,8));

			BOOST_CHECK_EQUAL((UINT)parent.data().size(), 64);

			const Vec3u part_size(2,2,2);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1, -1, -1)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i( 0,  0,  0)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i( 1,  1,  1)).size(), part_size);
			BOOST_CHECK_EQUAL(parent(Vec3i( 2,  2,  2)).size(), part_size);

			BOOST_CHECK_EQUAL(
				parent(Vec3i(-1, -1, -1)).offset(), Vec3i(-3, -3, -3)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 0,  0,  0)).offset(), Vec3i(-1, -1, -1)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 1,  1,  1)).offset(), Vec3i( 1,  1,  1)
			);
			BOOST_CHECK_EQUAL(
				parent(Vec3i( 2,  2,  2)).offset(), Vec3i( 3,  3,  3)
			);
		}
	}

	/**
	 * Simple get and set values.
	 */
	BOOST_AUTO_TEST_CASE(get_and_set_simple)
	{
		// ==== Setup ====
		PartitionedGrid<FLOAT, 3> grid(Vec3u(4,4,4), Vec3i(-2, -2,-2));
		PartitionedGrid<FLOAT, 3>::BranchGrid& parent = grid.branch();

		// ==== Action ====
		grid.fill(-1.0f);

		// ==== Confirm ====
		for (INT x = -2; x <= 1; x++)
			for (INT y = -2; y <= 1; y++)
				for (INT z = -2; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), -1.0f);
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
		BOOST_CHECK_EQUAL(grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 2.0f);
		BOOST_CHECK_EQUAL(grid(pos3), 3.0f);
		BOOST_CHECK_EQUAL(grid(pos4), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos5), 5.0f);
		BOOST_CHECK_EQUAL((FLOAT)grad(0), 0.0f);
		BOOST_CHECK_EQUAL((FLOAT)grad(1), 3.5f);
		BOOST_CHECK_EQUAL((FLOAT)grad(2), 0.0f);
	}

	/**
	 * @name LookupPartitionedGrid
	 * @ref felt::LookupPartitionedGrid
	 */

	/**
	 * Simple lookup get and set values.
	 */
	BOOST_AUTO_TEST_CASE(partitioned_lookup)
	{
		typedef LookupPartitionedGrid<3, 3> Grid_t;
		typedef LookupPartitionedGrid<3, 3>::BranchGrid BranchGrid;
		typedef LookupPartitionedGrid<3, 3>::BranchGrid::Lookup LookupGrid;

		Grid_t grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		BranchGrid& branch = grid.branch();
		LookupGrid& lookup = branch.lookup();


		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), Grid_t::NULL_IDX_TUPLE);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), Grid_t::NULL_IDX_TUPLE);
				}


		const Vec3i pos1(1, -4, -1);
		const Vec3i pos2(2, -3, -2);
		const Vec3i pos3(3, -2, -3);
		const Vec3i pos4(4, -1, -4);
		const Vec3i part1(0, -1,  0);
		const Vec3i part2_3(1, -1, -1);
		const Vec3i part4(1,  0, -1);

		grid.add(pos1);
		grid.add(pos2, 0);
		grid.add(pos3, 0);
		grid.add(pos4, 2);

		BOOST_CHECK_EQUAL((UINT)grid(pos1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2)(0), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos4)(2), 0);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)branch(part4)(pos4)(2), 0);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(0)[0], part1);
		BOOST_CHECK_EQUAL(branch.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(branch.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : branch.list(i))
				for (const Vec3i pos : branch(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(2);

		BOOST_CHECK_EQUAL(branch.list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos4), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL(lookup(part4), Grid_t::NULL_IDX_TUPLE);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid(pos2), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL(lookup(part1), Grid_t::NULL_IDX_TUPLE);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(branch.list(i).size(), 0);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), Grid_t::NULL_IDX_TUPLE);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), Grid_t::NULL_IDX_TUPLE);
					for (UINT i = 0; i < 3; i++)
						BOOST_CHECK_EQUAL(branch(pos).list(i).size(), 0);
				}
	}

	/**
	 * @name SharedLookupPartitionedGrid
	 * @ref felt::SharedLookupPartitionedGrid
	 */

	BOOST_AUTO_TEST_CASE(partitioned_shared_lookup)
	{
		typedef SharedLookupPartitionedGrid<3, 3> GridType;
		typedef GridType::BranchGrid BranchGrid;
		typedef BranchGrid::Lookup LookupGrid;
		const Vec3u& BRANCH_NULL_IDX = LookupGrid::NULL_IDX_TUPLE;
		const UINT& CHILD_NULL_IDX = GridType::NULL_IDX;

		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		BranchGrid& branch = grid.branch();
		LookupGrid& lookup = branch.lookup();


		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), CHILD_NULL_IDX);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), BRANCH_NULL_IDX);
				}


		const Vec3i pos1(1, -4, -1);
		const Vec3i pos2(2, -3, -2);
		const Vec3i pos3(3, -2, -3);
		const Vec3i pos4(4, -1, -4);
		const Vec3i part1(0, -1,  0);
		const Vec3i part2_3(1, -1, -1);
		const Vec3i part4(1,  0, -1);

		grid.add(pos1);
		grid.add(pos2, 0);
		grid.add(pos3, 0);
		grid.add(pos4, 2);

		BOOST_CHECK_EQUAL((UINT)grid(pos1), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), 0);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 1);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), 0);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)branch(part4)(pos4), 0);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(0)[0], part1);
		BOOST_CHECK_EQUAL(branch.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(branch.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : branch.list(i))
				for (const Vec3i pos : branch(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(2);

		BOOST_CHECK_EQUAL(branch.list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos4), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part4), BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid(pos2), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part1), BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(branch.list(i).size(), 0);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), CHILD_NULL_IDX);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), BRANCH_NULL_IDX);
					for (UINT i = 0; i < 3; i++)
						BOOST_CHECK_EQUAL(branch(pos).list(i).size(), 0);
				}
	}

	/**
	 * @name SharedTrackedPartitionedGrid
	 * @ref felt::SharedTrackedPartitionedGrid
	 */

	BOOST_AUTO_TEST_CASE(partitioned_shared_tracked)
	{
		typedef SharedTrackedPartitionedGrid<FLOAT, 3, 3> GridType;
		typedef GridType::BranchGrid BranchGrid;
		typedef BranchGrid::Lookup LookupGrid;
		const Vec3u& BRANCH_NULL_IDX = LookupGrid::NULL_IDX_TUPLE;
		const UINT& CHILD_NULL_IDX = GridType::Child::Lookup::NULL_IDX;

		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		BranchGrid& branch = grid.branch();
		LookupGrid& lookup = branch.lookup();

		grid.fill(-1.0f);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					const Vec3i pos(x, y, z);
					const Vec3i pos_child = grid.pos_child(pos);
					BOOST_CHECK_EQUAL(grid(pos), -1.0f);
					BOOST_CHECK_EQUAL(
						grid.child(pos_child).lookup()(pos), CHILD_NULL_IDX
					);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), BRANCH_NULL_IDX);
				}


		const Vec3i pos1(1, -4, -1);
		const Vec3i pos2(2, -3, -2);
		const Vec3i pos3(3, -2, -3);
		const Vec3i pos4(4, -1, -4);
		const Vec3i part1(0, -1,  0);
		const Vec3i part2_3(1, -1, -1);
		const Vec3i part4(1,  0, -1);

		grid.add(pos1, 1.0f, 0);
		grid.add(pos2, 2.0f, 0);
		grid.add(pos3, 3.0f, 0);
		grid.add(pos4, 4.0f, 2);

		BOOST_CHECK_EQUAL((UINT)grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL((UINT)grid(pos2), 2.0f);
		BOOST_CHECK_EQUAL((UINT)grid(pos3), 3.0f);
		BOOST_CHECK_EQUAL((UINT)grid(pos4), 4.0f);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)branch(part4)(pos4), 4.0f);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(0)[0], part1);
		BOOST_CHECK_EQUAL(branch.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(branch.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : branch.list(i))
				for (const Vec3i pos : branch(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(-2.0f, 2);

		BOOST_CHECK_EQUAL(grid(pos4), -2.0f);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part4).lookup()(pos4), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part4), BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(grid(pos2), 2.0f);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part2_3).lookup()(pos2), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(branch(part1).lookup()(pos1), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part1), BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(branch.list(i).size(), 0);

		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), BRANCH_NULL_IDX);
					for (UINT i = 0; i < 3; i++)
						BOOST_CHECK_EQUAL(branch(pos).list(i).size(), 0);
				}
	}

	/**
	 * @name PartitionedArray
	 * @ref felt::PartitionedArray
	 */

	BOOST_AUTO_TEST_CASE(partitioned_array)
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

		BOOST_CHECK_EQUAL(grid.branch().list().size(), 3);
		BOOST_CHECK_EQUAL(grid.child(part1).size(), 1);
		BOOST_CHECK_EQUAL(grid.child(part2_3).size(), 2);
		BOOST_CHECK_EQUAL(grid.child(part4).size(), 1);

		BOOST_CHECK_EQUAL(grid.child(part1)[0], 1.0f);
		BOOST_CHECK_EQUAL(grid.child(part2_3)[0], 2.0f);
		BOOST_CHECK_EQUAL(grid.child(part2_3)[1], 3.0f);
		BOOST_CHECK_EQUAL(grid.child(part4)[0], 4.0f);

		grid.reset();

		BOOST_CHECK_EQUAL(grid.branch().list().size(), 0);
		BOOST_CHECK_EQUAL(grid.child(part1).size(), 0);
		BOOST_CHECK_EQUAL(grid.child(part2_3).size(), 0);
		BOOST_CHECK_EQUAL(grid.child(part4).size(), 0);
	}

BOOST_AUTO_TEST_SUITE_END()

/** @} */ // End group Tests.

/**
 *  @class felt::PartitionedGrid
 *  @test see @ref Tests
 *  @class felt::LookupPartitionedGrid
 *  @test see @ref Tests
 *  @class felt::SharedLookupPartitionedGrid
 *  @test see @ref Tests
 *  @class felt::SharedTrackedPartitionedGrid
 *  @test see @ref Tests
 *  @class felt::PartitionedArray
 *  @test see @ref Tests
 */
