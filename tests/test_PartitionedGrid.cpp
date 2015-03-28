#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "PartitionedGrid.hpp"

using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_PartitionedGrid)

	BOOST_AUTO_TEST_CASE(init_simple)
	{
		{
			typedef PartitionedGrid<FLOAT, 3, 3> Grid_t;
			// Default constructor.
			Grid_t grid;
			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(0,0,0));
		}

		{
			typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;
			typedef Grid_t::BranchGrid BranchGrid;

			Grid_t grid(Vec3u(4,4,4), Vec3i(-2, -2,-2));
			BranchGrid& parent = grid.branch();

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
			typedef PartitionedGrid<FLOAT, 3, 3> Grid_t;
			typedef Grid_t::BranchGrid PartGrid_t;

			Grid_t grid(Vec3u(9,9,9), Vec3i(-4, -4, -4));
			PartGrid_t& parent = grid.branch();

			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(9,9,9));

			BOOST_CHECK_EQUAL(grid.data().size(), 0);
			BOOST_CHECK_EQUAL((UINT)parent.data().size(), 27u);

			const Vec3u part_dims(3,3,3);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1, -1, -1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i( 1,  1,  1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i( 0,  0,  0)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1,  0,  1)).dims(), part_dims);
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
			typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;

			Grid_t grid(Vec3u(8,8,8), Vec3i(-3,-3,-3));
			const Grid_t::BranchGrid& parent = grid.branch();

			BOOST_CHECK_EQUAL(grid.dims(), Vec3u(8,8,8));

			BOOST_CHECK_EQUAL((UINT)parent.data().size(), 64);

			const Vec3u part_dims(2,2,2);
			BOOST_CHECK_EQUAL(parent(Vec3i(-1, -1, -1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i( 0,  0,  0)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i( 1,  1,  1)).dims(), part_dims);
			BOOST_CHECK_EQUAL(parent(Vec3i( 2,  2,  2)).dims(), part_dims);

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


	BOOST_AUTO_TEST_CASE(get_and_set_simple)
	{
		typedef PartitionedGrid<FLOAT, 3, 2> Grid_t;
		typedef Grid_t::BranchGrid PartGrid_t;

		Grid_t grid(Vec3u(4,4,4), Vec3i(-2, -2,-2));
		PartGrid_t& parent = grid.branch();

		grid.fill(-1.0f);

		for (INT x = -2; x <= 1; x++)
			for (INT y = -2; y <= 1; y++)
				for (INT z = -2; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), -1.0f);
				}

		const Vec3i pos1(-2, -2, -2);
		const Vec3i pos2(-1, -1, -1);
		const Vec3i pos3( 0,  0,  0);
		const Vec3i pos4( 1,  1,  1);
		const Vec3i pos5(-2, -1,  1);
		const Vec3i pos6( 0,  1,  0);

		grid(pos1) = 1.0f;
		grid(pos2) = 2.0f;
		grid(pos3) = 3.0f;
		grid(pos4) = 4.0f;
		grid(pos5) = 5.0f;
		grid(pos6) = 6.0f;

		const Vec3f grad = grid.grad(pos3);

		BOOST_CHECK_EQUAL(grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL(grid(pos2), 2.0f);
		BOOST_CHECK_EQUAL(grid(pos3), 3.0f);
		BOOST_CHECK_EQUAL(grid(pos4), 4.0f);
		BOOST_CHECK_EQUAL(grid(pos5), 5.0f);
		BOOST_CHECK_EQUAL((FLOAT)grad(0), 0.0f);
		BOOST_CHECK_EQUAL((FLOAT)grad(1), 3.5f);
		BOOST_CHECK_EQUAL((FLOAT)grad(2), 0.0f);
	}


	BOOST_AUTO_TEST_CASE(array_mapped)
	{
		typedef MappedPartitionedGrid<FLOAT, 3, 3, 3> Grid_t;

		Grid_t grid(Vec3u(9,9,9), Vec3i(-4,-4,-4));
		Grid_t::BranchGrid& branch = grid.branch();

		BOOST_CHECK_EQUAL((UINT)branch.data().size(), 27);

		grid.fill(7.0f);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					const Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), 7.0f);
					const Vec3i& pos_child = grid.pos_child(pos);
					BOOST_CHECK_EQUAL(branch(pos_child).data().size(), 27);

					for (UINT l = 0; l < 3; l++)
						BOOST_CHECK_EQUAL(branch(pos_child).list(l).size(), 0);
				}

		const Vec3i pos1( 1, -4, -4);
		grid.add(pos1, 1.0f, 1);

		const Vec3i pos_child( 0, -1, -1);
		BOOST_CHECK_EQUAL(grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL(branch.list(1).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(1)[0], pos_child);
		BOOST_CHECK_EQUAL(branch(pos_child).list(1).size(), 1);
		BOOST_CHECK_EQUAL(branch(pos_child).list(1)[0], pos1);
		BOOST_CHECK_EQUAL(branch(pos_child)(pos1), 1.0f);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);


		const Vec3i pos2( 2, -4, -4);
		const Vec3i pos3( 3, -4, -4);
		const Vec3i pos4( 4, -4, -4);
		const Vec3i pos5(-1, -4, -4);

		const Vec3i part1_3_4( 1, -1, -1);
		const Vec3i part2_5( 0, -1, -1);

		grid.add(pos2, 2.0f, 1);
		grid.add(pos3, 3.0f, 0);
		grid.add(pos4, 4.0f);
		grid.add(pos5, 5.0f, 2);

		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(1).size(), 2);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(0)[0], part1_3_4);
		BOOST_CHECK_EQUAL(branch.list(1)[0], part2_5);
		BOOST_CHECK_EQUAL(branch.list(1)[1], part1_3_4);

		BOOST_CHECK_EQUAL(branch(part1_3_4).list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part2_5).list(1).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1_3_4).list(1).size(), 1);
		BOOST_CHECK_EQUAL(branch(part2_5).list(2).size(), 1);

		grid.reset(-1.0f, 1);

		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch.list(1).size(), 0);
		BOOST_CHECK_EQUAL(branch.list(2).size(), 1);

		BOOST_CHECK_EQUAL(branch(part1_3_4).list(0).size(), 2);
		BOOST_CHECK_EQUAL(branch(part2_5).list(1).size(), 0);
		BOOST_CHECK_EQUAL(branch(part1_3_4).list(1).size(), 0);

		BOOST_CHECK_EQUAL(grid(pos1), -1.0f);
		BOOST_CHECK_EQUAL(grid(pos2), -1.0f);
		BOOST_CHECK_EQUAL(grid(pos3),  3.0f);
		BOOST_CHECK_EQUAL(grid(pos4),  4.0f);
		BOOST_CHECK_EQUAL(grid(pos5),  5.0f);

		branch(part1_3_4).remove(0, 0);

		BOOST_CHECK_EQUAL(branch(part1_3_4).list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1_3_4).list(0)[0], pos4);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);

		branch(part1_3_4).remove(0, 0);

		BOOST_CHECK_EQUAL(branch(part1_3_4).list(0).size(), 0);
		BOOST_CHECK_EQUAL(branch.list(0).size(), 1);

		branch.remove(0, 2);

		BOOST_CHECK_EQUAL(branch.list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part2_5).list(2).size(), 1);
	}


	BOOST_AUTO_TEST_CASE(array_pos_mapped)
	{
		typedef PosMappedPartitionedGrid<3, 3, 3> Grid_t;
		typedef PosMappedPartitionedGrid<3, 3, 3>::BranchGrid BranchGrid;
		typedef PosMappedPartitionedGrid<3, 3, 3>::LookupGrid LookupGrid;

		Grid_t grid(Vec3u(9,9,9), Vec3i(-4,-4,-4));
		BranchGrid& branch = grid.branch();
		LookupGrid& lookup = grid.lookup();


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
		BOOST_CHECK_EQUAL(lookup.list(0).size(), 2);
		BOOST_CHECK_EQUAL(lookup.list(2).size(), 1);
		BOOST_CHECK_EQUAL(lookup.list(0)[0], part1);
		BOOST_CHECK_EQUAL(lookup.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(lookup.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : grid.lookup().list(i))
				for (const Vec3i pos : branch(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(2);

		BOOST_CHECK_EQUAL(lookup.list(2).size(), 0);
		BOOST_CHECK_EQUAL(branch(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos4), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL(lookup(part4), Grid_t::NULL_IDX_TUPLE);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(branch(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid(pos2), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(lookup.list(0).size(), 1);
		BOOST_CHECK_EQUAL(branch(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), Grid_t::NULL_IDX_TUPLE);
		BOOST_CHECK_EQUAL(lookup(part1), Grid_t::NULL_IDX_TUPLE);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(lookup.list(i).size(), 0);

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

BOOST_AUTO_TEST_SUITE_END()
