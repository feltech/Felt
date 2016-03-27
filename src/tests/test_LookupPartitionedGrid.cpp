#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include <Felt/LookupPartitionedGrid.hpp>

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_LookupPartitionedGrid)
	/**
	 * Simple lookup get and set values.
	 */
	BOOST_AUTO_TEST_CASE(initialise_and_populate)
	{
		typedef LookupPartitionedGrid<3, 3> GridType;
		typedef LookupPartitionedGrid<3, 3>::ChildrenGrid ChildrenGrid;
		typedef LookupPartitionedGrid<3, 3>::ChildrenGrid::Lookup LookupGrid;

		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		ChildrenGrid& children = grid.children();
		LookupGrid& lookup = children.lookup();


		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), GridType::MixinType::Traits::NULL_IDX_DATA);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), GridType::MixinType::Traits::NULL_IDX_DATA);
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
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)children(part4)(pos4)(2), 0);
		BOOST_CHECK_EQUAL(children.list(0).size(), 2);
		BOOST_CHECK_EQUAL(children.list(2).size(), 1);
		BOOST_CHECK_EQUAL(children.list(0)[0], part1);
		BOOST_CHECK_EQUAL(children.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(children.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : children.list(i))
				for (const Vec3i pos : children(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(2);

		BOOST_CHECK_EQUAL(children.list(2).size(), 0);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos4), GridType::MixinType::Traits::NULL_IDX_DATA);
		BOOST_CHECK_EQUAL(lookup(part4), GridType::MixinType::Traits::NULL_IDX_DATA);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid(pos2), GridType::MixinType::Traits::NULL_IDX_DATA);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(children.list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), GridType::MixinType::Traits::NULL_IDX_DATA);
		BOOST_CHECK_EQUAL(lookup(part1), GridType::MixinType::Traits::NULL_IDX_DATA);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(children.list(i).size(), 0);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(grid(pos), GridType::MixinType::Traits::NULL_IDX_DATA);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					BOOST_CHECK_EQUAL(lookup(pos), GridType::MixinType::Traits::NULL_IDX_DATA);
					for (UINT i = 0; i < 3; i++)
						BOOST_CHECK_EQUAL(children(pos).list(i).size(), 0);
				}
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_SharedLookupPartitionedGrid)
	BOOST_AUTO_TEST_CASE(initialise_and_populate)
	{
		typedef SharedLookupPartitionedGrid<3, 3> GridType;
		typedef GridType::ChildrenGrid ChildrenGrid;
		typedef ChildrenGrid::Lookup LookupGrid;
		const Vec3u& BRANCH_NULL_IDX = LookupGrid::Traits::NULL_IDX_DATA;
		const UINT CHILD_NULL_IDX = GridType::NULL_IDX;

		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		ChildrenGrid& children = grid.children();
		LookupGrid& lookup = children.lookup();


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
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)children(part4)(pos4), 0);
		BOOST_CHECK_EQUAL(children.list(0).size(), 2);
		BOOST_CHECK_EQUAL(children.list(2).size(), 1);
		BOOST_CHECK_EQUAL(children.list(0)[0], part1);
		BOOST_CHECK_EQUAL(children.list(0)[1], part2_3);
		BOOST_CHECK_EQUAL(children.list(2)[0], part4);
		BOOST_CHECK_EQUAL((UINT)lookup(part1)(0), 0);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);
		BOOST_CHECK_EQUAL((UINT)lookup(part4)(2), 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : children.list(i))
				for (const Vec3i pos : children(pos_child).list(i))
					apos.push_back(pos);

		BOOST_CHECK_EQUAL(apos[0], pos1);
		BOOST_CHECK_EQUAL(apos[1], pos2);
		BOOST_CHECK_EQUAL(apos[2], pos3);
		BOOST_CHECK_EQUAL(apos[3], pos4);

		grid.reset(2);

		BOOST_CHECK_EQUAL(children.list(2).size(), 0);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos4), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part4), BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid(pos2), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(children.list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid(pos1), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part1), BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(children.list(i).size(), 0);

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
						BOOST_CHECK_EQUAL(children(pos).list(i).size(), 0);
				}
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_LazySharedLookupPartitionedGrid)

	struct Fixture {
		const UINT NULL_IDX = LazySharedLookupPartitionedGrid<3, 3>::NULL_IDX;
		LazySharedLookupPartitionedGrid<3, 3> grid;
		Fixture()
		: grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3))
		{}
	};

	struct ResetFixture : Fixture {
		const UINT NULL_IDX = LazySharedLookupPartitionedGrid<3, 3>::NULL_IDX;
		PartitionedGrid<FLOAT, 3> grid_master;
		ResetFixture()
		: Fixture(),
		  grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3))
		{}
	};

	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazySharedLookupPartitionedGrid initialisation]
		// ==== Setup ====
		LazySharedLookupPartitionedGrid<3, 3> grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		const UINT NULL_IDX = LazySharedLookupPartitionedGrid<3, 3>::NULL_IDX;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL((bool)grid.children().get(Vec3i(1,1,1)).is_active(), false);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).data().size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).background(), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySharedLookupPartitionedGrid initialisation]
	}

	BOOST_FIXTURE_TEST_CASE(add_should_activate_once, Fixture)
	{
		// ==== Setup ====
		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos_child(-1, -1, -1);

		// ==== Action ====
		grid.add(pos1, 0);
		grid.add(pos2, 1);

		// ==== Confirm ====
		BOOST_CHECK(grid.children().get(pos_child).is_active());
		BOOST_CHECK_EQUAL(grid.children().list().size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos1), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos2), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(1).size(), 1);
	}

	BOOST_FIXTURE_TEST_CASE(remove_should_deactivate_when_child_is_inactive, Fixture)
	{
		// ==== Setup ====
		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos_child(-1, -1, -1);
		grid.add(pos1, 0);
		grid.add(pos2, 1);

		// ==== Action ====
		grid.remove(pos1, 0);

		// ==== Confirm ====
		BOOST_CHECK(grid.children().get(pos_child).is_active());
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos1), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos2), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(1).size(), 1);

		// ==== Action ====
		grid.remove(pos2, 1);

		// ==== Confirm ====
		BOOST_CHECK(!grid.children().get(pos_child).is_active());
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().list(1).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos1), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos2), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(1).size(), 0);
	}

	BOOST_FIXTURE_TEST_CASE(reset_shouldnt_deactivate_when_other_list_still_active, ResetFixture)
	{
		// ==== Setup ====
		const Vec3i pos_child(-1, -1, -1);
		const Vec3i pos(-4, -4, -4);
		grid.add(pos, 0);

		// ==== Action ====
		grid.reset(grid_master, 1);

		// ==== Confirm ====

		BOOST_CHECK_EQUAL(grid.get(pos), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).data().size(), 3*3*3);
	}

	BOOST_FIXTURE_TEST_CASE(reset_shouldnt_deactivate_when_master_grid_is_tracking, ResetFixture)
	{
		// ==== Setup ====
		const Vec3i pos_child(-1, -1, -1);
		const Vec3i pos(-4, -4, -4);

		grid_master.add_child(pos_child);
		grid.add(pos, 0);

		// ==== Action ====
		grid.reset(grid_master, 0);

		// ==== Confirm ====

		BOOST_CHECK_EQUAL(grid.get(pos), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).data().size(), 3*3*3);
	}

	BOOST_AUTO_TEST_CASE(reset_mixed_cases)
	{
		/// [LazySharedLookupPartitionedGrid reset_mixed_cases]
		// ==== Setup ====
		const UINT NULL_IDX = LazySharedLookupPartitionedGrid<3, 3>::NULL_IDX;
		PartitionedGrid<FLOAT, 3> grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		LazySharedLookupPartitionedGrid<3, 3> grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));

		const Vec3i pos_deactivated(0, 0, 0);
		const Vec3i pos_active_because_master(-4, 0, 4);
		const Vec3i pos_active_because_other_list(4, 0, 0);
		const Vec3i pos_child_deactivated(0, 0, 0);
		const Vec3i pos_child_active_because_master(-1, 0, 1);
		const Vec3i pos_child_active_because_other_list(1, 0, 0);

		grid_master.add_child(pos_child_active_because_master);
		grid.add(pos_active_because_master);
		grid.add(pos_deactivated);
		grid.add(pos_active_because_other_list, 1);

		// ==== Action ====
		grid.reset(grid_master, 0);

		// ==== Confirm ====

		// Behaves like standard lookup grid at the child level,
		BOOST_CHECK_EQUAL(grid.get(pos_active_because_master), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.get(pos_deactivated), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.get(pos_active_because_other_list), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child_active_because_master).list().size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child_deactivated).list().size(), 0);
		BOOST_CHECK_EQUAL(
			grid.children().get(pos_child_active_because_other_list).list(1).size(), 1
		);

		// but destroys inactive partitions,
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(pos_child_deactivated)).is_active(), false);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(pos_child_deactivated)).data().size(), 0);

		// except for partitions being tracked by the master grid,
		BOOST_CHECK_EQUAL(grid.children().get(pos_child_active_because_master).is_active(), true);
		BOOST_CHECK_EQUAL(
			grid.children().get(pos_child_active_because_master).data().size(), 3*3*3
		);
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 1);

		// and except for partitions that still have active lists.
		BOOST_CHECK_EQUAL(grid.children().list(1).size(), 1);
		BOOST_CHECK_EQUAL(
			grid.children().get(pos_child_active_because_other_list).is_active(), true
		);
		BOOST_CHECK_EQUAL(
			grid.children().get(pos_child_active_because_other_list).data().size(), 3*3*3
		);
		/// [LazySharedLookupPartitionedGrid reset_mixed_cases]
	}
BOOST_AUTO_TEST_SUITE_END()
