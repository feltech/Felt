#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include <Felt/LookupPartitionedGrid.hpp>

using namespace felt;


SCENARIO("LazySingleLookupPartitionedGrid")
{
	WHEN("initialise_and_populate")
	{
		// ==== Setup ====
		typedef LookupPartitionedGrid<3, 3> GridType;
		typedef GridType::ChildrenGrid ChildrenGrid;
		typedef ChildrenGrid::Lookup LookupGrid;
		const Vec3u& BRANCH_NULL_IDX = LookupGrid::Traits::NULL_IDX_DATA;
		const UINT CHILD_NULL_IDX = GridType::NULL_IDX;

		// ==== Action ====
		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		ChildrenGrid& children = grid.children();
		LookupGrid& lookup = children.lookup();

		// ==== Confirm ====
		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					CHECK(grid(pos) == CHILD_NULL_IDX);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					CHECK(lookup(pos) == BRANCH_NULL_IDX);
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

		CHECK((UINT)grid(pos1) == 0);
		CHECK((UINT)grid(pos2) == 0);
		CHECK((UINT)grid(pos3) == 1);
		CHECK((UINT)grid(pos4) == 0);
		CHECK(children(part1).list(0).size() == 1);
		CHECK(children(part2_3).list(0).size() == 2);
		CHECK(children(part4).list(2).size() == 1);
		CHECK((UINT)children(part4)(pos4) == 0);
		CHECK(children.list(0).size() == 2);
		CHECK(children.list(2).size() == 1);
		CHECK(children.list(0)[0] == part1);
		CHECK(children.list(0)[1] == part2_3);
		CHECK(children.list(2)[0] == part4);
		CHECK((UINT)lookup(part1)(0) == 0);
		CHECK((UINT)lookup(part2_3)(0) == 1);
		CHECK((UINT)lookup(part4)(2) == 0);

		std::vector<Vec3i> apos;
		for (UINT i = 0; i < 3; i++)
			for (const Vec3i& pos_child : children.list(i))
				for (const Vec3i pos : children(pos_child).list(i))
					apos.push_back(pos);

		CHECK(apos[0] == pos1);
		CHECK(apos[1] == pos2);
		CHECK(apos[2] == pos3);
		CHECK(apos[3] == pos4);

		grid.reset(2);

		CHECK(children.list(2).size() == 0);
		CHECK(children(part4).list(2).size() == 0);
		CHECK(grid(pos4) == CHILD_NULL_IDX);
		CHECK(lookup(part4) == BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		CHECK(children(part2_3).list(0).size() == 1);
		CHECK(grid(pos2) == CHILD_NULL_IDX);
		CHECK((UINT)lookup(part2_3)(0) == 1);

		grid.remove(pos1, 0);

		CHECK(children.list(0).size() == 1);
		CHECK(children(part1).list(0).size() == 0);
		CHECK(grid(pos1) == CHILD_NULL_IDX);
		CHECK(lookup(part1) == BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			CHECK(children.list(i).size() == 0);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					Vec3i pos(x, y, z);
					CHECK(grid(pos) == CHILD_NULL_IDX);
				}
		for (INT x = -1; x <= 1; x++)
			for (INT y = -1; y <= 1; y++)
				for (INT z = -1; z <= 1; z++)
				{
					Vec3i pos(x, y, z);
					CHECK(lookup(pos) == BRANCH_NULL_IDX);
					for (UINT i = 0; i < 3; i++)
						CHECK(children(pos).list(i).size() == 0);
				}
	}

	WHEN("initialisation")
	{
		/// [LazySingleLookupPartitionedGrid initialisation]
		// ==== Setup ====
		LookupPartitionedGrid<3, 3> grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));
		const UINT NULL_IDX = LookupPartitionedGrid<3, 3>::NULL_IDX;

		// ==== Confirm ====
		CHECK((bool)grid.children().get(Vec3i(1,1,1)).is_active() == false);
		CHECK(grid.children().get(Vec3i(1,1,1)).data().size() == 0);
		CHECK(grid.children().get(Vec3i(1,1,1)).background() == NULL_IDX);
		CHECK(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)) == NULL_IDX);
		/// [LazySingleLookupPartitionedGrid initialisation]
	}

	WHEN("reset_mixed_cases")
	{
		/// [LazySingleLookupPartitionedGrid reset_mixed_cases]
		// ==== Setup ====
		const UINT NULL_IDX = LookupPartitionedGrid<3, 3>::NULL_IDX;
		PartitionedGrid<FLOAT, 3> grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3));
		LookupPartitionedGrid<3, 3> grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));

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
		CHECK(grid.get(pos_active_because_master) == NULL_IDX);
		CHECK(grid.get(pos_deactivated) == NULL_IDX);
		CHECK(grid.get(pos_active_because_other_list) == 0);
		CHECK(grid.children().get(pos_child_active_because_master).list(0).size() == 0);
		CHECK(grid.children().get(pos_child_deactivated).list(0).size() == 0);
		CHECK(
			grid.children().get(pos_child_active_because_other_list).list(1).size() == 1
		);

		// but destroys inactive partitions,
		CHECK(grid.children().get(Vec3i(pos_child_deactivated)).is_active() == false);
		CHECK(grid.children().get(Vec3i(pos_child_deactivated)).data().size() == 0);

		// except for partitions being tracked by the master grid,
		CHECK(grid.children().get(pos_child_active_because_master).is_active() == true);
		CHECK(
			grid.children().get(pos_child_active_because_master).data().size() == 3*3*3
		);
		CHECK(grid.children().list(0).size() == 0);

		// and except for partitions that still have active lists.
		CHECK(grid.children().list(1).size() == 1);
		CHECK(
			grid.children().get(pos_child_active_because_other_list).is_active() == true
		);
		CHECK(
			grid.children().get(pos_child_active_because_other_list).data().size() == 3*3*3
		);
		/// [LazySingleLookupPartitionedGrid reset_mixed_cases]
	}

	WHEN("reset_all")
	{
		/// [LazySingleLookupPartitionedGrid reset_all]
		// ==== Setup ====

		const UINT NULL_IDX = LookupPartitionedGrid<3, 3>::NULL_IDX;
		PartitionedGrid<FLOAT, 3> grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3));
		LookupPartitionedGrid<3, 3> grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3));

		const Vec3i pos_list_0(0, 0, 0);
		const Vec3i pos_active_because_master(-4, 0, 4);
		const Vec3i pos_list_1(4, 0, 0);
		const Vec3i pos_child_list_0(0, 0, 0);
		const Vec3i pos_child_active_because_master(-1, 0, 1);
		const Vec3i pos_child_list_1(1, 0, 0);

		grid_master.add_child(pos_child_active_because_master);
		grid.add(pos_active_because_master);
		grid.add(pos_list_0);
		grid.add(pos_list_1, 1);

		// ==== Action ====

		grid.reset_all(grid_master);

		// ==== Confirm ====

		// Resets all children.
		CHECK(grid.get(pos_active_because_master) == NULL_IDX);
		CHECK(grid.get(pos_list_0) == NULL_IDX);
		CHECK(grid.get(pos_list_1) == NULL_IDX);
		CHECK(grid.children().get(pos_child_active_because_master).list().size() == 0);
		CHECK(grid.children().get(pos_child_list_0).list().size() == 0);
		CHECK(grid.children().get(pos_child_list_1).list(1).size() == 0);

		// Deactivates partitions not being tracked by the master.
		CHECK(grid.children().get(Vec3i(pos_child_list_0)).is_active() == false);
		CHECK(grid.children().get(Vec3i(pos_child_list_0)).data().size() == 0);
		CHECK(grid.children().get(pos_child_list_1).is_active() == false);
		CHECK(grid.children().get(pos_child_list_1).data().size() == 0);
		CHECK(grid.children().list(1).size() == 0);

		// Leaves active those partitions being tracked by the master grid.
		CHECK(grid.children().get(pos_child_active_because_master).is_active() == true);
		CHECK(
			grid.children().get(pos_child_active_because_master).data().size() == 3*3*3
		);
		CHECK(grid.children().list(0).size() == 0);
		/// [LazySingleLookupPartitionedGrid reset_all]
	}

}

struct LazySingleLookupPartitionedGridFixture {
	const UINT NULL_IDX = LookupPartitionedGrid<3, 3>::NULL_IDX;
	LookupPartitionedGrid<3, 3> grid;
	LazySingleLookupPartitionedGridFixture()
			: grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), Vec3u(3, 3, 3))
	{}
};

struct LazySingleLookupPartitionedGridResetFixture : LazySingleLookupPartitionedGridFixture {
	const UINT NULL_IDX = LookupPartitionedGrid<3, 3>::NULL_IDX;
	PartitionedGrid<FLOAT, 3> grid_master;
	LazySingleLookupPartitionedGridResetFixture()
			: LazySingleLookupPartitionedGridFixture(),
			  grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3))
	{}
};


TEST_CASE_METHOD(LazySingleLookupPartitionedGridFixture, "LazySingleLookupPartitionedGrid: add_should_activate_once", "[LazySingleLookupPartitionedGrid]")
{
	// ==== Setup ====
	const Vec3i pos1(-4, -4, -4);
	const Vec3i pos2(-3, -4, -4);
	const Vec3i pos_child(-1, -1, -1);

	// ==== Action ====
	grid.add(pos1, 0);
	grid.add(pos2, 1);

	// ==== Confirm ====
	CHECK(grid.children().get(pos_child).is_active());
	CHECK(grid.children().list(0).size() == 1);
	CHECK(grid.children().get(pos_child).get(pos1) == 0);
	CHECK(grid.children().get(pos_child).get(pos2) == 0);
	CHECK(grid.children().get(pos_child).list(0).size() == 1);
	CHECK(grid.children().get(pos_child).list(1).size() == 1);
}

TEST_CASE_METHOD(LazySingleLookupPartitionedGridFixture, "LazySingleLookupPartitionedGrid: remove_should_deactivate_when_child_is_inactive", "[LazySingleLookupPartitionedGrid]")
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
	CHECK(grid.children().get(pos_child).is_active());
	CHECK(grid.children().list(0).size() == 0);
	CHECK(grid.children().list(1).size() == 1);
	CHECK(grid.children().get(pos_child).get(pos1) == NULL_IDX);
	CHECK(grid.children().get(pos_child).get(pos2) == 0);
	CHECK(grid.children().get(pos_child).list(0).size() == 0);
	CHECK(grid.children().get(pos_child).list(1).size() == 1);

	// ==== Action ====
	grid.remove(pos2, 1);

	// ==== Confirm ====
	CHECK(!grid.children().get(pos_child).is_active());
	CHECK(grid.children().list(0).size() == 0);
	CHECK(grid.children().list(1).size() == 0);
	CHECK(grid.children().get(pos_child).get(pos1) == NULL_IDX);
	CHECK(grid.children().get(pos_child).get(pos2) == NULL_IDX);
	CHECK(grid.children().get(pos_child).list(0).size() == 0);
	CHECK(grid.children().get(pos_child).list(1).size() == 0);
}

TEST_CASE_METHOD(LazySingleLookupPartitionedGridResetFixture, "LazySingleLookupPartitionedGrid: reset_shouldnt_deactivate_when_other_list_still_active", "[LazySingleLookupPartitionedGrid]")
{
	// ==== Setup ====
	const Vec3i pos_child(-1, -1, -1);
	const Vec3i pos(-4, -4, -4);
	grid.add(pos, 0);

	// ==== Action ====
	grid.reset(grid_master, 1);

	// ==== Confirm ====

	CHECK(grid.get(pos) == 0);
	CHECK(grid.children().get(pos_child).list(0).size() == 1);
	CHECK(grid.children().get(pos_child).is_active() == true);
	CHECK(grid.children().get(pos_child).data().size() == 3*3*3);
}

TEST_CASE_METHOD(LazySingleLookupPartitionedGridResetFixture, "LazySingleLookupPartitionedGrid: reset_shouldnt_deactivate_when_master_grid_is_tracking", "[LazySingleLookupPartitionedGrid]")
{
	// ==== Setup ====
	const Vec3i pos_child(-1, -1, -1);
	const Vec3i pos(-4, -4, -4);

	grid_master.add_child(pos_child);
	grid.add(pos, 0);

	// ==== Action ====
	grid.reset(grid_master, 0);

	// ==== Confirm ====

	CHECK(grid.get(pos) == NULL_IDX);
	CHECK(grid.children().list(0).size() == 0);
	CHECK(grid.children().get(pos_child).list(0).size() == 0);
	CHECK(grid.children().get(pos_child).is_active() == true);
	CHECK(grid.children().get(pos_child).data().size() == 3*3*3);
}

