#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include <Felt/TrackedPartitionedGrid.hpp>

using namespace felt;

SCENARIO("SingleTrackedPartitionedGrid")
{
	WHEN("initialise_and_populate")
	{
		typedef SingleTrackedPartitionedGrid<FLOAT, 3, 3> GridType;
		typedef GridType::ChildrenGrid ChildrenGrid;
		typedef ChildrenGrid::Lookup LookupGrid;
		const Vec3u& BRANCH_NULL_IDX = LookupGrid::Traits::NULL_IDX_DATA;
		const UINT CHILD_NULL_IDX = GridType::Child::Lookup::NULL_IDX;

		GridType grid(Vec3u(9,9,9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3));
		ChildrenGrid& children = grid.children();
		LookupGrid& lookup = children.lookup();

		grid.fill(-1.0f);

		for (INT x = -4; x <= 4; x++)
			for (INT y = -4; y <= 4; y++)
				for (INT z = -4; z <= 4; z++)
				{
					const Vec3i pos(x, y, z);
					const Vec3i pos_child = grid.pos_child(pos);
					CHECK(grid(pos) == -1.0f);
					CHECK(
						grid.children().get(pos_child).lookup()(pos) == CHILD_NULL_IDX
					);
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

		grid.add(pos1, 1.0f, 0);
		grid.add(pos2, 2.0f, 0);
		grid.add(pos3, 3.0f, 0);
		grid.add(pos4, 4.0f, 2);

		CHECK((UINT)grid(pos1) == 1.0f);
		CHECK((UINT)grid(pos2) == 2.0f);
		CHECK((UINT)grid(pos3) == 3.0f);
		CHECK((UINT)grid(pos4) == 4.0f);
		CHECK(children(part1).list(0).size() == 1);
		CHECK(children(part2_3).list(0).size() == 2);
		CHECK(children(part4).list(2).size() == 1);
		CHECK((UINT)children(part4)(pos4) == 4.0f);
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

		grid.reset(-2.0f, 2);

		CHECK(grid(pos4) == -2.0f);
		CHECK(children.list(2).size() == 0);
		CHECK(children(part4).list(2).size() == 0);
		CHECK(children(part4).lookup()(pos4) == CHILD_NULL_IDX);
		CHECK(lookup(part4) == BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		CHECK(grid(pos2) == 2.0f);
		CHECK(children.list(0).size() == 2);
		CHECK(children(part2_3).list(0).size() == 1);
		CHECK(children(part2_3).lookup()(pos2) == CHILD_NULL_IDX);
		CHECK((UINT)lookup(part2_3)(0) == 1);

		grid.remove(pos1, 0);

		CHECK(grid(pos1) == 1.0f);
		CHECK(children.list(0).size() == 1);
		CHECK(children(part1).list(0).size() == 0);
		CHECK(children(part1).lookup()(pos1) == CHILD_NULL_IDX);
		CHECK(lookup(part1) == BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			CHECK(children.list(i).size() == 0);

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
}


SCENARIO("LazySingleTrackedPartitionedGrid")
{

	WHEN("initialisation")
	{
		/// [LazySingleTrackedPartitionedGrid initialisation]
		// ==== Setup ====
		using GridType = LazySingleTrackedPartitionedGrid<FLOAT, 3, 3>;
		const UINT NULL_IDX = GridType::Child::Lookup::NULL_IDX;

		// ==== Action ====
		GridType grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 7.0f, Vec3u(3, 3, 3));

		// ==== Confirm ====
		CHECK((bool)grid.children().get(Vec3i(1,1,1)).is_active() == false);
		CHECK(grid.children().get(Vec3i(1,1,1)).background() == 7);
		CHECK(grid.children().get(Vec3i(1,1,1)).data().size() == 0);
		CHECK(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)) == 7);
		CHECK(grid.children().get(Vec3i(1,1,1)).lookup().data().size() == 0);
		CHECK(grid.children().get(Vec3i(1,1,1)).lookup().get(Vec3i(1,1,1)) == NULL_IDX);
		/// [LazySingleTrackedPartitionedGrid initialisation]
	}

}


struct LazySingleTrackedPartitionedGridFixture
{
	using GridType = LazySingleTrackedPartitionedGrid<FLOAT, 3, 3>;
	const UINT NULL_IDX = GridType::Lookup::NULL_IDX;
	GridType grid;
	LazySingleTrackedPartitionedGridFixture()
			: grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 7, Vec3u(3, 3, 3))
	{}
};


struct LazySingleTrackedPartitionedGridResetFixture : LazySingleTrackedPartitionedGridFixture
{
	PartitionedGrid<FLOAT, 3> grid_master;
	LazySingleTrackedPartitionedGridResetFixture()
			: LazySingleTrackedPartitionedGridFixture(),
			  grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3))
	{}
};


TEST_CASE_METHOD(LazySingleTrackedPartitionedGridFixture, "LazySingleTrackedPartitionedGrid: add_should_activate_once", "[LazySingleTrackedPartitionedGrid]")
{
	// ==== Setup ====
	const Vec3i pos1(-4, -4, -4);
	const Vec3i pos2(-3, -4, -4);
	const Vec3i pos_child(-1, -1, -1);

	// ==== Action ====
	grid.add(pos1, 3, 0);
	grid.add(pos2, 4, 1);

	// ==== Confirm ====
	CHECK(grid.children().get(pos_child).is_active());
	CHECK(grid.children().list().size() == 1);
	CHECK(grid.children().get(pos_child).get(pos1) == 3);
	CHECK(grid.children().get(pos_child).get(pos2) == 4);
	CHECK(grid.children().get(pos_child).lookup().list(0).size() == 1);
	CHECK(grid.children().get(pos_child).lookup().list(1).size() == 1);
	CHECK(grid.children().get(pos_child).lookup().get(pos1) == 0);
	CHECK(grid.children().get(pos_child).lookup().get(pos2) == 0);
	CHECK(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)) == 7);
	CHECK(grid.children().get(Vec3i(1,1,1)).lookup().get(Vec3i(1,1,1)) == NULL_IDX);
}

TEST_CASE_METHOD(LazySingleTrackedPartitionedGridFixture, "LazySingleTrackedPartitionedGrid: remove_should_deactivate_when_child_is_inactive", "[LazySingleTrackedPartitionedGrid]")
{
	// ==== Setup ====
	const Vec3i pos1(-4, -4, -4);
	const Vec3i pos2(-3, -4, -4);
	const Vec3i pos_child(-1, -1, -1);
	grid.add(pos1, 3, 0);
	grid.add(pos2, 4, 1);

	// ==== Action ====
	grid.remove(pos1, 0);

	// ==== Confirm ====
	CHECK(grid.children().get(pos_child).is_active());
	CHECK(grid.children().get(pos_child).lookup().is_active());
	CHECK(grid.children().list(0).size() == 0);
	CHECK(grid.children().list(1).size() == 1);
	CHECK(grid.children().get(pos_child).lookup().get(pos1) == NULL_IDX);
	CHECK(grid.children().get(pos_child).lookup().get(pos2) == 0);
	CHECK(grid.children().get(pos_child).list(0).size() == 0);
	CHECK(grid.children().get(pos_child).list(1).size() == 1);

	// ==== Action ====
	grid.remove(pos2, 1);

	// ==== Confirm ====
	CHECK(!grid.children().get(pos_child).is_active());
	CHECK(!grid.children().get(pos_child).lookup().is_active());
	CHECK(grid.children().list(0).size() == 0);
	CHECK(grid.children().list(1).size() == 0);
	CHECK(grid.children().get(pos_child).lookup().get(pos1) == NULL_IDX);
	CHECK(grid.children().get(pos_child).lookup().get(pos2) == NULL_IDX);
	CHECK(grid.children().get(pos_child).list(0).size() == 0);
	CHECK(grid.children().get(pos_child).list(1).size() == 0);
}

TEST_CASE_METHOD(LazySingleTrackedPartitionedGridResetFixture, "LazySingleTrackedPartitionedGrid: reset_should_deactivate", "[LazySingleTrackedPartitionedGrid]")
{
	// ==== Setup ====

	const Vec3i pos_child(-1, -1, -1);
	const Vec3i pos(-4, -4, -4);
	grid.add(pos, 4, 0);

	// ==== Action ====

	grid.reset(grid_master, 0);

	// ==== Confirm ====

	// Value reset.
	CHECK(grid.get(pos) == 7);
	// Child still tracked.
	CHECK(
		grid.children().lookup().get(pos_child) == Vec3u(NULL_IDX, NULL_IDX, NULL_IDX)
	);
	// Child inactive.
	CHECK(grid.children().get(pos_child).is_active() == false);
	CHECK(grid.children().get(pos_child).data().size() == 0);
	// Child lookup inactive.
	CHECK(grid.children().get(pos_child).lookup().is_active() == false);
	CHECK(grid.children().get(pos_child).lookup().data().size() == 0);
	// Position no longer tracked in child.
	CHECK(grid.children().get(pos_child).lookup().get(pos) == NULL_IDX);
	CHECK(grid.children().get(pos_child).lookup().list(0).size() == 0);
}

TEST_CASE_METHOD(LazySingleTrackedPartitionedGridResetFixture, "LazySingleTrackedPartitionedGrid: reset_shouldnt_deactivate_when_other_list_still_active", "[LazySingleTrackedPartitionedGrid]")
{
	// ==== Setup ====

	const Vec3i pos_child(-1, -1, -1);
	const Vec3i pos(-4, -4, -4);
	grid.add(pos, 4, 0);

	// ==== Action ====

	grid.reset(grid_master, 1);

	// ==== Confirm ====

	// Value unchanged.
	CHECK(grid.get(pos) == 4);
	// Child still tracked.
	CHECK(grid.children().lookup().get(pos_child) == Vec3u(0, NULL_IDX, NULL_IDX));
	// Child still active.
	CHECK(grid.children().get(pos_child).is_active() == true);
	CHECK(grid.children().get(pos_child).data().size() == 3*3*3);
	// Child lookup still active.
	CHECK(grid.children().get(pos_child).lookup().is_active() == true);
	CHECK(grid.children().get(pos_child).lookup().data().size() == 3*3*3);
	// Position still tracked in child.
	CHECK(grid.children().get(pos_child).lookup().get(pos) == 0);
	CHECK(grid.children().get(pos_child).lookup().list(0).size() == 1);
}

TEST_CASE_METHOD(LazySingleTrackedPartitionedGridResetFixture, "LazySingleTrackedPartitionedGrid: reset_shouldnt_deactivate_when_master_grid_is_tracking", "[LazySingleTrackedPartitionedGrid]")
{
	// ==== Setup ====

	const Vec3i pos_child(-1, -1, -1);
	const Vec3i pos(-4, -4, -4);

	grid_master.add_child(pos_child);
	grid.add(pos, 4, 0);

	// ==== Action ====

	grid.reset(grid_master, 0);

	// ==== Confirm ====

	// Value reset.
	CHECK(grid.get(pos) == 7);
	// Child no longer tracked.
	CHECK(grid.children().lookup().get(pos_child) == Vec3u(NULL_IDX, NULL_IDX, NULL_IDX));
	// Child still active.
	CHECK(grid.children().get(pos_child).is_active() == true);
	CHECK(grid.children().get(pos_child).data().size() == 3*3*3);
	// Child lookup still active.
	CHECK(grid.children().get(pos_child).lookup().is_active() == true);
	CHECK(grid.children().get(pos_child).lookup().data().size() == 3*3*3);
	// Position no longer tracked in child.
	CHECK(grid.children().get(pos_child).lookup().get(pos) == NULL_IDX);
	CHECK(grid.children().get(pos_child).lookup().list(0).size() == 0);
}

