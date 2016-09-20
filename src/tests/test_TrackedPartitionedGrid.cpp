#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include <Felt/SingleTrackedPartitionedGrid.hpp>

using namespace felt;

SCENARIO("LazySingleTrackedPartitionedGrid")
{
	WHEN("initialisation")
	{
		/// [LazySingleTrackedPartitionedGrid initialisation]
		// ==== Setup ====
		using GridType = SingleTrackedPartitionedGrid<FLOAT, 3, 3>;
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
	using GridType = SingleTrackedPartitionedGrid<FLOAT, 3, 3>;
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
	grid.remove(pos1, 0, 7);

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
	grid.remove(pos2, 1, 7);

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

