#include "catch.hpp"

#define _TESTING

#include <Felt/SingleTrackedGrid.hpp>

using namespace felt;

SCENARIO("LazySingleTrackedGrid")
{
	WHEN("initialisation")
	{
		// ==== Setup ====
		LazySingleTrackedGrid<FLOAT, 3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3);
		const UINT NULL_IDX = LazySingleTrackedGrid<FLOAT, 3, 3>::Lookup::NULL_IDX;

		// ==== Confirm ====
		CHECK(grid.is_active() == false);
		CHECK(grid.data().size() == 0);
		CHECK(grid.background() == 3);
		CHECK(grid.get(Vec3i(1,1,1)) == 3);
		CHECK(grid.lookup().is_active() == false);
		CHECK(grid.lookup().data().size() == 0);
		CHECK(grid.lookup().background() == NULL_IDX);
		CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
	}
}


struct LazySingleTrackedGridFixture {
	const UINT NULL_IDX = LazySingleTrackedGrid<FLOAT, 3, 3>::Lookup::NULL_IDX;
	LazySingleTrackedGrid<FLOAT, 3, 3> grid;
	LazySingleTrackedGridFixture()
			: grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 3)
	{}
};

TEST_CASE_METHOD(LazySingleTrackedGridFixture, "LazySingleTrackedGrid: activate_should_activate_lookup", "[LazySingleTrackedGrid]")
{
	/// [LazySingleTrackedGrid activate]
	// ==== Action ====
	grid.activate();

	// ==== Confirm ====
	CHECK(grid.is_active() == true);
	CHECK(grid.data().size() == 3*3*3);
	CHECK(grid.get(Vec3i(1,1,1)) == 3);
	CHECK(grid.lookup().is_active() == true);
	CHECK(grid.lookup().data().size() == 3*3*3);
	CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
	/// [LazySingleTrackedGrid activate]
}

TEST_CASE_METHOD(LazySingleTrackedGridFixture, "LazySingleTrackedGrid: deactivate_should_deactivate_lookup", "[LazySingleTrackedGrid]")
{
	/// [LazySingleTrackedGrid deactivate]
	// ==== Action ====
	grid.activate();
	grid.deactivate();

	// ==== Confirm ====
	CHECK(grid.is_active() == false);
	CHECK(grid.data().size() == 0);
	CHECK(grid.background() == 3);
	CHECK(grid.get(Vec3i(1,1,1)) == 3);
	CHECK(grid.lookup().is_active() == false);
	CHECK(grid.lookup().data().size() == 0);
	CHECK(grid.lookup().background() == NULL_IDX);
	CHECK(grid.lookup().get(Vec3i(1,1,1)) == NULL_IDX);
	/// [LazySingleTrackedGrid deactivate]
}

