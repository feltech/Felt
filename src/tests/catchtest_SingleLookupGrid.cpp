#include "catch.hpp"

#define _TESTING

#include "Felt/SingleLookupGrid.hpp"

using namespace felt;

SCENARIO("SingleLookupGrid")
{
	WHEN("intialise_and_populate")
	{
		using GridType = SingleLookupGrid<3, 3>;
		GridType grid(Vec3u(10,10,10), Vec3i(0, -5, -5));

		const Vec3i pos1(1, 0, -1);
		const Vec3i pos2(2, 1, 0);
		const Vec3i pos3(3, -1, 0);
		const Vec3i pos4(4, -1, 2);
		const Vec3i pos5(5, -2, 1);
		const Vec3i pos6(6, -2, 2);

		// Add the positions to the array and set index lookup values.
		grid.add(pos1, 0);
		grid.add(pos2, 1);
		grid.add(pos3, 1);
		grid.add(pos4, 2);

		CHECK(grid.list(0).size() == 1);
		CHECK(grid.list(1).size() == 2);
		CHECK(grid.list(2).size() == 1);
		CHECK(grid.list(0)[0] == pos1);
		CHECK(grid.list(1)[0] == pos2);
		CHECK(grid.list(1)[1] == pos3);
		CHECK(grid.list(2)[0] == pos4);
		CHECK((UINT)grid(pos1) == 0);
		CHECK((UINT)grid(pos2) == 0);
		CHECK((UINT)grid(pos3) == 1);
		CHECK((UINT)grid(pos4) == 0);

		grid.remove(pos2, 1);

		CHECK(grid.list(0).size() == 1);
		CHECK(grid.list(1).size() == 1);
		CHECK(grid.list(2).size() == 1);
		CHECK(grid.list(0)[0] == pos1);
		CHECK(grid.list(1)[0] == pos3);
		CHECK(grid.list(2)[0] == pos4);
		CHECK((UINT)grid(pos1) == 0);
		CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos3) == 0);
		CHECK((UINT)grid(pos4) == 0);

		grid.add(pos5, 2);
		grid.add(pos6, 2);

		CHECK(grid.list(0).size() == 1);
		CHECK(grid.list(1).size() == 1);
		CHECK(grid.list(2).size() == 3);
		CHECK(grid.list(0)[0] == pos1);
		CHECK(grid.list(1)[0] == pos3);
		CHECK(grid.list(2)[0] == pos4);
		CHECK(grid.list(2)[1] == pos5);
		CHECK(grid.list(2)[2] == pos6);
		CHECK((UINT)grid(pos1) == 0);
		CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos3) == 0);
		CHECK((UINT)grid(pos4) == 0);
		CHECK((UINT)grid(pos5) == 1);
		CHECK((UINT)grid(pos6) == 2);

		grid.remove(pos4, 2);
		grid.remove(0, 0);

		CHECK(grid.list(0).size() == 0);
		CHECK(grid.list(1).size() == 1);
		CHECK(grid.list(2).size() == 2);
		CHECK(grid.list(0)[0] == pos1);
		CHECK(grid.list(1)[0] == pos3);
		CHECK(grid.list(2)[1] == pos5);
		CHECK(grid.list(2)[0] == pos6);
		CHECK((UINT)grid(pos1) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos3) == 0);
		CHECK((UINT)grid(pos4) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos5) == 1);
		CHECK((UINT)grid(pos6) == 0);

		grid.reset(2);

		CHECK(grid.list(0).size() == 0);
		CHECK(grid.list(1).size() == 1);
		CHECK(grid.list(2).size() == 0);
		CHECK(grid.list(1)[0] == pos3);
		CHECK((UINT)grid(pos1) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos2) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos3) == 0);
		CHECK((UINT)grid(pos4) == GridType::NULL_IDX);
		CHECK((UINT)grid(pos5) == GridType::NULL_IDX);

	}
}


SCENARIO("LazySingleLookupGrid")
{
	WHEN("initialisation")
	{
		/// [LazySingleLookupGrid initialisation]
		// ==== Setup ====
		LazySingleLookupGrid<3, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
		const UINT NULL_IDX_DATA = LazySingleLookupGrid<3, 3>::NULL_IDX;

		// ==== Confirm ====
		CHECK(grid.is_active() == false);
		CHECK(grid.data().size() == 0);
		CHECK(grid.background() == NULL_IDX_DATA);
		CHECK(grid.get(Vec3i(1,1,1)) == NULL_IDX_DATA);
		/// [LazySingleLookupGrid initialisation]
	}
}
