#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include <Felt/TrackedPartitionedGrid.hpp>

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_SingleTrackedPartitionedGrid)
	BOOST_AUTO_TEST_CASE(initialise_and_populate)
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
					BOOST_CHECK_EQUAL(grid(pos), -1.0f);
					BOOST_CHECK_EQUAL(
						grid.children().get(pos_child).lookup()(pos), CHILD_NULL_IDX
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
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 2);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 1);
		BOOST_CHECK_EQUAL((UINT)children(part4)(pos4), 4.0f);
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

		grid.reset(-2.0f, 2);

		BOOST_CHECK_EQUAL(grid(pos4), -2.0f);
		BOOST_CHECK_EQUAL(children.list(2).size(), 0);
		BOOST_CHECK_EQUAL(children(part4).list(2).size(), 0);
		BOOST_CHECK_EQUAL(children(part4).lookup()(pos4), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part4), BRANCH_NULL_IDX);

		grid.remove(pos2, 0);

		BOOST_CHECK_EQUAL(grid(pos2), 2.0f);
		BOOST_CHECK_EQUAL(children.list(0).size(), 2);
		BOOST_CHECK_EQUAL(children(part2_3).list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part2_3).lookup()(pos2), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL((UINT)lookup(part2_3)(0), 1);

		grid.remove(pos1, 0);

		BOOST_CHECK_EQUAL(grid(pos1), 1.0f);
		BOOST_CHECK_EQUAL(children.list(0).size(), 1);
		BOOST_CHECK_EQUAL(children(part1).list(0).size(), 0);
		BOOST_CHECK_EQUAL(children(part1).lookup()(pos1), CHILD_NULL_IDX);
		BOOST_CHECK_EQUAL(lookup(part1), BRANCH_NULL_IDX);

		grid.remove(pos3, 0);

		for (UINT i = 0; i < 3; i++)
			BOOST_CHECK_EQUAL(children.list(i).size(), 0);

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


BOOST_AUTO_TEST_SUITE(test_LazySingleTrackedPartitionedGrid)

	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazySingleTrackedPartitionedGrid initialisation]
		// ==== Setup ====
		using GridType = LazySingleTrackedPartitionedGrid<FLOAT, 3, 3>;
		const UINT NULL_IDX = GridType::Child::Lookup::NULL_IDX;

		// ==== Action ====
		GridType grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 7.0f, Vec3u(3, 3, 3));

		// ==== Confirm ====
		BOOST_CHECK_EQUAL((bool)grid.children().get(Vec3i(1,1,1)).is_active(), false);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).background(), 7);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).data().size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)), 7);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).lookup().data().size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).lookup().get(Vec3i(1,1,1)), NULL_IDX);
		/// [LazySingleTrackedPartitionedGrid initialisation]
	}

	struct Fixture
	{
		using GridType = LazySingleTrackedPartitionedGrid<FLOAT, 3, 3>;
		const UINT NULL_IDX = GridType::Lookup::NULL_IDX;
		GridType grid;
		Fixture()
			: grid(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 7, Vec3u(3, 3, 3))
		{}
	};

	BOOST_FIXTURE_TEST_CASE(add_should_activate_once, Fixture)
	{
		// ==== Setup ====
		const Vec3i pos1(-4, -4, -4);
		const Vec3i pos2(-3, -4, -4);
		const Vec3i pos_child(-1, -1, -1);

		// ==== Action ====
		grid.add(pos1, 3, 0);
		grid.add(pos2, 4, 1);

		// ==== Confirm ====
		BOOST_CHECK(grid.children().get(pos_child).is_active());
		BOOST_CHECK_EQUAL(grid.children().list().size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos1), 3);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).get(pos2), 4);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().list(0).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos1), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos2), 0);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).get(Vec3i(1,1,1)), 7);
		BOOST_CHECK_EQUAL(grid.children().get(Vec3i(1,1,1)).lookup().get(Vec3i(1,1,1)), NULL_IDX);
	}

	BOOST_FIXTURE_TEST_CASE(remove_should_deactivate_when_child_is_inactive, Fixture)
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
		BOOST_CHECK(grid.children().get(pos_child).is_active());
		BOOST_CHECK(grid.children().get(pos_child).lookup().is_active());
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().list(1).size(), 1);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos1), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos2), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(1).size(), 1);

		// ==== Action ====
		grid.remove(pos2, 1);

		// ==== Confirm ====
		BOOST_CHECK(!grid.children().get(pos_child).is_active());
		BOOST_CHECK(!grid.children().get(pos_child).lookup().is_active());
		BOOST_CHECK_EQUAL(grid.children().list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().list(1).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos1), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos2), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(0).size(), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).list(1).size(), 0);
	}


	struct ResetFixture : Fixture {
		PartitionedGrid<FLOAT, 3> grid_master;
		ResetFixture()
		: Fixture(),
		  grid_master(Vec3u(9, 9, 9), Vec3i(-4,-4,-4), 0, Vec3u(3, 3, 3))
		{}
	};

	BOOST_FIXTURE_TEST_CASE(reset_should_deactivate, ResetFixture)
	{
		// ==== Setup ====

		const Vec3i pos_child(-1, -1, -1);
		const Vec3i pos(-4, -4, -4);
		grid.add(pos, 4, 0);

		// ==== Action ====

		grid.reset(grid_master, 0);

		// ==== Confirm ====

		// Value reset.
		BOOST_CHECK_EQUAL(grid.get(pos), 7);
		// Child still tracked.
		BOOST_CHECK_EQUAL(
			grid.children().lookup().get(pos_child), Vec3u(NULL_IDX, NULL_IDX, NULL_IDX)
		);
		// Child inactive.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).is_active(), false);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).data().size(), 0);
		// Child lookup inactive.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().is_active(), false);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().data().size(), 0);
		// Position no longer tracked in child.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().list(0).size(), 0);
	}

	BOOST_FIXTURE_TEST_CASE(reset_shouldnt_deactivate_when_other_list_still_active, ResetFixture)
	{
		// ==== Setup ====

		const Vec3i pos_child(-1, -1, -1);
		const Vec3i pos(-4, -4, -4);
		grid.add(pos, 4, 0);

		// ==== Action ====

		grid.reset(grid_master, 1);

		// ==== Confirm ====

		// Value unchanged.
		BOOST_CHECK_EQUAL(grid.get(pos), 4);
		// Child still tracked.
		BOOST_CHECK_EQUAL(grid.children().lookup().get(pos_child), Vec3u(0, NULL_IDX, NULL_IDX));
		// Child still active.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).data().size(), 3*3*3);
		// Child lookup still active.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().data().size(), 3*3*3);
		// Position still tracked in child.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos), 0);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().list(0).size(), 1);
	}

	BOOST_FIXTURE_TEST_CASE(reset_shouldnt_deactivate_when_master_grid_is_tracking, ResetFixture)
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
		BOOST_CHECK_EQUAL(grid.get(pos), 7);
		// Child still tracked.
		BOOST_CHECK_EQUAL(grid.children().lookup().get(pos_child), Vec3u(0, NULL_IDX, NULL_IDX));
		// Child still active.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).data().size(), 3*3*3);
		// Child lookup still active.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().is_active(), true);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().data().size(), 3*3*3);
		// Position no longer tracked in child.
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().get(pos), NULL_IDX);
		BOOST_CHECK_EQUAL(grid.children().get(pos_child).lookup().list(0).size(), 0);
	}
BOOST_AUTO_TEST_SUITE_END()

