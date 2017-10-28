#include <Felt/Impl/Poly.hpp>
#include <Felt/Polys.hpp>
#include <Felt/Surface.hpp>

#include "Utils.hpp"
#include "catch.hpp"

using namespace Felt;


SCENARIO("Impl::Poly::Single")
{

GIVEN("an empty 3D polygonisaton and a 9x9x9 3-layer surface with 3x3x3 partitions")
{
	using Surface = Surface<3, 3>;
	using IsoGrid = typename Surface::IsoGrid;
	using IsoChild = typename IsoGrid::Child;
	using Poly = Impl::Poly::Single<IsoGrid>;

	Surface surface{Vec3i{9,9,9}, Vec3i{3,3,3}};
	Poly poly{surface.isogrid()};

	THEN("poly is initially empty")
	{
		CHECK(poly.vtxs().size() == 0);
		CHECK(poly.spxs().size() == 0);
	}

	WHEN("poly is sized to cover central partition and activated")
	{
		const IsoChild& isochild = surface.isogrid().children().get(Vec3i{0,0,0});

		poly.resize(isochild.size(), isochild.offset());
		poly.bind(isochild.lookup());
		poly.activate();

		AND_WHEN("attempting to polygonise when no surface has been constructed")
		{
			poly.march();

			THEN("poly is still empty")
			{
				CHECK(poly.vtxs().size() == 0);
				CHECK(poly.spxs().size() == 0);
			}
		}

		AND_WHEN("surface is seeded and expanded slightly")
		{
			surface.seed(Vec3i(0,0,0));
			surface.update([](const auto& pos_, const auto& grid_) {
				(void)pos_; (void)grid_;
				return -0.4f;
			});

			THEN("poly is still empty")
			{
				CHECK(poly.vtxs().size() == 0);
				CHECK(poly.spxs().size() == 0);
			}

			AND_WHEN("partial isogrid is polygonised")
			{
				poly.march();

				THEN("number of vertices is as expected")
				{
					CHECK(poly.vtxs().size() == 6);
				}

				THEN("number of simplices is as expected")
				{
					CHECK(poly.spxs().size() == 8);
				}

				THEN("vertices are correct")
				{
					CHECK(poly.vtxs()[0].pos == ApproxVec(Vec3f{0.4f, 0.0f, 0.0f}));
					CHECK(poly.vtxs()[1].pos == ApproxVec(Vec3f{0.0f, 0.0f, 0.4f}));
					CHECK(poly.vtxs()[2].pos == ApproxVec(Vec3f{0.0f, 0.4f, 0.0f}));
					CHECK(poly.vtxs()[3].pos == ApproxVec(Vec3f{-0.4f, 0.0f, 0.0f}));
					CHECK(poly.vtxs()[4].pos == ApproxVec(Vec3f{0.0f, -0.4f, 0.0f}));
					CHECK(poly.vtxs()[5].pos == ApproxVec(Vec3f{0.0f, 0.0f, -0.4f}));

					CHECK(poly.vtxs()[0].norm == ApproxVec(Vec3f{1, 0, 0}));
					CHECK(poly.vtxs()[1].norm == ApproxVec(Vec3f{0, 0, 1}));
					CHECK(poly.vtxs()[2].norm == ApproxVec(Vec3f{0, 1, 0}));
					CHECK(poly.vtxs()[3].norm == ApproxVec(Vec3f{-1, 0, 0}));
					CHECK(poly.vtxs()[4].norm == ApproxVec(Vec3f{0, -1, 0}));
					CHECK(poly.vtxs()[5].norm == ApproxVec(Vec3f{0, 0, -1}));
				}

				THEN("simplices are correct")
				{
					CHECK(poly.spxs()[0].idxs == (Vec3u{1,0,2}));
					CHECK(poly.spxs()[1].idxs == (Vec3u{1,2,3}));
					CHECK(poly.spxs()[2].idxs == (Vec3u{1,4,0}));
					CHECK(poly.spxs()[3].idxs == (Vec3u{0,5,2}));
					CHECK(poly.spxs()[4].idxs == (Vec3u{4,1,3}));
					CHECK(poly.spxs()[5].idxs == (Vec3u{3,2,5}));
					CHECK(poly.spxs()[6].idxs == (Vec3u{0,4,5}));
					CHECK(poly.spxs()[7].idxs == (Vec3u{4,3,5}));
				}

				AND_WHEN("poly is deactivated")
				{
					poly.deactivate();

					THEN("poly is empty and deallocated")
					{
						CHECK(poly.vtxs().capacity() == 0);
						CHECK(poly.spxs().capacity() == 0);
					}
				}

				AND_WHEN("poly is reset")
				{
					poly.reset();

					THEN("poly is empty but not deallocated")
					{
						CHECK(poly.vtxs().size() == 0);
						CHECK(poly.vtxs().size() == 0);
						CHECK(poly.spxs().capacity() > 0);
						CHECK(poly.spxs().capacity() > 0);
					}
				}
			}
		}
	}
}

GIVEN("an empty 2D polygonisaton and a 9x9 3-layer surface with 3x3 partitions")
{
	using Surface = Surface<2, 3>;
	using IsoGrid = typename Surface::IsoGrid;
	using IsoChild = typename IsoGrid::Child;
	using Poly = Impl::Poly::Single<IsoGrid>;

	Surface surface{Vec2i{9,9}, Vec2i{3,3}};
	Poly poly{surface.isogrid()};

	THEN("poly is initially empty")
	{
		CHECK(poly.vtxs().size() == 0);
		CHECK(poly.spxs().size() == 0);
	}

	WHEN("poly is sized to cover central partition and activated")
	{
		const IsoChild& isochild = surface.isogrid().children().get(Vec2i{0,0});

		poly.resize(isochild.size(), isochild.offset());
		poly.bind(isochild.lookup());
		poly.activate();

		AND_WHEN("attempting to polygonise when no surface has been constructed")
		{
			poly.march();

			THEN("poly is still empty")
			{
				CHECK(poly.vtxs().size() == 0);
				CHECK(poly.spxs().size() == 0);
			}
		}

		AND_WHEN("surface is seeded and expanded slightly")
		{
			surface.seed(Vec2i{0,0});
			surface.update([](const auto& pos_, const auto& grid_) {
				(void)pos_; (void)grid_;
				return -0.4f;
			});

			THEN("poly is still empty")
			{
				CHECK(poly.vtxs().size() == 0);
				CHECK(poly.spxs().size() == 0);
			}

			AND_WHEN("partial isogrid is polygonised")
			{
				poly.march();

				THEN("number of vertices is as expected")
				{
					CHECK(poly.vtxs().size() == 4);
				}

				THEN("number of simplices is as expected")
				{
					CHECK(poly.spxs().size() == 4);
				}

				THEN("vertices are correct")
				{
					CHECK(poly.vtxs()[0].pos == ApproxVec(Vec2f{0.4f, 0.0f}));
					CHECK(poly.vtxs()[1].pos == ApproxVec(Vec2f{0.0f, 0.4f}));
					CHECK(poly.vtxs()[2].pos == ApproxVec(Vec2f{-0.4f, 0.0f}));
					CHECK(poly.vtxs()[3].pos == ApproxVec(Vec2f{0.0f, -0.4f}));
				}

				THEN("simplices are correct")
				{
					CHECK(poly.spxs()[0].idxs == (Vec2u{0,1}));
					CHECK(poly.spxs()[1].idxs == (Vec2u{1,2}));
					CHECK(poly.spxs()[2].idxs == (Vec2u{3,0}));
					CHECK(poly.spxs()[3].idxs == (Vec2u{2,3}));
				}

				AND_WHEN("poly is deactivated")
				{
					poly.deactivate();

					THEN("poly is empty and deallocated")
					{
						CHECK(poly.is_active() == false);
						CHECK(poly.vtxs().capacity() == 0);
						CHECK(poly.spxs().capacity() == 0);
					}
				}

				AND_WHEN("poly is reset")
				{
					poly.reset();

					THEN("poly is empty but not deallocated")
					{
						CHECK(poly.is_active() == true);
						CHECK(poly.vtxs().size() == 0);
						CHECK(poly.vtxs().size() == 0);
						CHECK(poly.spxs().capacity() > 0);
						CHECK(poly.spxs().capacity() > 0);
					}
				}
			}
		}
	}
}
}

/**
 * Utility: assert Poly::Grid matches simple Poly::Single polygonisation.
 *
 * Forward declaration.
 */
template <class Surface>
ListIdx assert_partitioned_matches_baseline (
	const Polys<Surface>& polys_,
	const Impl::Poly::Single<typename Surface::IsoGrid>& poly_
);

/**
 * Utility: construct a baseline Poly::Single for testing Poly::Grid against.
 * @param surface
 * @return
 */
template <class TSurface>
Impl::Poly::Single<typename TSurface::IsoGrid> baseline_poly(const TSurface& surface_);


SCENARIO("Polys")
{

GIVEN("an empty 3D polygonisaton and a 15x15x15 3-layer surface with 5x5x5 partitions")
{
	using Surface = Surface<3, 3>;
	using PolyGrid = Polys<Surface>;
	using IsoGrid = typename Surface::IsoGrid;
	using Poly = Impl::Poly::Single<IsoGrid>;

	// Surface to polygonise.
	Surface surface{Vec3i{15,15,15}, Vec3i{5,5,5}};

	// The Poly::Grid to test.
	PolyGrid polys{surface};

	THEN("grid has a matching number of children polys to the isogrid")
	{
		CHECK(polys.children().data().size() == surface.isogrid().children().data().size());
	}

	THEN("child poly size is one greater than the isogrid child size")
	{
		const Vec3i& one = Vec3i::Constant(1);
		const Vec3i& two = Vec3i::Constant(2);

		CHECK(
			polys.children().get(Vec3i(0,0,0)).size() == (
				surface.isogrid().children().get(Vec3i(0,0,0)).size() + two
			)
		);
		CHECK(
			polys.children().get(Vec3i(0,0,0)).offset() == (
				surface.isogrid().children().get(Vec3i(0,0,0)).offset() - one
			)
		);

		CHECK(
			polys.children().get(Vec3i(-1,-1,-1)).size() == (
				surface.isogrid().children().get(Vec3i(-1,-1,-1)).size() + two
			)
		);
		CHECK(
			polys.children().get(Vec3i(-1,-1,-1)).offset() == (
				surface.isogrid().children().get(Vec3i(-1,-1,-1)).offset() - one
			)
		);
	}

	THEN("child polys are inactive")
	{
		CHECK(polys.children().get(Vec3i(0,0,0)).is_active() == false);
	}

	THEN("child polys are bound to the correct isogrid child lookup")
	{
		CHECK(
			polys.children().get(Vec3i(0,0,0)).bind() ==
				&surface.isogrid().children().get(Vec3i(0,0,0)).lookup()
		);
		CHECK(
			polys.children().get(Vec3i(-1,-1,-1)).bind() ==
				&surface.isogrid().children().get(Vec3i(-1,-1,-1)).lookup()
		);
	}

	WHEN("surface is seeded and expanded")
	{
		surface.seed(Vec3i(0,0,0));
		surface.update([](const auto&, const auto&){ return -0.5f; });
		INFO(stringify_grid_slice(surface.isogrid()));
		surface.update([](const auto&, const auto&){ return -0.5f; });
		INFO(stringify_grid_slice(surface.isogrid()));

		AND_WHEN("surface is polygonised")
		{
			polys.notify();
			polys.march();

			THEN("central partition has correct number of vertices and simplices")
			{
				CHECK(polys.children().get(Vec3i(0,0,0)).vtxs().size() == 30);
				CHECK(polys.children().get(Vec3i(0,0,0)).spxs().size() == 56);
			}

			AND_WHEN("surface is contracted and expanded back to how it was then polygonised")
			{
				surface.update([](const auto&, const auto&){ return 1.0f; });
				polys.notify();
				surface.update([](const auto&, const auto&){ return -1.0f; });
				polys.notify();
				polys.march();

				THEN("central partition has correct number of vertices and simplices")
				{
					CHECK(polys.children().get(Vec3i(0,0,0)).vtxs().size() == 30);
					CHECK(polys.children().get(Vec3i(0,0,0)).spxs().size() == 56);
				}
			}

			THEN("we can get a list of the updated partitions")
			{
				const std::set<PosIdx> pos_idxs_expected {
					polys.children().index(Vec3i{ 0, 0, -1}),
					polys.children().index(Vec3i{ 0,-1, 0}),
					polys.children().index(Vec3i{-1, 0, 0}),
					polys.children().index(Vec3i{ 0, 0, 0}),
					polys.children().index(Vec3i{ 0, 0, 1}),
					polys.children().index(Vec3i{ 0, 1, 0}),
					polys.children().index(Vec3i{ 1, 0, 0})
				};
				const std::set<PosIdx> pos_idxs_changed{
					polys.changes().begin(), polys.changes().end()
				};

				CHECK(pos_idxs_changed == pos_idxs_expected);
			}

			AND_WHEN("one point is modified and poly is notified")
			{
				surface.update_start();
				surface.delta(Vec3i(0,1,0), -0.3f);
				surface.update_end();
				polys.notify();

				THEN("list of the updated partitions still hasn't changed")
				{
					const std::set<PosIdx> pos_idxs_expected {
						polys.children().index(Vec3i{ 0, 0,-1}),
						polys.children().index(Vec3i{ 0,-1, 0}),
						polys.children().index(Vec3i{-1, 0, 0}),
						polys.children().index(Vec3i{ 0, 0, 0}),
						polys.children().index(Vec3i{ 0, 0, 1}),
						polys.children().index(Vec3i{ 0, 1, 0}),
						polys.children().index(Vec3i{ 1, 0, 0})
					};
					const std::set<PosIdx> pos_idxs_changed{
						polys.changes().begin(), polys.changes().end()
					};

					CHECK(pos_idxs_changed.size() == pos_idxs_expected.size());
					CHECK(pos_idxs_changed == pos_idxs_expected);
				}

				AND_WHEN("surface is polygonised")
				{
					polys.march();

//					Impl::Grid::Simple<bool, 3> grid_changed{
//						surface.isogrid().children().size(),
//						surface.isogrid().children().offset(),
//						false
//					};
//					std::stringstream ss;
//					for (const PosIdx pos_idx : polys.changes())
//					{
//						grid_changed.set(pos_idx, true);
//						ss << polys.children().index(pos_idx) << std::endl;
//					}
//					INFO(ss.str());
//					INFO(stringify_grid_slice(grid_changed));

					THEN("list of updated partitions has now changed")
					{
						const std::set<PosIdx> pos_idxs_expected {
							polys.children().index(Vec3i{-1, 0, 0}),
							polys.children().index(Vec3i{ 0, 0, 0}),
							polys.children().index(Vec3i{ 0, 1, 0}),
							polys.children().index(Vec3i{ 1, 0, 0}),
							polys.children().index(Vec3i{ 0, 0,-1}),
							polys.children().index(Vec3i{ 0, 0, 1})
						};

						const std::set<PosIdx> pos_idxs_changed{
							polys.changes().begin(), polys.changes().end()
						};

						CHECK(pos_idxs_changed.size() == pos_idxs_expected.size());
						CHECK(pos_idxs_changed == pos_idxs_expected);
					}
				}
			}

			AND_WHEN("surface is expanded and polygonised")
			{
				surface.update([](const auto&, const auto&){ return -0.5f; });
				polys.notify();
				surface.update([](const auto&, const auto&){ return -0.5f; });
				polys.notify();
				polys.march();

				THEN("poly grid matches single poly of whole surface")
				{
					const Poly& poly = baseline_poly(surface);

					const ListIdx total_vtx = assert_partitioned_matches_baseline(polys, poly);

					// Total vertices will have duplicates at the border of the spatial
					// partitions.
					// The 'tip' of the shape at the three lowest corners (5 vertices making
					// up a pyramid) will be outside the central partition. The central
					// partition will thus have three points missing, one at each extremity,
					// since they fall entirely outside the partition.  Thus 4x4 = 12
					// vertices are duplicates of another 12 across the partition lines.
					// So, 12 duplicates + 3 end points - 3 cut from the central partition.
					CHECK(total_vtx == poly.vtxs().size() + 12 + 3 - 3);
					// As mentioned above, each lower extremity non-central partition has 5
					// vertices, making up the endpoint pyramids at those extremities.
					CHECK(polys.children().get(Vec3i(-1,0,0)).vtxs().size() == 5);
					CHECK(polys.children().get(Vec3i(0,-1,0)).vtxs().size() == 5);
					CHECK(polys.children().get(Vec3i(0,0,-1)).vtxs().size() == 5);
				}


				AND_WHEN("surface is contracted and polygonised")
				{
					surface.update([](const auto&, const auto&){ return 1.0f; });
					polys.notify();
					polys.march();

					THEN("poly grid matches single poly of whole surface")
					{
						const ListIdx total_vtx =
							assert_partitioned_matches_baseline(polys, baseline_poly(surface));
						CHECK(total_vtx == 30);
						CHECK(polys.children().get(Vec3i(0,0,0)).vtxs().size() == 30);
						CHECK(polys.children().get(Vec3i(0,0,0)).spxs().size() == 56);
					}

					THEN("poly has the same childs active as the isogrid")
					{
						for (
							PosIdx pos_idx_child = 0;
							pos_idx_child < polys.children().data().size(); pos_idx_child++
						) {
							const Vec3i& pos_child = polys.children().index(pos_idx_child);

							INFO("Check if child " + Felt::format(pos_child) + " should be active");
							CHECK(
								polys.children().get(pos_idx_child).is_active() ==
									surface.isogrid().children().get(pos_idx_child).is_active()
							);
						}
					}

					AND_WHEN("surface is contracted to destruction and polygonised")
					{
						surface.update([](const auto&, const auto&){ return 1.0f; });
						polys.notify();
						surface.update([](const auto&, const auto&){ return 1.0f; });
						polys.notify();
						polys.march();

						THEN("poly grid matches single poly of whole surface")
						{
							const ListIdx total_vtx =
								assert_partitioned_matches_baseline(polys, baseline_poly(surface));
							CHECK(total_vtx == 0);
							CHECK(polys.children().get(Vec3i(0,0,0)).vtxs().size() == 0);
							CHECK(polys.children().get(Vec3i(0,0,0)).spxs().size() == 0);
						}

						THEN("poly has the same childs active as the isogrid (i.e. none)")
						{
							for (
								PosIdx pos_idx_child = 0;
								pos_idx_child < polys.children().data().size(); pos_idx_child++
							) {
								const Vec3i& pos_child = polys.children().index(pos_idx_child);

								INFO(
									"Check if child " + Felt::format(pos_child) +
									" should be active"
								);
								CHECK(
									polys.children().get(pos_idx_child).is_active() ==
										surface.isogrid().children().get(pos_idx_child).is_active()
								);
							}
						}
					}
				} // End WHEN surface is contracted and polygonised

//				AND_WHEN("poly is reset")
//				{
//					poly.reset();
//
//					THEN("child polys are deactivated")
//					{
//						for (const Poly& child : poly.children().data())
//						{
//							CHECK(child.vtx().size() == 0);
//							CHECK(child.spx().size() == 0);
//							CHECK(child.is_active() == false);
//						}
//
//						CHECK(poly.changes().list().size() == 0);
//					}
//				}
			} // End WHEN surface is expanded and polygonised.

			AND_WHEN(
				"surface is expanded with one 'tip' pushed back into central partition, then"
				" polygonised"
			) {
				surface.update([](const auto&, const auto&){ return -1.0f; });
				polys.notify();
				surface.update([](const auto&, const auto&){ return -0.3f; });
				polys.notify();
				surface.update_start();
				surface.delta(Vec3i(0,-2,0), 1.0f);
				surface.update_end();
				polys.notify();

				polys.march();

				THEN("poly grid matches single poly of whole surface")
				{
					const Poly& poly = baseline_poly(surface);
					const ListIdx total_vtx =
						assert_partitioned_matches_baseline(polys, poly);

					// One of the 'tips' have been pushed back into the central partition,
					// So, now just 8 duplicates + 2 endpoints - 2 cut from the central
					// partition.
					CHECK(total_vtx == poly.vtxs().size() + 8 + 2 - 2);
				}
			}
		}

		// Failed originally because of std::vector reinitialisation invalidating references
		// during march (in v1).
		AND_WHEN("notify + expand + notify + march")
		{
			polys.notify();
			surface.update([](const auto&, const auto&){ return -1.0f; });
			polys.notify();

			polys.march();

			THEN("poly grid matches single poly of whole surface")
			{
				const Poly& poly = baseline_poly(surface);
				const ListIdx total_vtx =
					assert_partitioned_matches_baseline(polys, poly);

				CHECK(total_vtx == poly.vtxs().size() + 12);
			}
		}

		AND_WHEN("a point is modified without notifying poly")
		{
			surface.update_start();
			surface.delta(Vec3i(0,-1,0), -1.0f);
			surface.update_end();
			surface.update_start();
			surface.delta(Vec3i(0,-2,0), -1.0f);
			surface.update_end();
			surface.update_start();
			surface.delta(Vec3i(0,-3,0), 0.3f);
			surface.update_end();

			AND_WHEN("poly is notified and marched")
			{
				polys.notify();
				polys.march();

				THEN("polygonisation has only been done on most recently modified partitions")
				{
					const std::set<PosIdx> pos_idxs_expected {
						polys.children().index(Vec3i{ 0,-1, 0}),
						polys.children().index(Vec3i{ 0,-1,-1}),
						polys.children().index(Vec3i{ 0,-1, 1}),
						polys.children().index(Vec3i{-1,-1, 0}),
						polys.children().index(Vec3i{ 1,-1, 0})
					};
					const std::set<PosIdx> pos_idxs_changed{
						polys.changes().begin(), polys.changes().end()
					};

//					Impl::Grid::Simple<bool, 3> grid_changed{
//						surface.isogrid().children().size(), surface.isogrid().children().offset(),
//						false
//					};
//					for (const PosIdx pos_idx : polys.changes())
//						grid_changed.set(pos_idx, true);
//					INFO(stringify_grid_slice(grid_changed));
/*
      1,    0,    0,
      1,    0,    0,
      1,    0,    0,
 */
					CHECK(pos_idxs_changed.size() == pos_idxs_expected.size());
					CHECK(pos_idxs_changed == pos_idxs_expected);

					Poly poly{surface.isogrid()};
					poly.resize(surface.isogrid().size(), surface.isogrid().offset());
					poly.activate();

					for (const PosIdx pos_idx_child : pos_idxs_expected)
					{
						poly.bind(surface.isogrid().children().get(pos_idx_child).lookup());
						poly.march();
					}

					assert_partitioned_matches_baseline(polys, poly);
					// Should fail:
//					assert_partitioned_matches_baseline(polys, baseline_poly(surface));
				}
			}

			AND_WHEN("poly is invalidated and polygonised")
			{
				polys.invalidate();
				polys.march();

				THEN("polygonisation has been done over whole surface")
				{
					const Poly& poly = baseline_poly(surface);
					const ListIdx total_vtx =
						assert_partitioned_matches_baseline(polys, poly);

					CHECK(total_vtx == poly.vtxs().size() + 12);
				}
			}
		}
	} // End WHEN surface is seeded and expanded.
} // End GIVEN 15x15x15 3-layer surface with 5x5x5 partitions.


GIVEN(
	"a polygonisation of a surface with 16x16x16 isogrid in two 16x8x16 partitions, with two seeds"
	" in separate partitions"
) {
	using Surface = Surface<3, 3>;
	using PolyGrid = Polys<Surface>;

	// Surface to polygonise.
	Surface surface{Vec3i{16,16,16}, Vec3i{16,8,16}};
	// The Poly::Grid to test.
	PolyGrid polys{surface};

	surface.seed(Vec3i{0,-4,0});
	surface.seed(Vec3i{0,2,0});
	surface.update([](const auto&, const auto&) {
		return -1.0f;
	});

	polys.notify();
	polys.march();

	INFO(stringify_grid_slice(surface.isogrid()));
	// Record count of simplices in each spatial partition.
	Impl::Grid::Simple<ListIdx, 3> grid_spx_count{
		surface.isogrid().children().size(), surface.isogrid().children().offset(), 0
	};
	for (PosIdx pos_idx_child = 0; pos_idx_child < polys.children().data().size(); pos_idx_child++)
		grid_spx_count.set(pos_idx_child, polys.children().get(pos_idx_child).spxs().size());


	WHEN("surface is expanded and contracted across partitions, polygonising along the way")
	{
		// Expand - expanding across to other partition.
		surface.update_start();
		surface.delta(Vec3i(0,1,0), -1.0f);
		surface.update_end();
		polys.notify();
		polys.march();

		// Contract.
		surface.update_start();
		surface.delta(Vec3i(0,0,0), 1.0f);
		surface.delta(Vec3i(-1,1,0), 1.0f);
		surface.delta(Vec3i(1,1,0), 1.0f);
		surface.delta(Vec3i(0,1,-1), 1.0f);
		surface.delta(Vec3i(0,1,1), 1.0f);
		surface.update_end();
		polys.notify();
		polys.march();

		THEN("poly grid matches baseline poly")
		{
			assert_partitioned_matches_baseline(polys, baseline_poly(surface));
		}

		THEN("number of simplices is unchanged from before expand/contract")
		{
			for (
				PosIdx pos_idx_child = 0; pos_idx_child < polys.children().data().size();
				pos_idx_child++
			) {
				CHECK(
					polys.children().get(pos_idx_child).spxs().size() ==
						grid_spx_count.get(pos_idx_child)
				);
			}
		}
	}
}

}


/**
 * Utility: construct a baseline Poly::Single for testing Poly::Grid against.
 * @param surface
 * @return
 */
template <class TSurface>
Impl::Poly::Single<typename TSurface::IsoGrid> baseline_poly(const TSurface& surface_)
{
	using Surface = TSurface;
	using IsoGrid = typename Surface::IsoGrid;
	using IsoChild = typename IsoGrid::Child;
	using Poly = Impl::Poly::Single<IsoGrid>;

	// Create a Poly::Single to encompass whole isogrid, for checking Poly::Grid against.
	Poly poly{surface_.isogrid()};
	poly.resize(surface_.isogrid().size(), surface_.isogrid().offset());
	poly.activate();

	// Loop every child spatial partition, bind child and polygonise,
	// resulting in one big polygonisation.
	for (const IsoChild& isochild : surface_.isogrid().children().data())
	{
		poly.bind(isochild.lookup());
		poly.march();
	}

	return poly;
}


/**
 * Utility: assert PolyGrid matches simple Poly polygonisation.
 */
template <class TSurface>
ListIdx assert_partitioned_matches_baseline (
	const Polys<TSurface>& polys_,
	const Impl::Poly::Single<typename TSurface::IsoGrid>& poly_
) {
	using IsoGrid = typename TSurface::IsoGrid;
	using Poly = Impl::Poly::Single<IsoGrid>;
	using Simplex = typename Poly::Simplex;
	static const Dim dims = Impl::Traits<IsoGrid>::t_dims;


	ListIdx total_vtx = 0;
	ListIdx total_spx = 0;
	for (const Poly& child : polys_.children().data())
	{
		total_vtx += child.vtxs().size();
		total_spx += child.spxs().size();
		const VecDi<dims>& pos_part_start = child.offset();
		const VecDi<dims>& pos_part_end = child.offset() + child.size();

		if (child.vtxs().size() > 0)
			INFO(
				"Partition "
				+ Felt::format(pos_part_start) + "-" + Felt::format(pos_part_end)
				+ " vtxs = "
				+ std::to_string(child.vtxs().size())
				+ ", spxs = "
				+ std::to_string(child.spxs().size())
			);

		for (const Simplex& polys_spx : child.spxs())
		{
			Vec3f polys_vtxs[3];
			polys_vtxs[0] = child.vtxs()[polys_spx.idxs(0)].pos;
			polys_vtxs[1] = child.vtxs()[polys_spx.idxs(1)].pos;
			polys_vtxs[2] = child.vtxs()[polys_spx.idxs(2)].pos;
			auto it = std::find_if(
				poly_.spxs().begin(), poly_.spxs().end(),
				[&](const Simplex& poly_spx) {
					Vec3f poly_vtxs[3];
					poly_vtxs[0] = poly_.vtxs()[poly_spx.idxs(0)].pos;
					poly_vtxs[1] = poly_.vtxs()[poly_spx.idxs(1)].pos;
					poly_vtxs[2] = poly_.vtxs()[poly_spx.idxs(2)].pos;
					return (
						poly_vtxs[0] == polys_vtxs[0] &&
						poly_vtxs[1] == polys_vtxs[1] &&
						poly_vtxs[2] == polys_vtxs[2]
					);
				}
			);

			INFO(
				"Simplex from partition "
				+ Felt::format(polys_vtxs[0]) + "-"
				+ Felt::format(polys_vtxs[1]) + "-"
				+ Felt::format(polys_vtxs[2])
				+ " found in baseline"
			);
			CHECK(it != child.spxs().end());
		}
	}


	for (const Simplex& poly_spx : poly_.spxs())
	{
		Vec3f poly_vtxs[3];
		poly_vtxs[0] = poly_.vtxs()[poly_spx.idxs(0)].pos;
		poly_vtxs[1] = poly_.vtxs()[poly_spx.idxs(1)].pos;
		poly_vtxs[2] = poly_.vtxs()[poly_spx.idxs(2)].pos;

		bool found_match = false;

		for (const Poly& child : polys_.children().data())
		{
			for (const Simplex& polys_spx : child.spxs())
			{
				Vec3f polys_vtxs[3];
				polys_vtxs[0] = child.vtxs()[polys_spx.idxs(0)].pos;
				polys_vtxs[1] = child.vtxs()[polys_spx.idxs(1)].pos;
				polys_vtxs[2] = child.vtxs()[polys_spx.idxs(2)].pos;

				found_match = (
					poly_vtxs[0] == polys_vtxs[0] &&
					poly_vtxs[1] == polys_vtxs[1] &&
					poly_vtxs[2] == polys_vtxs[2]
				);
				if (found_match)
					break;
			}
			if (found_match)
				break;
		}

		INFO(
			"Simplex from baseline "
			+ Felt::format(poly_vtxs[0]) + "-"
			+ Felt::format(poly_vtxs[1]) + "-"
			+ Felt::format(poly_vtxs[2])
			+ " found in partition"
		);
		CHECK(found_match);
	}

	INFO("Total: " + std::to_string(total_spx) + " spxs");
	INFO("Total: " + std::to_string(total_vtx) + " vtxs");

	CHECK(total_spx == poly_.spxs().size());

	return total_vtx;
}

