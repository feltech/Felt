#include <Felt/Impl/Poly.hpp>

#include "../catch.hpp"

#include "Utils.hpp"
#include "Felt/Impl/Surface.hpp"


using namespace Felt;

SCENARIO("Impl::Poly::Single")
{

GIVEN("an empty 3D polygonisaton and a 9x9x9 3-layer surface with 3x3x3 partitions")
{
	using SurfaceType = Surface<3, 3>;
	using IsoGridType = typename SurfaceType::IsoGrid;
	using IsoChildType = typename IsoGridType::ChildType;
	using PolyType = Impl::Poly::Single<IsoGridType>;

	SurfaceType surface{Vec3i{9,9,9}, Vec3i{3,3,3}};
	PolyType poly{surface.isogrid()};

	THEN("poly is initially empty")
	{
		CHECK(poly.vtxs().size() == 0);
		CHECK(poly.spxs().size() == 0);
	}

	WHEN("poly is sized to cover central partition and activated")
	{
		const IsoChildType& isochild = surface.isogrid().children().get(Vec3i{0,0,0});

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
	using SurfaceType = Surface<2, 3>;
	using IsoGridType = typename SurfaceType::IsoGrid;
	using IsoChildType = typename IsoGridType::ChildType;
	using PolyType = Impl::Poly::Single<IsoGridType>;

	SurfaceType surface{Vec2i{9,9}, Vec2i{3,3}};
	PolyType poly{surface.isogrid()};

	THEN("poly is initially empty")
	{
		CHECK(poly.vtxs().size() == 0);
		CHECK(poly.spxs().size() == 0);
	}

	WHEN("poly is sized to cover central partition and activated")
	{
		const IsoChildType& isochild = surface.isogrid().children().get(Vec2i{0,0});

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
}


SCENARIO("Impl::Poly::Grid")
{
GIVEN("an empty 3D polygonisaton and a 9x9x9 3-layer surface with 3x3x3 partitions")
{
	using SurfaceType = Surface<3, 3>;
	using PolyGridType = Impl::Poly::Grid<SurfaceType>;

	SurfaceType surface{Vec3i{9,9,9}, Vec3i{3,3,3}};

	PolyGridType poly{surface};

	THEN("the grid has a matching number of children polys to the isogrid.")
	{
		CHECK(poly.children().data().size() == 27);
	}
}
}
