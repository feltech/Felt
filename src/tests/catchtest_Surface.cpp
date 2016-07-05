#include <boost/numeric/ublas/io.hpp>
#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include "Felt/Surface.hpp"
#include "Utils.hpp"

using namespace felt;

template <typename... Args>
auto s(Args&&... args) -> decltype(std::to_string(std::forward<Args>(args)...)) {
	return std::to_string(std::forward<Args>(args)...);
}


SCENARIO("Surface")
{
/*
 * Basic initialisation.
 */
WHEN("init")
{
	// ==== Setup ====
	// Basic initialisation of 2D surface with 2 layers in a 5x5 embedding.
	Surface<2, 2> surface(Vec2u(7, 7), Vec2u(3, 3));

	// ==== Confirm ====
	CHECK(surface.isogrid().size() == Vec2u(7, 7));
	CHECK(surface.isogrid().children().data().size() == 9);
	CHECK(surface.isogrid().children().get(Vec2i(0, 0)).size() == Vec2u(3, 3));
	CHECK(surface.isogrid().children().get(Vec2i(0, 0)).data().size() == 0);
	CHECK(surface.isogrid().size() == Vec2u(7, 7));
	// Grid is initialised to all points 'outside' the surface (since there is no surface yet).
	CHECK(surface.isogrid().get(Vec2i(0, 0)) == 3);
}

/*
 * Narrow band layers.
 */
WHEN("layers")
{
	// 3D surface with default (=2) number of layers.
	Surface<3> surface(Vec3u(7, 7, 7));
	Surface<3>::IsoGrid::ChildrenGrid children;
	Vec3i pos = Vec3i(0, 0, 0);

	CHECK(children.list(surface.layer_idx(-2)).size() == 0);
	CHECK(children.list(surface.layer_idx(-1)).size() == 0);
	CHECK(children.list(surface.layer_idx(0)).size() == 0);
	CHECK(children.list(surface.layer_idx(1)).size() == 0);
	CHECK(children.list(surface.layer_idx(2)).size() == 0);

	// Add a single zero-layer point.
	surface.isogrid(pos) = 0;
	surface.layer_add(pos, 0);

	// Check zero-layer array has registered point.
	CHECK(surface.layer(0).size() == 1);
	CHECK((*surface.layer(0).begin() - pos) == Vec3i::Zero());

	// Check layer calculation from value.
	// -- zero-layer point just added.
	CHECK(surface.layer_id(pos) == 0);

	// Add three arbitrary points to layer -1.
	surface.layer_add(Vec3i(0, 0, 1), -1);
	surface.layer_add(Vec3i(0, 0, 2), -1);
	surface.layer_add(Vec3i(0, 0, 3), -1);

	// Remove two points from layer -1
	surface.layer_remove(Vec3i(0, 0, 1), -1);
	CHECK(surface.layer(-1).size() == 2);
	surface.layer_remove(Vec3i(0, 0, 3), -1);
	CHECK(surface.layer(-1).size() == 1);

	// Move a point from layer 0 to layer -1
	surface.layer_move(pos, 0, -1);
	CHECK(surface.layer(-1).size() == 2);

	// Check lists updated.
	CHECK(surface.layer(0).size() == 0);
	CHECK(surface.layer(-1).size() == 2);
}

/*
 * Placing a single singularity point.
 */
GIVEN("a surface with 5x5 isogrid")
{
	Surface<2, 2> surface(Vec2u(5, 5));

	WHEN("we add a singularity seed at the centre")
	{
		surface.seed(Vec2i(0, 0));

		// Trivially check centre of seed is indeed a zero-level point (i.e. point
		// on the surface).
		THEN("the value at the centre of the grid is 0")
		{
			const FLOAT val_centre = surface.isogrid().get(Vec2i(0, 0));
			CHECK(val_centre == 0);
		}

		THEN("the data in the grid matches expected data")
		{
			// A 2D 2-layer singularity (seed) point should look like the following.

			Grid<FLOAT, 2> isogrid_check(Vec2u(5, 5), Vec2i::Zero(), 0);
			isogrid_check.data() = {
				3, 3, 2, 3, 3,	// |
				3, 2, 1, 2, 3,	// -
				2, 1, 0, 1, 2,	// x
				3, 2, 1, 2, 3,	// +
				3, 3, 2, 3, 3	// |
			//	|____ - y + ____|
			};
		//	std::cerr << isogrid.data() << std::endl << std::endl;
		//	std::cerr << isogrid_check.data() << std::endl << std::endl;

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			INFO(stringifyGridSlice(surface.isogrid().snapshot()));
			INFO(stringifyGridSlice(isogrid_check));

			const FLOAT diff = isogrid_check.vdata().sum();

			CHECK(diff == 0);
		}

		THEN("the surface layers track the correct number of points")
		{
			// Check appropriate points have been added to narrow band layers.
			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 0);
			CHECK(surface.layer(0).size() == 1);
			CHECK(surface.layer(1).size() == 4);
			CHECK(surface.layer(2).size() == 8);
		}
	}
}

/*
 * Given a grid point, find neighbouring point closest to zero-curve.
 */
WHEN("next_closest_grid_point")
{
	// Create seed point, as above, and navigate to centre.

	Surface<2, 2> surface(Vec2u(5, 5));
	Surface<2, 2>::IsoGrid& isogrid = surface.isogrid();

	surface.seed(Vec2i(0, 0));

	Vec2i pos_next = Vec2i(-1, -2);
	CHECK(isogrid(pos_next) == 3);

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == 2);

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == 1);

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == 0);

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == 0);

	// Ensure it also works with negative distances.
	// NOTE: row-major (y,x) element ordering...
	surface.isogrid().snapshot().data() = {
		2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -2,
		-2, -2, -2, -2 };
	surface.isogrid().flush_snapshot();
	// NOTE: ...but accessed as (x,y)
	pos_next = Vec2i(2, 0);

	CHECK(isogrid(pos_next) == -2);
	CHECK(pos_next == Vec2i(2, 0));

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == -1);
	CHECK(pos_next == Vec2i(1, 0));

	pos_next = surface.next_closest(pos_next);

	CHECK(isogrid(pos_next) == 0);
	CHECK(pos_next == Vec2i(0, 0));
}

/*
 * Using delta isogrid grid/list.
 */
WHEN("delta_isogrid_clamping")
{
	//! [Delta isogrid clamping]
	// ==== Setup ====
	// Basic non-threaded check.
	Surface<3, 2> surface(Vec3u(5, 5, 5));
	Vec3i pos(0, 0, 0);

	// ==== Action ====
	// Apply a delta to the surface.
	surface.delta(pos, -2);

	// ==== Confirm ====
	// Check delta was stored in underlying grid - will be clamped to -1.
	CHECK(surface.delta(pos) == -1);
	// Check position vector of point in surface grid that delta was
	// applied to is stored in a corresponding list to be iterated over.
	CHECK(surface.delta().leafs(surface.layer_idx(0)).size() == 1);
	//! [Delta isogrid clamping]
}

/*
 * Update isogrid with delta isogrid.
 */
WHEN("delta_isogrid_update")
{
	//! [Simple delta isogrid update]
	// ==== Setup ====
	Surface<3, 2> surface(Vec3u(5, 5, 5));
	Surface<3, 2>::DeltaIsoGrid& delta = surface.delta();

	// ==== Action ====
	// Put in 'dirty' state, to check update_start is doing it's job.
	surface.delta(Vec3i(0, 0, 0), 0.5f);

	// ==== Confirm ====
	CHECK(surface.delta().children().list(surface.layer_idx(0)).size() == 1);
	CHECK(surface.delta().get(Vec3i(0, 0, 0)) == 0.5f);

	// ==== Action ====
	// Clear delta isogrid.
	surface.update_start();

	// ==== Confirm ====
	// Check update_start cleared the above surface.delta changes.
	CHECK(surface.delta().children().list(surface.layer_idx(0)).size() == 0);
	CHECK(surface.delta().get(Vec3i(0, 0, 0)) == 0.0f);

	// ==== Action ====
	// Add a zero-layer point.
	surface.layer_add(0.0f, Vec3i(0, 0, 0));

	// Clear delta isogrid.
	surface.update_start();
	// Do nothing.
	surface.delta(Vec3i(0, 0, 0), 0);
	// Apply delta isogrid.
	surface.update_end();

	// ==== Confirm ====
	// Ensure nothing was changed.  Every point in 5x5x5 grid == 3, except centre which == 0.
	CHECK(surface.isogrid().snapshot().vdata().sum() == 3 * 5 * 5 * 5 - 3);
	// Delta isogrid position vector list should still contain one point.
	CHECK(surface.delta().leafs(surface.layer_idx(0)).size() == 1);
	// Delta isogrid grid itself should have reset back to zero.
	CHECK(delta(Vec3i(0, 0, 0)) == 0);

	// ==== Action ====
	// Clear delta isogrid.
	surface.update_start();
	// Apply small update.
	surface.delta(Vec3i(0, 0, 0), 0.4f);
	// Apply delta isogrid.
	surface.update_end();

	// ==== Confirm ====
	// Ensure change applied.  Every point in grid == 3, except centre which == 0.4.
	CHECK(surface.isogrid().snapshot().vdata().sum() == 3 * 5 * 5 * 5 - 3 + 0.4f);
	CHECK(surface.isogrid().get(Vec3i(0, 0, 0)) == 0.4f);
	//! [Simple delta isogrid update]
}

/*
 * Update signed distance transform of outer layer points.
 */
WHEN("distance_transform")
{
	// Check distance calculation for a single point.
	{
		typedef Surface<3, 2> SurfaceT;
		SurfaceT surface(Vec3u(5, 5, 5));
		SurfaceT::IsoGrid& isogrid = surface.isogrid();

		surface.seed(Vec3i(0, 0, 0));

		// Basic distance calculation.
		isogrid(Vec3i(0, 0, 0)) = -0.6f;
		const FLOAT dist = surface.distance(Vec3i(-1, 0, 0), 1);
		CHECK(dist == Approx(0.4f).epsilon(0.0001f));
	}
	// Update seed point by less than |0.5| and check outer layer
	// distances are updated.
	{
		typedef Surface<2, 2> SurfaceT;
		SurfaceT surface(Vec2u(5, 5), Vec2u(5, 5));

		surface.seed(Vec2i(0, 0));

		Grid<FLOAT, 2> isogrid_check(Vec2u(5, 5), Vec2i::Zero(), 0);
		isogrid_check.data() = {
			3, 3, 1.6, 3, 3, 3, 1.6, 0.6, 1.6, 3, 1.6, 0.6, -0.4, 0.6, 1.6, 3,
			1.6, 0.6, 1.6, 3, 3, 3, 1.6, 3, 3};

		surface.update_start();
		{
			surface.delta(Vec2i(0, 0), -0.4f);
		}
		surface.update_end();

		surface.update_start();
		{
			// Check update_start cleared the above surface.delta changes.
			for (const Vec2i& pos_child : surface.delta().children())
				for (const Vec2i& pos : surface.delta().children().get(pos_child))
					CHECK(surface.delta().get(pos) == 0);

		}
		surface.update_end();


		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();

		CHECK(diff == 0);
	}
}

/*
 * Update layers.
 */
WHEN("layer_update")
{
	typedef Surface<2, 2> Surface_t;
	Surface_t surface(Vec2u(9, 9));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(Vec2u(9, 9), Vec2i::Zero(), 0);

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update_start();
	{
		surface.delta(Vec2i(0, 0), -0.6f);
	}
	surface.update_end();

	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3, 3, 3, 3, 3, 3,
			2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 0.4, 1.4, 2.4, 3, 3, 3,
			2.4, 1.4, 0.4, -0.6, 0.4, 1.4, 2.4, 3, 3, 3, 2.4, 1.4, 0.4, 1.4,
			2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));

		CHECK(surface.layer(0).size() == 4);
		CHECK(surface.layer(-1).size() == 1);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 8);
		CHECK(surface.layer(2).size() == 12);
	}

	// Cycle new zero-layer points and move back to original signed distance.
//	surface.update_start();
//	{
//		for (const Vec2i& pos : surface.layer(0))
//			surface.delta(pos, 0.6f);
//	}
//	surface.update_end();

	// Update using lambda.
	surface.update([](const Vec2i& pos, const Surface_t::IsoGrid& isogrid) {
		return 0.6f;
	});


	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 2, 1, 0, 1, 2, 3,
			3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}

	// Collapse the seed completely, leaving no zero-layer, only outer layers.
	surface.update_start();
	{
		for (const Vec2i& pos : surface.layer(0))
			surface.delta(pos, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3,
			3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		CHECK(surface.layer(0).size() == 0);
		CHECK(surface.layer(-1).size() == 0);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 1);
		CHECK(surface.layer(2).size() == 4);

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}

	// Collapse still further, so there is only the outermost layer.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (const Vec2i& pos : surface.layer(0))
			surface.delta(pos, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		CHECK(surface.layer(0).size() == 0);
		CHECK(surface.layer(-1).size() == 0);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 0);
		CHECK(surface.layer(2).size() == 1);

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}

	// Final collapse leaves the whole grid as 'outside' points.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (const Vec2i& pos : surface.layer(0))
			surface.delta(pos, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		CHECK(surface.layer(0).size() == 0);
		CHECK(surface.layer(-1).size() == 0);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 0);
		CHECK(surface.layer(2).size() == 0);

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}

	// Further updates have no effect.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (const Vec2i& pos : surface.layer(0))
			surface.delta(pos, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.isogrid().data() << std::endl << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

		CHECK(surface.layer(0).size() == 0);
		CHECK(surface.layer(-1).size() == 0);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 0);
		CHECK(surface.layer(2).size() == 0);

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}
}

/*
 * Iterating the zero-layer
 */
WHEN("iterate_layers")
{
	Surface<3> surface(Vec3u(9, 9, 9));

	// Values to compute for testing against.
	int counter;
	Vec3i pos_sum;

	// Create seed point and expand the narrow band.
	surface.seed(Vec3i(0, 0, 0));
	surface.update_start();
	{
		surface.delta(Vec3i(0, 0, 0), -1.0f);
	}
	surface.update_end();

	CHECK(surface.layer(0).size() == 6);

	// Iterate over surface, using partitioned grid.
	// Only version that can be parallelised easily using OpenMP.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	const UINT zeroLayerIdx = surface.layer_idx(0);

	#pragma omp parallel for
	for (UINT part_idx = 0; part_idx < surface.parts().size(); part_idx++)
	{
		const Vec3i& pos_part = surface.parts()[part_idx];
		for (const Vec3i& pos : surface.layer(pos_part, 0))
		{
			FLOAT val = surface(pos);
			#pragma omp critical
			{
				CHECK(val == 0);
				counter++;
				pos_sum += pos;
			}
		}
	};
	CHECK(counter == 6);
	CHECK(pos_sum == Vec3i(0, 0, 0));


	counter = 0;
	pos_sum = Vec3i(0, 0, 0);

	for (
		INT layer_id = surface.LAYER_MIN; layer_id <= surface.LAYER_MAX;
		layer_id++
	)
		for (Vec3i& part : surface.parts(layer_id))
			for (const Vec3i& pos : surface.layer(part, layer_id))
			{
				FLOAT val = surface(pos);
				#pragma omp critical
				{
					CHECK(val == layer_id);
					counter++;
					pos_sum += pos;
				}
			}

	CHECK(counter == 63);
	CHECK(pos_sum == Vec3i(0, 0, 0));

	// Iterate over zero-layer using STL for_each and lambda callback function.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	std::for_each(
		surface.layer(0).begin(), surface.layer(0).end(),
		[&](const Vec3i& pos)
		{
			pos_sum += pos;
			counter++;
		}
	);
	CHECK(counter == 6);
	CHECK(pos_sum == Vec3i(0, 0, 0));

	// Iterate over zero-layer using range based for loop.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	for (auto pos : surface.layer(0))
	{
		pos_sum += pos;
		counter++;
	}

	CHECK(counter == 6);
	CHECK(pos_sum == Vec3i(0, 0, 0));

}

GIVEN("a 9x9 surface initialised with a singulaty point in the centre")
{
	typedef Surface<2, 2> SurfaceT;
	SurfaceT surface(Vec2u(9, 9));
	surface.seed(Vec2i(0, 0));

	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(Vec2u(9, 9), Vec2i::Zero(), 0);

	WHEN("we expand the surface one unit outwards")
	{
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});

		THEN("the grid data matches a surface of radius 1")
		{
			isogrid_check.data() = {
			      3,    3,    3,    3,    3,    3,    3,    3,    3,
			      3,    3,    3,    3,    2,    3,    3,    3,    3,
			      3,    3,    3,    2,    1,    2,    3,    3,    3,
			      3,    3,    2,    1,    0,    1,    2,    3,    3,
			      3,    2,    1,    0,   -1,    0,    1,    2,    3,
			      3,    3,    2,    1,    0,    1,    2,    3,    3,
			      3,    3,    3,    2,    1,    2,    3,    3,    3,
			      3,    3,    3,    3,    2,    3,    3,    3,    3,
			      3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();

			INFO(stringifyGridSlice(surface.isogrid()));
			CHECK(diff == Approx(0));
		}

		AND_WHEN("we expand by one unit again")
		{
			surface.update([](auto& pos, auto& isogrid) {
				return -1.0f;
			});

			THEN("the grid data matches a surface of radius 2")
			{
				isogrid_check.data() = {
				      3,    3,    3,    3,    2,    3,    3,    3,    3,
				      3,    3,    3,    2,    1,    2,    3,    3,    3,
				      3,    3,    2,    1,    0,    1,    2,    3,    3,
				      3,    2,    1,    0,   -1,    0,    1,    2,    3,
				      2,    1,    0,   -1,   -2,   -1,    0,    1,    2,
				      3,    2,    1,    0,   -1,    0,    1,    2,    3,
				      3,    3,    2,    1,    0,    1,    2,    3,    3,
				      3,    3,    3,    2,    1,    2,    3,    3,    3,
				      3,    3,    3,    3,    2,    3,    3,    3,    3
				};

				isogrid_check.vdata() = (
					isogrid_check.vdata() - surface.isogrid().snapshot().vdata()
				);
				const FLOAT diff = isogrid_check.vdata().sum();

				INFO(stringifyGridSlice(surface.isogrid()));
				CHECK(diff == Approx(0));
			}

			AND_WHEN("we expand by one unit 9 more times")
			{
				for (UINT i = 0; i < 9; i++)
					surface.update([](auto& pos, auto& isogrid) {
						return -1.0f;
					});

				THEN("the grid data matches an area completely consumed by the surface")
				{
					isogrid_check.data() = {
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
					     -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3
					};

					isogrid_check.vdata() = (
						isogrid_check.vdata() - surface.isogrid().snapshot().vdata()
					);
					const FLOAT diff = isogrid_check.vdata().sum();

					INFO(stringifyGridSlice(surface.isogrid()));
					CHECK(diff == Approx(0));
				}

				THEN("the layer lists are empty")
				{
					CHECK(surface.layer(0).size() == 0);
					CHECK(surface.layer(-1).size() == 0);
					CHECK(surface.layer(-2).size() == 0);
					CHECK(surface.layer(1).size() == 0);
					CHECK(surface.layer(2).size() == 0);
				}
			}
		}
	}
}


WHEN("affected_outer_layers")
{
	//![Calculate affected outer layers for localised narrow band updates]
	// ==== Setup ====
	using PosArray = Surface<2, 2>::IsoGrid::PosArray;
	Surface<2, 2> surface(Vec2u(9, 9));
	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update_start();
	{
		for (auto pos : surface.layer(0))
		{
			surface.delta(pos, -1.0f);
		}
	}
	surface.update_end();
	// Clean up from previous update.
	surface.update_start();
	// Add a couple of points that could affect the narrow band.
	surface.delta(Vec2i(0, 1), 0.3f);
	surface.delta(Vec2i(1, 0), 0.3f);
	//		3.0,	3.0,	3.0,	 2.0,	3.0,	3.0,	3.0,
	//		3.0,	3.0,	2.0,	 1.0,	2.0,	3.0,	3.0,
	//		3.0,	2.0,	1.0,	 0.0,	1.0,	2.0,	3.0,
	//		2.0,	1.0,	0.0,	-1.0,	0.3,	1.0,	2.0,
	//		3.0,	2.0,	1.0,	 0.3,	1.0,	2.0,	3.0,
	//		3.0,	3.0,	2.0,	 1.0,	2.0,	3.0,	3.0,
	//		3.0,	3.0,	3.0,	 2.0,	3.0,	3.0,	3.0;

	// ==== Action ====

	surface.calc_affected();

	// ==== Confirm ====

	PosArray check_layers_pos[5];
	check_layers_pos[2 + -2] = PosArray();
	check_layers_pos[2 + -1] = PosArray({
		Vec2i(0, 0)
	});
	check_layers_pos[2 + 0] = PosArray({
	// We don't care for now about zero-layer points.
//		Vec2i(0,1),
//		Vec2i(1,0)
	});
	check_layers_pos[2 + 1] = PosArray({
		// For (0,1):
		Vec2i(-1, 1), Vec2i(1, 1), Vec2i(0, 2),
		// For (1,0):
		Vec2i(2, 0), Vec2i(1, -1)
	});

	check_layers_pos[2 + 2] = PosArray({
		// For (0,1):
		Vec2i(-2, 1), Vec2i(2, 1),
		Vec2i(-1, 2), Vec2i(1, 2),
		Vec2i(0, 3),
		// For (1,0):
		Vec2i(3, 0), Vec2i(1, -2), Vec2i(2, -1)
	});

	for (INT layer_id = -2; layer_id <= 2; layer_id++)
	{
		if (layer_id == 0)
			continue;

		const INT layer_idx = 2 + layer_id;

		INFO(
			"Layer " + s(layer_id) + " at index " + s(layer_idx) +
			" number of leafs " +
			s(surface.affected().leafs(layer_idx).size()) + " == " +
			s(check_layers_pos[layer_idx].size())
		);
		CHECK(
			surface.affected().leafs(layer_idx).size() == check_layers_pos[layer_idx].size());

		for (auto pos : check_layers_pos[layer_idx])
		{
			auto iter = std::find(
				surface.affected().leafs(layer_idx).begin(),
				surface.affected().leafs(layer_idx).end(), pos
			);

			INFO(
				"Affected grid layer " + s(layer_id) + " at index " + s(layer_idx) +
				" should contain (" + s(pos(0)) + "," + s(pos(1)) + ")"
			);
			CHECK(
				iter != surface.affected().leafs(layer_idx).end());
		}

		for (auto pos : surface.affected().leafs(layer_idx))
		{
			auto iter = std::find(
				check_layers_pos[layer_idx].begin(),
				check_layers_pos[layer_idx].end(), pos
			);

			INFO(
				"Checking list layer " + s(layer_id) + " at index " + s(layer_idx) +
				" should contain (" + s(pos(0)) + "," + s(pos(1)) + ")"
			);
			CHECK(
				iter != check_layers_pos[layer_idx].end());
		}
	}
	//![Calculate affected outer layers for localised narrow band updates]
}

/*
 * Localised update.
 */
WHEN("local_update")
{
	typedef Surface<2, 2> SurfaceT;
	SurfaceT surface(Vec2u(9, 9));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(Vec2u(9, 9), Vec2i::Zero(), 0);

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
//	3,	3,	3,	3,	3,	3,	3,
//	3,	3,	3,	2,	3,	3,	3,
//	3,	3,	2,	1,	2,	3,	3,
//	3,	2,	1,	0,	1,	2,	3,
//	3,	3,	2,	1,	2,	3,	3,
//	3,	3,	3,	2,	3,	3,	3,
//	3,	3,	3,	3,	3,	3,	3;
	surface.update_start();
	{
		surface.delta(Vec2i(0, 0), -0.6f);
	}
	// Using localised update, which will only update outer layers that are
	// affected by changes to the modified zero layer points.  In this test
	// case, all outer layer points are affected, same as a global update.
	surface.update_end_local();

	{
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3, 3, 3, 3, 3, 3,
			2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 0.4, 1.4, 2.4, 3, 3, 3,
			2.4, 1.4, 0.4, -0.6, 0.4, 1.4, 2.4, 3, 3, 3, 2.4, 1.4, 0.4, 1.4,
			2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));

		CHECK(surface.layer(0).size() == 4);
		CHECK(surface.layer(-1).size() == 1);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(1).size() == 8);
		CHECK(surface.layer(2).size() == 12);
	}

	// Cycle new zero-layer points and move back to original signed distance.
	surface.update_start();
	{
		for (const Vec2i& pos : surface.layer(0))
			surface.delta(pos, 0.6f);
	}
	surface.update_end_local();

	{
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 2, 1, 0, 1, 2, 3,
			3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0).epsilon(0.000001f));
	}
}


WHEN("deactivates_with_inside_background_value")
{
	// ==== Setup ====
	using SurfaceT = Surface<2, 2>;
	Vec2u size(21, 21);
	SurfaceT surface(size, Vec2u(2, 2));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(size, Vec2i::Zero(), 0);

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));

	// ==== Action ====
	for (UINT i = 0; i < 5; i++)
		surface.update([](auto& pos_, auto& grid_) {
			return -1.0f;
		});

	// ==== Confirm ====

	INFO(stringifyGridSlice(surface.isogrid()));
	isogrid_check.data() = {
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,
		3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,
		3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
		3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
	};

	isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
	const FLOAT diff = isogrid_check.vdata().sum();
	CHECK(diff == Approx(0).epsilon(0.000001f));
}


/**
 * Test walking the zero layer out to a given distance.
 */
WHEN("walk")
{
	{
		// ==== Setup ====
		Surface<2, 2> surface(Vec2u(16, 16));

		// Create seed point and expand the narrow band.
		surface.seed(Vec2i(0, 0));
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});

		// ==== Action ====
		SingleLookupGrid<2, surface.NUM_LAYERS> lookup = surface.walk_band<2>(
			Vec2i(-3,0)
		);

		// ==== Contirm ====
		CHECK(lookup.list(surface.layer_idx(-2)).size() == 1);
		CHECK(lookup.list(surface.layer_idx(-1)).size() == 1);
		CHECK(lookup.list(surface.layer_idx(0)).size() == 3);
		CHECK(lookup.list(surface.layer_idx(1)).size() == 3);
		CHECK(lookup.list(surface.layer_idx(2)).size() == 5);

		CHECK(lookup.list(surface.layer_idx(0))[0] == Vec2i(-3, 0));
		CHECK(lookup.list(surface.layer_idx(0))[1] == Vec2i(-2, -1));
		CHECK(lookup.list(surface.layer_idx(0))[2] == Vec2i(-2, 1));
	}
	{
		// ==== Setup ====
		Surface<3, 2> surface(Vec3u(9, 9, 9));
		using MultiLookup = SingleLookupGrid<3, surface.NUM_LAYERS>;

		// Create seed point and expand the narrow band.
		surface.seed(Vec3i(0, 0, 0));

		// ==== Action ====
		MultiLookup& lookup1 = surface.walk_band<1>(Vec3i(0,0,0));

		// ==== Contirm ====
		CHECK(lookup1.list(surface.layer_idx(-2)).size() == 0);
		CHECK(lookup1.list(surface.layer_idx(-1)).size() == 0);
		CHECK(lookup1.list(surface.layer_idx(0)).size() == 1);
		CHECK(lookup1.list(surface.layer_idx(1)).size() == 6);
		CHECK(lookup1.list(surface.layer_idx(2)).size() == 0);

		// ==== Action ====
		MultiLookup& lookup2 = surface.walk_band<1>(Vec3i(0,0,0));

		// ==== Contirm ====
		CHECK(&lookup1 == &lookup2);

		// ==== Action ====
		MultiLookup& lookup3 = surface.walk_band<2>(Vec3i(0,0,0));

		// ==== Contirm ====
		CHECK(&lookup1 == &lookup2);
		CHECK(&lookup1 != &lookup3);
		CHECK(&lookup2 != &lookup3);
	}

	{
		// ==== Setup ====
		Surface<3, 2> surface(Vec3u(16, 16, 16));

		// Create seed point and expand the narrow band.
		surface.seed(Vec3i(0, 0, 0));
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});

//		INFO(stringifyGridSlice(surface.isogrid()));
/*
/home/dave/Dropbox/Workspace/Felt/src/tests/test_Surface.cpp(924): Message:
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -3 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/
		// ==== Action ====
		SingleLookupGrid<3, surface.NUM_LAYERS> lookup = surface.walk_band<1>(
			Vec3i(0,0,0)
		);

		// ==== Contirm ====
		CHECK(lookup.list(surface.layer_idx(-2)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(-1)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(0)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(1)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(2)).size() == 0);

		// ==== Action ====
		lookup = surface.walk_band<2>(Vec3i(-5,0,0));

		// ==== Contirm ====
		CHECK(lookup.list(surface.layer_idx(-2)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(-1)).size() == 0);
		CHECK(lookup.list(surface.layer_idx(0)).size() == 1);
		CHECK(lookup.list(surface.layer_idx(1)).size() == 1);
		CHECK(lookup.list(surface.layer_idx(2)).size() == 5);

		CHECK(lookup(Vec3i(-4, 0, 0)) == 0);
		CHECK(lookup(Vec3i(-3, 0, 0)) == 0);
		CHECK(
			[](const UINT idx) {
				return 0 <= idx <= 4;
			}((lookup(Vec3i(-5, 0, 0)))
		));
		CHECK(lookup(Vec3i(-6, 0, 0)) == lookup.NULL_IDX);

	}
}

/**
 * Test delta isogrid update spread using Gaussian distribution given a list of
 * positions to update.
 */
WHEN("gaussian_from_list")
{
	// ==== Setup ====
	Surface<2, 2> surface(Vec2u(16, 16));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});

	SingleLookupGrid<2, surface.NUM_LAYERS> lookup = surface.walk_band<2>(
		Vec2i(-3,0)
	);

	// ==== Action ====
	surface.update_start();
	surface.delta_gauss(
		lookup.list(surface.layer_idx(0)), Vec2f(-3.5f, 0), 0.5f, 0.2f
	);
	surface.update_end();


	// === Confirm ===

	CHECK(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1)) == Approx(0.5f).epsilon(0.0000001f));

	CHECK(
		surface.delta().get(Vec2i(-3, 0)) == Approx(0.3457f).epsilon(0.0001f));
	CHECK(
		surface.delta().get(Vec2i(-2, -1)) == Approx(0.07714f).epsilon(0.0001f));
	CHECK(
		surface.delta().get(Vec2i(-2, 1)) == Approx(0.07714f).epsilon(0.0001f));
}

/**
 * Test delta isogrid update spread using Gaussian distribution given point and
 * radius.
 */
WHEN("gaussian_from_dist")
{
	// ==== Setup ====
	Surface<2, 2> surface(Vec2u(16, 16));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});


	// ==== Action ====
	surface.update_start();
	surface.delta_gauss<2>(
		Vec2f(-3.0f, 0), 0.5f, 0.2f
	);
	surface.update_end();


	// === Confirm ===

	CHECK(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1)) == Approx(0.5f).epsilon(0.0000001f));

	CHECK(
		surface.delta().get(Vec2i(-3, 0)) == Approx(0.28805843f).epsilon(0.0001f));
	CHECK(
		surface.delta().get(Vec2i(-2, -1)) == Approx(0.105970778f).epsilon(0.0001f));
	CHECK(
		surface.delta().get(Vec2i(-2, 1)) == Approx(0.105970778f).epsilon(0.0001f));
}

/**
 * Test raycasting to zero curve.
 */
WHEN("ray")
{
	// ==== Setup ====
	Surface<3, 3> surface(Vec3u(32, 32, 32), Vec3u(5, 5, 5));
	Vec3f pos_hit;

	// Create seed point and expand the narrow band.
	surface.seed(Vec3i(0, 0, 0));
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});

//	INFO(stringifyGridSlice(surface.isogrid()));
/*
/home/dave/Dropbox/Workspace/Felt/src/tests/test_Surface.cpp(1,185): Message:
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -3 |   -2 |   -1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    1 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    2 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    3 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
|    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |    4 |
*/

	// ==== Action ====
	// Simplest "dead on" case - from outside grid.
	pos_hit = surface.ray(Vec3f(-100.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	CHECK(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm() <= 0.00001f
	);

	// ==== Action ====
	// Simplest "dead on" case - from inside grid.
	pos_hit = surface.ray(Vec3f(-6.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	CHECK(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm() <= 0.00001f
	);

	// ==== Action ====
	// Simplest "dead on" case - from inside surface.
	pos_hit = surface.ray(Vec3f(0.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	CHECK(pos_hit == surface.NULL_POS<FLOAT>());

	// ==== Action ====
	// Simplest "dead on" case - from zero layer.
	pos_hit = surface.ray(Vec3f(-3.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	CHECK(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm() <= 0.00001f
	);

	// ==== Setup ====
	surface.update([](auto& pos, auto& isogrid) {
		return -0.3f;
	});

	// ==== Action ====
	// Ray  interpolate to zero curve.
	pos_hit = surface.ray(Vec3f(-10.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	CHECK(
		(pos_hit - Vec3f(-3.3f, 0, 0)).squaredNorm() <= 0.00001f
	);

	// ==== Setup ====
	surface.update([](auto& pos, auto& isogrid) {
		return 0.3f;
	});
	INFO(stringifyGridSlice(surface.isogrid()));

//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(0, -1, 0)).offset()));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(0, -1, 0)).size()));
//	INFO(stringifyGridSlice(surface.isogrid().children().get(Vec3i(0, -1, 0))));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(-1, 0, 0)).offset()));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(-1, 0, 0)).size()));
//	INFO(stringifyGridSlice(surface.isogrid().children().get(Vec3i(-1, 0, 0))));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(-1, -1, 0)).offset()));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(-1, -1, 0)).size()));
//	INFO(stringifyGridSlice(surface.isogrid().children().get(Vec3i(-1, -1, 0))));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(0, 0, 0)).offset()));
//	INFO(stringifyVector(surface.isogrid().children().get(Vec3i(0, 0, 0)).size()));
//	INFO(stringifyGridSlice(surface.isogrid().children().get(Vec3i(0, 0, 0))));

	// ==== Action ====
	// Ray at an angle.
	pos_hit = surface.ray(
		Vec3f(-10.0f, -10.0f, 0.0f), Vec3f(1, 1, 0).normalized()
	);

	// ==== Confirm ====
	CHECK(
		(pos_hit - Vec3f(-1.5, -1.5, 0)).squaredNorm() <=
		0.00001f
	);

	pos_hit = surface.ray(Vec3f(10, 10, 10), Vec3f(-1, -1, -1).normalized());

	// ==== Confirm ====
	CHECK(pos_hit != surface.NULL_POS<FLOAT>());

	// ==== Action ====
	// Rotating ray.
	using MatrixType = Eigen::Matrix3f;
	MatrixType mat_rot(MatrixType::Identity());
	pos_hit = surface.ray(Vec3f(6.72, -6.55, -3.45), Vec3f(-0.672, 0.655, 0.345));

	CHECK(pos_hit != surface.NULL_POS<FLOAT>());

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f::UnitY()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10.0f);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		INFO(
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
		CHECK(
			pos_hit != surface.NULL_POS<FLOAT>());
	}

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f(1, 1, 1).normalized()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		INFO(
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
		CHECK(
			pos_hit != surface.NULL_POS<FLOAT>());
	}

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f(0, 1, 1).normalized()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		INFO(
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
		CHECK(
			pos_hit != surface.NULL_POS<FLOAT>());
	}
	// ==== Confirm ====
//	INFO(stringifyGridSlice(surface.isogrid()));
/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -3 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |

*/
}


WHEN("gaussian_from_ray")
{
	// ==== Setup ====
	Surface<2, 2> surface(Vec2u(16, 16), Vec2u(3, 3));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& isogrid) {
		return -1.0f;
	});
	FLOAT leftover;


	// ==== Action ====
	surface.update_start();
	leftover = surface.delta_gauss<2>(
		Vec2f(-2.4f, -10.0f), Vec2f(0, 1), 0.5f, 0.2f
	);
	surface.update_end();


	// === Confirm ===
//	INFO(stringifyGridSlice(surface.isogrid()));
/*
/home/dave/Dropbox/Workspace/Felt/src/tests/test_Surface.cpp(1,198): Message:
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -3 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |

*/
	CHECK(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1))
		+ surface.delta().get(Vec2i(-1, -2)) == Approx(0.5f).epsilon(0.000001f));

	CHECK(leftover <= 0.000001f);

//	CHECK(
//		surface.delta().get(Vec2i(-3, 0)) == Approx(0.152139202f).epsilon(0.0001f
//));
//	CHECK(
//		surface.delta().get(Vec2i(-2, -1)) == Approx(0.24048205f).epsilon(0.0001f
//));
//	CHECK(
//		surface.delta().get(Vec2i(-2, 1)) == Approx(0.0559686795).epsilon(0.0001f
//));
//	CHECK(
//		surface.delta().get(Vec2i(-1, -2)) == Approx(0.0514338501).epsilon(0.0001f
//));
}

}
