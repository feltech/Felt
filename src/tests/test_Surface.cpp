#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include "Felt/Surface.hpp"
#include "Utils.hpp"

using namespace felt;

/// Alias std::to_string for compactness.
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
	surface.isogrid().get(pos) = 0;
	surface.layer_add(pos, 0);

	// Check zero-layer array has registered point.
	CHECK(surface.layer(0).size() == 1);
	CHECK((*surface.layer(0).begin() - pos) == Vec3i::Zero());

	// Check layer calculation from value.
	// -- zero-layer point just added.
	CHECK(surface.layer_id(pos) == 0);

	// Move a point from layer 0 to layer -1
	surface.layer_move(pos, 0, -1);
	CHECK(surface.layer(-1).size() == 1);
	CHECK(surface.layer(0).size() == 0);
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

	pos_next = surface.next_closest(pos_next, 1);

	CHECK(isogrid(pos_next) == 2);

	pos_next = surface.next_closest(pos_next, 1);

	CHECK(isogrid(pos_next) == 1);

	pos_next = surface.next_closest(pos_next, 1);

	CHECK(isogrid(pos_next) == 0);

	pos_next = surface.next_closest(pos_next, 1);

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

	pos_next = surface.next_closest(pos_next, -1);

	CHECK(isogrid(pos_next) == -1);
	CHECK(pos_next == Vec2i(1, 0));

	pos_next = surface.next_closest(pos_next, -1);

	CHECK(isogrid(pos_next) == 0);
	CHECK(pos_next == Vec2i(0, 0));
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
	surface.layer_add(Vec3i(0, 0, 0), 0.0f);

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
			FLOAT val = surface.isogrid().get(pos);
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
				FLOAT val = surface.isogrid().get(pos);
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


GIVEN("a 9x9 2-layer surface with a singularity seed at the centre")
{
	// 2D surface with 2x narrow band layers, respectively.
	using Surface_t = Surface<2, 2>;
	// Construct the surface.
	Surface_t surface(Vec2u(9, 9));

	// Grid to use for checking data.
	Grid<FLOAT, 2> isogrid_check(Vec2u(5, 5), Vec2i::Zero(), 0);

	// Create seed point in the centre
	surface.seed(Vec2i(0, 0));

	INFO(stringify_grid_slice(surface.isogrid()));

	// Trivially check centre of seed is indeed a zero-level point (i.e. point
	// on the surface).
	THEN("the value at the centre of the grid is 0")
	{
		const FLOAT val_centre = surface.isogrid().get(Vec2i(0, 0));
		CHECK(val_centre == 0);
	}

	THEN("the surface data matches a singularity seed point")
	{
		// A 2D 2-layer singularity (seed) point should look like the following.

		isogrid_check.data() = {
			  3,    3,    3,    3,    3,    3,    3,    3,    3,
			  3,    3,    3,    3,    3,    3,    3,    3,    3,
			  3,    3,    3,    3,    2,    3,    3,    3,    3,
			  3,    3,    3,    2,    1,    2,    3,    3,    3,
			  3,    3,    2,    1,    0,    1,    2,    3,    3,
			  3,    3,    3,    2,    1,    2,    3,    3,    3,
			  3,    3,    3,    3,    2,    3,    3,    3,    3,
			  3,    3,    3,    3,    3,    3,    3,    3,    3,
			  3,    3,    3,    3,    3,    3,    3,    3,    3
		};

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();

		const FLOAT diff = isogrid_check.vdata().sum();

		CHECK(diff == 0);

		// Check appropriate points have been added to narrow band layers.
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(-1).size() == 0);
		CHECK(surface.layer(0).size() == 1);
		CHECK(surface.layer(1).size() == 4);
		CHECK(surface.layer(2).size() == 8);
	}

	AND_WHEN("we expand the surface one unit outwards")
	{
		surface.update([](auto& pos, auto& isogrid) {
			return -1.0f;
		});

		INFO(stringify_grid_slice(surface.isogrid()));

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

			isogrid_check.vdata() = isogrid_check.vdata() -
				surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();

			CHECK(diff == Approx(0));
		}

		AND_WHEN("we expand by one unit again")
		{
			surface.update([](auto& pos, auto& isogrid) {
				return -1.0f;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

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

				CHECK(diff == Approx(0));
			}

			AND_WHEN("we expand by one unit 9 more times")
			{
				for (UINT i = 0; i < 9; i++)
				{
					surface.update([](auto& pos, auto& isogrid) {
						return -1.0f;
					});
				}

				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("the surface data matches an area completely consumed by the surface")
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

					CHECK(diff == Approx(0));

					CHECK(surface.layer(0).size() == 0);
					CHECK(surface.layer(-1).size() == 0);
					CHECK(surface.layer(-2).size() == 0);
					CHECK(surface.layer(1).size() == 0);
					CHECK(surface.layer(2).size() == 0);
				}
			}
		}

		AND_WHEN("we contract the surface by 1 unit inwards")
		{
			surface.update_start();
			{
				for (const Vec2i pos : surface.layer(0))
					surface.delta(pos, 1.0f);
			}
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the surface has collapsed back to a singularity")
			{
				// A 2D 2-layer singularity (seed) point should look like the following.

				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.vdata() = isogrid_check.vdata() -
					surface.isogrid().snapshot().vdata();
				INFO(stringify_grid_slice(surface.isogrid().snapshot()));
				INFO(stringify_grid_slice(isogrid_check));

				const FLOAT diff = isogrid_check.vdata().sum();

				CHECK(diff == 0);

				CHECK(surface.layer(-2).size() == 0);
				CHECK(surface.layer(-1).size() == 0);
				CHECK(surface.layer(0).size() == 1);
				CHECK(surface.layer(1).size() == 4);
				CHECK(surface.layer(2).size() == 8);
			}

			THEN("iterating over layer 0 gives 1 point")
			{
				UINT total_iterations = 0;
				for (auto& pos : surface.layer(0))
					total_iterations++;

				CHECK(total_iterations == 1);
			}

			AND_WHEN("we contract the surface by 1 unit inwards again")
			{
				surface.update_start();
				{
					for (const Vec2i pos : surface.layer(0))
						surface.delta(pos, 1.0f);
				}
				surface.update_end();

				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("iterating over layer 0 gives 0 points")
				{
					UINT total_iterations = 0;
					for (auto& pos : surface.layer(0))
						total_iterations++;

					CHECK(total_iterations == 0);
				}

				THEN("the surface data matches an area completely outside the surface")
				{
					isogrid_check.data() = {
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3
					};

					isogrid_check.vdata() = isogrid_check.vdata() -
						surface.isogrid().snapshot().vdata();

					const FLOAT diff = isogrid_check.vdata().sum();

					CHECK(diff == 0);

					CHECK(surface.layer(-2).size() == 0);
					CHECK(surface.layer(-1).size() == 0);
					CHECK(surface.layer(0).size() == 0);
					CHECK(surface.layer(1).size() == 0);
					CHECK(surface.layer(2).size() == 0);
				}
			}
		} // End when contract by 1.
	} // End when expand by 1.

	WHEN("we expand by 0.6")
	{
		// Expand the narrow band.
		surface.update([](auto& pos, auto& isogrid) {
			return -0.6f;
		});

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the grid data and layers are as expected")
		{
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
				  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
				  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
				  3,  2.4,  1.4,  0.4, -0.6,  0.4,  1.4,  2.4,    3,
				  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
				  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
				  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			// snapshot() copies the spatially partitioned grid into a single
			// unpartitioned grid. vdata() returns an Eigen vector with the grid data, giving
			// access to arithmetic and other BLAS functions.
			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();
			CHECK(diff == Approx(0));

			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 1);
			CHECK(surface.layer(0).size() == 4);
			CHECK(surface.layer(1).size() == 8);
			CHECK(surface.layer(2).size() == 12);
		}

		AND_WHEN("we contract by 0.6")
		{
			// Update using lambda.
			surface.update([](const Vec2i& pos, const Surface_t::IsoGrid& isogrid) {
				return 0.6f;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the surface is once more a seed")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				CHECK(surface.layer(-2).size() == 0);
				CHECK(surface.layer(-1).size() == 0);
				CHECK(surface.layer(0).size() == 1);
				CHECK(surface.layer(1).size() == 4);
				CHECK(surface.layer(2).size() == 8);

				isogrid_check.vdata() =
					isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
				const FLOAT diff = isogrid_check.vdata().sum();
				CHECK(diff == Approx(0));
			}

			AND_WHEN("we contract by 0.6 again")
			{
				// Update using lambda.
				surface.update([](const Vec2i& pos, const Surface_t::IsoGrid& isogrid) {
					return 0.6f;
				});

				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("the surface has completely collapsed and all points are outside")
				{
					isogrid_check.data() = {
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3
					};

					CHECK(surface.layer(-2).size() == 0);
					CHECK(surface.layer(-1).size() == 0);
					CHECK(surface.layer(0).size() == 0);
					CHECK(surface.layer(1).size() == 0);
					CHECK(surface.layer(2).size() == 0);

					isogrid_check.vdata() =
						isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
					const FLOAT diff = isogrid_check.vdata().sum();
					CHECK(diff == Approx(0));
				}
			}
		}
	}
}


GIVEN("a 16x9 2-layer surface with two small regions side-by-side")
{
	Surface<2, 2> surface(Vec2u(16, 9), Vec2u::Constant(3));
	Grid<FLOAT, 2> isogrid_check(Vec2u(16, 9), Vec2i::Zero(), 0);
	// Create two seed points and expand the narrow band.
	surface.seed(Vec2i(-4, 0));
	surface.seed(Vec2i(4, 0));
	surface.update([](auto& pos, auto& grid) {
		return -1.0f;
	});

	INFO(stringify_grid_slice(surface.isogrid()));


	THEN("outermost layers in central partitions are as expected")
	{
		CHECK(surface.layer(Vec2i(0,0), 2).size() == 3);
		CHECK(surface.layer(Vec2i(1,0), 2).size() == 3);
	}

	THEN("the surface is in the expected state")
	{
		isogrid_check.data() = {
			  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
			  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
			  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
			  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
			  3,    2,    1,    0,   -1,    0,    1,    2,    3,    2,    1,    0,   -1,    0,    1,    2,
			  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
			  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
			  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
			  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
		};

		isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		CHECK(diff == Approx(0));

		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(-1).size() == 2);
		CHECK(surface.layer(0).size() == 8);
		CHECK(surface.layer(1).size() == 16);
		CHECK(surface.layer(2).size() == 24);
	}

	WHEN("we expand the subsurfaces towards one-another")
	{
		surface.update_start();
		surface.delta(Vec2i(-3, 0), -1.0f);
		surface.delta(Vec2i(3, 0), -1.0f);
		surface.update_end();

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the centremost partitions contain the expected number of outer layer points")
		{
			CHECK(surface.layer(Vec2i(0,0), 2).size() == 3);
			CHECK(surface.layer(Vec2i(1,0), 2).size() == 2);
		}

		THEN("the surface is in the expected state")
		{
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    2,    2,    3,    3,    3,    3,    3,    2,    2,    3,    3,    3,
				  3,    3,    3,    2,    1,    1,    2,    3,    3,    3,    2,    1,    1,    2,    3,    3,
				  3,    3,    2,    1,    0,    0,    1,    2,    3,    2,    1,    0,    0,    1,    2,    3,
				  3,    2,    1,    0,   -1,   -1,    0,    1,    2,    1,    0,   -1,   -1,    0,    1,    2,
				  3,    3,    2,    1,    0,    0,    1,    2,    3,    2,    1,    0,    0,    1,    2,    3,
				  3,    3,    3,    2,    1,    1,    2,    3,    3,    3,    2,    1,    1,    2,    3,    3,
				  3,    3,    3,    3,    2,    2,    3,    3,    3,    3,    3,    2,    2,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();
			CHECK(diff == Approx(0));

			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 4);
			CHECK(surface.layer(0).size() == 12);
			CHECK(surface.layer(1).size() == 20);
			CHECK(surface.layer(2).size() == 27);
		}
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

	INFO(stringify_grid_slice(surface.isogrid()));
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

}


SCENARIO("Local updating")
{
GIVEN("a 9x9 2-layer surface with a seed point in the centre")
{
	typedef Surface<2, 2> SurfaceT;
	SurfaceT surface(Vec2u(9, 9), Vec2u(2, 2));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(Vec2u(9, 9), Vec2i::Zero(), 0);

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));

	WHEN("we contract the surface by 1 unit inwards")
	{
		surface.update_start();
		surface.delta(Vec2i(0, 0), 1.0f);
		surface.update_end_local();
		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the surface data matches an area completely outside the surface")
		{
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.vdata() = isogrid_check.vdata() -
				surface.isogrid().snapshot().vdata();

			const FLOAT diff = isogrid_check.vdata().sum();

			CHECK(diff == 0);

			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 0);
			CHECK(surface.layer(0).size() == 0);
			CHECK(surface.layer(1).size() == 0);
			CHECK(surface.layer(2).size() == 0);
		}
	}

	WHEN("we expand by 1 unit")
	{
		//![Calculate affected outer layers for localised narrow band updates]
		surface.update_start();
		{
			for (auto pos : surface.layer(0))
			{
				surface.delta(pos, -1.0f);
			}
		}
		surface.update_end();

		AND_WHEN("we modify a couple of points and calculate the affected narrow band points")
		{
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

			THEN("the affected narrow band points are as expected")
			{
				using PosArray = SurfaceT::PosArray;
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
						surface.affected().leafs(layer_idx).size()
						==
						check_layers_pos[layer_idx].size()
					);

					for (auto pos : check_layers_pos[layer_idx])
					{
						auto iter = std::find(
							surface.affected().leafs(layer_idx).begin(),
							surface.affected().leafs(layer_idx).end(), pos
						);

						INFO(
							"Affected grid layer " + s(layer_id) + " at index " +
							s(layer_idx) + " should contain (" + s(pos(0)) + "," +
							s(pos(1)) + ")"
						);
						CHECK(iter != surface.affected().leafs(layer_idx).end());
					}

					for (auto pos : surface.affected().leafs(layer_idx))
					{
						auto iter = std::find(
							check_layers_pos[layer_idx].begin(),
							check_layers_pos[layer_idx].end(), pos
						);

						INFO(
							"Checking list layer " + s(layer_id) + " at index " +
							s(layer_idx) + " should contain (" + s(pos(0)) + "," +
							s(pos(1)) + ")"
						);
						CHECK(iter != check_layers_pos[layer_idx].end());
					}
				}
			}
		}
		//![Calculate affected outer layers for localised narrow band updates]

		AND_WHEN("we cycle a square region partially containing the surface")
		{
			UINT num_visited = 0;
			Vec2i pos_visited;
			surface.update(
				Vec2i(1, 0), Vec2i(3, 3),
				[&num_visited, &pos_visited](const Vec2i& pos, const SurfaceT::IsoGrid& grid) {
					num_visited++;
					pos_visited = pos;
					return 0;
				}
			);

			THEN("we only visit the points in the region")
			{
				CHECK(num_visited == 1);
				CHECK(pos_visited == Vec2i(1, 0));
			}
		}

		AND_WHEN("we cycle a square region completely containing the surface")
		{
			UINT num_visited = 0;
			Vec2i pos_visited;
			surface.update(
				Vec2i(-100, -100), Vec2i(100, 100),
				[&num_visited, &pos_visited](const Vec2i& pos, const SurfaceT::IsoGrid& grid) {
					num_visited++;
					return 0;
				}
			);

			THEN("we only visit the valid points")
			{
				CHECK(num_visited == 4);
			}
		}

	} // End WHEN("we expand by 1 unit")

	WHEN("we expand the centre point")
	{
		surface.update_start();
		surface.delta(Vec2i(0, 0), -0.6f);
		// Using localised update, which will only update outer layers that are
		// affected by changes to the modified zero layer points.  In this test
		// case, all outer layer points are affected, same as a global update.
		surface.update_end_local();

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the grid data is as expected")
		{
			isogrid_check.data() = {
				3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3, 3, 3, 3, 3, 3,
				2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 0.4, 1.4, 2.4, 3, 3, 3,
				2.4, 1.4, 0.4, -0.6, 0.4, 1.4, 2.4, 3, 3, 3, 2.4, 1.4, 0.4, 1.4,
				2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3,
				3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();
			CHECK(diff == Approx(0));

			CHECK(surface.layer(0).size() == 4);
			CHECK(surface.layer(-1).size() == 1);
			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(1).size() == 8);
			CHECK(surface.layer(2).size() == 12);
		}

		AND_WHEN("we contract the centre point by the same amount using a local update")
		{

			// Cycle new zero-layer points and move back to original signed distance.
			surface.update_start();
			{
				for (const Vec2i& pos : surface.layer(0))
					surface.delta(pos, 0.6f);
			}
			surface.update_end_local();

			THEN("the grid data is back to how it was")
			{
				isogrid_check.data() = {
					3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
					2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 2, 1, 0, 1, 2, 3,
					3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
					3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

				isogrid_check.vdata() =
					isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
				const FLOAT diff = isogrid_check.vdata().sum();

				CHECK(diff == Approx(0).epsilon(0.000001f));
			}
		}
	}
}

GIVEN("a 16x9 2-layer surface with two small regions side-by-side")
{
	Surface<2, 2> surface(Vec2u(16, 9), Vec2u::Constant(3));
	Grid<FLOAT, 2> isogrid_check(Vec2u(16, 9), Vec2i::Zero(), 0);
	// Create two seed points and expand the narrow band.
	surface.seed(Vec2i(-4, 0));
	surface.seed(Vec2i(4, 0));
	surface.update([](auto& pos, auto& grid) {
		return -1.0f;
	});

	INFO(stringify_grid_slice(surface.isogrid()));

	WHEN("we expand/contract the subsurfaces towards one-another")
	{
		for (UINT i = 0; i < 10; i++)
		{
			surface.update_start();
			surface.delta(Vec2i(-3, 0), -1.0f);
			surface.delta(Vec2i(3, 0), -1.0f);
			surface.update_end_local();

			surface.update_start();
			surface.delta(Vec2i(-3, 1), 1.0f);
			surface.delta(Vec2i(-2, 0), 1.0f);
			surface.delta(Vec2i(-3, -1), 1.0f);

			surface.delta(Vec2i(3, 1), 1.0f);
			surface.delta(Vec2i(2, 0), 1.0f);
			surface.delta(Vec2i(3, -1), 1.0f);
			surface.update_end_local();
		}

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("outermost layers in central partitions are as expected")
		{
			CHECK(surface.layer(Vec2i(0,0), 2).size() == 3);
			CHECK(surface.layer(Vec2i(1,0), 2).size() == 3);
		}

		THEN("the surface has returned to its inital state")
		{
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
				  3,    2,    1,    0,   -1,    0,    1,    2,    3,    2,    1,    0,   -1,    0,    1,    2,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();
			CHECK(diff == Approx(0));

			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 2);
			CHECK(surface.layer(0).size() == 8);
			CHECK(surface.layer(1).size() == 16);
			CHECK(surface.layer(2).size() == 24);
		}
	}
}


GIVEN("a 11x11x11 3-layer surface with 3x3x3 partitions initialised with a 1 unit radius surface")
{
	Surface<3,3> surface(Vec3u(11,11,11), Vec3u(3,3,3));
	surface.seed(Vec3i(0,0,0));
	surface.update([](auto& pos, auto& grid) {
		return -1.0f;
	});

	INFO(stringify_grid_slice(surface.isogrid()));

	THEN("the layer lists have the expected size")
	{
		CHECK(surface.layer(-3).size() == 0);
		CHECK(surface.layer(-2).size() == 0);
		CHECK(surface.layer(-1).size() == 1);
		CHECK(surface.layer(0).size() == 6);
		CHECK(surface.layer(1).size() == 18);
		CHECK(surface.layer(2).size() == 38);
		CHECK(surface.layer(3).size() == 66);
	}

	AND_WHEN("we expand and contract two nearby points using local update")
	{
		surface.update_start();
		surface.delta(Vec3i(0,1,0), -1);
		surface.update_end_local();

		surface.update_start();
		surface.delta(Vec3i(0,2,0), 1);
		surface.delta(Vec3i(1,1,0), 1);
		surface.delta(Vec3i(-1,1,0), 1);
		surface.delta(Vec3i(0,1,1), 1);
		surface.delta(Vec3i(0,1,-1), 1);
		surface.update_end_local();

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the layer lists have the same size as before")
		{
			CHECK(surface.layer(-3).size() == 0);
			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 1);
			CHECK(surface.layer(0).size() == 6);
			CHECK(surface.layer(1).size() == 18);
			CHECK(surface.layer(2).size() == 38);
			CHECK(surface.layer(3).size() == 66);
		}
	}
}
}


SCENARIO("Complex layer interactions")
{
	GIVEN(
		"a 12x12 3-layer grid with two seeds diagonally opposite, and the left seed expanded"
		" towards the right"
	) {
		// ==== Setup ====
		using SurfaceT = Surface<2, 3>;
		Vec2u size(12, 12);
		SurfaceT surface(size, Vec2u(2, 2));
		// Grid to set values of manually, for checking against.
		Grid<FLOAT, 2> isogrid_check(size, Vec2i::Zero(), 0);

		// Create seed point and expand the narrow band.
		surface.seed(Vec2i(-2, -2));
		surface.seed(Vec2i(2, 2));

		INFO(stringify_grid_slice(surface.isogrid()));

//    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
//    4,    4,    4,    4,    3,    4,    4,    4,    4,    4,    4,    4,
//    4,    4,    4,    3,    2,    3,    4,    4,    4,    4,    4,    4,
//    4,    4,    3,    2,    1,    2,    3,    4,    4,    4,    4,    4,
//    4,    3,    2,    1,    0,    1,    2,    3,    4,    4,    4,    4,
//    4,    4,    3,    2,    1,    2,    3,    4,    3,    4,    4,    4,
//    4,    4,    4,    3,    2,    3,    4,    3,    2,    3,    4,    4,
//    4,    4,    4,    4,    3,    4,    3,    2,    1,    2,    3,    4,
//    4,    4,    4,    4,    4,    3,    2,    1,    0,    1,    2,    3,
//    4,    4,    4,    4,    4,    4,    3,    2,    1,    2,    3,    4,
//    4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//    4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4

		surface.update_start();
		surface.delta(Vec2i(-2, -2), -1);
		surface.update_end();
		surface.update_start();
		surface.delta(Vec2i(-1, -2), -1);
		surface.update_end();
		surface.update_start();
		surface.delta(Vec2i(0, -2), -1);
		surface.update_end();

		INFO(stringify_grid_slice(surface.isogrid()));

//	      4,    4,    4,    4,    3,    3,    3,    4,    4,    4,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    3,    4,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    2,    3,    4,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    1,    2,    3,    4,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,    0,    1,    2,    3,    4,
//	      4,    3,    2,    1,    0,    0,    0,    1,    2,    3,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    2,    2,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    3,    3,    2,    1,    0,    1,    2,    3,
//	      4,    4,    4,    4,    4,    4,    3,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4


//		Bad:
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    3,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    3,    3,    3,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4

//		Also bad:
//	      4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    4,    4,    3,    3,    3,    3,    3,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4

		WHEN(
			"we simultaneously expand the left seed and contract the right, then expand the "
			"left again, using local updates"
		) {

			surface.update_start();
			surface.delta(Vec2i(1, -2), -1);
			surface.delta(Vec2i(2, 2), 1);
			surface.update_end_local();

			INFO(stringify_grid_slice(surface.isogrid()));

			surface.update_start();
			surface.delta(Vec2i(1, -1), -1);
			surface.update_end_local();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data is as expected")
			{

				isogrid_check.data() = {
					  4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
					  4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
					  4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
					  4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
					  3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
					  4,    3,    2,    1,    0,    0,    0,   -1,    0,    1,    2,    3,
					  4,    4,    3,    2,    1,    1,    1,    0,    1,    2,    3,    4,
					  4,    4,    4,    3,    2,    2,    2,    1,    2,    3,    4,    4,
					  4,    4,    4,    4,    3,    3,    3,    2,    3,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4
				};

				isogrid_check.vdata() =
					isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
				const FLOAT diff = isogrid_check.vdata().sum();
				CHECK(diff == Approx(0).epsilon(0.000001f));

				CHECK(surface.layer(-3).size() == 0);
				CHECK(surface.layer(-2).size() == 0);
				CHECK(surface.layer(-1).size() == 5);
				CHECK(surface.layer(0).size() == 11);
				CHECK(surface.layer(1).size() == 15);
				CHECK(surface.layer(2).size() == 19);
				CHECK(surface.layer(3).size() == 23);
			}
		}

		WHEN(
			"we simultaneously expand the left seed and contract the right, then expand the "
			"left again, using global updates"
		) {

			surface.update_start();
			surface.delta(Vec2i(1, -2), -1);
			surface.delta(Vec2i(2, 2), 1);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			surface.update_start();
			surface.delta(Vec2i(1, -1), -1);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data is as expected")
			{
				isogrid_check.data() = {
					  4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
					  4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
					  4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
					  4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
					  3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
					  4,    3,    2,    1,    0,    0,    0,   -1,    0,    1,    2,    3,
					  4,    4,    3,    2,    1,    1,    1,    0,    1,    2,    3,    4,
					  4,    4,    4,    3,    2,    2,    2,    1,    2,    3,    4,    4,
					  4,    4,    4,    4,    3,    3,    3,    2,    3,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4
				};

				isogrid_check.vdata() =
					isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
				const FLOAT diff = isogrid_check.vdata().sum();
				CHECK(diff == Approx(0).epsilon(0.000001f));

				CHECK(surface.layer(-3).size() == 0);
				CHECK(surface.layer(-2).size() == 0);
				CHECK(surface.layer(-1).size() == 5);
				CHECK(surface.layer(0).size() == 11);
				CHECK(surface.layer(1).size() == 15);
				CHECK(surface.layer(2).size() == 19);
				CHECK(surface.layer(3).size() == 23);
			}
		}
	}
}


SCENARIO("Raycasting")
{
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

//	INFO(stringify_grid_slice(surface.isogrid()));
/*
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
	pos_hit = surface.ray(Vec3f(-35.0f, 0, 0), Vec3f(1, 0, 0));

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
	pos_hit = surface.ray(Vec3f(0, 0, 0), Vec3f(1, 0, 0));

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
	INFO(stringify_grid_slice(surface.isogrid()));

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
			"Ray hit from " + felt::format(origin) + " in direction "
			+ felt::format(dir) + " should not be NULL_POS"
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
			"Ray hit from " + felt::format(origin) + " in direction "
			+ felt::format(dir) + " should not be NULL_POS"
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
			"Ray hit from " + felt::format(origin) + " in direction "
			+ felt::format(dir) + " should not be NULL_POS"
		);
		CHECK(
			pos_hit != surface.NULL_POS<FLOAT>());
	}
	// ==== Confirm ====
//	INFO(stringify_grid_slice(surface.isogrid()));
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

GIVEN("a 3-layer flat surface in an 20x20x20 grid with 16x16x16 partitions")
{
	Surface<3, 3> surface(Vec3u(20, 20, 20), Vec3u(16, 16, 16));
	surface.seed(Vec3i(0,0,0));
	surface.update([](auto& pos, auto& phi)->FLOAT {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& phi)->FLOAT {
		return -1.0f;
	});
	for (UINT i = 0; i < 10; i++)
		surface.update([](auto& pos, auto& grid) {
			if (std::abs(pos(1)) > 1)
				return 0;
			else
				return -1;
		});
//	INFO(stringify_grid_slice(surface.isogrid()));
/*
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
	  4,    4,    4,    4,    4,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    4,    4,    4,    4,
*/
	WHEN("we cast a ray diagonally downward from outside the isogrid")
	{
		const Vec3f& pos_hit = surface.ray(
			Vec3f(-5.45783, 44.8901, -57.4607), Vec3f(0.134944, -0.616392, 0.77579).normalized()
		);

		//pos + 69.5*dir = (3.9205,2.051,-3.5433)

		THEN("the surface is hit")
		{
			CHECK(pos_hit != surface.NULL_POS<FLOAT>());
		}
	}
}

//Casting: (-1.29043  49.6148 -66.8919) => 0.0725882 -0.660291  0.747493
GIVEN("a 3-layer flat periodic surface in an 50x50x50 grid with 16x16x16 partitions")
{
	Surface<3, 3> surface(Vec3u(50, 50, 50), Vec3u(16, 16, 16));
	surface.seed(Vec3i(0,0,0));
	surface.update([](auto& pos, auto& phi)->FLOAT {
		return -1.0f;
	});
	surface.update([](auto& pos, auto& phi)->FLOAT {
		return -1.0f;
	});
	for (UINT i = 0; i < 20; i++)
		surface.update([](auto& pos, auto& grid) {
			if (std::abs(pos(1)) > 1)
				return 0;
			else
				return -1;
		});
//	INFO(stringify_grid_slice(surface.isogrid()));


	WHEN("we cast rays diagonally downward from outside the isogrid")
	{
		// | -25 -- -9 -- 7 -- 23 -- 50
		//pos + 69.5*dir = (3.9205,2.051,-3.5433)
		const Vec3f& pos_hit1 = surface.ray(
			Vec3f(-1.29043, 49.6148, -66.8919),
			Vec3f(0.0725882, -0.660291, 0.747493).normalized()
		);
		//pos + 32.5*dir = (-3.73342,1.94405,-18.64452)
		const Vec3f& pos_hit2 = surface.ray(
			Vec3f(-0.0219189, 18.1713, -46.5578),
			Vec3f(-0.114205,-0.499295, 0.858872).normalized()
		);
		//pos + 34.7*dir = (-1.33501,2.01918,-15.87545)
		const Vec3f& pos_hit3 = surface.ray(
			Vec3f(-0.0139845, 18.1755, -46.5565),
			Vec3f(-0.0380706, -0.465599, 0.884177).normalized()
		);

		THEN("the surface is hit")
		{
			CHECK(pos_hit1 != surface.NULL_POS<FLOAT>());
			CHECK(pos_hit2 != surface.NULL_POS<FLOAT>());
			CHECK(pos_hit3 != surface.NULL_POS<FLOAT>());
		}
	}
}


GIVEN("a 2-layer surface of radius 3 in a 16x16 grid with 3x3 partitions")
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

	INFO(stringify_grid_slice(surface.isogrid()));


	WHEN("we cast a ray upward from just below the grid")
	{
		const Vec2f& pos_hit = surface.ray(Vec2f(-2.4f, -10.0f), Vec2f(0, 1));

		THEN("the surface is hit where expected")
		{
			CHECK(pos_hit == ApproxVec(Vec2f(-2.21609, -0.78391)).epsilon(0.1));
		}
	}
}
}
