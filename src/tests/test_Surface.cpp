#include <boost/numeric/ublas/io.hpp>
#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "Felt/Surface.hpp"
#include "Utils.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_Surface)
/*
 * Basic initialisation.
 */
BOOST_AUTO_TEST_CASE(init)
{
	// ==== Setup ====
	// Basic initialisation of 2D surface with 2 layers in a 5x5 embedding.
	Surface<2, 2> surface(Vec2u(7, 7));

	// ==== Confirm ====
	BOOST_CHECK_EQUAL(surface.isogrid().size(), Vec2u(7, 7));
	BOOST_CHECK_EQUAL(surface.isogrid().children().data().size(), 1);
	BOOST_CHECK_EQUAL(surface.isogrid().children().get(Vec2i(0, 0)).size(), Vec2u(7, 7));
	BOOST_CHECK_EQUAL(surface.isogrid().children().get(Vec2i(0, 0)).data().size(), 7*7);
	// Usable isogrid should have size=size-layers and be offset by half grid width.
	// In this test case, only the centre point is actually usable, all other points are
	// reserved for outer layers.
	BOOST_CHECK_EQUAL(surface.pos_min(), Vec2i(0, 0));
	BOOST_CHECK_EQUAL(surface.pos_min(), Vec2i(0, 0));
	// But actual isogrid isogrid should have size equal to size.
	BOOST_CHECK_EQUAL(surface.isogrid().size(), Vec2u(7, 7));
	// Grid is initialised to all points 'outside' the surface (since there is no surface yet).
	BOOST_CHECK_EQUAL(surface.isogrid().get(Vec2i(0, 0)), 3);
}

/*
 * Narrow band layers.
 */
BOOST_AUTO_TEST_CASE(layers)
{
	// 3D surface with default (=2) number of layers.
	Surface<3> surface(Vec3u(7, 7, 7));
	Surface<3>::IsoGrid::ChildrenGrid children;
	Vec3i pos = Vec3i(0, 0, 0);

	BOOST_CHECK_EQUAL(children.list(surface.layer_idx(-2)).size(), 0);
	BOOST_CHECK_EQUAL(children.list(surface.layer_idx(-1)).size(), 0);
	BOOST_CHECK_EQUAL(children.list(surface.layer_idx(0)).size(), 0);
	BOOST_CHECK_EQUAL(children.list(surface.layer_idx(1)).size(), 0);
	BOOST_CHECK_EQUAL(children.list(surface.layer_idx(2)).size(), 0);

	// Add a single zero-layer point.
	surface.isogrid(pos) = 0;
	surface.layer_add(pos, 0);

	// Check zero-layer array has registered point.
	BOOST_CHECK_EQUAL(surface.layer(0).size(), 1);
	BOOST_CHECK_EQUAL((*surface.layer(0).begin() - pos), Vec3i::Zero());

	// Check layer calculation from value.
	// -- zero-layer point just added.
	BOOST_CHECK_EQUAL(surface.layer_id(pos), 0);

	// Add three arbitrary points to layer -1.
	surface.layer_add(Vec3i(0, 0, 1), -1);
	surface.layer_add(Vec3i(0, 0, 2), -1);
	surface.layer_add(Vec3i(0, 0, 3), -1);

	// Remove two points from layer -1
	surface.layer_remove(Vec3i(0, 0, 1), -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 2);
	surface.layer_remove(Vec3i(0, 0, 3), -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);

	// Move a point from layer 0 to layer -1
	surface.layer_move(pos, 0, -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 2);

	// Check lists updated.
	BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 2);
}

/*
 * Placing a single singularity point.
 */
BOOST_AUTO_TEST_CASE(seed)
{
	Surface<2, 2> surface(Vec2u(5, 5));
	Surface<2, 2>::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec2i(0, 0));

	surface.seed(Vec2i(0, 0));

	// Trivially check centre of seed is indeed a zero-level point (i.e. point
	// on the surface).

	const FLOAT val_centre = isogrid(Vec2i(0, 0));
	BOOST_CHECK_EQUAL(val_centre, 0);

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

	isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();

	const FLOAT diff = isogrid_check.vdata().sum();

	BOOST_CHECK_EQUAL(diff, 0);

	// Check appropriate points have been added to narrow band layers.
	BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(0).size(), 1);
	BOOST_CHECK_EQUAL(surface.layer(1).size(), 4);
	BOOST_CHECK_EQUAL(surface.layer(2).size(), 8);
}

/*
 * Given a grid point, find neighbouring point closest to zero-curve.
 */
BOOST_AUTO_TEST_CASE(next_closest_grid_point)
{
	// Create seed point, as above, and navigate to centre.

	Surface<2, 2> surface(Vec2u(5, 5));
	Surface<2, 2>::IsoGrid& isogrid = surface.isogrid();

	surface.seed(Vec2i(0, 0));

	Vec2i pos_next = Vec2i(-1, -2);
	BOOST_CHECK_EQUAL(isogrid(pos_next), 3);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), 2);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), 1);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), 0);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), 0);

	// Ensure it also works with negative distances.
	// NOTE: row-major (y,x) element ordering...
	surface.isogrid().data() = {
		2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -2,
		-2, -2, -2, -2 };
	surface.isogrid().flush_snapshot();
	// NOTE: ...but accessed as (x,y)
	pos_next = Vec2i(2, 0);

	BOOST_CHECK_EQUAL(isogrid(pos_next), -2);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(2, 0));

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), -1);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(1, 0));

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(isogrid(pos_next), 0);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(0, 0));
}

/*
 * Using delta isogrid grid/list.
 */
BOOST_AUTO_TEST_CASE(delta_isogrid_clamping)
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
	BOOST_CHECK_EQUAL(surface.delta(pos), -1);
	// Check position vector of point in surface grid that delta was
	// applied to is stored in a corresponding list to be iterated over.
	BOOST_CHECK_EQUAL(surface.delta().leafs(surface.layer_idx(0)).size(), 1);
	//! [Delta isogrid clamping]
}

/*
 * Update isogrid with delta isogrid.
 */
BOOST_AUTO_TEST_CASE(delta_isogrid_update)
{
	//! [Simple delta isogrid update]
	// ==== Setup ====
	Surface<3, 2> surface(Vec3u(5, 5, 5));
	Surface<3, 2>::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec3i(0,0,0));
	Surface<3, 2>::DeltaIsoGrid& delta = surface.delta();

	// ==== Action ====
	// Put in 'dirty' state, to check update_start is doing it's job.
	surface.delta(Vec3i(0, 0, 0), 0.5f);

	// ==== Confirm ====
	BOOST_CHECK_EQUAL(surface.delta().children().list(surface.layer_idx(0)).size(), 1);
	BOOST_CHECK_EQUAL(surface.delta().get(Vec3i(0, 0, 0)), 0.5f);

	// ==== Action ====
	// Clear delta isogrid.
	surface.update_start();

	// ==== Confirm ====
	// Check update_start cleared the above surface.delta changes.
	BOOST_CHECK_EQUAL(surface.delta().children().list(surface.layer_idx(0)).size(), 0);
	BOOST_CHECK_EQUAL(surface.delta().get(Vec3i(0, 0, 0)), 0.0f);

	// ==== Action ====
	// Add a zero-layer point.
	surface.isogrid(Vec3i(0, 0, 0), 0.0f);

	// Clear delta isogrid.
	surface.update_start();
	// Do nothing.
	surface.delta(Vec3i(0, 0, 0), 0);
	// Apply delta isogrid.
	surface.update_end();

	// ==== Confirm ====
	// Ensure nothing was changed.  Every point in 5x5x5 grid == 3, except centre which == 0.
	BOOST_CHECK_EQUAL(isogrid.vdata().sum(), 3 * 5 * 5 * 5 - 3);
	// Delta isogrid position vector list should still contain one point.
	BOOST_CHECK_EQUAL(surface.delta().leafs(surface.layer_idx(0)).size(), 1);
	// Delta isogrid grid itself should have reset back to zero.
	BOOST_CHECK_EQUAL(delta(Vec3i(0, 0, 0)), 0);

	// ==== Action ====
	// Clear delta isogrid.
	surface.update_start();
	// Apply small update.
	surface.delta(Vec3i(0, 0, 0), 0.4f);
	// Apply delta isogrid.
	surface.update_end();

	// ==== Confirm ====
	// Ensure change applied.  Every point in grid == 3, except centre which == 0.4.
	BOOST_CHECK_EQUAL(isogrid.vdata().sum(), 3 * 5 * 5 * 5 - 3 + 0.4f);
	BOOST_CHECK_EQUAL(isogrid.get(Vec3i(0, 0, 0)), 0.4f);
	//! [Simple delta isogrid update]
}

/*
 * Update signed distance transform of outer layer points.
 */
BOOST_AUTO_TEST_CASE(distance_transform)
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
		BOOST_CHECK_CLOSE(dist, 0.4f, 0.0001f);
	}
	// Update seed point by less than |0.5| and check outer layer
	// distances are updated.
	{
		typedef Surface<2, 2> SurfaceT;
		SurfaceT surface(Vec2u(5, 5));
		SurfaceT::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec2i(0, 0));

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
					BOOST_CHECK_EQUAL(surface.delta().get(pos), 0);

		}
		surface.update_end();


		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();

		BOOST_CHECK_EQUAL(diff, 0);
	}
}

/*
 * Update layers.
 */
BOOST_AUTO_TEST_CASE(layer_update)
{
	typedef Surface<2, 2> Surface_t;
	Surface_t surface(Vec2u(9, 9));
	Surface_t::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec2i(0, 0));
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

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 4);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 8);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 12);
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

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
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

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 4);

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
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

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 1);

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
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

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 0);

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
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

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 0);

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}
}

/*
 * Iterating the zero-layer
 */
BOOST_AUTO_TEST_CASE(iterate_layers)
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

	BOOST_CHECK_EQUAL(surface.layer(0).size(), 6);

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
				BOOST_CHECK_EQUAL(val, 0);
				counter++;
				pos_sum += pos;
			}
		}
	};
	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));


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
					BOOST_CHECK_EQUAL(val, layer_id);
					counter++;
					pos_sum += pos;
				}
			}

	BOOST_CHECK_EQUAL(counter, 63);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

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
	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

	// Iterate over zero-layer using range based for loop.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	for (auto pos : surface.layer(0))
	{
		pos_sum += pos;
		counter++;
	}

	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

}

/*
 * Check that isogrid grid is bounded, that is, we cannot cause the surface to
 * attempt to leave the isogrid embedding.
 */
BOOST_AUTO_TEST_CASE(check_bounded)
{
	typedef Surface<2, 2> SurfaceT;
	SurfaceT surface(Vec2u(9, 9));
	SurfaceT::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec2i(0, 0));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> isogrid_check(Vec2u(9, 9), Vec2i::Zero(), 0);

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

	// Attempt to expand to outside the grid.
	// delta isogrid should be modified from -1.0 to approx -0.5.
	surface.update_start();
	{
		for (auto pos : surface.layer(0))
		{
			surface.delta(pos, -1.0f);
		}
	}
	surface.update_end();

	// Try to expand again.
	// delta isogrid should be modified from -1.0 to 0.
	surface.update_start();
	{
		for (auto pos : surface.layer(0))
		{
			surface.delta(pos, -1.0f);
		}
	}
	surface.update_end();

	// Test it.
	{
//		std::cerr << surface.isogrid().data() << std::endl;
		isogrid_check.data() = {
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1.5, 3, 3, 3, 3, 3, 3, 3,
			1.5, 0.5, 1.5, 3, 3, 3, 3, 3, 1.5, 0.5, -0.5, 0.5, 1.5, 3, 3, 3,
			1.5, 0.5, -0.5, -1.5, -0.5, 0.5, 1.5, 3, 3, 3, 1.5, 0.5, -0.5, 0.5,
			1.5, 3, 3, 3, 3, 3, 1.5, 0.5, 1.5, 3, 3, 3, 3, 3, 3, 3, 1.5, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		// isogrid_check uses 'whole' 0.5s, but internally, to prevent rounding,
		// max isogrid at grid boundary is 0.5-epsilon*2.
		BOOST_CHECK_SMALL(
			diff, std::numeric_limits<FLOAT>::epsilon() * 7 * 7 * 2
		);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 4);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 8);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 12);
	}
}

template <typename... Args>
auto s(Args&&... args) -> decltype(std::to_string(std::forward<Args>(args)...)) {
  return std::to_string(std::forward<Args>(args)...);
}

BOOST_AUTO_TEST_CASE(affected_outer_layers)
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

		BOOST_CHECK_MESSAGE(
			surface.affected().leafs(layer_idx).size() == check_layers_pos[layer_idx].size(),
			"Layer " + s(layer_id) + " at index " + s(layer_idx) +
			" number of leafs " +
			s(surface.affected().leafs(layer_idx).size()) + " == " +
			s(check_layers_pos[layer_idx].size())
		);

		for (auto pos : check_layers_pos[layer_idx])
		{
			auto iter = std::find(
				surface.affected().leafs(layer_idx).begin(),
				surface.affected().leafs(layer_idx).end(), pos
			);

			BOOST_CHECK_MESSAGE(
				iter != surface.affected().leafs(layer_idx).end(),
				"Affected grid layer " + s(layer_id) + " at index " + s(layer_idx) +
				" should contain (" + s(pos(0)) + "," + s(pos(1)) + ")"
			);
		}

		for (auto pos : surface.affected().leafs(layer_idx))
		{
			auto iter = std::find(
				check_layers_pos[layer_idx].begin(),
				check_layers_pos[layer_idx].end(), pos
			);

			BOOST_CHECK_MESSAGE(
				iter != check_layers_pos[layer_idx].end(),
				"Checking list layer " + s(layer_id) + " at index " + s(layer_idx) +
				" should contain (" + s(pos(0)) + "," + s(pos(1)) + ")"
			);
		}
	}
	//![Calculate affected outer layers for localised narrow band updates]
}

/*
 * Localised update.
 */
BOOST_AUTO_TEST_CASE(local_update)
{
	typedef Surface<2, 2> SurfaceT;
	SurfaceT surface(Vec2u(9, 9));
	SurfaceT::IsoGrid::Child& isogrid = surface.isogrid().children().get(Vec2i(0, 0));
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

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 4);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 8);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 12);
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

		isogrid_check.vdata() = isogrid_check.vdata() - isogrid.vdata();
		const FLOAT diff = isogrid_check.vdata().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}
}

/**
 * Test walking the zero layer out to a given distance.
 */
BOOST_AUTO_TEST_CASE(walk)
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
		SharedLookupGrid<2, surface.NUM_LAYERS> lookup = surface.walk_band<2>(
			Vec2i(-3,0)
		);

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-2)).size(), 1);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-1)).size(), 1);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0)).size(), 3);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(1)).size(), 3);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(2)).size(), 5);

		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0))[0], Vec2i(-3, 0));
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0))[1], Vec2i(-2, -1));
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0))[2], Vec2i(-2, 1));
	}
	{
		// ==== Setup ====
		Surface<3, 2> surface(Vec3u(9, 9, 9));
		using Lookup = SharedLookupGrid<3, surface.NUM_LAYERS>;

		// Create seed point and expand the narrow band.
		surface.seed(Vec3i(0, 0, 0));

		// ==== Action ====
		Lookup& lookup1 = surface.walk_band<1>(Vec3i(0,0,0));

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(lookup1.list(surface.layer_idx(-2)).size(), 0);
		BOOST_CHECK_EQUAL(lookup1.list(surface.layer_idx(-1)).size(), 0);
		BOOST_CHECK_EQUAL(lookup1.list(surface.layer_idx(0)).size(), 1);
		BOOST_CHECK_EQUAL(lookup1.list(surface.layer_idx(1)).size(), 6);
		BOOST_CHECK_EQUAL(lookup1.list(surface.layer_idx(2)).size(), 0);

		// ==== Action ====
		Lookup& lookup2 = surface.walk_band<1>(Vec3i(0,0,0));

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(&lookup1, &lookup2);

		// ==== Action ====
		Lookup& lookup3 = surface.walk_band<2>(Vec3i(0,0,0));

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(&lookup1, &lookup2);
		BOOST_CHECK_NE(&lookup1, &lookup3);
		BOOST_CHECK_NE(&lookup2, &lookup3);
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

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
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
		SharedLookupGrid<3, surface.NUM_LAYERS> lookup = surface.walk_band<1>(
			Vec3i(0,0,0)
		);

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-2)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-1)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(1)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(2)).size(), 0);

		// ==== Action ====
		lookup = surface.walk_band<2>(Vec3i(-5,0,0));

		// ==== Contirm ====
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-2)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(-1)).size(), 0);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(0)).size(), 1);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(1)).size(), 1);
		BOOST_CHECK_EQUAL(lookup.list(surface.layer_idx(2)).size(), 5);

		BOOST_CHECK_EQUAL(lookup(Vec3i(-4, 0, 0)), 0);
		BOOST_CHECK_EQUAL(lookup(Vec3i(-3, 0, 0)), 0);
		BOOST_CHECK_PREDICATE(
			[](const UINT idx) {
				return 0 <= idx <= 4;
			}, (lookup(Vec3i(-5, 0, 0)))
		);
		BOOST_CHECK_EQUAL(lookup(Vec3i(-6, 0, 0)), lookup.NULL_IDX);

	}
}
	
/**
 * Test delta isogrid update spread using Gaussian distribution given a list of
 * positions to update.
 */
BOOST_AUTO_TEST_CASE(gaussian_from_list)
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

	SharedLookupGrid<2, surface.NUM_LAYERS> lookup = surface.walk_band<2>(
		Vec2i(-3,0)
	);

	// ==== Action ====
	surface.update_start();
	surface.delta_gauss(
		lookup.list(surface.layer_idx(0)), Vec2f(-3.5f, 0), 0.5f, 0.2f
	);
	surface.update_end();
		
	
	// === Confirm ===

	BOOST_CHECK_CLOSE(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1)),
		0.5f, 0.0000001f
	);

	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-3, 0)), 0.3457f, 0.0001f
	);
	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-2, -1)), 0.07714f, 0.0001f
	);
	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-2, 1)), 0.07714f, 0.0001f
	);
}

/**
 * Test delta isogrid update spread using Gaussian distribution given point and
 * radius.
 */
BOOST_AUTO_TEST_CASE(gaussian_from_dist)
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

	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1)),
		0.5f, 0.0000001f
	);

	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-3, 0)), 0.28805843f, 0.0001f
	);
	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-2, -1)), 0.105970778f, 0.0001f
	);
	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-2, 1)), 0.105970778f, 0.0001f
	);
}

/**
 * Test raycasting to zero curve.
 */
BOOST_AUTO_TEST_CASE(ray)
{
	// ==== Setup ====
	Surface<3, 3> surface(Vec3u(16, 16, 16), Vec3u(5, 5, 5));
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

//	BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
/*
/home/dave/Dropbox/Workspace/Felt/src/tests/test_Surface.cpp(1,115): Message:
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
	// Simplest "dead on" case - from outside grid.
	pos_hit = surface.ray(Vec3f(-100.0f, 0, 0), Vec3f(1, 0, 0));
	
	// ==== Confirm ====
	BOOST_CHECK_LE(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm(), 0.00001f
	);

	// ==== Action ====
	// Simplest "dead on" case - from inside grid.
	pos_hit = surface.ray(Vec3f(-6.0f, 0, 0), Vec3f(1, 0, 0));
	
	// ==== Confirm ====
	BOOST_CHECK_LE(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm(), 0.00001f
	);

	// ==== Action ====
	// Simplest "dead on" case - from inside surface.
	pos_hit = surface.ray(Vec3f(0.0f, 0, 0), Vec3f(1, 0, 0));
	
	// ==== Confirm ====
	BOOST_CHECK_EQUAL(pos_hit, surface.NULL_POS<FLOAT>());
	
	// ==== Action ====
	// Simplest "dead on" case - from zero layer.
	pos_hit = surface.ray(Vec3f(-3.0f, 0, 0), Vec3f(1, 0, 0));
	
	// ==== Confirm ====
	BOOST_CHECK_LE(
		(pos_hit - Vec3f(-3.0f, 0, 0)).squaredNorm(), 0.00001f
	);
	
	// ==== Setup ====
	surface.update([](auto& pos, auto& isogrid) {
		return -0.3f;
	});

	// ==== Action ====
	// Ray  interpolate to zero curve.
	pos_hit = surface.ray(Vec3f(-10.0f, 0, 0), Vec3f(1, 0, 0));

	// ==== Confirm ====
	BOOST_CHECK_LE(
		(pos_hit - Vec3f(-3.3f, 0, 0)).squaredNorm(), 0.00001f
	);
	
	// ==== Setup ====
	surface.update([](auto& pos, auto& isogrid) {
		return 0.3f;
	});

	// ==== Action ====
	// Ray at an angle.
	pos_hit = surface.ray(
		Vec3f(-10.0f, -10.0f, 0.0f), Vec3f(1, 1, 0).normalized()
	);

	// ==== Confirm ====
	BOOST_CHECK_LE(
		(pos_hit - Vec3f(-1.5, -1.5, 0)).squaredNorm(),
		0.00001f
	);

	pos_hit = surface.ray(Vec3f(10, 10, 10), Vec3f(-1, -1, -1).normalized());

	// ==== Confirm ====
	BOOST_CHECK_NE(pos_hit, surface.NULL_POS<FLOAT>());

	// ==== Action ====
	// Rotating ray.

	pos_hit = surface.ray(Vec3f(6.72, -6.55, -3.45), Vec3f(-0.672, 0.655, 0.345));

	BOOST_CHECK_NE(pos_hit, surface.NULL_POS<FLOAT>());

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		Eigen::Matrix3f mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f::UnitY()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10.0f);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		BOOST_CHECK_MESSAGE(
			pos_hit != surface.NULL_POS<FLOAT>(),
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
	}

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		Eigen::Matrix3f mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f(1, 1, 1).normalized()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		BOOST_CHECK_MESSAGE(
			pos_hit != surface.NULL_POS<FLOAT>(),
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
	}

	for (FLOAT rot_mult = 0; rot_mult < 2.0f; rot_mult += 0.1f)
	{
		// ==== Setup ====
		Eigen::Matrix3f mat_rot = Eigen::AngleAxisf(
			rot_mult * M_PI, Vec3f(0, 1, 1).normalized()
		).matrix();
		const Vec3f origin = mat_rot*Vec3f(0, 0, -10);
		const Vec3f dir = (mat_rot*Vec3f(0, 0, 1)).normalized();

		// ==== Action ====
		pos_hit = surface.ray(origin, dir);

		// ==== Confirm ====
		BOOST_CHECK_MESSAGE(
			pos_hit != surface.NULL_POS<FLOAT>(),
			"Ray hit from " + stringifyVector(origin) + " in direction "
			+ stringifyVector(dir) + " should not be NULL_POS"
		);
	}
	// ==== Confirm ====
//	BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
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


BOOST_AUTO_TEST_CASE(gaussian_from_ray)
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
	FLOAT leftover;


	// ==== Action ====
	surface.update_start();
	leftover = surface.delta_gauss<2>(
		Vec2f(-2.4f, -10.0f), Vec2f(0, 1), 0.5f, 0.2f
	);
	surface.update_end();


	// === Confirm ===
//	BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
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
	BOOST_CHECK_CLOSE_FRACTION(
		surface.delta().get(Vec2i(-3, 0))
		+ surface.delta().get(Vec2i(-2, 1))
		+ surface.delta().get(Vec2i(-2, -1))
		+ surface.delta().get(Vec2i(-1, -2)),
		0.5f, 0.000001f
	);

	BOOST_CHECK_LE(leftover, 0.000001f);

//	BOOST_CHECK_CLOSE_FRACTION(
//		surface.delta().get(Vec2i(-3, 0)), 0.152139202f, 0.0001f
//	);
//	BOOST_CHECK_CLOSE_FRACTION(
//		surface.delta().get(Vec2i(-2, -1)), 0.24048205f, 0.0001f
//	);
//	BOOST_CHECK_CLOSE_FRACTION(
//		surface.delta().get(Vec2i(-2, 1)), 0.0559686795, 0.0001f
//	);
//	BOOST_CHECK_CLOSE_FRACTION(
//		surface.delta().get(Vec2i(-1, -2)), 0.0514338501, 0.0001f
//	);
}

BOOST_AUTO_TEST_SUITE_END()
