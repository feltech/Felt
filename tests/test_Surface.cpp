#include <boost/numeric/ublas/io.hpp>
#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "Felt.hpp"
#include "Surface.hpp"

using namespace felt;

BOOST_AUTO_TEST_SUITE(test_Surface)
/*
 * Basic initialisation.
 */
BOOST_AUTO_TEST_CASE(init)
{
	// Basic initialisation of 2D surface with 2 layers in a 5x5 embedding.

	Surface<2, 2> surface(Vec2u(7, 7));
	{
		const Vec2u vec_dims = surface.dims();

		BOOST_CHECK_EQUAL((UINT)vec_dims(0), 7);
		BOOST_CHECK_EQUAL((UINT)vec_dims(1), 7);

		// Usable isogrid should have size=dims-layers and be offset by half
		// grid width.
		// In this test case, only the centre point is actually usable, all
		// other points are reserved for outer layers.

		const Vec2i pos_min = surface.pos_min();
		const Vec2i pos_max = surface.pos_max();

		BOOST_CHECK_EQUAL(pos_min, Vec2i(0, 0));
		BOOST_CHECK_EQUAL(pos_max, Vec2i(0, 0));

		// But actual phi isogrid should have size equal to dims.

		Grid<FLOAT, 2>& phi = surface.phi();

		const Vec2u phi_dims = phi.dims();
		BOOST_CHECK_EQUAL(phi_dims, Vec2u(7, 7));

		// Grid is initialised to all points 'outside' the surface (since there
		// is no surface yet).

		const FLOAT val_centre = phi(Vec2i(0, 0));
		BOOST_CHECK_EQUAL(val_centre, 3);
	}

	{
		// Check OpenMP thread support.

		const UINT num_threads = omp_get_max_threads();
		BOOST_WARN_GT(omp_get_max_threads(), 1);
		BOOST_TEST_MESSAGE("Num OpenMP threads: " << num_threads);

		BOOST_CHECK_EQUAL(omp_get_max_threads(), omp_get_num_procs());

		for (UINT threadIdx = 0; threadIdx < num_threads; threadIdx++)
		{
			std::vector<Vec2i>* pary_omp_dPhi = &surface.dphi(threadIdx);
			BOOST_CHECK_EQUAL(pary_omp_dPhi->size(), 0u);
		}

		// Check delta phi grid.

		Grid<FLOAT, 2>& dphi = surface.dphi();

		const Vec2u dphi_dims = dphi.dims();
		BOOST_CHECK_EQUAL(dphi_dims, Vec2u(7, 7));

		// Initialised to zero.
		const FLOAT val_centre = dphi(Vec2i(0, 0));
		BOOST_CHECK_EQUAL(val_centre, 0);
	}
}

/*
 * Narrow band layers.
 */
BOOST_AUTO_TEST_CASE(layers)
{
	// 3D surface with default (=2) number of layers.
	Surface<3> surface(Vec3u(7, 7, 7));
	Grid<FLOAT, 3>& phi = surface.phi();
	Grid<UINT, 3>& idx = surface.idx();

	Vec3i pos = Vec3i(0, 0, 0);

	BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
	BOOST_CHECK_EQUAL(surface.layer(2).size(), 0);

	// Check layer index lookup initialisation.
	BOOST_CHECK_EQUAL(idx(pos), surface.null_idx());

	// Add a single zero-layer point.
	phi(pos) = 0;
	surface.layer_add(0, pos);

	// Check zero-layer array has registered point.
	BOOST_CHECK_EQUAL(surface.layer(0).size(), 1);
	BOOST_CHECK_EQUAL((surface.layer(0)[0] - pos), Vec3i::Zero());

	// Check layer calculation from value.
	// -- zero-layer point just added.
	BOOST_CHECK_EQUAL(surface.layerID(pos), 0);

	// Check index grid has registered new zero-layer point.
	BOOST_CHECK_EQUAL(idx(pos), 0);

	// Add three arbitrary points to layer -1.
	surface.layer_add(-1, Vec3i(0, 0, 1));
	surface.layer_add(-1, Vec3i(0, 0, 2));
	surface.layer_add(-1, Vec3i(0, 0, 3));

	// Remove two points from layer -1
	surface.layer_remove(Vec3i(0, 0, 1), -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 2);
	surface.layer_remove(Vec3i(0, 0, 3), -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);

	// Move a point from layer 0 to layer -1
	surface.layer_move(pos, 0, -1);
	BOOST_CHECK_EQUAL(surface.layer(-1).size(), 2);

	// Arbitrary point @ 0, so moved point @ 1.
	BOOST_CHECK_EQUAL(idx(pos), 1);
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
	Grid<FLOAT, 2>& phi = surface.phi();

	surface.seed(Vec2i(0, 0));

	// Trivially check centre of seed is indeed a zero-level point (i.e. point
	// on the surface).

	const FLOAT val_centre = phi(Vec2i(0, 0));
	BOOST_CHECK_EQUAL(val_centre, 0);

	// A 2D 2-layer singularity (seed) point should look like the following.

	Grid<FLOAT, 2> phi_check(Vec2u(5, 5));
	phi_check.data() << 3, 3, 2, 3, 3,	// |
	3, 2, 1, 2, 3,	// -
	2, 1, 0, 1, 2,	// x
	3, 2, 1, 2, 3,	// +
	3, 3, 2, 3, 3;	// |
	//	|____ - y + ____|
//	std::cerr << phi.data() << std::endl << std::endl;
//	std::cerr << phi_check.data() << std::endl << std::endl;

	phi_check.data() = phi_check.data() - phi.data();
	const FLOAT diff = phi_check.data().sum();

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
	Grid<FLOAT, 2>& phi = surface.phi();

	surface.seed(Vec2i(0, 0));

	Vec2i pos_next = Vec2i(-1, -2);
	BOOST_CHECK_EQUAL(phi(pos_next), 3);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), 2);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), 1);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), 0);

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), 0);

	// Ensure it also works with negative distances.
	// NOTE: row-major (y,x) element ordering...
	surface.phi().data() << (
		2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -2,
		-2, -2, -2, -2
	);
	// NOTE: ...but accessed as (x,y)
	pos_next = Vec2i(2, 0);

	BOOST_CHECK_EQUAL(phi(pos_next), -2);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(2, 0));

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), -1);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(1, 0));

	pos_next = surface.next_closest(pos_next);

	BOOST_CHECK_EQUAL(phi(pos_next), 0);
	BOOST_CHECK_EQUAL(pos_next, Vec2i(0, 0));
}

/*
 * Using delta phi grid/list.
 */
BOOST_AUTO_TEST_CASE(delta_phi)
{
	// Basic non-threaded check.
	{
		Surface<3, 2> surface(Vec3u(5, 5, 5));
		Grid<FLOAT, 3>& dphi = surface.dphi();

		Vec3i pos = Vec3i(0, 0, 0);
		// Apply a delta to the surface.
		surface.dphi(pos, -2);
		// Check delta was stored in underlying grid.
		BOOST_CHECK_EQUAL(dphi(pos), -2);
		// Check position vector of point in surface grid that delta was
		// applied to is stored in a corresponding list to be iterated over.
		UINT sum = 0;
		for (UINT threadIdx = 0; threadIdx < surface.num_threads(); threadIdx++)
			sum += surface.dphi(threadIdx).size();
		BOOST_CHECK_EQUAL(sum, 1);
	}

	// Multi-threaded check.
	{
		Surface<3, 2> surface(Vec3u(5, 5, 5));
		BOOST_WARN_MESSAGE(
			omp_get_max_threads() >= 1,
			(
				"only " << omp_get_max_threads() <<
				" OpenMP thread available, not a good test of OpenMP."
			)
		);
		surface.num_threads(omp_get_max_threads());
		Grid<FLOAT, 3>& dphi = surface.dphi();

		#pragma omp parallel for// num_threads(4)
		for (UINT threadIdx = 0; threadIdx < surface.num_threads(); threadIdx++)
		{
			#pragma omp critical
			BOOST_REQUIRE_EQUAL(omp_get_thread_num(), threadIdx);

			Vec3i pos = Vec3i(0, 0, threadIdx);
			surface.dphi(pos, threadIdx + 1);
		}

		for (UINT threadIdx = 0; threadIdx < surface.num_threads(); threadIdx++)
		{
			Vec3i pos = Vec3i(0, 0, threadIdx);
			BOOST_CHECK_EQUAL(dphi(pos), threadIdx + 1);

			std::vector<Vec3i>& apos_thread = surface.dphi(threadIdx);

			UINT uThreadSize = apos_thread.size();
			BOOST_REQUIRE_EQUAL(uThreadSize, 1);

			Vec3i pos_thread = apos_thread[0];

			INT dist = (pos - pos_thread).sum();
			BOOST_CHECK_EQUAL(dist, 0);
		}
	}
}

/*
 * Update phi with delta phi.
 */
BOOST_AUTO_TEST_CASE(delta_phi_update)
{
	UINT sum;
	Surface<3, 2> surface(Vec3u(5, 5, 5));
	omp_set_dynamic(0);
	omp_set_num_threads(4);
	surface.num_threads(4);

	Grid<FLOAT, 3>& phi = surface.phi();
	Grid<FLOAT, 3>& dphi = surface.dphi();

	// Put in 'dirty' state, to check update_start is doing it's job.
	surface.dphi(Vec3i(0, 0, 0), 1.0f);

	// Clear delta phi.
	surface.update_start();
	{
		// Check update_start cleared the above surface.dphi changes.
		for (UINT threadIdx = 0; threadIdx < surface.num_threads(); threadIdx++)
			BOOST_CHECK_EQUAL(surface.dphi(threadIdx).size(), 0u);
		BOOST_CHECK_EQUAL(dphi(Vec3i(0, 0, 0)), 0);
	}
	// Apply delta phi.
	surface.update_end();

	// Add a zero-layer point.
	surface.phi(Vec3i(0, 0, 0), 0.0f);

	// Clear delta phi.
	surface.update_start();
	{
		// Do nothing.
		surface.dphi(Vec3i(0, 0, 0), 0);
	}
	// Apply delta phi.
	surface.update_end();

	// Ensure nothing was changed.  Every point in 5x5x5 grid == 3, except
	// centre which == 0.
	BOOST_CHECK_EQUAL(phi.data().sum(), 3 * 5 * 5 * 5 - 3);
	// Delta phi position vector list should still contain one point.
	sum = 0;
	for (UINT threadIdx = 0; threadIdx < surface.num_threads(); threadIdx++)
		sum += surface.dphi(threadIdx).size();
	BOOST_CHECK_EQUAL(sum, 1u);
	// Delta phi grid itself should have reset back to zero.
	BOOST_CHECK_EQUAL(dphi(Vec3i(0, 0, 0)), 0);

	// Clear delta phi.
	surface.update_start();
	{
		// Apply small update.
		surface.dphi(Vec3i(0, 0, 0), 0.4f);
	}
	// Apply delta phi.
	surface.update_end();

	// Ensure change applied.  Every point in grid == 3, except centre which
	// == 0.4.
	BOOST_CHECK_EQUAL(phi.data().sum(), 3 * 5 * 5 * 5 - 3 + 0.4f);

	omp_set_dynamic(1);
}

/*
 * Update signed distance transform of outer layer points.
 */
BOOST_AUTO_TEST_CASE(distance_transform)
{
	// Check distance calculation for a single point.
	{
		Surface<3, 2> surface(Vec3u(5, 5, 5));
		Grid<FLOAT, 3>& phi = surface.phi();

		surface.seed(Vec3i(0, 0, 0));

		// Basic distance calculation.
		phi(Vec3i(0, 0, 0)) = -0.6f;
		const FLOAT dist = surface.distance(Vec3i(-1, 0, 0), 1);
		BOOST_CHECK_CLOSE(dist, 0.4f, 0.0001f);
	}
	// Update seed point by less than |0.5| and check outer layer
	// distances are updated.
	{
		Surface<2, 2> surface(Vec2u(5, 5));
		Grid<FLOAT, 2>& phi = surface.phi();

		surface.seed(Vec2i(0, 0));

		Grid<FLOAT, 2> phi_check(Vec2u(5, 5));
		phi_check.data() << (
			3, 3, 1.6, 3, 3, 3, 1.6, 0.6, 1.6, 3, 1.6, 0.6, -0.4, 0.6, 1.6, 3,
			1.6, 0.6, 1.6, 3, 3, 3, 1.6, 3, 3
		);

		surface.update_start();
		{
			surface.dphi(Vec2i(0, 0), -0.4f);
		}
		surface.update_end();

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();

		BOOST_CHECK_EQUAL(diff, 0);
	}
}

/*
 * Update layers.
 */
BOOST_AUTO_TEST_CASE(layer_update)
{
	Surface<2, 2> surface(Vec2u(9, 9));
	Grid<FLOAT, 2>& phi = surface.phi();
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> phi_check(Vec2u(9, 9));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update_start();
	{
		surface.dphi(Vec2i(0, 0), -0.6f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3, 3, 3, 3, 3, 3,
			2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 0.4, 1.4, 2.4, 3, 3, 3,
			2.4, 1.4, 0.4, -0.6, 0.4, 1.4, 2.4, 3, 3, 3, 2.4, 1.4, 0.4, 1.4,
			2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
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
		for (UINT posIdx = 0; posIdx < surface.layer().size(); posIdx++)
			surface.dphi(posIdx, 0.6f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 2, 1, 0, 1, 2, 3,
			3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}

	// Collapse the seed completely, leaving no zero-layer, only outer layers.
	surface.update_start();
	{
		for (UINT posIdx = 0; posIdx < surface.size(); posIdx++)
			surface.dphi(posIdx, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3,
			3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 4);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}

	// Collapse still further, so there is only the outermost layer.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (UINT posIdx = 0; posIdx < surface.size(); posIdx++)
			surface.dphi(posIdx, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 1);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}

	// Final collapse leaves the whole grid as 'outside' points.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (UINT posIdx = 0; posIdx < surface.layer().size(); posIdx++)
			surface.dphi(posIdx, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 0);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}

	// Further updates have no effect.
	surface.update_start();
	{
		// Has no effect, since zero-layer is gone size is 0.
		for (UINT posIdx = 0; posIdx < surface.size(); posIdx++)
			surface.dphi(posIdx, 1.0f);
	}
	surface.update_end();

	{
//		std::cerr << surface.phi().data() << std::endl << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 0);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}
}

/*
 * Iterating the zero-layer
 */
BOOST_AUTO_TEST_CASE(iterate_zero_layer)
{
	Surface<3> surface(Vec3u(9, 9, 9));

	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 3> phi_check(Vec3u(9, 9, 9));

	// Values to compute for testing against.
	int counter;
	Vec3i pos_sum;

	// Create seed point and expand the narrow band.
	surface.seed(Vec3i(0, 0, 0));
	surface.update_start();
	{
		surface.dphi(Vec3i(0, 0, 0), -1.0f);
	}
	surface.update_end();

	BOOST_CHECK_EQUAL(surface.size(), 6);

	// Iterate over surface, using parameterised index.
	// Only version that can be parallelised easily using OpenMP.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
#pragma omp parallel for
	for (UINT i = 0; i < surface.size(); i++)
	{
		Vec3i pos = surface[i];
		FLOAT val = surface(pos);
#pragma omp critical
		{
			BOOST_CHECK_EQUAL(val, 0);
			counter++;
			pos_sum += pos;
		}
	};
	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

	// Iterate over zero-layer using STL for_each and lambda callback function.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	std::for_each(surface.begin(), surface.end(), [&](Vec3i& pos)
	{
		pos_sum += pos;
		counter++;
	});
	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

	// Iterate over zero-layer using wrapped for_each.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	surface.each([&](Vec3i pos)
	{
		pos_sum += pos;
		counter++;
	});

	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

	// Iterate over zero-layer using range based for loop.
	counter = 0;
	pos_sum = Vec3i(0, 0, 0);
	for (auto pos : surface)
	{
		pos_sum += pos;
		counter++;
	}

	BOOST_CHECK_EQUAL(counter, 6);
	BOOST_CHECK_EQUAL(pos_sum, Vec3i(0, 0, 0));

}

/*
 * Check that phi grid is bounded, that is, we cannot cause the surface to
 * attempt to leave the phi embedding.
 */
BOOST_AUTO_TEST_CASE(check_bounded)
{
	Surface<2, 2> surface(Vec2u(9, 9));
	Grid<FLOAT, 2>& phi = surface.phi();
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> phi_check(Vec2u(9, 9));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update_start();
	{
		for (auto pos : surface)
		{
			surface.dphi(pos, -1.0f);
		}
	}
	surface.update_end();

	// Attempt to expand to outside the grid.
	// delta phi should be modified from -1.0 to approx -0.5.
	surface.update_start();
	{
		for (auto pos : surface)
		{
			surface.dphi(pos, -1.0f);
		}
	}
	surface.update_end();

	// Try to expand again.
	// delta phi should be modified from -1.0 to 0.
	surface.update_start();
	{
		for (auto pos : surface)
		{
			surface.dphi(pos, -1.0f);
		}
	}
	surface.update_end();

	// Test it.
	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1.5, 3, 3, 3, 3, 3, 3, 3,
			1.5, 0.5, 1.5, 3, 3, 3, 3, 3, 1.5, 0.5, -0.5, 0.5, 1.5, 3, 3, 3,
			1.5, 0.5, -0.5, -1.5, -0.5, 0.5, 1.5, 3, 3, 3, 1.5, 0.5, -0.5, 0.5,
			1.5, 3, 3, 3, 3, 3, 1.5, 0.5, 1.5, 3, 3, 3, 3, 3, 3, 3, 1.5, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		// phi_check uses 'whole' 0.5s, but internally, to prevent rounding,
		// max phi at grid boundary is 0.5-epsilon*2.
		BOOST_CHECK_SMALL(diff,
				std::numeric_limits<FLOAT>::epsilon() * 7 * 7 * 2);

		BOOST_CHECK_EQUAL(surface.layer(0).size(), 4);
		BOOST_CHECK_EQUAL(surface.layer(-1).size(), 1);
		BOOST_CHECK_EQUAL(surface.layer(-2).size(), 0);
		BOOST_CHECK_EQUAL(surface.layer(1).size(), 8);
		BOOST_CHECK_EQUAL(surface.layer(2).size(), 12);
	}
}

BOOST_AUTO_TEST_CASE(affected_outer_layers)
{
	Surface<2, 2> surface(Vec2u(9, 9));
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> phi_check(Vec2u(9, 9));

	// Create seed point and expand the narrow band.
	surface.seed(Vec2i(0, 0));
	surface.update_start();
	{
		for (auto pos : surface)
		{
			surface.dphi(pos, -1.0f);
		}
	}
	surface.update_end();

	surface.update_start();
	{
		surface.dphi(Vec2i(0, 1), 0.3f);
		surface.dphi(Vec2i(1, 0), 0.3f);

		std::vector<Vec2i> aAffected[5];
		surface.affected(aAffected);

//		3.0,	3.0,	3.0,	 2.0,	3.0,	3.0,	3.0,
//		3.0,	3.0,	2.0,	 1.0,	2.0,	3.0,	3.0,
//		3.0,	2.0,	1.0,	 0.0,	1.0,	2.0,	3.0,
//		2.0,	1.0,	0.0,	-1.0,	0.3,	1.0,	2.0,
//		3.0,	2.0,	1.0,	 0.3,	1.0,	2.0,	3.0,
//		3.0,	3.0,	2.0,	 1.0,	2.0,	3.0,	3.0,
//		3.0,	3.0,	3.0,	 2.0,	3.0,	3.0,	3.0;

		std::vector<Vec2i> aposCheck[5];
		aposCheck[2 + -2] = std::vector<Vec2i>();
		aposCheck[2 + -1] = std::vector<Vec2i>(
		{ Vec2i(0, 0), });
		aposCheck[2 + 0] = std::vector<Vec2i>(
		{
		// We don't care for now abouve zero-layer points.
//			Vec2i(0,1),
//			Vec2i(1,0)
				});
		aposCheck[2 + 1] = std::vector<Vec2i>(
		{

		// For (0,1):
				Vec2i(-1, 1), Vec2i(1, 1), Vec2i(0, 2),

				// For (1,0):
				Vec2i(2, 0), Vec2i(1, -1) });

		aposCheck[2 + 2] = std::vector<Vec2i>(
		{

		// For (0,1):
				Vec2i(-2, 1), Vec2i(2, 1),

				Vec2i(-1, 2), Vec2i(1, 2),

				Vec2i(0, 3),

				// For (1,0):

				Vec2i(3, 0), Vec2i(1, -2), Vec2i(2, -1) });

		for (INT layerID = -2; layerID <= 2; layerID++)
		{
			if (layerID == 0)
				continue;

			const INT layerIdx = 2 + layerID;
//			std::cerr << "Affected points: layer " << layerID << std::endl;
			BOOST_CHECK_EQUAL(aAffected[layerIdx].size(),
					aposCheck[layerIdx].size());

			for (auto pos : aAffected[layerIdx])
			{
				auto iter = std::find(aposCheck[layerIdx].begin(),
						aposCheck[layerIdx].end(), pos);
				BOOST_CHECK(iter != aposCheck[layerIdx].end());
			}
		}
	}
	surface.update_end();
}

/*
 * Localised update.
 */
BOOST_AUTO_TEST_CASE(local_update)
{
	Surface<2, 2> surface(Vec2u(9, 9));
	Grid<FLOAT, 2>& phi = surface.phi();
	// Grid to set values of manually, for checking against.
	Grid<FLOAT, 2> phi_check(Vec2u(9, 9));

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
		surface.dphi(Vec2i(0, 0), -0.6f);
	}
	// Using localised update, which will only update outer layers that are
	// affected by changes to the modified zero layer points.  In this test
	// case, all outer layer points are affected, same as a global update.
	surface.update_end_local();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3, 3, 3, 3, 3, 3,
			2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 0.4, 1.4, 2.4, 3, 3, 3,
			2.4, 1.4, 0.4, -0.6, 0.4, 1.4, 2.4, 3, 3, 3, 2.4, 1.4, 0.4, 1.4,
			2.4, 3, 3, 3, 3, 3, 2.4, 1.4, 2.4, 3, 3, 3, 3, 3, 3, 3, 2.4, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
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
		for (UINT posIdx = 0; posIdx < surface.size(); posIdx++)
			surface.dphi(posIdx, 0.6f);
	}
	surface.update_end_local();

	{
//		std::cerr << surface.phi().data() << std::endl;
		phi_check.data() << (
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 2, 1, 0, 1, 2, 3,
			3, 3, 3, 3, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3,
			3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
		);

		phi_check.data() = phi_check.data() - phi.data();
		const FLOAT diff = phi_check.data().sum();
		BOOST_CHECK_SMALL(diff, 0.000001f);
	}
}

BOOST_AUTO_TEST_SUITE_END()
