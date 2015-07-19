#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include <array>
#include <string>

#define _TESTING

#include "Felt/Surface.hpp"
#include "Felt/Poly.hpp"
#include "Felt/PolyGrid.hpp"

#include "Utils.hpp"

using namespace felt;

/**
 * Tests for the PolyGrid class.
 */
BOOST_AUTO_TEST_SUITE(test_PolyGrid)

	/**
	 * Test basic initialisation.
	 */
	BOOST_AUTO_TEST_CASE(initialise)
	{
		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 2> Surface_t;
		// Initialise a surface.
		Surface<3, 2> surface(Vec3u(9,9,9), Vec3u(3, 3, 3));
		PolyGrid<3> poly(surface);

		BOOST_CHECK_EQUAL((UINT)poly.data().size(), 27);
//		BOOST_CHECK_EQUAL(poly.lookup().data().size(), 27);

//		BOOST_CHECK_EQUAL((UINT)poly.changes().branch().data().size(), 27);
	}


	/**
	 * Utility function to assert tracked changes are as expected.
	 *
	 * @param poly
	 * @param apos_expected
	 */
	template <UINT N>
	void assert_expected_changes_tracked (
		const PolyGrid<3>& poly,  const std::array<Vec3i, N>& apos_expected
	) {
		typedef PolyGrid<3> PolyGrid_t;

		BOOST_CHECK_EQUAL(poly.changes().branch().list().size(), N);

		{
			auto begin_iter = poly.changes().branch().list().begin();
			auto end_iter = poly.changes().branch().list().end();
			for (auto pos : apos_expected)
				BOOST_CHECK_MESSAGE(
					std::find(begin_iter, end_iter, pos) != end_iter,
					stringifyVector(pos) << " was expected but not found."
				);
		}
		{
			auto begin_iter = apos_expected.begin();
			auto end_iter = apos_expected.end();
			for (auto pos : poly.changes().branch().list())
				BOOST_CHECK_MESSAGE(
					std::find(begin_iter, end_iter, pos) != end_iter,
					stringifyVector(pos) << " was found but unexpected."
				);
		}
	}

	/**
	 * Test that changes to the underlying surface are tracked as expected.
	 */
	BOOST_AUTO_TEST_CASE(changes_expand)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 2> Surface_t;
		// Initialise a surface.
		Surface<3, 2> surface(Vec3u(15,15,15), Vec3u(5, 5, 5));
		PolyGrid<3> poly(surface);

		// Initialise a seed.
		surface.seed(Vec3i(0,0,0));
		// TODO: dummy zero-layer update require to initialise dphi grid.
		surface.update_start();
		surface.dphi(Vec3i(0,0,0), 0);
		surface.update_end();
/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |

*/
		// ==== Action ====

		poly.notify(surface);


		// ==== Confirm ====

		assert_expected_changes_tracked(
			poly, std::array<Vec3i, 1>{ Vec3i(0,0,0) }
		);


		// ==== Action ====

		// Expand the surface outward.
		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1);
		surface.update_end();
		poly.notify(surface);
		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1);
		surface.update_end();
		poly.notify(surface);
		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1);
		surface.update_end();
		// Notify will add new tracking points using status change list rather
		// than delta phi tracking list (delta phi is still within central
		// partition).
		poly.notify(surface);


//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.phi()));

/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|  * 3 |    3 |    3 |    3 |    3 |  * 2 |    1 |    0 |    1 |    2 |  * 3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -3 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|  * 3 |    3 |    3 |    3 |    2 |  * 1 |    0 |   -1 |    0 |    1 |  * 2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|  * 3 |    3 |    3 |    3 |    3 |  * 3 |    3 |    3 |    3 |    3 |  * 3 |    3 |    3 |    3 |    3 |

*/
		// ==== Confirm ====

		assert_expected_changes_tracked(
			poly, std::array<Vec3i, 7>{
				Vec3i(0,0,0), Vec3i(0,0,-1), Vec3i(0,0,1), Vec3i(0,-1,0),
				Vec3i(0,1,0), Vec3i(-1,0,0), Vec3i(1,0,0)
			}
		);


		// ==== Action ====

		surface.update_start();
		surface.dphi(Vec3i(0,-3,0), 1);
		surface.update_end();
		poly.notify(surface);

/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |

*/

		// ==== Confirm ====

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.phi()));

		assert_expected_changes_tracked(
			poly, std::array<Vec3i, 6>{
				Vec3i(0,0,0), Vec3i(0,0,-1), Vec3i(0,0,1),
				Vec3i(0,1,0), Vec3i(-1,0,0), Vec3i(1,0,0)
			}
		);
	}


	/**
	 * Test (re-)polygonisations based on tracked changes.
	 */
	BOOST_AUTO_TEST_CASE(poly_cubes)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 2> Surface_t;
		// Initialise a surface.
		Surface<3, 2> surface(Vec3u(15,15,15), Vec3u(5, 5, 5));
		PolyGrid<3> poly(surface);
		Poly<3> poly_single(surface.phi().dims(), surface.phi().offset());

		// Initialise a seed.
		surface.seed(Vec3i(0,0,0));
		surface.update_start();
		surface.dphi(Vec3i(0,0,0), -1.0f);
		surface.update_end();

		// ==== Action ====

		poly.notify(surface);
		poly.poly_cubes(surface);

		// ==== Confirm ====

		poly_single.surf(surface);

		BOOST_CHECK_EQUAL(poly.get(Vec3i(0,0,0)).vtx().size(), 30);
		BOOST_CHECK_EQUAL(poly.get(Vec3i(0,0,0)).spx().size(), 56);
		BOOST_CHECK_EQUAL(
			poly.get(Vec3i(0,0,0)).vtx().size(),
			poly_single.vtx().size()
		);
		BOOST_CHECK_EQUAL(
			poly.get(Vec3i(0,0,0)).spx().size(),
			poly_single.spx().size()
		);

		// ==== Action ====

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1);
		surface.update_end();

		poly.notify(surface);
		poly.poly_cubes(surface);

/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/
		// ==== Confirm ====

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.phi()));
		poly_single.reset();
		poly_single.surf(surface);

		UINT total_vtx = 0;
		UINT total_spx = 0;
		for (const Vec3i& pos_child : poly)
		{
			total_vtx += poly.get(pos_child).vtx().size();
			if (poly.get(pos_child).vtx().size() > 0)
				BOOST_TEST_MESSAGE(
					stringifyVector(pos_child)
					+ " "
					+ std::to_string(poly.get(pos_child).vtx().size())
				);
			total_spx += poly.get(pos_child).spx().size();
		}

		// Total simplices should be the same.
		BOOST_CHECK_EQUAL(total_spx, poly_single.spx().size());

		// Total vertices will have duplicates at the border of the spatial
		// partitions.
		// The 'tip' of the shape at the three lowest corners (5 vertices making
		// up a pyramid) will be outside the central partition. The central
		// partition will thus have three points missing, one at each extremity,
		// since they fall entirely outside the partition.  Thus 4x4 = 12
		// vertices are duplicates of another 12 across the partition lines.
		// So, 12 duplicates + 3 end points - 3 cut from the central partition.
		BOOST_CHECK_EQUAL(total_vtx, poly_single.vtx().size() + 12 + 3 - 3);
		// As mentioned above, each lower extremity non-central partition has 5
		// vertices, making up the endpoint pyramids at those extremities.
		BOOST_CHECK_EQUAL(poly.get(Vec3i(-1,0,0)).vtx().size(), 5);
		BOOST_CHECK_EQUAL(poly.get(Vec3i(0,-1,0)).vtx().size(), 5);
		BOOST_CHECK_EQUAL(poly.get(Vec3i(0,0,-1)).vtx().size(), 5);


		// ==== Action ====

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -0.3);
		surface.update_end();
		poly.notify(surface);
		surface.update_start();
		surface.dphi(Vec3i(0,-2,0), 1);
		surface.update_end();
		poly.notify(surface);

		poly.poly_cubes(surface);
/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 | -1.3 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |  2.7 |  1.7 |  0.7 | -0.3 | -1.3 | -1.3 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 | -1.3 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/
		// ==== Confirm ====

		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.phi()));
		poly_single.reset();
		poly_single.surf(surface);

		total_vtx = 0;
		total_spx = 0;
		for (const Vec3i& pos_child : poly)
		{
			total_vtx += poly.get(pos_child).vtx().size();
			if (poly.get(pos_child).vtx().size() > 0)
				BOOST_TEST_MESSAGE(
					stringifyVector(pos_child)
					+ " "
					+ std::to_string(poly.get(pos_child).vtx().size())
				);
			total_spx += poly.get(pos_child).spx().size();
		}

		// Total simplices should be the same.
		BOOST_CHECK_EQUAL(total_spx, poly_single.spx().size());
		// One of the 'tips' have been pushed back into the central partition,
		// So, now just 8 duplicates + 2 endpoints - 2 cut from the central
		// partition.
		BOOST_CHECK_EQUAL(total_vtx, poly_single.vtx().size() + 8 + 2 - 2);


		// ==== Action ====
		poly.reset();


		// ==== Confirm ====
		total_vtx = 0;
		total_spx = 0;
		for (const Vec3i& pos_child : poly)
		{
			total_vtx += poly.get(pos_child).vtx().size();
			BOOST_TEST_MESSAGE(
				stringifyVector(pos_child)
				+ " "
				+ std::to_string(poly.get(pos_child).vtx().size())
			);
			total_spx += poly.get(pos_child).spx().size();
		}
		BOOST_CHECK_EQUAL(total_vtx, 0);
		BOOST_CHECK_EQUAL(total_spx, 0);
		BOOST_CHECK_EQUAL(poly.changes().leafs().size(), 0);
	}

	/**
	 * Test (re-)polygonisations based on tracked changes.
	 *
	 * Test case similar to above, but failed originally because of std::vector
	 * reinitialisation invalidating references during poly_cubes.
	 */
	BOOST_AUTO_TEST_CASE(poly_cubes_2)
	{
		// ==== Setup ====
		Surface<3> surface(Vec3u(13,13,13), Vec3u(4,4,4));

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1.0);
		surface.update_end();

		PolyGrid<3> polys(surface);
		Poly<3> poly(surface.phi().dims(), surface.phi().offset());

		surface.seed(Vec3i(0,0,0));

		// ==== Action ====

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1.0);
		surface.update_end();

		polys.notify(surface);

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.dphi(pos, -1.0);
		surface.update_end();

		polys.notify(surface);
		polys.poly_cubes(surface);


		// ==== Confirm ====
		poly.surf(surface);

		UINT total_vtx = 0;
		UINT total_spx = 0;
		for (const Vec3i& pos_child : polys)
		{
			total_vtx += polys.get(pos_child).vtx().size();
			if (polys.get(pos_child).vtx().size() > 0)
				BOOST_TEST_MESSAGE(
					stringifyVector(pos_child)
					+ " "
					+ std::to_string(polys.get(pos_child).vtx().size())
				);
			total_spx += polys.get(pos_child).spx().size();

			for (const Poly<3>::Simplex& polys_spx : polys.get(pos_child).spx())
			{
				Vec3f polys_vtxs[3];
				polys_vtxs[0] = polys.get(pos_child).vtx()[
					polys_spx.idxs(0)
				].pos;
				polys_vtxs[1] = polys.get(pos_child).vtx()[
					polys_spx.idxs(1)
				].pos;
				polys_vtxs[2] = polys.get(pos_child).vtx()[
					polys_spx.idxs(2)
				].pos;
				auto it = std::find_if(
					poly.spx().begin(), poly.spx().end(),
					[&](const Poly<3>::Simplex& poly_spx) {
						Vec3f poly_vtxs[3];
						poly_vtxs[0] = poly.vtx()[poly_spx.idxs(0)].pos;
						poly_vtxs[1] = poly.vtx()[poly_spx.idxs(1)].pos;
						poly_vtxs[2] = poly.vtx()[poly_spx.idxs(2)].pos;
						return (
							poly_vtxs[0] == polys_vtxs[0]
							&& poly_vtxs[1] == polys_vtxs[1]
							&& poly_vtxs[2] == polys_vtxs[2]
						);
					}
				);

				BOOST_CHECK_MESSAGE(
					it != polys.get(pos_child).spx().end(),
					(
						"Simplex "
						+ stringifyVector(polys_vtxs[0]) + "-"
						+ stringifyVector(polys_vtxs[1]) + "-"
						+ stringifyVector(polys_vtxs[2])
						+ " found in partition"
					)
				);
			}
		}

		BOOST_TEST_MESSAGE(std::to_string(total_spx) + " spxs");
		BOOST_TEST_MESSAGE(std::to_string(total_vtx) + " vtxs");

		BOOST_CHECK_EQUAL(total_spx, poly.spx().size());
	}
BOOST_AUTO_TEST_SUITE_END()
