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
 * @ingroup Tests
 * @defgroup PolygonisationTests
 * @{
 * 	@name PolyGrid
 * 	@ref felt::PolyGrid
 */
BOOST_AUTO_TEST_SUITE(test_PolyGrid)

	/**
	 * Test basic initialisation.
	 */
	BOOST_AUTO_TEST_CASE(initialise)
	{
		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a surface.
		Surface_t surface(Vec3u(9,9,9), Vec3u(3, 3, 3));
		PolyGrid<3> poly(surface);

		BOOST_CHECK_EQUAL((UINT)poly.data().size(), 27);
	}

	/**
	 * Utility: assert PolyGrid matches simple Poly polygonisation.
	 *
	 * @param polys
	 * @param poly
	 * @return
	 */
	template <UINT D>
	UINT assert_partitioned_matches_baseline (
		const PolyGrid<D>& polys, const Poly<D>& poly
	) {

		UINT total_vtx = 0;
		UINT total_spx = 0;
		for (const Vec3i& pos_child : polys)
		{
			total_vtx += polys.get(pos_child).vtx().size();
			total_spx += polys.get(pos_child).spx().size();

			if (polys.get(pos_child).vtx().size() > 0)
				BOOST_TEST_MESSAGE(
					"Partition "
					+ stringifyVector(pos_child)
					+ " vtxs = "
					+ std::to_string(polys.get(pos_child).vtx().size())
					+ ", spxs = "
					+ std::to_string(polys.get(pos_child).spx().size())
				);

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
						"Simplex from partition "
						+ stringifyVector(polys_vtxs[0]) + "-"
						+ stringifyVector(polys_vtxs[1]) + "-"
						+ stringifyVector(polys_vtxs[2])
						+ " found in baseline"
					)
				);
			}
		}


		for (const Poly<3>::Simplex& poly_spx : poly.spx())
		{
			Vec3f poly_vtxs[3];
			poly_vtxs[0] = poly.vtx()[poly_spx.idxs(0)].pos;
			poly_vtxs[1] = poly.vtx()[poly_spx.idxs(1)].pos;
			poly_vtxs[2] = poly.vtx()[poly_spx.idxs(2)].pos;

			bool found_match = false;

			for (const Vec3i& pos_child : polys)
			{
				for (
					const Poly<3>::Simplex& polys_spx
					: polys.get(pos_child).spx()
				) {
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

					found_match = (
						poly_vtxs[0] == polys_vtxs[0]
						&& poly_vtxs[1] == polys_vtxs[1]
						&& poly_vtxs[2] == polys_vtxs[2]
					);
					if (found_match)
						break;
				}
				if (found_match)
					break;
			}

			BOOST_CHECK_MESSAGE(
				found_match,
				(
					"Simplex from baseline "
					+ stringifyVector(poly_vtxs[0]) + "-"
					+ stringifyVector(poly_vtxs[1]) + "-"
					+ stringifyVector(poly_vtxs[2])
					+ " found in partition"
				)
			);
		}

		BOOST_TEST_MESSAGE("Total: " + std::to_string(total_spx) + " spxs");
		BOOST_TEST_MESSAGE("Total: " + std::to_string(total_vtx) + " vtxs");

		BOOST_CHECK_EQUAL(total_spx, poly.spx().size());

		return total_vtx;
	}


	/**
	 * Test (re-)polygonisations based on tracked changes.
	 */
	BOOST_AUTO_TEST_CASE(poly_cubes)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a surface.
		Surface_t surface(Vec3u(15,15,15), Vec3u(5, 5, 5));
		PolyGrid<3> poly(surface);
		Poly<3> poly_single(surface.isogrid().size(), surface.isogrid().offset());

		// Initialise a seed.
		surface.seed(Vec3i(0,0,0));
		surface.update_start();
		surface.disogrid(Vec3i(0,0,0), -1.0f);
		surface.update_end();

		// ==== Action ====

		poly.notify(surface);
		poly.poly_cubes(surface);
		poly.update_end();

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
			surface.disogrid(pos, -1);
		surface.update_end();

		poly.notify(surface);
		poly.poly_cubes(surface);
		poly.update_end();

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

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
		poly_single.reset();
		poly_single.surf(surface);

		UINT total_vtx = assert_partitioned_matches_baseline(poly, poly_single);

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
			surface.disogrid(pos, -0.3);
		surface.update_end();
		poly.notify(surface);
		surface.update_start();
		surface.disogrid(Vec3i(0,-2,0), 1);
		surface.update_end();
		poly.notify(surface);

		poly.poly_cubes(surface);
		poly.update_end();
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

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));
		poly_single.reset();
		poly_single.surf(surface);

		total_vtx = assert_partitioned_matches_baseline(poly, poly_single);

		// One of the 'tips' have been pushed back into the central partition,
		// So, now just 8 duplicates + 2 endpoints - 2 cut from the central
		// partition.
		BOOST_CHECK_EQUAL(total_vtx, poly_single.vtx().size() + 8 + 2 - 2);


		// ==== Action ====
		poly.reset();

		// ==== Confirm ====
		total_vtx = 0;
		UINT total_spx = 0;
		for (const Vec3i& pos_child : poly)
		{
			total_vtx += poly.get(pos_child).vtx().size();
			total_spx += poly.get(pos_child).spx().size();
		}
		BOOST_CHECK_EQUAL(total_vtx, 0);
		BOOST_CHECK_EQUAL(total_spx, 0);
		BOOST_CHECK_EQUAL(poly.changes().list().size(), 0);
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
		Surface<3, 3> surface(Vec3u(13,13,13), Vec3u(4,4,4));
		PolyGrid<3> polys(surface);
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());

		surface.seed(Vec3i(0,0,0));

		// ==== Action ====

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.disogrid(pos, -1.0);
		surface.update_end();

		polys.notify(surface);

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.disogrid(pos, -1.0);
		surface.update_end();

		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/

		poly.surf(surface);

		assert_partitioned_matches_baseline(polys, poly);
	}


	BOOST_AUTO_TEST_CASE(poly_cubes_twice)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a surface with a a single partition.
		Surface_t surface(Vec3u(16,16,16), Vec3u(16, 16, 16));
		PolyGrid<3> polys(surface);

		// Initialise a seed.
		surface.seed(Vec3i(0,0,0));

		surface.update([](const auto& pos, const auto& grid) {
			return -1.0f;
		});

		// ==== Action ====

		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

		surface.update([](const auto& pos, const auto& grid) {
			return -0.01f;
		});

		UINT num_spxs_before = 0;
		UINT num_vtxs_before = 0;
		for (const Vec3i& pos_child : polys)
		{
			num_spxs_before += polys.get(pos_child).spx().size();
			num_vtxs_before += polys.get(pos_child).vtx().size();
		}

		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

		BOOST_CHECK_EQUAL(surface.isogrid().children().data().size(), 1u);

		UINT num_spxs_after = 0;
		UINT num_vtxs_after = 0;
		for (const Vec3i& pos_child : polys)
		{
			num_vtxs_after += polys.get(pos_child).vtx().size();
			num_spxs_after += polys.get(pos_child).spx().size();
		}

		BOOST_CHECK_EQUAL(num_spxs_after, num_spxs_before);
		BOOST_CHECK_EQUAL(num_vtxs_after, num_vtxs_before);
	}

	BOOST_AUTO_TEST_CASE(poly_cubes_expand_contract_single_partition)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a surface with a single partition.
		Surface_t surface(Vec3u(16,16,16), Vec3u(16, 16, 16));
		PolyGrid<3> polys(surface);

		// Initialise a seed.
		surface.seed(Vec3i(0,0,0));
		
		surface.update([](const auto& pos, const auto& grid) {
			return -1.0f;
		});

		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		UINT num_spxs_before = 0;
		UINT num_vtxs_before = 0;
		for (const Vec3i& pos_child : polys)
		{
			num_spxs_before += polys.get(pos_child).spx().size();
			num_vtxs_before += polys.get(pos_child).vtx().size();
		}

		// ==== Action ====

		// Expand.
		surface.update_start();
		surface.disogrid(Vec3i(-1,0,0), -1.0f);
		surface.update_end();
		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// Contract
		surface.update_start();
		surface.disogrid(Vec3i(-2,0,0), 1.0f);
		surface.disogrid(Vec3i(-1,-1,0), 1.0f);
		surface.disogrid(Vec3i(-1,1,0), 1.0f);
		surface.disogrid(Vec3i(-1,0,-1), 1.0f);
		surface.disogrid(Vec3i(-1,0,1), 1.0f);
		surface.update_end();
		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

		BOOST_CHECK_EQUAL(surface.isogrid().children().data().size(), 1);

		UINT num_spxs_after = 0;
		UINT num_vtxs_after = 0;
		for (const Vec3i& pos_child : polys)
		{
			num_vtxs_after += polys.get(pos_child).vtx().size();
			num_spxs_after += polys.get(pos_child).spx().size();
		}

		BOOST_CHECK_EQUAL(num_spxs_after, num_spxs_before);
		BOOST_CHECK_EQUAL(num_vtxs_after, num_vtxs_before);
	}


	BOOST_AUTO_TEST_CASE(poly_cubes_expand_contract_across_partition)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a 16x16x16 surface with two 16x8x16 partitions.
		Surface_t surface(Vec3u(16,16,16), Vec3u(16, 8, 16));
		PolyGrid<3> polys(surface);
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
		Grid<UINT, 3> grid_spxs_before(polys.size(), polys.offset());

		// Initialise a seed.
		surface.seed(Vec3i(0,-4,0));
		surface.seed(Vec3i(0,2,0));

		surface.update([](const auto& pos, const auto& grid) {
			return -1.0f;
		});

		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		UINT num_spxs_before = 0;
		for (const Vec3i& pos_child : polys)
		{
			grid_spxs_before(pos_child) = polys.get(pos_child).spx().size();
			num_spxs_before += polys.get(pos_child).spx().size();
		}

		// ==== Action ====

		// Expand - expanding across to other partition.
		surface.update_start();
		surface.disogrid(Vec3i(0,1,0), -1.0f);
		surface.update_end();
		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// Contract.
		surface.update_start();
		surface.disogrid(Vec3i(0,0,0), 1.0f);
		surface.disogrid(Vec3i(-1,1,0), 1.0f);
		surface.disogrid(Vec3i(1,1,0), 1.0f);
		surface.disogrid(Vec3i(0,1,-1), 1.0f);
		surface.disogrid(Vec3i(0,1,1), 1.0f);
		surface.update_end();
		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

		poly.surf(surface);

		assert_partitioned_matches_baseline(polys, poly);
		for (const Vec3i& pos_child : polys)
		{
			BOOST_CHECK_MESSAGE(
				grid_spxs_before(pos_child) == polys.get(pos_child).spx().size(),
				"Partition " + stringifyVector(pos_child) + " "
				+ std::to_string(polys.get(pos_child).spx().size())
				+ " spxs now == "
				+ std::to_string(grid_spxs_before(pos_child)) + " spxs before"
			);
		}

		UINT num_spxs_after = 0;
		for (const Vec3i& pos_child : polys)
			num_spxs_after += polys.get(pos_child).spx().size();

		BOOST_CHECK_EQUAL(num_spxs_after, num_spxs_before);
	}


	BOOST_AUTO_TEST_CASE(poly_cubes_in_partition_hosting_neighbours_polys)
	{
		// ==== Setup ====

		typedef PolyGrid<3> PolyGrid_t;
		typedef Surface<3, 3> Surface_t;
		// Initialise a 24x24x24 surface with two 24x8x24 partitions.
		Surface_t surface(Vec3u(24,24,24), Vec3u(24, 12, 24));
		PolyGrid<3> polys(surface);
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
		Grid<UINT, 3> grid_spxs_before(polys.size(), polys.offset());

		// Initialise a seed.
		surface.seed(Vec3i(0,-5,0));
		surface.seed(Vec3i(0,4,0));

		surface.update([](const auto& pos, const auto& grid) {
			return -1.0f;
		});
		polys.notify(surface);

		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// Expand - expanding across to other partition.
		surface.update_start();
		surface.disogrid(Vec3i(0,3,0), -1.0f);
		surface.update_end();
		polys.notify(surface);
		polys.poly_cubes(surface);
		polys.update_end();

		poly.surf(surface);
		assert_partitioned_matches_baseline(polys, poly);
		poly.reset();

		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		UINT num_spxs_before = 0;
		for (const Vec3i& pos_child : polys)
		{
			grid_spxs_before(pos_child) = polys.get(pos_child).spx().size();
			num_spxs_before += polys.get(pos_child).spx().size();
		}

		// ==== Action ====

		surface.update_start();
		surface.disogrid(Vec3i(0,5,0), 0.0f);
		surface.update_end();

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

		poly.surf(surface);
		assert_partitioned_matches_baseline(polys, poly);
		poly.reset();

		for (const Vec3i& pos_child : polys)
		{
			BOOST_CHECK_MESSAGE(
				grid_spxs_before(pos_child) == polys.get(pos_child).spx().size(),
				"Partition " + stringifyVector(pos_child) + " "
				+ std::to_string(polys.get(pos_child).spx().size())
				+ " spxs now == "
				+ std::to_string(grid_spxs_before(pos_child)) + " spxs before"
			);
		}

		UINT num_spxs_after = 0;
		for (const Vec3i& pos_child : polys)
			num_spxs_after += polys.get(pos_child).spx().size();

		BOOST_CHECK_EQUAL(num_spxs_after, num_spxs_before);
	}


	/**
	 * Test polygonisation of entire grid.
	 *
	 */
	BOOST_AUTO_TEST_CASE(poly_all)
	{
		// ==== Setup ====
		Surface<3, 3> surface(Vec3u(13,13,13), Vec3u(4,4,4));
		PolyGrid<3> polys(surface);
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());

		surface.seed(Vec3i(0,0,0));

		// ==== Action ====

		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.disogrid(pos, -1.0);
		surface.update_end();


		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.disogrid(pos, -1.0);
		surface.update_end();

		polys.surf(surface);

//		BOOST_TEST_MESSAGE(stringifyGridSlice(surface.isogrid()));

		// ==== Confirm ====

/*
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    2 |    1 |    0 |   -1 |   -2 |   -1 |    0 |    1 |    2 |    3 |    3 |
|    3 |    3 |    3 |    2 |    1 |    0 |   -1 |    0 |    1 |    2 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    2 |    1 |    0 |    1 |    2 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    2 |    1 |    2 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    2 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/

		poly.surf(surface);
		assert_partitioned_matches_baseline(polys, poly);
	}

BOOST_AUTO_TEST_SUITE_END()
/** @} */ // End group Tests.

/**
 *  @class felt::PolyGrid
 *  @test see @ref Tests
 */
