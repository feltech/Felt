#include "catch.hpp"
#include <boost/test/output_test_stream.hpp>

#define _TESTING

#include "Felt/Poly.hpp"

#include "Utils.hpp"

using namespace felt;
/**
 * @defgroup Tests
 * @defgroup PolygonisationTests Polygonisation Tests
 *
 * Tests for polygonisation of scalar field within narrow band of felt::Surface.
 *
 * @{
 * 	@name Poly
 * 	@ref felt::Poly
 */
SCENARIO("Poly")
{
	/**
	 * Initialsing.
	 */
	WHEN("init")
	{
		Surface<2> surface2D(Vec2u(9,9));
		Surface<3> surface3D(Vec3u(9,9,9));
		// Create a 2D polygonisation in a 9x9 embedding, offset by (-4,-4)
		// so that (0,0) in coordinate space translates to (5,5) in grid space.
		Poly<2> poly2D(surface2D.isogrid().size(), surface2D.isogrid().offset());
		// Similarly, create a 3D polygonisation in a 9x9x9 embedding.
		Poly<3> poly3D(surface3D.isogrid().size(), surface3D.isogrid().offset());

		// Create a 2D vertex, consisting simply of position.
		Poly<2>::Vertex vertex2D;
		vertex2D.pos(0) = 1;
		vertex2D.pos(1) = 1;

		// Create a 3D vertex, consisting of position and normal.
		Poly<3>::Vertex vertex3D;
		vertex3D.pos(0) = 1;
		vertex3D.pos(1) = 1;
		vertex3D.pos(2) = 1;
		vertex3D.norm(0) = 1;
		vertex3D.norm(1) = 1;
		vertex3D.norm(2) = 1;

		// Create an (uninitialised) 3D simplex (i.e. triangle).
		Poly<3>::Simplex triangle;

		CHECK(poly2D.vtx().size() == 0);

		CHECK(poly3D.vtx().size() == 0);

		// Add dummy vertex and simplex to the polygonisation object.
		poly3D.vtx().push_back(std::move(vertex3D));
		poly3D.spx().push_back(triangle);
		CHECK(poly3D.spx().size() == 1);

		// Reset the polygonisation.
		poly3D.reset();
		// Ensure vertices and simplices are destroyed.
		CHECK(poly3D.vtx().size() == 0);
		CHECK(poly3D.spx().size() == 0);
	}

	GIVEN("a polygonisation attached to a 7x7 surface of 0.4 units radius")
	{
		Surface<2> surface(Vec2u(7,7));
		surface.seed(Vec2i(0,0));
		surface.update([](auto& pos, auto& grid) {
			return -0.4f;
		});

		INFO(stringifyGridSlice(surface.isogrid()));

		Poly<2> poly(surface.isogrid().size(), surface.isogrid().offset());

		WHEN("we calculate the vertex position between (0,0) and (1,0)")
		{
			// Index in vertex array of vertex along edge from centre to +x.
			UINT idx = poly.idx(Vec2i(0,0), 0, surface.isogrid());

			THEN("the returned index is 0")
			{
				CHECK(idx == 0);
			}

			AND_WHEN("we retrieve the vertex")
			{
				const Poly<2>::Vertex& vertex = poly.vtx(idx);

				THEN("the vertex position is interpolated away from the centre")
				{
					CHECK(vertex.pos == Vec2f(0.4f, 0));
				}
			}
		}
	}
	
	GIVEN("a polygonisation attached to a 7x7x7 surface of 0.4 units radius")
	{
		Surface<3> surface(Vec3u(7,7,7));
		surface.seed(Vec3i(0,0,0));
		surface.update([](auto& pos, auto& grid) {
			return -0.4f;
		});

		INFO(stringifyGridSlice(surface.isogrid()));

		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
		// Reserve space to prevent reallocation when testing for equality of memory addresses.
		poly.vtx().reserve(2);

		WHEN("we request the vertex between (0,0,0) and (0,0,1)")
		{
			// Index in vertex array of vertex along edge from centre to +z.
			UINT idx = poly.idx(Vec3i(0,0,0), 2, surface.isogrid());

			THEN("the returned index is 0")
			{
				CHECK(idx == 0);
			}

			AND_WHEN("we retrieve the vertex")
			{
				const Poly<3>::Vertex& vertex = poly.vtx(idx);

				THEN("the vertex position is interpolated away from the centre")
				{
					CHECK(vertex.pos == Vec3f(0, 0, 0.4f));
				}
				
				THEN("the vertex normal is in the +z direction")
				{
					CHECK(vertex.norm == Vec3f(0, 0, 1.0f));
				}
				
				AND_WHEN("we retrieve another vertex")
				{
					UINT idx = poly.idx(Vec3i(0,0,-1), 2, surface.isogrid());
					
					THEN("the new vertex's index is 1")
					{
						CHECK(idx == 1);
					}
					
					AND_WHEN("we request the previous vertex again")
					{
						UINT idx = poly.idx(Vec3i(0,0,0), 2, surface.isogrid());
						const Poly<3>::Vertex& vertex_again = poly.vtx(idx);
						
						THEN("the index is unchanged and the vertex is the same object")
						{
							CHECK(idx == 0);
							CHECK(&vertex == &vertex_again);
						}
					}
				}
			}
		}
	}

	/**
	 * Test the cube corner inside/outside status bitmask.
	 */
	WHEN("mask_2D")
	{
		// Initialise a 2D grid for testing.
		Surface<2> surface(Vec2u(9,9), Vec2u(9,9));
		Poly<2> poly(surface.isogrid().size(), surface.isogrid().offset());
		surface.isogrid().add_child(Vec2i(0, 0));
		surface.isogrid().snapshot().data() = {
			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
			 2,	 1,	 0,	-1,	-2,	-1,	 0,  1,	 2,
			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3
		};

		surface.isogrid().flush_snapshot();

		unsigned short mask;
		mask = Poly<2>::mask(surface.isogrid(), Vec2i(-3,-3));
		// All outside = 1111.
		CHECK(mask == 15);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(0,0));
		// All inside = 0000
		CHECK(mask == 0);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(-1,-1));
		// 0000
		CHECK(mask == 0);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(1,-1));
		// 0010
		CHECK(mask == 2);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(2,1));
		// 1111
		CHECK(mask == 15);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(-2,0));
		// 1000
		CHECK(mask == 8);

		mask = Poly<2>::mask(surface.isogrid(), Vec2i(-1,-2));
		// 0001
		CHECK(mask == 1);
	}

	/**
	 * Test the cube corner inside/outside status bitmask.
	 */
	WHEN("mask_3D")
	{

		// 3D.
		{
			// Initialise a surface.
			Surface<3> surface(Vec3u(13,13,13));
			Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
			unsigned short mask;
			// At time of init, all points are "outside" the surface (there is
			// no surface).
			mask = poly.mask(surface.isogrid(), Vec3i(0,0,0));
			// All outside = 11111111.
			CHECK(mask == 255);

			// Initialise a seed and expand it.
			surface.seed(Vec3i(0,0,0));
			surface.update_start();
			surface.delta(Vec3i(0,0,0), -1);
			surface.update_end();

			// Relative position of corners in bitmask order (LSB first,
			// MSB last):
//			(0, 0, 0),
//			(1, 0, 0),
//			(1, 0,-1),
//			(0, 0,-1),
//			(0, 1, 0),
//			(1, 1, 0),
//			(1, 1,-1),
//			(0, 1,-1)

			// Cross section of surface now looks like this:
//			 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
//			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
//			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
//			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
//			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3;


			mask = poly.mask(surface.isogrid(), Vec3i(0,0,0));

			// The mask of cube starting at (0,0,0)
			CHECK(mask == 0b11100100);

			// Expand the surface outwards twice.
			surface.update_start();;
			for (auto pos : surface.layer(0))
				surface.delta(pos, -1);
			surface.update_end();
			surface.update_start();;
			for (auto pos : surface.layer(0))
				surface.delta(pos, -1);
			surface.update_end();

			// The central cube is now completely inside the surface.
			mask = Poly<3>::mask(surface.isogrid(), Vec3i(0,0,0));

//			for (unsigned bitIdx =0; bitIdx < 8; bitIdx++)
//				std::cerr << (1 & (mask >> (7-bitIdx)));
//			std::cerr << std::endl;

			// All inside.
			CHECK(mask == 0);
		}
	}

	/**
	 * Ensure corner bitmask translates to edge mask and vertex order lookup.
	 *
	 * Calculate vertices from edge mask and join them to make CCW ordered
	 * simplices using vertex ordering lookup. 2D.
	 */
	WHEN("edge_vertices_2D")
	{
		Surface<2> surface(Vec2u(9,9), Vec2u(9,9));
		Poly<2> poly(surface.isogrid().size(), surface.isogrid().offset());
		surface.isogrid().add_child(Vec2i(0, 0));
		surface.isogrid().snapshot().data() = {
			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
			 2,	 1,	 0,	-1,	-2,	-1,	 0,  1,	 2,
			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3
		};
		surface.isogrid().flush_snapshot();

		unsigned short mask, vtx_mask;


		// 0010
		mask = Poly<2>::mask(surface.isogrid(), Vec2i(1,-1));
		// 0, -1
		// 1,  0

		vtx_mask = Poly<2>::vtx_mask[mask];
		CHECK(vtx_mask == 0b0011);

		// Map of edge index to axis in {0,1} and offset in
		// {(0,0), (1,0), (0,1)}.
		CHECK(Poly<2>::edges[0].axis == 0);
		CHECK(Poly<2>::edges[0].offset == Vec2i(0,0));
		CHECK(Poly<2>::edges[1].axis == 1);
		CHECK(Poly<2>::edges[1].offset == Vec2i(1,0));

		// CCW ordering of edge vertices.
		const short* vtx_order = Poly<2>::vtx_order[mask];
		CHECK(vtx_order[0] == 0);
		CHECK(vtx_order[1] == 1);
		CHECK(vtx_order[2] == -1);
		CHECK(vtx_order[3] == -1);

		// Simplex (line) at given position.
		Poly<2>::SpxArray& spxs = poly.spx();
		poly.spx(Vec2i(1,-1), surface.isogrid());
		// Check only one simplex.
		REQUIRE(spxs.size() == 1);

		// Check ordering of indexes into vertices making up the simplex.
		CHECK((UINT)spxs[0].idxs(0) == 0);
		CHECK((UINT)spxs[0].idxs(1) == 1);

		// Check position of vertices at the endpoints of the simplex.
		Vec2f vtx1_pos = poly.vtx((UINT)spxs[0].idxs(0)).pos;
		Vec2f vtx2_pos = poly.vtx((UINT)spxs[0].idxs(1)).pos;
		CHECK((FLOAT)vtx1_pos(0) == 1);
		CHECK((FLOAT)vtx1_pos(1) == -1);
		CHECK((FLOAT)vtx2_pos(0) == 2);
		CHECK((FLOAT)vtx2_pos(1) == 0);

		// Check degenerate case: cube where corner is precisely zero.
		// 0,  1
		// 1,  2
		// TODO: doesn't work, see discussion in 3D test below.
//		poly.reset();
//		poly.spx(Vec2i(2,0));
//		CHECK(spxs.size() == 0);
	}


	/**
	 * Test corner bitmask translates to edge mask and vertex order lookup.
	 *
	 * Calculate vertices from edge mask and join them to make CCW ordered
	 * simplices using vertex ordering lookup. 3D.
	 */
	WHEN("edge_vertices_3D")
	{
		unsigned short mask, vtx_mask;

		// Initialise a surface.
		Surface<3> surface(Vec3u(13,13,13));
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
		Poly<3>::SpxArray& spxs = poly.spx();
		Poly<3>::VtxArray& vtxs = poly.vtx();

		// At time of init, all points are "outside" the surface
		// (there is no surface).
		mask = poly.mask(surface.isogrid(), Vec3i(0,0,0));
		// All outside = 11111111.
		vtx_mask = Poly<3>::vtx_mask[mask];

		CHECK(vtx_mask == 0b0000);

		surface.isogrid().fill(-1);
		// All inside = 00000000.
		vtx_mask = Poly<3>::vtx_mask[mask];

		// Reset back to 'all outside' status.
		surface.isogrid().fill(3);

		CHECK(vtx_mask == 0b0000);
		// Initialise a seed and expand it.
		surface.seed(Vec3i(0,0,0));

		// Attempt to generate triangle mesh for cube at (0,0,0).
		poly.spx(Vec3i(0,0,0), surface.isogrid());

		// TODO: Currently, we have a degenerate case -- corners that are at
		// precisely zero (i.e. points or lines rather than triangles),
		// so no simplices should be created.
		// 3x edges of the cube are cut, but interpolation yields all 3
		// cut points come from the same corner, the singularity seed point.
		// Should find a way to strip simplices/vertices of degenerate
		// triangles.
		CHECK(vtxs.size() == 3);
		CHECK(spxs.size() == 1);

		// Expand the surface outward.
		surface.update_start();
		surface.delta(Vec3i(0,0,0), -1);
		surface.update_end();


		// Relative position of corners in bitmask order
		// (LSB first, MSB last):
//			(0, 0, 0),
//			(1, 0, 0),
//			(1, 0,-1),
//			(0, 0,-1),
//			(0, 1, 0),
//			(1, 1, 0),
//			(1, 1,-1),
//			(0, 1,-1)

		// Cross section of surface now looks like this:
//			 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
//			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
//			 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
//			 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
//			 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
//			 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3,	 3;


		mask = Poly<3>::mask(surface.isogrid(), Vec3i(0,0,0));
/*
		== 0b11100100 (see test 'mask').
		(0, 0, 0) == inside
		(1, 0, 0) == inside
		(1, 0,-1) == outside
		(0, 0,-1) == inside
		(0, 1, 0) == inside
		(1, 1, 0) == outside
		(1, 1,-1) == outside
		(0, 1,-1) == outside
*/

		vtx_mask = Poly<3>::vtx_mask[mask];
/*
		( 1,  0,  0 ) --- ( 1,  0, -1 ) == e1
		( 1,  0, -1 ) --- ( 0,  0, -1 ) == e2
		( 0,  1,  0 ) --- ( 1,  1,  0 ) == e4
		( 0,  1,  0 ) --- ( 0,  1, -1 ) == e7
		( 1,  0,  0 ) --- ( 1,  1,  0 ) == e9
		( 0,  0, -1 ) --- ( 0,  1, -1 ) == e11
*/

		INFO(
			std::to_string(mask) +
			" = " + stringifyBitmask(mask, 8) +
			" => " + stringifyBitmask(vtx_mask, 12)
		);
		CHECK(vtx_mask == 0b101010010110);

		// Map of edge index to axis and offset.
		CHECK(Poly<3>::edges[1].axis == 2);
		CHECK(Poly<3>::edges[1].offset == Vec3i(1,0,-1));
		CHECK(Poly<3>::edges[7].axis == 2);
		CHECK(Poly<3>::edges[7].offset == Vec3i(0,1,-1));
		CHECK(Poly<3>::edges[9].axis == 1);
		CHECK(Poly<3>::edges[9].offset == Vec3i(1,0,0));

		// CCW ordering of edge vertices.
		const short* vtx_order = Poly<3>::vtx_order[mask];
		// Triangle 1.
		CHECK(vtx_order[0] == 4);
		CHECK(vtx_order[1] == 11);
		CHECK(vtx_order[2] == 7);
		// Triangle 2.
		CHECK(vtx_order[3] == 9);
		CHECK(vtx_order[4] == 11);
		CHECK(vtx_order[5] == 4);
		// Triangle 3.
		CHECK(vtx_order[6] == 9);
		CHECK(vtx_order[7] == 2);
		CHECK(vtx_order[8] == 11);
		// Triangle 4.
		CHECK(vtx_order[9] == 9);
		CHECK(vtx_order[10] == 1);
		CHECK(vtx_order[11] == 2);
		// No triangle.
		CHECK(vtx_order[12] == -1);
		CHECK(vtx_order[13] == -1);
		CHECK(vtx_order[14] == -1);
		CHECK(vtx_order[15] == -1);

		// Check that edge bitmask matches vertex order array.
		for (UINT idx = 0; idx < 16; idx++)
			if (vtx_order[idx] >= 0)
			{
				INFO(
					stringifyBitmask(vtx_mask, 12) + " >> " +
					std::to_string(vtx_order[idx])
				);
				CHECK(
					(bool) ((vtx_mask >> vtx_order[idx]) & 1)
				);
			}

/*
----+y
|
|
+x

|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |  1.7 |  0.7 | -0.3 | -1.3 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |  1.7 |  0.7 | -0.3 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |  1.7 |  0.7 |  1.7 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |  1.7 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
|    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |    3 |
*/

		// Check that the corner inside/outside status mask is indeed still
		// the same.
		REQUIRE(
			Poly<3>::mask(surface.isogrid(), Vec3i(0,0,0)) == mask
		);

		// Reset the polygonisation.
		poly.reset();

		// Recalculate the polygonisation (triangle mesh) for the updated
		// isogrid grid.
		poly.spx(Vec3i(0,0,0), surface.isogrid());

		// Check 4 triangles are now created from 6 vertices.
		CHECK(vtxs.size() == 6);
		CHECK(spxs.size() == 4);


		// Expand the surface a bit, but not enough to change the edges
		// that cross the zero curve. This will mean that interpolation
		// gives a vertex along the cube  edge, rather than precisely at the
		// corner, so no degenerate triangles.
		surface.update_start();
		for (auto pos : surface.layer(0))
			surface.delta(pos, -0.3);
		surface.update_end();

		INFO(stringifyGridSlice(surface.isogrid()));

		// Check that the corner inside/outside status mask is indeed still
		// the same.
		REQUIRE(
			Poly<3>::mask(surface.isogrid(), Vec3i(0,0,0)) == mask
		);

		// Reset the polygonisation.
		poly.reset();

		// Recalculate the polygonisation (triangle mesh) for the updated
		// isogrid grid.
		poly.spx(Vec3i(0,0,0), surface.isogrid());

		// Check 4 triangles are now created from 6 vertices.
		CHECK(vtxs.size() == 6);
		CHECK(spxs.size() == 4);

	}


	WHEN("poly_whole_surface")
	{
		// Initialise a surface.
		Surface<3> surface(Vec3u(13,13,13));
		Poly<3> poly(surface.isogrid().size(), surface.isogrid().offset());
		// Initialise a seed and expand it.
		surface.seed(Vec3i(0,0,0));
		surface.update_start();
		surface.delta(Vec3i(0,0,0), -1.3f);
		surface.update_end();

		// Polygonise zero-layer.
		poly.surf(surface);

		CHECK(poly.spx().size() == 56);
		CHECK(poly.vtx().size() == 30);
	}
}


/** @} */ // End group Tests.

/**
 *  @class felt::Poly
 *  @test see @ref Tests
 */