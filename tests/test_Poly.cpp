#include <iostream>
#include <fstream>
#include <boost/test/unit_test.hpp>

#define _TESTING

#include "Surface.hpp"
#include "Poly.hpp"

using namespace felt;
/*
 * Test the Poly library.
 */
BOOST_AUTO_TEST_SUITE(test_Poly)

	std::string stringifyBitmask(long mask, short length = 8)
	{
		std::string str;
		for (unsigned bitIdx = 0; bitIdx < length; bitIdx++)
			str += std::to_string(1 & (mask >> (length-1-bitIdx)));
		return str;
	}


	/*
	 * Initialsing.
	 */
	BOOST_AUTO_TEST_CASE(init)
	{
		Poly<2> poly2D(Vec2u(9,9), Vec2i(-4,-4));


		Poly<3> poly3D(Vec3u(9,9,9), Vec3i(-4,-4,-4));

		Poly<2>::Vertex vertex2D;
		vertex2D.pos(0) = 1;
		vertex2D.pos(1) = 1;

		Poly<3>::Vertex vertex3D;
		vertex3D.pos(0) = 1;
		vertex3D.pos(1) = 1;
		vertex3D.pos(2) = 1;
		vertex3D.norm(0) = 1;
		vertex3D.norm(1) = 1;
		vertex3D.norm(2) = 1;

		Poly<3>::Simplex triangle;

		Grid<Vec2u, 2>& gridVert2D = poly2D.idx();
		Grid<Vec3u, 3>& gridVert3D = poly3D.idx();

		// Check offset parameter has been applied to the underlying grid.
		BOOST_CHECK_EQUAL(gridVert2D.offset(), Vec2i(-4,-4));
		// Check the grid has been initialised with "null" values.
		BOOST_CHECK_EQUAL(gridVert2D(Vec2i(1,1)), Vec2u::Constant(std::numeric_limits<UINT>::max()));
		BOOST_CHECK_EQUAL(gridVert3D(Vec3i(1,1,1)), Vec3u::Constant(std::numeric_limits<UINT>::max()));


		// Add dummy vertex and simplex to arrays.
		poly3D.vtx().push_back(vertex3D);
		poly3D.spx().push_back(triangle);
		// Set edge vertex indices at (0,0,0) to dummy values.
		gridVert3D(Vec3i(0,0,0)) = Vec3u(1,2,3);
		// Ensure vertex and simplex were added.
		BOOST_CHECK_EQUAL(poly3D.vtx().size(), 1);
		BOOST_CHECK_EQUAL(poly3D.spx().size(), 1);
		// Ensure index lookup grid is updated.
		BOOST_CHECK_EQUAL((gridVert3D(0,0,0) - Vec3u(1,2,3)).sum(), 0);
		// Reset the polygonisation.
		poly3D.reset();
		// Ensure vertices and simplices are destroyed.
		BOOST_CHECK_EQUAL(poly3D.vtx().size(), 0);
		BOOST_CHECK_EQUAL(poly3D.spx().size(), 0);
		// Ensure grid is now back to null.
		BOOST_CHECK_EQUAL(gridVert3D(Vec3i(0,0,0)), Vec3u::Constant(std::numeric_limits<UINT>::max()));

	}

	/**
	 * Test calculation of vertices to eventually be joined to make triangles.
	 */
	BOOST_AUTO_TEST_CASE(lerp)
	{
		Surface<2> surface2D(Vec2u(7,7));
		Surface<3> surface3D(Vec3u(7,7,7));

		Poly<2>::Vertex vertex2D;
		Poly<3>::Vertex vertex3D;
		Poly<2> poly2D(surface2D.dims(), surface2D.offset());
		Poly<3> poly3D(surface3D.dims(), surface3D.offset());

		// Text extremities of grid, ensure no segmentation fault errors.
		poly2D.idx(surface2D.phi(), surface2D.pos_min(), 0);
		poly2D.idx(surface2D.phi(), surface2D.pos_max(), 0);
		poly2D.idx(surface2D.phi(), surface2D.pos_min(), 1);
		poly2D.idx(surface2D.phi(), surface2D.pos_max(), 1);

		poly3D.idx(surface3D.phi(), surface3D.pos_min(), 0);
		poly3D.idx(surface3D.phi(), surface3D.pos_max(), 0);
		poly3D.idx(surface3D.phi(), surface3D.pos_min(), 1);
		poly3D.idx(surface3D.phi(), surface3D.pos_max(), 1);
		poly3D.idx(surface3D.phi(), surface3D.pos_min(), 2);
		poly3D.idx(surface3D.phi(), surface3D.pos_max(), 2);

		// Reset vertex cache.
		poly2D.reset();
		poly3D.reset();

		// Create seed and expand outwards.
		// NOTE: will immediately hit edge of grid where max val is 0.5, so centre will be -0.5
		// and each neighbour will be +0.5.
		surface2D.seed(Vec2i(0,0));
		surface3D.seed(Vec3i(0,0,0));
		surface2D.update_start();
		surface2D.dphi(Vec2i(0,0), -1);
		surface2D.update_end();
		surface3D.update_start();
		surface3D.dphi(Vec3i(0,0,0), -1);
		surface3D.update_end();

//		std::cerr << surface3D.phi().data();


		// Index in vertex array of vertex along edge from centre to +x.
		UINT idx2D = poly2D.idx(surface2D.phi(), Vec2i(0,0), 0);
		// Index in vertex array of vertex along edge from centre to +z.
		UINT idx3D = poly3D.idx(surface3D.phi(), Vec3i(0,0,0), 2);
		// Vertex along these edges should be the first in the list.
		BOOST_CHECK_EQUAL(idx2D, 0);
		BOOST_CHECK_EQUAL(idx3D, 0);

		// Get the vertex at this index.
		vertex2D = poly2D.vtx(idx2D);
		vertex3D = poly3D.vtx(idx3D);
		// Ensure vertex is positioned correctly.
		BOOST_CHECK_SMALL((vertex2D.pos - Vec2f(0.5,0)).sum(), 0.00001f);
		BOOST_CHECK_SMALL((vertex3D.pos - Vec3f(0,0,0.5)).sum(), 0.00001f);
		// Ensure vertex normal is in correct direction (3D only).
		BOOST_CHECK_SMALL((vertex3D.norm - Vec3f(0,0,1)).sum(), 0.00001f);


		// Test cache is used for subsequent fetches:

		// First calculate another vertex.
		idx3D = poly3D.idx(surface3D.phi(), Vec3i(0,0,-1), 2);
		vertex3D = poly3D.vtx(idx3D);
		// This new vertex should be appended to array (index=1).
		BOOST_CHECK_EQUAL(idx3D, 1);
		// Check vertex position and normal is correct.
		BOOST_CHECK_SMALL((vertex3D.pos - Vec3f(0,0,-0.5)).sum(), 0.00001f);
		BOOST_CHECK_SMALL((vertex3D.norm - Vec3f(0,0,-1)).sum(), 0.00001f);

		// Now cache should be used for previous vertex, such that idx == 0, not 2.
		idx3D = poly3D.idx(surface3D.phi(), Vec3i(0,0,0), 2);
		vertex3D = poly3D.vtx(idx3D);
		BOOST_CHECK_EQUAL(idx3D, 0);
		// Check it's still at the correct position with the correct normal.
		BOOST_CHECK_SMALL((vertex3D.pos - Vec3f(0,0,0.5)).sum(), 0.00001f);
		BOOST_CHECK_SMALL((vertex3D.norm - Vec3f(0,0,1)).sum(), 0.00001f);
	}


	/**
	 * Test the cube corner inside/outside status bitmask.
	 */
	BOOST_AUTO_TEST_CASE(mask)
	{
		// 2D.
		{
			// Initialise a 2D grid for testing.
			Surface<2> surface(Vec2u(9,9));
			Poly<2> poly(surface.dims(), surface.offset());
			surface.phi().data() <<
				 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
				 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
				 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
				 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
				 2,	 1,	 0,	-1,	-2,	-1,	 0,  1,	 2,
				 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
				 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
				 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
				 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3;

			unsigned short mask;
			mask = Poly<2>::mask(surface.phi(), Vec2i(-3,-3));
			// All outside = 1111.
			BOOST_CHECK_EQUAL(mask, 15);

			mask = Poly<2>::mask(surface.phi(), Vec2i(0,0));
			// All inside = 0000
			BOOST_CHECK_EQUAL(mask, 0);

			mask = Poly<2>::mask(surface.phi(), Vec2i(-1,-1));
			// 0000
			BOOST_CHECK_EQUAL(mask, 0);

			mask = Poly<2>::mask(surface.phi(), Vec2i(1,-1));
			// 0010
			BOOST_CHECK_EQUAL(mask, 2);

			mask = Poly<2>::mask(surface.phi(), Vec2i(2,1));
			// 1111
			BOOST_CHECK_EQUAL(mask, 15);

			mask = Poly<2>::mask(surface.phi(), Vec2i(-2,0));
			// 1000
			BOOST_CHECK_EQUAL(mask, 8);

			mask = Poly<2>::mask(surface.phi(), Vec2i(-1,-2));
			// 0001
			BOOST_CHECK_EQUAL(mask, 1);
		}

		// 3D.
		{
			// Initialise a surface.
			Surface<3> surface(Vec3u(13,13,13));
			Poly<3> poly(surface.dims(), surface.offset());
			unsigned short mask;
			// At time of init, all points are "outside" the surface (there is no surface).
			mask = poly.mask(surface.phi(), Vec3i(0,0,0));
			// All outside = 11111111.
			BOOST_CHECK_EQUAL(mask, 255);

			// Initialise a seed and expand it.
			surface.seed(Vec3i(0,0,0));
			surface.update_start();
			surface.dphi(Vec3i(0,0,0), -1);
			surface.update_end();

			// Relative position of corners in bitmask order (LSB first, MSB last):
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


			mask = poly.mask(surface.phi(), Vec3i(0,0,0));

			// The mask of cube starting at (0,0,0)
			BOOST_CHECK_EQUAL(mask, 0b11100100);

			// Expand the surface outwards twice.
			surface.update_start();;
			for (auto pos : surface)
				surface.dphi(pos, -1);
			surface.update_end();
			surface.update_start();;
			for (auto pos : surface)
				surface.dphi(pos, -1);
			surface.update_end();

			// The central cube is now completely inside the surface.
			mask = Poly<3>::mask(surface.phi(), Vec3i(0,0,0));

//			for (unsigned bitIdx =0; bitIdx < 8; bitIdx++)
//				std::cerr << (1 & (mask >> (7-bitIdx)));
//			std::cerr << std::endl;

			// All inside.
			BOOST_CHECK_EQUAL(mask, 0);
		}
	}

	/**
	 * Test corner bitmask translates to edge mask and vertex order lookup.
	 * Calculate vertices from edge mask and join them to make CCW ordered
	 * simplices using vertex ordering lookup.
	 */
	BOOST_AUTO_TEST_CASE(edge_vertices)
	{
		// 2D.
		{
			Surface<2> surface(Vec2u(9,9));
			Poly<2> poly(surface.dims(), surface.offset());
			surface.phi().data() <<
				 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3,
				 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
				 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
				 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
				 2,	 1,	 0,	-1,	-2,	-1,	 0,  1,	 2,
				 3,	 2,	 1,	 0,	-1,	 0,	 1,	 2,	 3,
				 3,	 3,	 2,	 1,	 0,	 1,	 2,	 3,	 3,
				 3,	 3,	 3,	 2,	 1,	 2,	 3,	 3,	 3,
				 3,	 3,	 3,	 3,	 2,	 3,	 3,	 3,	 3;

			unsigned short mask, vtx_mask;

		
			// 0010
			mask = Poly<2>::mask(surface.phi(), Vec2i(1,-1));
			// 0, -1
			// 1,  0

			vtx_mask = Poly<2>::vtx_mask[mask];
			BOOST_CHECK_EQUAL(vtx_mask, 0b0011);

			// Map of edge index to axis in {0,1} and offset in {(0,0), (1,0), (0,1)}.
			BOOST_CHECK_EQUAL(Poly<2>::edges[0].axis, 0);
			BOOST_CHECK_EQUAL(Poly<2>::edges[0].offset, Vec2i(0,0));
			BOOST_CHECK_EQUAL(Poly<2>::edges[1].axis, 1);
			BOOST_CHECK_EQUAL(Poly<2>::edges[1].offset, Vec2i(1,0));

			// CCW ordering of edge vertices.
			const short* vtx_order = Poly<2>::vtx_order[mask];
			BOOST_CHECK_EQUAL(vtx_order[0], 0);
			BOOST_CHECK_EQUAL(vtx_order[1], 1);
			BOOST_CHECK_EQUAL(vtx_order[2], -1);
			BOOST_CHECK_EQUAL(vtx_order[3], -1);

			// Simplex (line) at given position.
			std::vector<Poly<2>::Simplex>& spxs = poly.spx();
			poly.spx(surface.phi(), Vec2i(1,-1));
			// Check only one simplex.
			BOOST_REQUIRE_EQUAL(spxs.size(), 1);

			// Check ordering of indexes into vertices making up the simplex.
			BOOST_CHECK_EQUAL((UINT)spxs[0].idxs(0), 0);
			BOOST_CHECK_EQUAL((UINT)spxs[0].idxs(1), 1);

			// Check position of vertices at the endpoints of the simplex.
			Vec2f vtx1_pos = poly.vtx((UINT)spxs[0].idxs(0)).pos;
			Vec2f vtx2_pos = poly.vtx((UINT)spxs[0].idxs(1)).pos;
			BOOST_CHECK_EQUAL((FLOAT)vtx1_pos(0), 1);
			BOOST_CHECK_EQUAL((FLOAT)vtx1_pos(1), -1);
			BOOST_CHECK_EQUAL((FLOAT)vtx2_pos(0), 2);
			BOOST_CHECK_EQUAL((FLOAT)vtx2_pos(1), 0);

			// Check degenerate case: cube where corner is precisely zero.
			// 0,  1
			// 1,  2
			poly.reset();
			poly.spx(surface.phi(), Vec2i(2,0));
			BOOST_CHECK_EQUAL(spxs.size(), 0);
		}

		// 3D.
		{
			// Initialise a surface.
			Surface<3> surface(Vec3u(13,13,13));
			Poly<3> poly(surface.dims(), surface.offset());
			unsigned short mask, vtx_mask;
			// At time of init, all points are "outside" the surface (there is no surface).
			mask = poly.mask(surface.phi(), Vec3i(0,0,0));
			// All outside = 11111111.
			vtx_mask = Poly<3>::vtx_mask[mask];

			BOOST_CHECK_EQUAL(vtx_mask, 0b0000);

			surface.phi().fill(-1);
			// All inside = 00000000.
			vtx_mask = Poly<3>::vtx_mask[mask];

			// Reset back to 'all outside' status.
			surface.phi().fill(3);

			BOOST_CHECK_EQUAL(vtx_mask, 0b0000);
			// Initialise a seed and expand it.
			surface.seed(Vec3i(0,0,0));
			surface.update_start();
			surface.dphi(Vec3i(0,0,0), -1);
			surface.update_end();


			// Relative position of corners in bitmask order (LSB first, MSB last):
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


			mask = Poly<3>::mask(surface.phi(), Vec3i(0,0,0));
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

			vtx_mask = PolyBase<3>::vtx_mask[mask];
/*
			( 1,  0,  0 ) --- ( 1,  0, -1 ) == e1
			( 1,  0, -1 ) --- ( 0,  0, -1 ) == e2
			( 0,  1,  0 ) --- ( 1,  1,  0 ) == e4
			( 0,  1,  0 ) --- ( 0,  1, -1 ) == e7
			( 1,  0,  0 ) --- ( 1,  1,  0 ) == e9
			( 0,  0, -1 ) --- ( 0,  1, -1 ) == e11
*/

			BOOST_TEST_MESSAGE(std::to_string(mask) + " = " + stringifyBitmask(mask, 8) + " => " + stringifyBitmask(vtx_mask, 12));
			BOOST_CHECK_EQUAL(vtx_mask, 0b101010010110);

			// Map of edge index to axis and offset.
			BOOST_CHECK_EQUAL(Poly<3>::edges[1].axis, 2);
			BOOST_CHECK_EQUAL(Poly<3>::edges[1].offset, Vec3i(1,0,-1));
			BOOST_CHECK_EQUAL(Poly<3>::edges[7].axis, 2);
			BOOST_CHECK_EQUAL(Poly<3>::edges[7].offset, Vec3i(0,1,-1));
			BOOST_CHECK_EQUAL(Poly<3>::edges[9].axis, 1);
			BOOST_CHECK_EQUAL(Poly<3>::edges[9].offset, Vec3i(1,0,0));

			// CCW ordering of edge vertices.
			const short* vtx_order = Poly<3>::vtx_order[mask];
			// Triangle 1.
			BOOST_CHECK_EQUAL(vtx_order[0], 4);
			BOOST_CHECK_EQUAL(vtx_order[1], 11);
			BOOST_CHECK_EQUAL(vtx_order[2], 7);
			// Triangle 2.
			BOOST_CHECK_EQUAL(vtx_order[3], 9);
			BOOST_CHECK_EQUAL(vtx_order[4], 11);
			BOOST_CHECK_EQUAL(vtx_order[5], 4);
			// Triangle 3.
			BOOST_CHECK_EQUAL(vtx_order[6], 9);
			BOOST_CHECK_EQUAL(vtx_order[7], 2);
			BOOST_CHECK_EQUAL(vtx_order[8], 11);
			// Triangle 4.
			BOOST_CHECK_EQUAL(vtx_order[9], 9);
			BOOST_CHECK_EQUAL(vtx_order[10], 1);
			BOOST_CHECK_EQUAL(vtx_order[11], 2);
			// No triangle.
			BOOST_CHECK_EQUAL(vtx_order[12], -1);
			BOOST_CHECK_EQUAL(vtx_order[13], -1);
			BOOST_CHECK_EQUAL(vtx_order[14], -1);
			BOOST_CHECK_EQUAL(vtx_order[15], -1);

			// Check that edge bitmask matches vertex order array.
			for (UINT idx = 0; idx < 16; idx++)
				if (vtx_order[idx] >= 0)
					BOOST_CHECK_MESSAGE((vtx_mask >> vtx_order[idx]) & 1, stringifyBitmask(vtx_mask, 12) + " >> " + std::to_string(vtx_order[idx]));

			// Simplices (triangles).
			std::vector<Poly<3>::Simplex>& spxs = poly.spx();
			std::vector<Poly<3>::Vertex>& vtxs = poly.vtx();
			// Attempt to generate triangle mesh for cube at (0,0,0).
			poly.spx(surface.phi(), Vec3i(0,0,0));

			// Currently, we have a degenerate case -- corners that are at precisely zero
			// (i.e. points or lines rather than triangles), so no simplices are created.
			// Check 0 triangles are created, but still 6 vertices.
			BOOST_CHECK_EQUAL(vtxs.size(), 6);
			BOOST_CHECK_EQUAL(spxs.size(), 0);

			// Expand the surface a bit, but not enough to change the edges that cross the zero curve.
			// This will mean that interpolation gives a vertex along the cube edge, rather than precisely
			// at the corner, so no degenerate triangles.
			surface.update_start();;
			for (auto pos : surface)
				surface.dphi(pos, -0.3);
			surface.update_end();

			// Check that the corner inside/outside status mask is indeed still the same.
			BOOST_REQUIRE_EQUAL(Poly<3>::mask(surface.phi(), Vec3i(0,0,0)), mask);

			// Reset the polygonisation.
			poly.reset();

			// Recalculate the polygonisation (triangle mesh) for the updated phi grid.
			poly.spx(surface.phi(), Vec3i(0,0,0));

			// Check 4 triangles are now created from 6 vertices.
			BOOST_CHECK_EQUAL(vtxs.size(), 6);
			BOOST_CHECK_EQUAL(spxs.size(), 4);

		}
	}


	BOOST_AUTO_TEST_CASE(local_reset)
	{
		// Initialise a surface.
		Surface<3> surface(Vec3u(13,13,13));
		Poly<3> poly(surface.dims(), surface.offset());
		Poly<3> poly_null(surface.dims(), surface.offset());

		// Initialise a seed and expand it.
		surface.seed(Vec3i(0,0,0));
		surface.update_start();
		surface.dphi(Vec3i(0,0,0), -1);
		surface.update_end();
		surface.update_start();;
		for (auto pos : surface)
			surface.dphi(pos, -0.3);
		surface.update_end();

		const Grid<FLOAT, 3>& phi = surface.phi();
		const Vec3u& dims = phi.dims();
		const Vec3i& offset = phi.offset();
		std::stringstream strGrid;
		
		INT z = 0;
		for (INT x = offset(0); x < (INT)dims(0) + offset(0); x++)
		{
			strGrid << std::endl << "|";
			for (INT y = offset(1); y < (INT)dims(1) + offset(1); y++)
				strGrid << std::setw(5) << (FLOAT)phi(x,y,z) << " |";
		}
		strGrid << std::endl;
		std::cerr << strGrid.str();

		// Simplices (triangles).
		std::vector<Poly<3>::Simplex>& spxs = poly.spx();
		std::vector<Poly<3>::Vertex>& vtxs = poly.vtx();
		// Attempt to generate triangle mesh for cube at (0,0,0).
		poly.spx(surface.phi(), Vec3i(0,0,0));

		BOOST_CHECK_EQUAL(spxs.size(), 4);
		BOOST_CHECK_EQUAL(vtxs.size(), 6);

		// Reset single point.
		poly.reset(Vec3i(0,0,0));

		// Check there's no vertices or triangles and that the index grid is clear.
		BOOST_CHECK_EQUAL(spxs.size(), 0);
		BOOST_CHECK_EQUAL(vtxs.size(), 0);
		for (unsigned i = 0; i < poly.idx().data().size(); i++)
			BOOST_CHECK((bool)(poly.idx().data()(i) == Poly<3>::null_idxs()));



	}

//	BOOST_AUTO_TEST_CASE(export_svg)
//	{
//		{
//			const UINT GRID_SIZE = 100;
//
//			Surface<2> surface(Vec2u(GRID_SIZE,GRID_SIZE));
//			Poly<2> poly(surface.dims(), surface.offset());
//			surface.seed(Vec2i(0,0));
//			surface.update_start();
//			surface.dphi(Vec2i(0,0), -1.0f);
//			surface.update_end();
//
//			for (UINT i = 0; i < 50; i++)
//			{
//				surface.update_start();
//				{
//					for (auto pos : surface)
//					{
//						surface.dphi(pos, surface.phi().gradE(pos).norm() * 0.5f*(-1.0f + 1.0f*surface.phi().curv(pos)));
//					}
//				}
//				surface.update_end();
//			}
//
//			std::ofstream fileHTML;
//			fileHTML.open("Poly.html", std::ios::out | std::ios::trunc);
//
//
//			fileHTML << "<!DOCTYPE html><html><body><svg height='500' width='500'>" << std::endl;
//			for (INT layerID = -2; layerID <= 2; layerID++)
//				for (auto pos : surface.layer(layerID))
//				{
//					const unsigned short mask = Poly<2>::mask(surface.phi(), pos);
//					std::vector<Poly<2>::Simplex>& spxs = poly.spx();
//					poly.spx(surface.phi(), pos);
//					for (UINT i = 0; i < spxs.size(); i++)
//					{
//						const Poly<2>::Simplex& spx = spxs[i];
//						const UINT idx1 = spx.idxs(0);
//						const UINT idx2 = spx.idxs(1);
//						const Poly<2>::Vertex& vtx1 = poly.vtx(idx1);
//						const Poly<2>::Vertex& vtx2 = poly.vtx(idx2);
//						const Vec2f pos1 = vtx1.pos;
//						const Vec2f pos2 = vtx2.pos;
//
//						const Vec2f svg_vtx1 = pos1 * 500/GRID_SIZE + Vec2f(250,250);
//						const Vec2f svg_vtx2 = pos2 * 500/GRID_SIZE + Vec2f(250,250);
//
//						fileHTML << "<line x1='" << svg_vtx1(0) << "' y1='" << svg_vtx1(1) << "' x2='" << svg_vtx2(0) << "' y2='" << svg_vtx2(1) << "' style='stroke:rgb(0,0,0)' />";
//					};
//				}
//			fileHTML << "</svg></html></body>";
//			fileHTML.close();
//		}
//	}
BOOST_AUTO_TEST_SUITE_END()

