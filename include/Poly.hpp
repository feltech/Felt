#ifndef POLYBASE_H
#define POLYBASE_H
#include <eigen3/Eigen/Dense>
#include <omp.h>
#include <vector>
#include <eigen3/Eigen/StdVector>

#include "Grid.hpp"


namespace felt {

	/**
	 * Forward declaration.
	 */
	template <UINT D>
	class PolyBase
	{
	};


	/**
	 * 2D-specific definitions.
	 */
	template <>
	class PolyBase<2> {
	protected:
		typedef Eigen::Matrix<UINT, 2, 1> SpxTuple;
	public:

		/**
		 * A 2D vertex (position only).
		 */
		struct Vertex {

			/**
			 * Position of vertex.
			 */
			Vec2f pos;

			/**
			 * Create a new vertex for a Grid at position pos.
			 * NOTE: Grid is unused for 2D vertex construction.
			 *
			 * @param <unused>
			 * @param pos
			 */
			template <typename PosType>
			Vertex(const Grid<FLOAT, 2>&, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
			}

			/**
			 * Create an uninitialised vertex.
			 */
			Vertex()
			{
			}
		};

		/**
		 * A 2D simplex (a line with 2 endpoints).
		 */
		struct Simplex {
			Vec2u idxs;
		};

		/**
		 * A 2D square edge (an offset from the bottom-left corner and an axis
		 * in {0,1})
		 */
		struct Edge {
			Vec2i offset;
			INT axis;
		};

#ifndef _TESTING
	protected:
#endif
		/**
		 * Number of edges on a square.
		 */
		static const short num_edges = 4;

		/**
		 * A lookup from inside/outside status bitmask to vertex ordering to
		 * create representative simplices (lines).
		 */
		static const short vtx_order [][4];


	protected:

		PolyBase<2>(const Vec2u& dims, const Vec2i& offset)
		{
		}

		~PolyBase<2>()
		{
		}
	};


	/**
	 * 3D-specific definitions.
	 */
	template <>
	class PolyBase<3> {
	protected:
		typedef Eigen::Matrix<UINT, 5, 1> SpxTuple;
	public:

		/**
		 * A 3D vertex (position and normal).
		 */
		struct Vertex {
			/**
			 * Position of vertex.
			 */
			Vec3f pos;
			/**
			 * Normal of vertex.
			 */
			Vec3f norm;

			/**
			 * Create a vertex for grid at position pos, calculating norm
			 * from the gradient in grid at pos.
			 *
			 * @param grid
			 * @param pos
			 */
			template <typename PosType>
			Vertex(const Grid<FLOAT, 3>& grid, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
				this->norm = grid.gradC(pos);
				this->norm.normalize();

			}

			/**
			 * Create an uninitialised vertex.
			 */
			Vertex()
			{
			}
		};

		/**
		 * A 3D simplex (a triangle with 3 endpoints).
		 */
		struct Simplex {
			Vec3u idxs;
		};

		/**
		 * A 3D cube edge (an offset from the left-bottom-forward corner and an
		 * axis in {0,1})
		 */
		struct Edge {
			Vec3i offset;
			INT axis;
		};

		/**
		 * Number of edges on a square.
		 */
		static const short num_edges = 12;

		/**
		 * A lookup from inside/outside status bitmask to vertex ordering to
		 * create representative simplices (triangles).
		 */
		static const short vtx_order [][16];

	protected:

		PolyBase<3>(const Vec3u& dims, const Vec3i& offset)
		{
		}

		~PolyBase<3>()
		{
		}
	};


	/**
	 * General polygonisation class.
	 */
	template <UINT D>
	class Poly : public PolyBase<D> {
	private:
		// Create typedefs of Eigen types for this D-dimensional polygonisation.
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;
	public:
		// Handy typedefs for storying vertices and references to them.

		/// D-dimensional vertex type.
		typedef typename PolyBase<D>::Vertex Vertex;
		/// Vertex tuple type (for spatial lookup grid).
		typedef VecDu VtxTuple;
		/// Vertex array type for primary vertex storage.
		typedef std::vector<Vertex, Eigen::aligned_allocator<Vertex> > VtxArray;
		/// Vertex spatial lookup grid type.
		typedef Grid<VtxTuple, D> VtxGrid;
		/// Simplex type (line or triangle).
		typedef typename PolyBase<D>::Simplex Simplex;
		/// Simplex tuple type (for spatial lookup grid).
		typedef typename PolyBase<D>::SpxTuple SpxTuple;
		/// Simplex spatial lookup grid type.
		typedef Grid<SpxTuple, D> SpxGrid;
		/// Simplex array type for primary simplex (line or triangle) storage.
		typedef std::vector<Simplex, Eigen::aligned_allocator<Simplex> >
			SpxArray;

		/**
		 * A lookup of offsets from start position to corners of a cube.
		 */
		static const VecDi corners [];

		/**
		 * A lookup of cube edges defined by offset and axis
		 */
		static const typename PolyBase<D>::Edge edges [];

		/**
		 * A lookup from corner inside/outside status bitmask to cut edge status
		 * bitmask.
		 */
		static const short vtx_mask [];

		/// Null index for flagging lack of reference into an array.
		static const UINT NullIdx;
		/// Null vertex tuple value for flagging no vertices at a grid position.
		static const VtxTuple NullVtxTuple;
		/// Null simplex tuple value for flagging no simplices at a grid
		/// position.
		static const SpxTuple NullSpxTuple;

		static const VecDi SpxGridPosOffset;
	protected:
		/**
		 * List of interpolated vertices.
		 */
		VtxArray m_a_vtx;
		/**
		 * List of simplexes (i.e. lines for 2D or triangles for 3D).
		 */
		SpxArray m_a_spx;
		/**
		 * Indices along each axis of phi grid of interpolated vertex.
		 */
		VtxGrid m_grid_vtx;
		/**
		 * Spatial lookup of grid pos to simplex (line) indices making up the
		 * polygonisation of this square.
		 */
		SpxGrid m_grid_spx;

	public:
		/**
		 * Defines a value to indicate small enough to consider as zero.
		 *
		 * @return
		 */
		static inline const FLOAT epsilon()
		{
			return std::numeric_limits<FLOAT>::epsilon();
		}

		/**
		 * Calculate corner mask of cube at pos, based on inside-outside status
		 * of corners in phi.
		 *
         * @param phi
         * @param pos
         * @return
         */
		static unsigned short mask(const Grid<FLOAT, D>& phi, const VecDi& pos)
		{
			// Num corners == 2^D.  That is, 4 for 2D, 8 for 3D.
			unsigned short mask = 0;
			const UINT num_corners = (1 << D);
			for (UINT idx = 0; idx < num_corners; idx++) {
				const VecDi corner = pos + Poly<D>::corners[idx];
				const FLOAT val = phi(corner);
				mask |= (val > 0) << idx;
			}
			return mask;
		}

		/**
		 * Construct a new polygonisation enclosing dims grid points with the
		 * origin shifted by offset.
		 *
         * @param dims
         * @param offset
          */
		Poly(const VecDu& dims, const VecDi& offset) :
		PolyBase<D>(dims, offset),
		m_grid_vtx(dims, offset),
		m_grid_spx(dims, offset)
		{
			this->reset();
		}

		/**
		 * Get the vertex array.
		 *
         * @return
         */
		VtxArray& vtx()
		{
			return m_a_vtx;
		}

		/**
		 * Get the vertex at idx.
		 *
         * @return
         */
		const Vertex& vtx(const UINT& idx)
		{
			return this->vtx()[idx];
		}

		/**
		 * Lookup or calculate the vertex at the zero-crossing in phi at pos_a
		 * along axis.
		 *
         * @param phi
         * @param pos_a
         * @param axis
         * @return
         */
		const Vertex& vtx(
			const Grid<FLOAT, D>& phi, const VecDi& pos_a, const UINT& axis
		)
		{
			return this->vtx(this->idx(phi, pos_a, axis));
		}

		/**
		 * Get the vertex lookup grid.
		 *
         * @return
         */
		VtxGrid& grid_vtx()
		{
			return m_grid_vtx;
		}

		SpxGrid& grid_spx()
		{
			return m_grid_spx;
		}


		/**
		 * Lookup, or calculate then store, and return the index into the vertex
		 * array of a vertex at the zero-crossing of phi at pos_a along axis.
		 *
         * @return
         */
		UINT idx(
			const Grid<FLOAT, D>& phi, const VecDi& pos_a, const UINT& axis
		)
		{
			// Check lookup to see if vertex has already been calculated.
			const UINT& idx_lookup = this->grid_vtx()(pos_a) (axis);
			if (idx_lookup != NullIdx) {
				return idx_lookup;
			}

			// Position of opposite endpoint.
			VecDi pos_b(pos_a);
			pos_b(axis) += 1;

			// Arbitrary small value, below which consider the vertex to be
			// precisely at one endpoint.
			const FLOAT val_small = Poly<D>::epsilon();

			// Value of phi at each endpoint of this edge.
			const FLOAT val_a = phi(pos_a);
			const FLOAT val_b = phi(pos_b);

			// The newly created vertex.
			Vertex vtx;

			// Check if lies very close to an endpoint or midpoint,
			// if so then no need (and possibly dangerous) to interpolate.
			if (std::abs(val_a) <= val_small) {
				vtx = Vertex(phi, pos_a);
			} else if (std::abs(val_b) <= val_small) {
				vtx = Vertex(phi, pos_b);
			} else {
				FLOAT mu;

				// If close to midpoint then put at midpoint.
				if (std::abs(val_a - val_b) <= val_small) {
					mu = (FLOAT) 0.5;
				} else
				// Otherwise interpolate between endpoints.
				{
					mu = val_a / (val_a - val_b);
				}

				const VecDf vec_a = pos_a.template cast<FLOAT>();
				const VecDf vec_b = pos_b.template cast<FLOAT>();
				const VecDf vec_c = vec_a + (vec_b - vec_a) * mu;

				vtx = Vertex(phi, vec_c);
			}

			// Append newly created vertex to the cache and return a reference
			// to it.
			const UINT idx = this->vtx().size();
			this->vtx().push_back(vtx);
			this->grid_vtx()(pos_a) (axis) = idx;
			return idx;
		}

		/**
		 * Get the array of simplices.
		 *
		 * @return
         */
		SpxArray& spx()
		{
			return m_a_spx;
		}


		/**
		 * Get the simplex at index idx in the array.
		 *
         * @return
         */
		const Simplex& spx(const UINT& idx)
		{
			return this->spx()[idx];
		}



		/**
		 * Generate simplex(es) for phi grid at position pos.
		 *
         * @param phi
         * @param pos
         * @param mask
         * @param spxs
         */
		void spx(const Grid<FLOAT, D>& phi, const VecDi& pos)
		{
			typedef typename PolyBase<D>::Edge Edge;

			// TODO: this is required for 3D polygonisation, since the marching
			// cubes implementation marches in the negative z-axis, but
			// positive x and y axes.  Each node of the simplex lookup grid
			// expects to be the bottom-back-left -most corner, so that
			// neighbourhood queries are more natural.  Hence an offset is
			// required so that the negative z-axis marching is compensated by
			// shifting the calculation in the +z direction by one grid node.
			const VecDi pos_calc = pos - SpxGridPosOffset;

			// Get corner inside-outside bitmask at this position.
			const unsigned short mask = Poly<D>::mask(phi, pos_calc);
			// Get a reference to the simplex array.
			SpxArray& spxs = this->spx();
			// Get a reference to the spatial simplex lookup grid.
			typename Poly<D>::SpxGrid& grid_spxs = this->grid_spx();

			// Array of indices of zero-crossing vertices along each axis from
			// this corner.
			UINT vtx_idxs[PolyBase<D>::num_edges];
			// Lookup the edges that are crossed from the corner mask.
			unsigned short vtx_mask = Poly<D>::vtx_mask[mask];
			// Loop over each crossed edge in the cube, looking up
			// (or calculating, if unavailable) the vertices at the
			// zero-crossing.
			for (UINT edge_idx = 0; edge_idx < PolyBase<D>::num_edges;
				edge_idx++
			) {
				// Check if current edge is crossed by the zero curve.
				if ((vtx_mask >> edge_idx) & 1) {
					Edge edge = Poly<D>::edges[edge_idx];
					// Edges are defined as an axis and an offset.
					// Lookup index of vertex along current edge.
					vtx_idxs[edge_idx] = this->idx(
						phi, pos_calc+edge.offset, edge.axis
					);
				}
			}

			// Check for degenerates. Compare every calculated vertex to every
			// other to ensure they are not located on top of one-another.
			// E.g. corners that lie at precisely zero will have D vertices that
			// all lie on that corner.
			for (UINT edge_idx1 = 0; edge_idx1 < PolyBase<D>::num_edges - 1;
				edge_idx1++
			) {
				for (UINT edge_idx2 = edge_idx1+1;
					edge_idx2 < PolyBase<D>::num_edges; edge_idx2++)
				{
					// Check both edges are bisected by the zero-curve.
					if (((vtx_mask >> edge_idx1) & 1)
						&& ((vtx_mask >> edge_idx2) & 1))
					{
						// Get the position vector component of the vertex
						// information for both edges.
						const VecDf& pos1 = this->vtx(vtx_idxs[edge_idx1]).pos;
						const VecDf& pos2 = this->vtx(vtx_idxs[edge_idx2]).pos;
						const FLOAT dist = (pos1 - pos2).squaredNorm();
						// If they are essentially the same vertex,
						// then there is no simplex for this cube.
						if (dist <= Poly<D>::epsilon())
							return;
					}
				}
			}

			// Store a count of current index in the array at this position
			// in the spatial lookup simplex grid.
			UINT grid_idx = 0;

			// Join the vertices along each edge that the surface crosses to
			// make a simplex (or simplices).
			// The vtx_order lookup translates corner in-out mask to CCW
			// vertex ordering. We take D elements at a time from the lookup,
			// with each successive subset of D elements forming the next
			// simplex.
			const short* vtx_order = PolyBase<D>::vtx_order[mask];
			for (UINT order_idx = 0; vtx_order[order_idx] != -1;
				order_idx += D
			) {
				Simplex simplex;
				// A simplex for number of dimensions D has D vertices,
				// i.e. D endpoints.
				for (UINT endpoint = 0; endpoint < D; endpoint++)
				{
					// Each vertex of the simplex is stored as an index
					// reference into the 'global' vertex array.
					simplex.idxs(endpoint) = vtx_idxs[
						vtx_order[order_idx+endpoint]
					];
				}
				// Append the simplex to the list of simplices that make up the
				// polygonisation of this grid location.
				spxs.push_back(simplex);
				grid_spxs(pos)[grid_idx] = spxs.size() - 1;
				grid_idx++;
			}
		} // End spx()


		/**
		 * Destroy vertices and fill the lookup grid with nulls.
		 *
         * @return
         */
		void reset()
		{
			// Fill simplex grid with null values.
			this->grid_spx().fill(NullSpxTuple);
			// Fill vertex grid with null values.
			this->grid_vtx().fill(NullVtxTuple);
			// Clear vertex and simplex arrays.
			this->vtx().resize(0);
			this->spx().resize(0);
		}


		void reset(const VecDi& pos)
		{
			VtxGrid& grid_vtx = this->grid_vtx();
			SpxGrid& grid_spx = this->grid_spx();

			const VecDu& dims = grid_spx.dims();

			// Store all 2^d bottom-left cube corners of cubes surrounding node
			// at pos in an array.
			std::vector<VecDi, Eigen::aligned_allocator<VecDi> >
			apos_corners(1 << dims.size());

			// Loop and create all cube corner positions.
			for (
				UINT corner_idx = 0; corner_idx < apos_corners.size();
				corner_idx++
			)
			{
				// 0 = 00 => (x,y)
				// 1 = 01 => (x-1,y)
				// 2 = 10 => (x,y-1)
				// 3 = 11 => (x-1,y-1)

				// Store position of current corner.
				VecDi pos_corner(dims.size());
				// Calculate position of current corner:
				// Loop each dimension, 2D (x,y) or 3D (x,y,z).
				for (INT dim = 0; dim < pos_corner.size(); dim++)
				{
					// Get the element of pos along this axis.
					INT pos_axis = pos(dim);
					// Use the current index in the corner array as a bitmask
					// hash to calculate this axis's offset in {0,1}.
					const INT dir = (corner_idx >> dim) & 1;
					// Offset by 0 or -1.
					pos_axis -= dir;
					// Update the element of the corner's position along this
					// axis to the (potentially) offset value of pos.
					pos_corner(dim) = pos_axis;
				}
				// Set the corner position in the array to the calculated
				// offset position.
				apos_corners[corner_idx] = pos_corner;
			}

			// Loop the cubes containing positions to invalidate.
			for (
				UINT corner_idx = 0; corner_idx < apos_corners.size();
				corner_idx++
			)
			{
				const VecDi& pos_corner = apos_corners[corner_idx];
				grid_spx(pos_corner) = NullSpxTuple;
			}
		}
	};


	/*
	 * Initialise null values.
	 */
	template <UINT D>
	const UINT Poly<D>::NullIdx = std::numeric_limits<UINT>::max();
	template <UINT D>
	const typename Poly<D>::VtxTuple Poly<D>::NullVtxTuple =
		Poly<D>::VtxTuple::Constant(Poly<D>::NullIdx);
	template <UINT D>
	const typename Poly<D>::SpxTuple Poly<D>::NullSpxTuple =
		Poly<D>::SpxTuple::Constant(Poly<D>::NullIdx);


	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	// 2D lookups.
	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////


	/*
		v = vertex
		e = edge
		s = simplex

		(0,0) = v0001
		(1,0) = v0010
		(1,1) = v0100
		(0,1) = v1000

			   e0100
		  v1000-----v0100
			|        |
	  e1000 |	     | e0010
			|	     |
		  v0001-----v0010
			   e0001

		0000

		v0 = inside
		v1 = outside
		____________________
		 v	 |	e	| s (CCW)
		--------------------
		0000 | 0000	|
		0001 | 1001	| 3,0
		0010 | 0011	| 0,1
		0011 | 1010	| 3,1
		0100 | 0110	| 1,2
		0101 | 1111	| 3,0 1,2
		0110 | 0101	| 0,2
		0111 | 1100	| 3,2
		1000 | 1100	| 2,3
		1001 | 0101 | 2,0
		1010 | 1111	| 2,1 0,3
		1011 | 0110	| 2,1
		1100 | 1010	| 3,1
		1101 | 0011	| 1,0
		1110 | 1001	| 0,3
		1111 | 0000	|
	*/

		/**
		 * Relative position of corners in CCW order.
		 */
		template<>
		const Vec2i Poly<2>::corners [] ={
			Vec2i(0, 0),
			Vec2i(1, 0),
			Vec2i(1, 1),
			Vec2i(0, 1)
		};
		/**
		 * Array of edge definitions (offset, direction) matching ::corners.
		 */
		template<>
		const Poly<2>::Edge Poly<2>::edges [] = {
			// (x,y,axis)
			{ Vec2i(0, 0), 0 },
			{ Vec2i(1, 0), 1 },
			{ Vec2i(0, 1), 0 },
			{ Vec2i(0, 0), 1 }
		};
		template<>
		const Vec2i Poly<2>::SpxGridPosOffset(0,0);

		/**
		 * Lookup from corner mask to edge mask.
		 */
		template<>
		const short Poly<2>::vtx_mask [] ={
			0b0000,
			0b1001,
			0b0011,
			0b1010,
			0b0110,
			0b1111,
			0b0101,
			0b1100,
			0b1100,
			0b0101,
			0b1111,
			0b0110,
			0b1010,
			0b0011,
			0b1001,
			0b0000
		};

		/**
		 * Ordering of vertices to build simplex(es).
		 */
		const short PolyBase<2>::vtx_order [][4] ={
			{ -1, -1, -1, -1 },
			{  3,  0, -1, -1 },
			{  0,  1, -1, -1 },
			{  3,  1, -1, -1 },
			{  1,  2, -1, -1 },
			{  3,  0,  1,  2 },
			{  0,  2, -1, -1 },
			{  3,  2, -1, -1 },
			{  2,  3, -1, -1 },
			{  2,  0, -1, -1 }, //
			{  2,  1,  0,  3 },
			{  2,  1, -1, -1 },
			{  3,  1, -1, -1 },
			{  1,  0, -1, -1 },
			{  0,  3, -1, -1 },
			{ -1, -1, -1, -1 }
		};



		///////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////
		// 3D lookups.
		///////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////


		/**
		 * Lookup from corner mask to edge mask.
		 */
		template<>
		const short Poly<3>::vtx_mask [] = {
				0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
				0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
				0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
				0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
				0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
				0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
				0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
				0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
				0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
				0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
				0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
				0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
				0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
				0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
				0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
				0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
				0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
				0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
				0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
				0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
				0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
				0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
				0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
				0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
				0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
				0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
				0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
				0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
				0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
				0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
				0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
				0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

		/**
		 * Given node, march cube toward back, up and right (0,0,0)->(1,1,-1).
		 */
		template<>
		const Vec3i Poly<3>::corners [] ={
			Vec3i( 0,  0,  0),	// c0
			Vec3i( 1,  0,  0),	// c1
			Vec3i( 1,  0, -1),	// c2
			Vec3i( 0,  0, -1),	// c3
			Vec3i( 0,  1,  0),	// c4
			Vec3i( 1,  1,  0),	// c5
			Vec3i( 1,  1, -1),	// c6
			Vec3i( 0,  1, -1)	// c7
		};

		/**
		 * Array of edge definitions (offset, direction) matching ::corners.
		 */
		template<>
		const Poly<3>::Edge Poly<3>::edges [] = {
			// (x,y,z, axis)
			{ Vec3i( 0,  0,  0), 	0 },	// e0
			{ Vec3i( 1,  0, -1), 	2 },	// e1
			{ Vec3i( 0,  0, -1), 	0 },	// e2
			{ Vec3i( 0,  0, -1), 	2 },	// e3
			{ Vec3i( 0,  1,  0), 	0 },	// e4
			{ Vec3i( 1,  1, -1), 	2 },	// e5
			{ Vec3i( 0,  1, -1), 	0 },	// e6
			{ Vec3i( 0,  1, -1), 	2 },	// e7
			{ Vec3i( 0,  0,  0), 	1 },	// e8
			{ Vec3i( 1,  0,  0), 	1 },	// e9
			{ Vec3i( 1,  0, -1), 	1 },	// e10
			{ Vec3i( 0,  0, -1), 	1 }		// e11
		};
		template<>
		const Vec3i Poly<3>::SpxGridPosOffset(0,0,-1);


		/**
		 * Ordering of vertices to build simplex(es).
		 */
		const short PolyBase<3>::vtx_order [][16] = {
				{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
				{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
				{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
				{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
				{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
				{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
				{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
				{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
				{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
				{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
				{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
				{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
				{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
				{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
				{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
				{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
				{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
				{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
				{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
				{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
				{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
				{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
				{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
				{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
				{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
				{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
				{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
				{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
				{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
				{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
				{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
				{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
				{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
				{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
				{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
				{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
				{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
				{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
				{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
				{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
				{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
				{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
				{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
				{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
				{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
				{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
				{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
				{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
				{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
				{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
				{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
				{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
				{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
				{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
				{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
				{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
				{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
				{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
				{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
				{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
				{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
				{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
				{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
				{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
				{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
				{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
				{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
				{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
				{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
				{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
				{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
				{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
				{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
				{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
				{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
				{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
				{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
				{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
				{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
				{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
				{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
				{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
				{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
				{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
				{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
				{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
				{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
				{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
				{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
				{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
				{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
				{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
				{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
				{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
				{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
				{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
				{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
				{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
				{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
				{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
				{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
				{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
				{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
				{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
				{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
				{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
				{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
				{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
				{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
				{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
				{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
				{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
				{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
				{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
				{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
				{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
				{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
				{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
				{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
				{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
				{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
				{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
				{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
				{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
				{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
				{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
				{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
				{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
				{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
				{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
				{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
				{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
				{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
				{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
				{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
				{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
				{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
				{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
				{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
				{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
				{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
				{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
				{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
				{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
				{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
				{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
				{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
				{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
				{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
				{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
				{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
				{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
				{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
				{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
				{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
				{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
				{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
				{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
				{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
				{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
				{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
				{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
				{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
				{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
				{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},	// [228]
				{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
				{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
				{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
				{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
				{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
				{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
				{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
				{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
				{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
				{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
				{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
				{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
				{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
				{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
				{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
				{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
		};

} // End namepsace felt.






#endif // POLYBASE_H
