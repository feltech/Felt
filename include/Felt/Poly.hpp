#ifndef POLY_H
#define POLY_H
#include <eigen3/Eigen/Dense>
#include <omp.h>
#include <array>
#include <vector>
#include <eigen3/Eigen/StdVector>

#include "Grid.hpp"
#include "MappedGrid.hpp"
#include "Surface.hpp"


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

		PolyBase<2>()
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

		PolyBase<3>()
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
	public:
		// Create typedefs of Eigen types for this D-dimensional polygonisation.
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;

		/// Signed distance grid type.
		typedef Grid<FLOAT, D>	PhiGrid;

		// Handy typedefs for storying vertices and references to them.

		/// D-dimensional vertex type.
		typedef typename PolyBase<D>::Vertex Vertex;
		/// Vertex tuple type (for spatial lookup grid).
		typedef VecDu VtxTuple;
		/// Vertex array type for primary vertex storage.
		typedef std::vector<Vertex> VtxArray;
		/// Vertex spatial look type.
		typedef TrackedGrid<VtxTuple, D> VtxGrid;
		/// Simplex type (line or triangle).
		typedef typename PolyBase<D>::Simplex Simplex;
		/// Simplex array type for primary simplex (line or triangle) storage.
		typedef std::vector<Simplex>
			SpxArray;

		/**
		 * A lookup of offsets from start position to corners of a cube.
		 */
		static const std::array<VecDi, 1 << D> corners;

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
		static const UINT NULL_IDX;
		/// Null vertex tuple value for flagging no vertices at a grid position.
		static const VtxTuple NULL_VTX_TUPLE;

		static const VecDi SpxGridPosOffset;
	protected:
		typedef LookupGrid<D>	LookupGrid_t;

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
		static unsigned short mask (const PhiGrid& phi, const VecDi& pos)
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
		 * Default constructor.
		 */
		Poly () :
			m_grid_vtx()
		{}


		/**
		 * Construct a new polygonisation enclosing a signed distance grid.
		 *
         * @param phi the signed distance grid to polygonise.
          */
		Poly (const VecDu& dims_, const VecDi& offset_) : m_grid_vtx()
		{
			this->init(dims_, offset_);
		}

		/**
		 * Initialise with given dimensions and offset.
		 *
		 * @param dims_
		 * @param offset_
		 */
		void init (const VecDu& dims_, const VecDi& offset_)
		{
			m_grid_vtx.init(dims_, offset_);
			m_grid_vtx.fill(NULL_VTX_TUPLE);
		}

		/**
		 * Get the vertex array.
		 *
         * @return
         */
		VtxArray& vtx ()
		{
			return m_a_vtx;
		}

		/**
		 * Get the vertex array.
		 *
         * @return
         */
		const VtxArray& vtx() const
		{
			return m_a_vtx;
		}

		/**
		 * Get the vertex at idx.
		 *
         * @return
         */
		const Vertex& vtx(const UINT& idx) const
		{
			return this->vtx()[idx];
		}


		/**
		 * Lookup, or calculate then store, and return the index into the vertex
		 * array of a vertex at the zero-crossing of phi at pos_a along axis.
		 *
         * @return
         */
		UINT idx(
			const VecDi& pos_a, const UINT& axis, const PhiGrid& grid_phi
		)
		{
			// Check lookup to see if vertex has already been calculated.
			const UINT& idx_lookup = m_grid_vtx(pos_a)(axis);
			if (idx_lookup != NULL_IDX) {
				return idx_lookup;
			}

			// Position of opposite endpoint.
			VecDi pos_b(pos_a);
			pos_b(axis) += 1;

			// Arbitrary small value, below which consider the vertex to be
			// precisely at one endpoint.
			const FLOAT val_small = Poly<D>::epsilon();

			// Value of phi at each endpoint of this edge.
			const FLOAT val_a = grid_phi(pos_a);
			const FLOAT val_b = grid_phi(pos_b);

			// The newly created vertex.
			Vertex vtx;

			// Check if lies very close to an endpoint or midpoint,
			// if so then no need (and possibly dangerous) to interpolate.
			if (std::abs(val_a) <= val_small) {
				vtx = Vertex(grid_phi, pos_a);
			} else if (std::abs(val_b) <= val_small) {
				vtx = Vertex(grid_phi, pos_b);
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

				vtx = Vertex(grid_phi, vec_c);
			}

			// Append newly created vertex to the cache and return a reference
			// to it.
			const UINT idx = this->vtx().size();
			this->vtx().push_back(vtx);
			m_grid_vtx(pos_a)(axis) = idx;
			m_grid_vtx.add(pos_a);
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
		 * Get the array of simplices.
		 *
		 * @return
         */
		const SpxArray& spx() const
		{
			return m_a_spx;
		}

		/**
		 * Get the simplex at index idx in the array.
		 *
         * @return
         */
		const Simplex& spx(const UINT& idx) const
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
		void spx(const VecDi& pos, const PhiGrid& grid_phi)
		{
			typedef typename PolyBase<D>::Edge Edge;

			// TODO: this is here for consistency only, since the marching
			// cubes implementation marches in the negative z-axis, but
			// positive x and y axes. Hence an offset is required so that the
			// negative z-axis marching is compensated by shifting the
			// calculation in the +z direction by one grid node.
			// (NOTE: has no effect for 2D).
			const VecDi pos_calc = pos - SpxGridPosOffset;

			// Get corner inside-outside bitmask at this position.
			const unsigned short mask = Poly<D>::mask(grid_phi, pos_calc);
			// Get a reference to the simplex array.
			SpxArray& spxs = this->spx();
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
						pos_calc + edge.offset, edge.axis, grid_phi
					);
				}
			}

			// Check for degenerates. Compare every calculated vertex to every
			// other to ensure they are not located on top of one-another.
			// E.g. corners that lie at precisely zero will have D vertices that
			// all lie on that corner.
			// TODO: can't just throw away all triangles like this - others may
			// be valid - must do a per-simplex degenerate check. Maybe not
			// worth it, though.
//			for (UINT edge_idx1 = 0; edge_idx1 < PolyBase<D>::num_edges - 1;
//				edge_idx1++
//			) {
//				for (UINT edge_idx2 = edge_idx1+1;
//					edge_idx2 < PolyBase<D>::num_edges; edge_idx2++)
//				{
//					// Check both edges are bisected by the zero-curve.
//					if (!(((vtx_mask >> edge_idx1) & 1)
//						&& ((vtx_mask >> edge_idx2) & 1)))
//						continue;
//
//					// Get the position vector component of the vertex
//					// information for both edges.
//					const VecDf& pos1 = this->vtx(vtx_idxs[edge_idx1]).pos;
//					const VecDf& pos2 = this->vtx(vtx_idxs[edge_idx2]).pos;
//					const FLOAT dist = (pos1 - pos2).squaredNorm();
//					// If they are essentially the same vertex,
//					// then there is no simplex for this cube.
//					if (dist <= Poly<D>::epsilon())
//						return;
//				}
//			}

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
			}
		} // End spx()


		/**
		 * Polygonise the zero-layer of a Surface.
		 *
		 * @param surface
		 */
		template <UINT L>
		void surf(const Surface<D, L>& surface)
		{
			LookupGrid_t grid_dupe(m_grid_vtx.dims(), m_grid_vtx.offset());
			for (VecDi pos_centre : surface.layer(0))
				for (UINT idx = 0; idx < corners.size(); idx++)
				{
					const VecDi pos_corner = (
						pos_centre - (corners[idx] - SpxGridPosOffset)
					);
					if (grid_dupe.add(pos_corner))
					{
						this->spx(pos_corner, surface.phi());
					}
				}
		}

		/**
		 * Destroy vertices and fill the lookup grid with nulls.
		 *
         * @return
         */
		void reset()
		{
			// Fill vertex grid with null values.
			m_grid_vtx.reset(NULL_VTX_TUPLE);
			// Clear vertex and simplex arrays.
			this->vtx().resize(0);
			this->spx().resize(0);
		}
	};


	/*
	 * Initialise null values.
	 */
	template <UINT D>
	const UINT Poly<D>::NULL_IDX = std::numeric_limits<UINT>::max();
	template <UINT D>
	const typename Poly<D>::VtxTuple Poly<D>::NULL_VTX_TUPLE =
		Poly<D>::VtxTuple::Constant(Poly<D>::NULL_IDX);


} // End namepsace felt.

#endif // POLYBASE_H
