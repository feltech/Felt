#ifndef POLY_H
#define POLY_H

#include <eigen3/Eigen/Dense>
#include <omp.h>
#include <array>
#include <vector>
#include <eigen3/Eigen/StdVector>

#include "PolyBase.hpp"
#include "Grid.hpp"
#include "MultiLookupGrid.hpp"
#include "MultiTrackedGrid.hpp"
#include "Surface.hpp"

namespace felt 
{
	/**
	 * General polygonisation class.
	 */
	template <UINT D>
	class Poly : public PolyBase<D, void> {
	public:
		/// Signed distance grid type.
		template <class Derived> using IsoGrid = GridBase<Derived>;

		// Create typedefs of Eigen types for this D-dimensional polygonisation.
		/**
		 * D-dimensional unsigned integer vector.
		 */
		using VecDu = felt::VecDu<D>;
		/**
		 * D-dimensional integer vector.
		 */
		using VecDi = felt::VecDi<D>;
		/**
		 * D-dimensional float vector.
		 */
		using VecDf = felt::VecDf<D>;

		/// D-dimensional vertex type.
		using typename PolyBase<D, void>::Vertex;
		using typename PolyBase<D, void>::Edge;

		/// Vertex tuple type (for spatial lookup grid).
		using VtxTuple = VecDu;
		/// Vertex array type for primary vertex storage.
		using VtxArray = std::vector<Vertex>;
		/// Vertex spatial lookup type.
		using VtxGrid = SingleTrackedGrid<VtxTuple, D>;
		/// Simplex type (line or triangle).
		using typename PolyBase<D, void>::Simplex;
		/// Simplex array type for primary simplex (line or triangle) storage.
		using SpxArray = std::vector<Simplex>;

		/// A lookup of offsets from start position to corners of a cube.
		using PolyBase<D, void>::corners;
		/// A lookup of cube edges defined by offset and axis
		using PolyBase<D, void>::edges;
		/// Defines simplices of cubes as D-tuples of vertex indices.
		using PolyBase<D, void>::vtx_order;		
		/**
		 * A lookup from corner inside/outside status bitmask to cut edge status
		 * bitmask.
		 */
		using PolyBase<D, void>::vtx_mask;		
		/// Offset to transform corners to more standard configuration.
		using PolyBase<D, void>::SpxGridPosOffset;
		
		/// Null index for flagging lack of reference into an array.
		static const UINT NULL_IDX;
		/// Null vertex tuple value for flagging no vertices at a grid position.
		static const VtxTuple NULL_VTX_TUPLE;

	protected:
		/// MultiLookup grid type used for avoiding duplicates.
		using MultiLookupGrid_t = MultiLookupGrid<D>;

		/**
		 * List of interpolated vertices.
		 */
		VtxArray m_a_vtx;
		/**
		 * List of simplexes (i.e. lines for 2D or triangles for 3D).
		 */
		SpxArray m_a_spx;
		/**
		 * Indices along each axis of isogrid grid of interpolated vertex.
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
		 * of corners in isogrid.
		 *
         * @param isogrid
         * @param pos
         * @return
         */
		template <class Derived>
		static unsigned short mask (const IsoGrid<Derived>& isogrid, const VecDi& pos)
		{
			// Num corners == 2^D.  That is, 4 for 2D, 8 for 3D.
			unsigned short mask = 0;
			const UINT num_corners = (1 << D);
			for (UINT idx = 0; idx < num_corners; idx++) {
				const VecDi corner = pos + Poly<D>::corners[idx];
				const FLOAT val = isogrid(corner);
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

		Poly (const Poly<D>& other)
		{
			m_a_spx = other.m_a_spx;
			m_a_vtx = other.m_a_vtx;
			m_grid_vtx = other.m_grid_vtx;
		}

		/**
		 * Construct a new polygonisation enclosing a signed distance grid.
		 *
         * @param isogrid the signed distance grid to polygonise.
          */
		Poly (const VecDu& size_, const VecDi& offset_) : m_grid_vtx()
		{
			this->init(size_, offset_);
		}

		/**
		 * Initialise with given dimensions and offset.
		 *
		 * @param size_
		 * @param offset_
		 */
		void init (const VecDu& size_, const VecDi& offset_)
		{
			m_grid_vtx.init(size_, offset_, NULL_VTX_TUPLE);
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
		const Vertex& vtx(const UINT idx) const
		{
			return this->vtx()[idx];
		}


		/**
		 * MultiLookup, or calculate then store, and return the index into the vertex
		 * array of a vertex at the zero-crossing of isogrid at pos_a along axis.
		 *
         * @return
         */
		template <class Derived>
		UINT idx(
			const VecDi& pos_a, const UINT axis, const IsoGrid<Derived>& grid_isogrid
		)
		{
			// Check lookup to see if vertex has already been calculated.
			const UINT idx_lookup = m_grid_vtx.get(pos_a)(axis);
			if (idx_lookup != NULL_IDX) {
				return idx_lookup;
			}

			// Position of opposite endpoint.
			VecDi pos_b(pos_a);
			pos_b(axis) += 1;

			// Arbitrary small value, below which consider the vertex to be
			// precisely at one endpoint.
			const FLOAT val_small = Poly<D>::epsilon();

			// Value of isogrid at each endpoint of this edge.
			const FLOAT val_a = grid_isogrid(pos_a);
			const FLOAT val_b = grid_isogrid(pos_b);

			// The newly created vertex.
			Vertex vtx;

			// Check if lies very close to an endpoint or midpoint,
			// if so then no need (and possibly dangerous) to interpolate.
			if (std::abs(val_a) <= val_small) {
				vtx = Vertex(grid_isogrid, pos_a);
			} else if (std::abs(val_b) <= val_small) {
				vtx = Vertex(grid_isogrid, pos_b);
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

				vtx = Vertex(grid_isogrid, vec_c);
			}

			// Append newly created vertex to the cache and return a reference
			// to it.
			const UINT idx = this->vtx().size();
			this->vtx().push_back(std::move(vtx));
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
		const Simplex& spx(const UINT idx) const
		{
			return this->spx()[idx];
		}



		/**
		 * Generate simplex(es) for isogrid grid at position pos.
		 *
         * @param isogrid
         * @param pos
         * @param mask
         * @param spxs
         */
		template <class Derived>
		void spx(const VecDi& pos, const IsoGrid<Derived>& grid_isogrid)
		{
			// TODO: this is here for consistency only, since the marching
			// cubes implementation marches in the negative z-axis, but
			// positive x and y axes. Hence an offset is required so that the
			// negative z-axis marching is compensated by shifting the
			// calculation in the +z direction by one grid node.
			// (NOTE: has no effect for 2D).
			const VecDi pos_calc = pos - SpxGridPosOffset;

			// Get corner inside-outside bitmask at this position.
			const unsigned short mask = this->mask(grid_isogrid, pos_calc);
			// Get a reference to the simplex array.
			SpxArray& spxs = this->spx();
			// Array of indices of zero-crossing vertices along each axis from
			// this corner.
			UINT vtx_idxs[this->num_edges];
			// MultiLookup the edges that are crossed from the corner mask.
			unsigned short vtx_mask = this->vtx_mask[mask];
			const short* vtx_order = this->vtx_order[mask];

			// Cube corners are all inside or all outside.
			if (vtx_order[0] == -1)
				return;

			// Loop over each crossed edge in the cube, looking up
			// (or calculating, if unavailable) the vertices at the
			// zero-crossing.
			for (UINT edge_idx = 0; edge_idx < this->num_edges;
				edge_idx++
			) {
				// Check if current edge is crossed by the zero curve.
				if ((vtx_mask >> edge_idx) & 1) {
					const Edge& edge = this->edges[edge_idx];
					// Edges are defined as an axis and an offset.
					// MultiLookup index of vertex along current edge.
					vtx_idxs[edge_idx] = this->idx(
						pos_calc + edge.offset, edge.axis, grid_isogrid
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
				spxs.push_back(std::move(simplex));
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
			MultiLookupGrid_t grid_dupe(m_grid_vtx.size(), m_grid_vtx.offset());
			for (VecDi pos_centre : surface.layer(0))
				for (const VecDi& pos_offset : corners)
				{
					const VecDi pos_corner = (
						pos_centre - (pos_offset - SpxGridPosOffset)
					);
					if (grid_dupe.add(pos_corner))
					{
						this->spx(pos_corner, surface.isogrid());
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
			m_grid_vtx.reset();

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
