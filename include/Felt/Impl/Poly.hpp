#ifndef INCLUDE_FELT_IMPL_POLY_HPP_
#define INCLUDE_FELT_IMPL_POLY_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/PolyMixin.hpp>
#include <Felt/Impl/Mixin/TrackedMixin.hpp>
#include <Felt/Surface.hpp>

namespace Felt
{
namespace Impl
{
namespace Poly
{

template <class TIsoGrid>
class Single :
	FELT_MIXINS(
		(Single<TIsoGrid>),
		(Grid::Access::ByRef)(Grid::Data)(Poly::Geom)(Tracked::Activate)(Tracked::SingleList::Reset)
		(Tracked::Resize)(Tracked::LookupInterface),
		(Grid::Activate)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using This = Single<TIsoGrid>;
	using Traits = Impl::Traits<This>;

	/// Small epsilon value within which we consider vertex position as "exact".
	static constexpr Distance epsilon = std::numeric_limits<Distance>::epsilon();
	/// Dimension of isogrid to polygonise.
	static constexpr Dim t_dims = Traits::t_dims;

	using ActivateImpl = Impl::Mixin::Tracked::Activate<This>;
	using GeomImpl = Impl::Mixin::Poly::Geom<This>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<This>;
	using ResetImpl = Impl::Mixin::Tracked::SingleList::Reset<This>;
	using ResizeImpl = Impl::Mixin::Tracked::Resize<This>;

	using Lookup = typename Traits::Lookup;
	/// Isogrid to (partially) polygonise.
	using IsoGrid = typename Traits::IsoGrid;
	/// Spatial partition type this poly will be responsible for.
	using IsoChild = typename IsoGrid::Child;
	/// Spatial partition type this poly will be responsible for.
	using IsoLookup = typename IsoChild::Lookup;
	/// Vertex index tuple type (for spatial lookup grid).
	using IdxTuple = typename Traits::Leaf;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	/// Float vector.
	using VecDf = Felt::VecDf<t_dims>;
	/// Cube edge type.
	using Edge = typename GeomImpl::Edge;
public:
	/// Vertex type.
	using Vertex = typename GeomImpl::Vertex;
	/// Simplex type.
	using Simplex = typename GeomImpl::Simplex;
private:
	/// Vertex array type for vertex storage.
	using VtxArray = std::vector<Vertex>;
	/// Simplex array type for simplex storage.
	using SpxArray = std::vector<Simplex>;


	/// Isogrid to (partially) polygonise.
	const IsoGrid*			m_pisogrid;
	const IsoLookup*		m_pisolookup;

	/// List of interpolated vertices.
	VtxArray	m_a_vtx;
	/// List of simplexes (i.e. lines for 2D or triangles for 3D).
	SpxArray	m_a_spx;

public:
	using ActivateImpl::is_active;
	using ActivateImpl::activate;
	using ResizeImpl::offset;
	using ResizeImpl::size;

	/**
	 * Construct a non-partitioned polygonisation of an isogrid.
	 *
	 * @param isogrid_ grid to be (partially) polygonised.
	 */
	Single(const IsoGrid& isogrid_) :
		ActivateImpl{IdxTuple::Constant(Felt::null_idx)}, LookupInterfaceImpl{Lookup{}},
		m_pisogrid{&isogrid_}, m_pisolookup{nullptr}
	{}

	/**
	 * Destroy the internal data array and lookup grid.
	 */
	void deactivate()
	{
		ActivateImpl::deactivate();
		m_a_vtx.resize(0);
		m_a_vtx.shrink_to_fit();
		m_a_spx.resize(0);
		m_a_spx.shrink_to_fit();
	}

	/**
	 * Reset without deallocating.
	 *
	 * Visit all vertices in lookup grid and set to null value then resize vertex and simplex lists.
	 */
	void reset()
	{
		ResetImpl::reset();
		m_a_vtx.resize(0);
		m_a_spx.resize(0);
	}

	/**
	 * Resize to fit size of isogrid child spatial partition.
	 *
	 * Will resize to one more than isochild size, since neighbouring Polys must overlap.
	 *
	 * @param size_ size of isogrid child partition.
	 * @param offset_ offset of isogrid child partition.
	 */
	void resize(const VecDi& size_, const VecDi& offset_)
	{
		static const VecDi one = VecDi::Constant(1);
		static const VecDi two = VecDi::Constant(2);

		const VecDi& size = size_ + two;
		const VecDi& offset = offset_ - one;

		ResizeImpl::resize(size, offset);
	}

	/**
	 * Update the polygonisation from the stored pointer to isogrid child lookup.
	 */
	void march()
	{
		for (TupleIdx list_idx = 0; list_idx < m_pisolookup->num_lists; list_idx++)
		{
			for (PosIdx pos_idx_leaf : m_pisolookup->list(list_idx))
			{
				spx(m_pisolookup->index(pos_idx_leaf));
			}
		}
	}

	/**
	 * Bind this Poly to the given Lookup grid giving positions to march over.
	 *
	 * I.e. a child spatial partition of the isogrid.
	 *
	 * @param pisolookup
	 */
	void bind(const IsoLookup& isolookup)
	{
		m_pisolookup = &isolookup;
	}

	/**
	 * Get a pointer to the isogrid child's lookup grid that gives points to polygonise.
	 *
	 * @return lookup grid to iterate over.
	 */
	IsoLookup const* bind() const
	{
		return m_pisolookup;
	}

	/**
	 * Get the vertex array.
	 *
	 * @return
	 */
	const VtxArray& vtxs() const
	{
		return m_a_vtx;
	}

	/**
	 * Get the array of simplices.
	 *
	 * @return
	 */
	const SpxArray& spxs() const
	{
		return m_a_spx;
	}

private:
	/**
	 * Generate simplex(es) for isogrid grid at position pos.
	 *
	 * @param pos
	 * @param m_isogrid
	 */
	void spx(const VecDi& pos)
	{
		// TODO: this is here for consistency only, since the marching
		// cubes implementation marches in the negative z-axis, but
		// positive x and y axes. Hence an offset is required so that the
		// negative z-axis marching is compensated by shifting the
		// calculation in the +z direction by one grid node.
		// (NOTE: has no effect for 2D).
		const VecDi pos_calc = pos - GeomImpl::SpxGridPosOffset;

		// Get corner inside-outside bitmask at this position.
		const unsigned short mask = this->mask(pos_calc);
		// Array of indices of zero-crossing vertices along each axis from
		// this corner.
		unsigned vtx_idxs[GeomImpl::num_edges];
		// Look up the edges that are crossed from the corner mask.
		unsigned short vtx_mask = GeomImpl::vtx_mask[mask];
		const short* vtx_order = GeomImpl::vtx_order[mask];

		// Cube corners are all inside or all outside.
		if (vtx_order[0] == -1)
			return;

		// Loop over each crossed edge in the cube, looking up
		// (or calculating, if unavailable) the vertices at the
		// zero-crossing.
		for (ListIdx edge_idx = 0; edge_idx < GeomImpl::num_edges; edge_idx++ )
		{
			// Check if current edge is crossed by the zero curve.
			if ((vtx_mask >> edge_idx) & 1)
			{
				const Edge& edge = GeomImpl::edges[edge_idx];
				// Edges are defined as an axis and an offset.
				// Look up index of vertex along current edge.
				vtx_idxs[edge_idx] = unsigned(idx(pos_calc + edge.offset, edge.axis));
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
//					if (dist <= This::epsilon())
//						return;
//				}
//			}

		// Join the vertices along each edge that the surface crosses to
		// make a simplex (or simplices).
		// The vtx_order lookup translates corner in-out mask to CCW
		// vertex ordering. We take D elements at a time from the lookup,
		// with each successive subset of D elements forming the next
		// simplex.
		for (ListIdx order_idx = 0; vtx_order[order_idx] != -1; order_idx += t_dims )
		{
			Simplex simplex;
			// A simplex for number of dimensions D has D vertices,
			// i.e. D endpoints.
			for (Dim endpoint = 0; endpoint < t_dims; endpoint++)
			{
				// Each vertex of the simplex is stored as an index
				// reference into the 'global' vertex array.
				simplex.idxs(endpoint) = vtx_idxs[ vtx_order[order_idx+ListIdx(endpoint)] ];
			}
			// Append the simplex to the list of simplices that make up the
			// polygonisation of this grid location.
			m_a_spx.push_back(std::move(simplex));
		}
	} // End spx.


	/**
	 * Lookup, or calculate then store, and return the index into the vertex
	 * array of a vertex at the zero-crossing of isogrid at pos_a along axis.
	 *
	 * @return
	 */
	ListIdx idx(const VecDi& pos_a, const Dim axis)
	{
		// Check lookup to see if vertex has already been calculated.
		const ListIdx idx_lookup = this->get(pos_a)(axis);
		if (idx_lookup != Felt::null_idx) {
			return idx_lookup;
		}

		// Position of opposite endpoint.
		VecDi pos_b(pos_a);
		pos_b(axis) += 1;

		// Value of isogrid at each endpoint of this edge.
		const Distance val_a = m_pisogrid->get(pos_a);
		const Distance val_b = m_pisogrid->get(pos_b);

		// The newly created vertex.
		Vertex vtx;

		// Check if lies very close to an endpoint or midpoint, if so then no need (and possibly
		// dangerous) to interpolate.
		if (std::abs(val_a) <= epsilon) {
			vtx = Vertex(m_pisogrid, pos_a);
		} else if (std::abs(val_b) <= epsilon) {
			vtx = Vertex(m_pisogrid, pos_b);
		} else {
			Distance mu;

			// If close to midpoint then put at midpoint.
			if (std::abs(val_a - val_b) <= epsilon) {
				mu = Distance(0.5f);
			} else
			// Otherwise interpolate between endpoints.
			{
				mu = val_a / (val_a - val_b);
			}

			const VecDf vec_a = pos_a.template cast<Distance>();
			const VecDf vec_b = pos_b.template cast<Distance>();
			const VecDf vec_c = vec_a + (vec_b - vec_a) * mu;

			vtx = Vertex(m_pisogrid, vec_c);
		}

		// Append newly created vertex to the cache and return a reference  to it.
		const ListIdx idx = m_a_vtx.size();
		m_a_vtx.push_back(std::move(vtx));
		this->get(pos_a)(axis) = idx;
		this->lookup().track(this->lookup().index(pos_a));
		return idx;
	}

	/**
	 * Calculate corner mask of cube at pos, based on inside-outside status
	 * of corners in isogrid.
	 *
	 * @param isogrid
	 * @param pos
	 * @return
	 */
	unsigned short mask (const VecDi& pos_) const
	{
		// Num corners == 2^D.  That is, 4 for 2D, 8 for 3D.
		unsigned short mask = 0;
		const ListIdx num_corners = (1 << t_dims);
		for (ListIdx idx = 0; idx < num_corners; idx++)
		{
			const VecDi corner = pos_ + GeomImpl::corners[idx];
			const Distance val = m_pisogrid->get(corner);
			mask = (unsigned short)(mask | ((val > 0) << idx));
		}
		return mask;
	}
};
} // Poly.
} // Impl.
} // Felt.


namespace Felt
{
namespace Impl
{
/**
 * Traits for Poly::Single.
 *
 * @tparam IsoGrid isogrid type to polygonise.
 */
template <class TIsoGrid>
struct Traits< Poly::Single<TIsoGrid> >
{
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<TIsoGrid>::t_dims;
	/// A vertex index for each positively directed edge stored at each grid node.
	using Leaf = VecDT<ListIdx, t_dims>;
	/// Type of lookup grid for tracking active positions.
	using Lookup = Lookup::LazySingleListSingleIdx<t_dims>;
	/// IsoGrid type that will be polygonised.
	using IsoGrid = TIsoGrid;
};

} // Impl.
} // Felt.

#endif /* INCLUDE_FELT_IMPL_POLY_HPP_ */
