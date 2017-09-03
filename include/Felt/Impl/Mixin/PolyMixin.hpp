#ifndef INCLUDE_FELT_IMPL_MIXIN_POLYMIXIN_HPP_
#define INCLUDE_FELT_IMPL_MIXIN_POLYMIXIN_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Mixin/PartitionedMixin.hpp>
#include <Felt/Impl/Mixin/TrackedMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Poly
{


template <class Derived>
class Activate : protected Tracked::Activate<Derived>
{
private:
	/// Base class.
	using Base = Felt::Impl::Mixin::Tracked::Activate<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Tuple of vertex indices.
	using IdxTuple = typename TraitsType::LeafType;

	// Do not expose from base.
	using Base::activate;
protected:
	Activate() : Base::Activate{IdxTuple::Constant(Felt::null_idx)}
	{}

	using Base::background;
	using Base::is_active;

	/**
	 * Allocate the internal data array and lookup grid.
	 */
	void activate()
	{
		Base::activate();
	}
	/**
	 * Destroy the internal data array and lookup grid.
	 */
	void deactivate()
	{
		Base::deactivate();
		pself->m_a_vtx.resize(0);
		pself->m_a_vtx.shrink_to_fit();
		pself->m_a_spx.resize(0);
		pself->m_a_spx.shrink_to_fit();
	}
};


template <class Derived>
class Children : protected Partitioned::Children<Derived>
{
private:
	using Base = Partitioned::Children<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Isogrid to (partially) polygonise.
	using IsoGridType = typename TraitsType::IsoGridType;
	/// Spatial partition type this poly will be responsible for.
	using IsoChildType = typename IsoGridType::ChildType;

protected:
	/**
	 * Construct and initialise children grid to hold child sub-grids.
	 */
	Children(
		const IsoGridType& isogrid_
	) : Base{isogrid_.size(), isogrid_.offset(), isogrid_.child_size(), ChildType(isogrid_)}
	{
		// Bind child polygonisation to isogrid child.
		for (PosIdx pos_idx = 0; pos_idx < this->children().data().size(); pos_idx++)
		{
			this->children().get(pos_idx).bind(isogrid_.children().get(pos_idx).lookup());
		}
	}
};


template <class Derived>
class Update
{
private:
	using Base = Partitioned::Children<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Surface type.
	using SurfaceType = typename TraitsType::SurfaceType;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Isogrid to (partially) polygonise.
	using IsoGridType = typename TraitsType::IsoGridType;
	/// Spatial partition type this poly will be responsible for.
	using IsoChildType = typename IsoGridType::ChildType;
	/// Lookup grid to track partitions containing zero-layer points.
	using ChangesGridType = Impl::Lookup::SingleListSingleIdx<Traits<IsoGridType>::t_dims>;

	SurfaceType const* m_psurface;
	std::unique_ptr<ChangesGridType> m_pgrid_update_pending;
	std::unique_ptr<ChangesGridType> m_pgrid_update_done;

protected:
	Update(const SurfaceType& surface_) :
		m_psurface{&surface_},
		m_pgrid_update_pending{
			std::make_unique<ChangesGridType>(
				surface_.isogrid().children().size(), surface_.isogrid().children().offset()
			)
		},
		m_pgrid_update_done{
			std::make_unique<ChangesGridType>(
				surface_.isogrid().children().size(), surface_.isogrid().children().offset()
			)
		}
	{}

	/**
	 * Notify of an update to the surface in order to track changes.
	 *
	 * This should be called whenever the surface is updated to ensure that eventual
	 * re-polygonisation only needs to update those spatial partitions that have actually changed.
	 */
	void notify()
	{
		static const TupleIdx num_lists = m_psurface->isogrid().children().lookup().num_lists;

		// Cycle outermost bands of delta update spatial partitions. We have three cases:
		// * Partition is currently polygonised, and since its in the delta grid it needs updating.
		// * Partition is not currently polygonised, but isogrid is tracking it, in which case it
		//   needs polygonising.
		// * Partition is not currently polygonised and isogrid is no longer tracking, so we no
		//   longer need to remember it.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1) {
			for (const PosIdx pos_idx_child : m_psurface->delta(layer_idx))
			{
				bool is_active = pself->children().get(pos_idx_child).is_active();

				for (
					TupleIdx layer_idx = 0; !is_active && layer_idx <= num_lists;
					layer_idx += num_lists - 1
				) {
					is_active = m_psurface->isogrid().children().lookup().is_tracked(
						pos_idx_child, layer_idx
					);
				}

				if (is_active)
					m_pgrid_update_pending->track(pos_idx_child);
				else
					m_pgrid_update_pending->untrack(pos_idx_child);
			}
		}

		// Cycle outermost status change lists, where a child may need resetting.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1)
		{
			for (
				const PosIdx pos_idx_child : m_psurface->status_change(layer_idx)
			) {
				if (pself->children().get(pos_idx_child).is_active())
					m_pgrid_update_pending->track(pos_idx_child);
			}
		}
	}

	/**
	 * Repolygonise partitions marked as changed since last polygonisation.
	 */
	void march()
	{
		#pragma omp parallel for
		for (ListIdx list_idx = 0; list_idx < m_pgrid_update_pending->list().size(); list_idx++)
		{
			const PosIdx pos_idx_child = m_pgrid_update_pending->list()[list_idx];

			ChildType& child = pself->children().get(pos_idx_child);

			// If isogrid child partition is active, then polygonise it.
			if (m_psurface->isogrid().children().get(pos_idx_child).is_active())
			{
				if (child.is_active())
					child.reset();
				else
					child.activate();
				child.march();

			} else {
				// Otherwise isogrid child partition has become inactive, so destroy the poly child.
				if (child.is_active())
					child.deactivate();
			}
		}

		std::swap(m_pgrid_update_pending, m_pgrid_update_done);
		m_pgrid_update_pending->reset();
	}


	/**
	 * Add all active poly childs and isogrid childs to change tracking for (re)polygonisation.
	 */
	void invalidate()
	{
		static const TupleIdx num_lists = m_psurface->isogrid().children().lookup().num_lists;

		// Remove pending changes, we're about to reconstruct the list.
		m_pgrid_update_pending->reset();
		// Flag curently active Poly::Single childs for re-polygonisation (or deactivation).
		for (const PosIdx pos_idx_child : pself->children().lookup().list())
			m_pgrid_update_pending->track(pos_idx_child);

		// Flag active outer-layer isogrid childs for re-polygonisation.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1)
		{
			for (
				const PosIdx pos_idx_child :
				m_psurface->isogrid().children().lookup().list(layer_idx)
			) {
				m_pgrid_update_pending->track(pos_idx_child);
			}
		}
	}

	/**
	 * Get list of partitions that were updated.
	 *
	 * @return list of position indices of partitions that were repolygonised in the last `march`.
	 */
	const PosArray& changes() const
	{
		return m_pgrid_update_done->list();
	}
};


template <class Derived>
class Reset : protected Tracked::SingleList::Reset<Derived>
{
private:
	using Base = Tracked::SingleList::Reset<Derived>;
	/// Dimension of the grid.
	static const Dim t_dims = Traits<Derived>::t_dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:

	/**
	 * Reset without deallocating.
	 *
	 * Visit all vertices in lookup grid and set to null value then resize vertex and simplex lists.
	 */
	void reset()
	{
		Base::reset();
		pself->m_a_vtx.resize(0);
		pself->m_a_spx.resize(0);
	}
};


template <class Derived>
class Resize : protected Tracked::Resize<Derived>
{
private:
	using Base = Tracked::Resize<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Isogrid child type.
	using IsoChild = typename TraitsType::IsoGridType::ChildType;
	/// Dimension of the grid.
	static const Dim t_dims = TraitsType::t_dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;
protected:

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

		Base::resize(size, offset);
	}
};


template <class Derived>
class March {
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Isogrid to (partially) polygonise.
	using IsoGridType = typename TraitsType::IsoGridType;
	/// Spatial partition type this poly will be responsible for.
	using IsoChildType = typename IsoGridType::ChildType;
	/// Spatial partition type this poly will be responsible for.
	using IsoLookupType = typename IsoChildType::LookupType;
	/// Vertex index tuple type (for spatial lookup grid).
	using IdxTuple = typename TraitsType::LeafType;
	/// Vertex type.
	using Vertex = typename TraitsType::Vertex;
	/// Simplex type.
	using Simplex = typename TraitsType::Simplex;
	/// Vertex array type for vertex storage.
	using VtxArray = std::vector<Vertex>;
	/// Simplex array type for simplex storage.
	using SpxArray = std::vector<Simplex>;
	/// Cube edge type.
	using Edge = typename TraitsType::Edge;
	/// Dimension of isogrid to polygonise.
	static constexpr Dim t_dims = TraitsType::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	/// Float vector.
	using VecDf = Felt::VecDf<t_dims>;

	/// Small epsilon value within which we consider vertex position as "exact".
	static constexpr Distance epsilon = std::numeric_limits<Distance>::epsilon();

	/// Isogrid to (partially) polygonise.
	const IsoGridType*		m_pisogrid;
	IsoLookupType const*	m_pisolookup;

protected:
	/// List of interpolated vertices.
	VtxArray	m_a_vtx;
	/// List of simplexes (i.e. lines for 2D or triangles for 3D).
	SpxArray	m_a_spx;

	March(const IsoGridType& isogrid_) : m_pisogrid{&isogrid_}, m_pisolookup{nullptr}
	{}

	/**
	 * Update the
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
	void bind(const IsoLookupType& isolookup)
	{
		m_pisolookup = &isolookup;
	}

	/**
	 * Get a pointer to the isogrid child's lookup grid that gives points to polygonise.
	 *
	 * @return lookup grid to iterate over.
	 */
	IsoLookupType const* bind() const
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
		const VecDi pos_calc = pos - TraitsType::SpxGridPosOffset;

		// Get corner inside-outside bitmask at this position.
		const unsigned short mask = this->mask(pos_calc);
		// Array of indices of zero-crossing vertices along each axis from
		// this corner.
		ListIdx vtx_idxs[TraitsType::num_edges];
		// Look up the edges that are crossed from the corner mask.
		unsigned short vtx_mask = TraitsType::vtx_mask[mask];
		const short* vtx_order = TraitsType::vtx_order[mask];

		// Cube corners are all inside or all outside.
		if (vtx_order[0] == -1)
			return;

		// Loop over each crossed edge in the cube, looking up
		// (or calculating, if unavailable) the vertices at the
		// zero-crossing.
		for (ListIdx edge_idx = 0; edge_idx < TraitsType::num_edges; edge_idx++ )
		{
			// Check if current edge is crossed by the zero curve.
			if ((vtx_mask >> edge_idx) & 1)
			{
				const Edge& edge = TraitsType::edges[edge_idx];
				// Edges are defined as an axis and an offset.
				// Look up index of vertex along current edge.
				vtx_idxs[edge_idx] = idx(pos_calc + edge.offset, edge.axis);
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
//					if (dist <= ThisType::epsilon())
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
		const ListIdx idx_lookup = pself->get(pos_a)(axis);
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

			const VecDf vec_a = pos_a.template cast<FLOAT>();
			const VecDf vec_b = pos_b.template cast<FLOAT>();
			const VecDf vec_c = vec_a + (vec_b - vec_a) * mu;

			vtx = Vertex(m_pisogrid, vec_c);
		}

		// Append newly created vertex to the cache and return a reference  to it.
		const ListIdx idx = m_a_vtx.size();
		m_a_vtx.push_back(std::move(vtx));
		pself->get(pos_a)(axis) = idx;
		pself->lookup().track(pself->lookup().index(pos_a));
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
			const VecDi corner = pos_ + TraitsType::corners[idx];
			const Distance val = m_pisogrid->get(corner);
			mask |= (unsigned short)((val > 0) << idx);
		}
		return mask;
	}
};


template <Dim D, typename Dummy>
class Traits
{
	static_assert(D == 2 || D == 3, "Poly only supports 2D or 3D polygonisations");
};


/**
 * 2D-specific definitions.
 */
template <typename Dummy>
struct Traits<2, Dummy> {
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
		template <typename PosType, class GridType>
		Vertex(const GridType* grid, const PosType& pos)
		{
			(void)grid;
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

	/**
	 * Number of edges on a square.
	 */
	static const short num_edges = 4;

	/**
	 * A lookup from inside/outside status bitmask to vertex ordering to
	 * create representative simplices (lines).
	 */
	static const short vtx_order [][4];

	/**
	 * A lookup of offsets from start position to corners of a cube.
	 */
	static const std::array<Vec2i, 1 << 2> corners;

	/**
	 * A lookup of cube edges defined by offset and axis
	 */
	static const Edge edges [];

	/**
	 * A lookup from corner inside/outside status bitmask to cut edge status
	 * bitmask.
	 */
	static const unsigned short vtx_mask [];

	/**
	 * Offset to normalise marching squares/cubes corner ordering.
	 *
	 * The marching cubes implementation marches in the negative z-axis, but positive x and y axes.
	 * Hence an offset is required so that the negative z-axis marching is compensated by shifting
	 * the calculation in the +z direction by one grid node.
	 * (NOTE: has no effect for 2D).
	 */
	static const Vec2i SpxGridPosOffset;
};

/**
 * 3D-specific definitions.
 */
template <typename Dummy>
struct Traits<3, Dummy> {
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
		template <typename PosType, class GridType>
		Vertex(const GridType* grid, const PosType& pos)
		{
			this->pos = pos.template cast<FLOAT>();
			this->norm = grid->grad(pos);
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

	/**
	 * A lookup of offsets from start position to corners of a cube.
	 */
	static const std::array<Vec3i, 1 << 3> corners;

	/**
	 * A lookup of cube edges defined by offset and axis
	 */
	static const Edge edges [];

	/**
	 * A lookup from corner inside/outside status bitmask to cut edge status
	 * bitmask.
	 */
	static const unsigned short vtx_mask [];

	/**
	 * @copydoc PolyBase<2,Derived>::SpxGridPosOffset
	 */
	static const Vec3i SpxGridPosOffset;
};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// 2D lookups.
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
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

/*
 * Relative position of corners in CCW order.
 */
template <typename Dummy>
const std::array<Vec2i, 4> Traits<2, Dummy>::corners = {
	Vec2i(0, 0),
	Vec2i(1, 0),
	Vec2i(1, 1),
	Vec2i(0, 1)
};

/*
 * Array of edge definitions (offset, direction) matching ::corners.
 */
template <typename Dummy>
const typename Traits<2, Dummy>::Edge Traits<2, Dummy>::edges [] = {
	// (x,y,axis)
	{ Vec2i(0, 0), 0 },
	{ Vec2i(1, 0), 1 },
	{ Vec2i(0, 1), 0 },
	{ Vec2i(0, 0), 1 }
};
template <typename Dummy>
const Vec2i Traits<2, Dummy>::SpxGridPosOffset(0,0);

/*
 * Lookup from corner mask to edge mask.
 */
template <typename Dummy>
const unsigned short Traits<2, Dummy>::vtx_mask [] ={
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

/*
 * A lookup from inside/outside status bitmask to vertex ordering to
 * create representative simplices (lines).
 */
template <typename Dummy>
const short Traits<2, Dummy>::vtx_order [][4] = {
	{ -1, -1, -1, -1 },
	{  3,  0, -1, -1 },
	{  0,  1, -1, -1 },
	{  3,  1, -1, -1 },
	{  1,  2, -1, -1 },
	{  3,  0,  1,  2 },
	{  0,  2, -1, -1 },
	{  3,  2, -1, -1 },
	{  2,  3, -1, -1 },
	{  2,  0, -1, -1 },
	{  2,  1,  0,  3 },
	{  2,  1, -1, -1 },
	{  3,  1, -1, -1 },
	{  1,  0, -1, -1 },
	{  0,  3, -1, -1 },
	{ -1, -1, -1, -1 }
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// 3D lookups.
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*
 * Lookup from corner mask to edge mask.
 */
template <typename Dummy>
const unsigned short Traits<3, Dummy>::vtx_mask [] = {
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
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

/**
 * Given node, march cube toward back, up and right (0,0,0)->(1,1,-1).
 */
template <typename Dummy>
const std::array<Vec3i, 8> Traits<3, Dummy>::corners = {
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
template <typename Dummy>
const typename Traits<3, Dummy>::Edge Traits<3, Dummy>::edges [] = {
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

template <typename Dummy>
const Vec3i Traits<3, Dummy>::SpxGridPosOffset(0,0,-1);

/**
 * Ordering of vertices to build simplex(es).
 */
template <typename Dummy>
const short Traits<3, Dummy>::vtx_order [][16] = {
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


} // Poly
} // Mixin
} // Impl
} // Felt
#endif /* INCLUDE_FELT_IMPL_MIXIN_POLYMIXIN_HPP_ */
