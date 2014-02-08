#ifndef POLYBASE_H
#define POLYBASE_H
#include <eigen3/Eigen/Dense>
#include <omp.h>
#include <vector>

#include "Grid.hpp"


namespace felt {

	template <UINT D>
	class PolyBase {
	public:
	};

	template <>
	class PolyBase<2> {
		typedef Eigen::Matrix<FLOAT, 2, 1> VecDf;
	public:

		struct Vertex {
			Vec2f pos;

			template <typename PosType>
			Vertex(const Grid<FLOAT, 2>&, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
			}

			Vertex()
			{
			}
		};

		struct Simplex {
			Vec2u idxs;
		};
		
		struct Edge {
			Vec2i offset;
			INT axis;
		};

		static const short num_edges = 4;
		static const Vec2i corners [];
		static const Edge edges [];
		static const short vtx_mask [];
		static const short vtx_order [][4];

	protected:

		~PolyBase<2>()
		{
		}
	};

	
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
	
	const Vec2i PolyBase<2>::corners [] ={
		Vec2i(0, 0),
		Vec2i(1, 0),
		Vec2i(1, 1),
		Vec2i(0, 1)
	};
	const PolyBase<2>::Edge PolyBase<2>::edges [] = {
		// (x,y,axis)
		{ Vec2i(0, 0), 0 },
		{ Vec2i(1, 0), 1 },
		{ Vec2i(0, 1), 0 },
		{ Vec2i(0, 0), 1 }
	};

	const short PolyBase<2>::vtx_mask [] ={
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

	template <>
	class PolyBase<3> {
		typedef Eigen::Matrix<FLOAT, 3, 1> VecDf;
	public:

		struct Vertex {
			Vec3f pos;
			Vec3f norm;

			template <typename PosType>
			Vertex(const Grid<FLOAT, 3>& grid, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
				this->norm = grid.gradC(pos);
				this->norm.normalize();

			}

			Vertex()
			{
			}
		};

		struct Simplex {
			Vec3u idxs;
		};

		static const Vec3i corners [];
	protected:

		~PolyBase<3>()
		{
		}
	};

	// Given node, march cube toward back, up and right (0,0,0)->(1,1,-1).
	const Vec3i PolyBase<3>::corners [] ={
		Vec3i(0, 0, 0),
		Vec3i(1, 0, 0),
		Vec3i(1, 0, -1),
		Vec3i(0, 0, -1),
		Vec3i(0, 1, 0),
		Vec3i(1, 1, 0),
		Vec3i(1, 1, -1),
		Vec3i(0, 1, -1)
	};


	template <UINT D>
	class Poly : public PolyBase<D> {
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;
	public:
		typedef typename PolyBase<D>::Vertex Vertex;
		typedef typename PolyBase<D>::Simplex Simplex;

	protected:
		// Indices along each axis of phi grid of interpolated vertex.
		Grid<VecDu, D> m_grid_idx;
		// List of interpolated vertices.
		std::vector<Vertex> m_a_vtx;
		// List of simplexes (i.e. lines for 2D or triangles for 3D).
		std::vector<Simplex> m_a_spx;

	public:

		static inline VecDu null_vtx()
		{
			return VecDu::Constant(null_idx());
		}

		static inline UINT null_idx()
		{
			return std::numeric_limits<UINT>::max();
		}

		static unsigned short mask(const Grid<FLOAT, D>& phi, const VecDi& pos)
		{
			// Num corners == 2^D.  That is, 4 for 2D, 8 for 3D.
			unsigned short mask = 0;
			const UINT num_corners = (1 << D);
			for (UINT idx = 0; idx < num_corners; idx++) {
				//				const VecDi offset = Grid<FLOAT,D>::index(idx, VecDu::Constant(2));
				const VecDi corner = pos + PolyBase<D>::corners[idx];
				const FLOAT val = phi(corner);
				mask |= (val > 0) << idx;
			}
			return mask;
		}

		Poly(const VecDu& dims, const VecDi& offset) :
		m_grid_idx(dims, offset)
		{
			this->reset();
		}

		void reset()
		{
			// Fill grid with null value.
			this->idx().fill(null_vtx());
			this->vtx().resize(0);
		}

		std::vector<Vertex>& vtx()
		{
			return m_a_vtx;
		}

		const Vertex& vtx(const UINT& idx)
		{
			return this->vtx()[idx];
		}

		const Vertex& vtx(const Grid<FLOAT, D>& phi, const VecDi& pos_a, const UINT& axis)
		{
			return this->vtx(this->idx(phi, pos_a, axis));
		}

		Grid<VecDu, D>& idx()
		{
			return m_grid_idx;
		}

		UINT idx(const Grid<FLOAT, D>& phi, const VecDi& pos_a, const UINT& axis)
		{
			// Check lookup to see if vertex has already been calculated.
			const UINT& idx_lookup = this->idx()(pos_a) (axis);
			if (idx_lookup != null_idx()) {
				return idx_lookup;
			}

			// Position of opposite endpoint.
			VecDi pos_b(pos_a);
			pos_b(axis) += 1;

			// Arbitrary small value, below which consider the vertex to be precisely at one endpoint.
			const FLOAT val_small = std::numeric_limits<FLOAT>::epsilon();

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

				if (std::abs(val_a - val_b) <= val_small) {
					mu = (FLOAT) 0.5;
				} else
					// Interpolate between endpoints.
				{
					mu = val_a / (val_a - val_b);
				}

				const VecDf vec_a = pos_a.template cast<FLOAT>();
				const VecDf vec_b = pos_b.template cast<FLOAT>();
				const VecDf vec_c = vec_a + (vec_b - vec_a) * mu;

				vtx = Vertex(phi, vec_c);
			}

			// Append newly created vertex to the cache and return a reference to it.
			const UINT idx = this->vtx().size();
			this->vtx().push_back(vtx);
			this->idx()(pos_a) (axis) = idx;
			return idx;
		}

		std::vector<Simplex>& spx()
		{
			return m_a_spx;
		}

		const Simplex& spx(const UINT& idx)
		{
			return this->spx()[idx];
		}

		void spx(const Grid<FLOAT, D>& phi, const VecDi& pos, const unsigned short mask, std::vector<Simplex>& spxs)
		{
			
			UINT vtx_idxs[PolyBase<D>::num_edges];
			unsigned short vtx_mask = PolyBase<D>::vtx_mask[mask];
			
			for (UINT edge_idx = 0; edge_idx < PolyBase<D>::num_edges; edge_idx++) {
				if ((vtx_mask >> edge_idx) & 1) {
					typename PolyBase<D>::Edge edge = PolyBase<D>::edges[edge_idx];
					vtx_idxs[edge_idx] = this->idx(phi, pos+edge.offset, edge.axis);					
				}
			}
			
			// Check for degenerates.
			for (UINT edge_idx1 = 0; edge_idx1 < PolyBase<D>::num_edges - 1; edge_idx1++) {
				for (UINT edge_idx2 = edge_idx1+1; edge_idx2 < PolyBase<D>::num_edges; edge_idx2++) {
					if ((vtx_mask >> edge_idx1) & 1 && (vtx_mask >> edge_idx2) & 1) {
						const VecDf& pos1 = this->vtx(vtx_idxs[edge_idx1]).pos;
						const VecDf& pos2 = this->vtx(vtx_idxs[edge_idx2]).pos;
						if ((pos1 - pos2).squaredNorm() <= std::numeric_limits<FLOAT>::epsilon())
							return;
					}
				}							
			}			
			
			for (UINT order_idx = 0; PolyBase<D>::vtx_order[mask][order_idx] != -1; order_idx += D) {
				Simplex simplex;
				for (UINT endpoint = 0; endpoint < D; endpoint++) {
					simplex.idxs(endpoint) = vtx_idxs[PolyBase<D>::vtx_order[mask][order_idx+endpoint]];
					if (simplex.idxs(endpoint)*simplex.idxs(endpoint) > 1000000) {
						return;
					}
				}
				spxs.push_back(simplex);
			}
//			spxs.push_back(this->idx(phi, pos, 0));
		}

		void spx_add(const VecDi& pos, const unsigned short& vtx_mask, const short* vtx_order, const Grid<FLOAT, D>& phi)
		{

		}

		void calc(const Grid<FLOAT, D>& phi, const VecDi& pos)
		{
			// Num corners == 2^D.  That is, 4 for 2D, 8 for 3D.
			std::vector<VecDi> apos(1 << (D - 1));
			for (UINT idx = 0; idx < apos.size(); idx++) {
				const VecDi offset = Grid<FLOAT, D>::index(idx, VecDu::Constant(2));
				apos[idx] = pos + offset;

			}
		}
	};
}

#endif // POLYBASE_H
