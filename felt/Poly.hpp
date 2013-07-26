#ifndef POLYBASE_H
#define POLYBASE_H
#include <eigen3/Eigen/Dense>
#include <omp.h>

#include "Grid.hpp"


namespace felt {

	template <UINT D>
	class PolyBase
	{
	};


	template <>
	class PolyBase<2>
	{
		typedef Eigen::Matrix<FLOAT, 2, 1> VecDf;
	public:
		struct Vertex {
			Vec2f pos;

			template <typename PosType>
			Vertex(const Grid<FLOAT,2>&, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
			}
			Vertex(){}
		};
	};


	template <>
	class PolyBase<3>
	{
		typedef Eigen::Matrix<FLOAT, 3, 1> VecDf;
	public:
		struct Vertex {
			Vec3f pos;
			Vec3f norm;

			template <typename PosType>
			Vertex(const Grid<FLOAT,3>& grid, const PosType& pos)
			{
				this->pos = pos.template cast<FLOAT>();
				this->norm = grid.gradC(pos);
				this->norm.normalize();

			}
			Vertex(){}

		};
	};


	template <UINT D>
	class Poly : public PolyBase<D>
	{
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;

	protected:
		Grid<UINT, D> m_grid_idx;

	public:
		Poly(const VecDu& dims, const VecDi& offset) :
		m_grid_idx(dims, offset)
		{
		}

		Grid<UINT, D>& idx()
		{
			return m_grid_idx;
		}

		static typename PolyBase<D>::Vertex Lerp(const Grid<FLOAT,D> phi, const VecDi& pos_a, const UINT& axis)
		{
			const FLOAT val_small = std::numeric_limits<FLOAT>::epsilon();
			VecDi pos_b(pos_a);
			pos_b(axis) += 1;

			const FLOAT val_a = phi(pos_a);
			const FLOAT val_b = phi(pos_b);


			// Check if lies very close to an endpoint or midpoint,
			// if so then no need (and possibly dangerous) to interpolate.
			if (std::abs(val_a) <= val_small)
			{
				return typename PolyBase<D>::Vertex(phi, pos_a);
			}
			else if (std::abs(val_b) <= val_small)
			{
				return typename PolyBase<D>::Vertex(phi, pos_b);
			}

			FLOAT mu;

			if (std::abs(val_a-val_b) <= val_small)
			{
				mu = (FLOAT)0.5;
			}
			else
			// Interpolate between endpoints.
			{
				mu = val_a / (val_a - val_b);
			}

			const VecDf vec_a = pos_a.template cast<FLOAT>();
			const VecDf vec_b = pos_b.template cast<FLOAT>();
			const VecDf vec_c = vec_a + (vec_b - vec_a) * mu;

			return typename PolyBase<D>::Vertex(phi, vec_c);
		}
	};
}

#endif // POLYBASE_H
