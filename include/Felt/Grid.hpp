#ifndef Grid_hpp
#define Grid_hpp

#include <inttypes.h>
#include <math.h>
#include <vector>
#include <functional>
#include <eigen3/Eigen/Dense>
#include <boost/iterator/iterator_facade.hpp>

namespace felt
{
	/**
	 * Use 32 bit float by default.
	 */
	typedef float FLOAT;

	/**
	 * Use 32 bit int by default.
	 */
	typedef int INT;

	/**
	 * Use 32 bit unsigned int by default.
	 */
	typedef size_t UINT;

	/**
	 * Shorthand for 2D float vector.
	 */
	typedef Eigen::Matrix<FLOAT, 2, 1> Vec2f;
	/**
	 * Shorthand for 2D unsigned integer vector.
	 */
	typedef Eigen::Matrix<UINT, 2, 1> Vec2u;
	/**
	 * Shorthand for 2D integer vector.
	 */
	typedef Eigen::Matrix<INT, 2, 1> Vec2i;
	/**
	 * Shorthand for 3D float vector.
	 */
	typedef Eigen::Matrix<FLOAT, 3, 1> Vec3f;
	/**
	 * Shorthand for 3D unsigned integer vector.
	 */
	typedef Eigen::Matrix<UINT, 3, 1> Vec3u;
	/**
	 * Shorthand for 3D integer vector.
	 */
	typedef Eigen::Matrix<INT, 3, 1> Vec3i;

	/**
	 * Get the sign of a value (+/-1).
	 *
	 * @param val
	 * @return
	 */
	template <typename T> int sgn(T val)
	{
		return (T(0) < val) - (val < T(0));
	}

	/**
	 * ASM optimised logarithm to base 2.
	 *
	 * @param x
	 * @return
	 */
	static inline uint32_t log2(const uint32_t x)
	{
	  uint32_t y;
	  asm ( "\tbsr %1, %0\n"
		  : "=r"(y)
		  : "r" (x)
	  );
	  return y;
	}


	template <UINT D, UINT N>
	class LookupGrid;

	/**
	 * Abstract base class for n-dimensional grid.
	 *
	 * Subclasses can override the get() method to return a datatype other
	 * than that stored in the grid, so that grid values can be mutated
	 * before e.g. using in gradient calculations.
	 *
	 * @tparam T the data type to store in the grid.
	 * @tparam D the number of dimensions of the grid.
	 * @tparam R the data type return by get().
	 */
	template <typename T, UINT D, typename R=T>
	class GridBase
	{
	public:
		/**
		 * D-dimensional unsigned integer vector.
		 */
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		/**
		 * D-dimensional integer vector.
		 */
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		/**
		 * D-dimensional float vector.
		 */
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;
		/**
		 * D-dimensional vector of type T.
		 */
		typedef Eigen::Matrix<T, D, 1> VecDT;
		/**
		 * Dynamic 1D vector (i.e. a resizeable array of data) of
		 * type T.
		 */
		typedef Eigen::Array<T, 1, Eigen::Dynamic> ArrayData;
		/**
		 * Resizeable array of VecDi (i.e. grid locations).
		 *
		 * Uses Eigen::aligned_allocator for optimal address alignment.
		 */
		typedef std::vector<
			VecDi, Eigen::aligned_allocator<VecDi>
		> PosArray;

		typedef GridBase<T, D, R>	GridBase_t;

		class iterator : public boost::iterator_facade<
			GridBase::iterator, const VecDi&, boost::forward_traversal_tag
		> {
		private:
			friend class boost::iterator_core_access;

			UINT			m_idx;
			VecDi			m_pos;
			const UINT&		m_num_elems;;
			const VecDu&	m_dims;
			const VecDi&	m_offset;
		public:
			iterator() : m_idx(0), m_num_elems(0),
			m_dims(VecDu::Constant(0)), m_offset(VecDi::Constant(0))
			{}

			iterator(const GridBase& grid, const UINT& startIdx = 0)
			: m_idx(startIdx), m_num_elems(grid.data().size()),
			  m_dims(grid.dims()), m_offset(grid.offset()),
			  m_pos(
				GridBase<T, D, R>::index(startIdx, grid.dims(), grid.offset())
			  )
			{}
		 private:

		    void increment() {
		    	m_idx++;
		    	m_pos = GridBase<T, D, R>::index(m_idx, m_dims, m_offset);
		    }

		    bool equal(iterator const& other) const
		    {
		        return this->m_idx == other.m_idx;
		    }

		    const VecDi& dereference() const {
		    	return m_pos;
		    }
		};


	protected:
		/**
		 * The translational offset of the grid's zero coordinate.
		 */
		VecDi m_offset;
		/**
		 * The dimensions (size) of the grid.
		 */
		VecDu m_dims;
		/**
		 * The physical size of a grid node (used for spatial derivatives).
		 */
		FLOAT m_dx;
		/**
		 * The actual grid data store.
		 */
		ArrayData m_data;

	public:
		/**
		 * Trivial destructor.
		 */
		virtual ~GridBase ()
		{}

		/**
		 * Initialise a zero-size grid.
		 *
		 * If custom initialisation is required, then use this constructor
		 * and call init() in the derived class.
		 */
		GridBase () :
		m_offset(VecDi::Zero()),
		m_dims(VecDu::Zero()),
		m_dx(1)
		{
		}

		/**
		 * Initialise a grid with given dimension, offset and delta x.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		GridBase (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) :
		m_dx(1)
		{
			this->init(dims, offset, delta);
		}

		/**
		 * Initialise the grid dimensions and offset.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		virtual void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			this->dx(delta);
			this->dims(dims);
			this->offset(offset);
		}

		/**
		 * Set grid offset.
		 *
		 * The offset is used to 'centre' the grid, so that e.g. negative
		 * grid positions can be used. It is equal to the spatial position
		 * of the zero coordinate.
		 *
		 * @param offset_new
		 */
		virtual void offset (const VecDi& offset_new)
		{
			m_offset = offset_new;
		}

		/**
		 * Get the grid offset parameter.
		 *
		 * @return
		 */
		const VecDi& offset () const
		{
			return m_offset;
		}

		/**
		 * Get grid's delta x.
		 * @return
		 */
		inline const FLOAT& dx () const
		{
			return m_dx;
		}

		/**
		 * Set grid's delta x.
		 * @param delta
		 */
		void dx (const FLOAT& delta)
		{
			m_dx = delta;
		}

		/**
		 * Shorthand to access grid values.
		 *
		 * @param pos
		 * @return
		 */
		R& operator() (const VecDi& pos)
		{
			return this->get(pos);
		}

		/**
		 * Shorthand to access grid values (const version).
		 * @param pos
		 * @return
		 */
		const R& operator() (const VecDi& pos) const
		{
			return this->get(pos);
		}

		/**
		 * Get interpolated grid value.
		 *
		 * Passing a floating point position vector will initiate a linear
		 * interpolation of the grid at that real-valued location.
		 *
		 * @param pos
		 * @return
		 */
		const T operator() (const VecDf& pos) const
		{
			return this->interp(pos);
		}

		/**
		 * Get interpolated grid value.
		 *
		 * Convenience operator to return linearly interpolated value given
		 * a real-valued location.
		 *
		 * @param pos
		 * @return
		 */
		T operator() (const VecDf& pos)
		{
			return this->interp(pos);
		}

		/**
		 * Abstract getter to be overriden by subclasses.
		 *
		 * This allows the subclass to mutate the value stored in the grid
		 * before it is used.
		 *
		 * @param pos
		 * @return
		 */
		virtual R& get (const VecDi& pos) = 0;


		/**
		 * Abstract getter to be overriden by subclasses (const version).
		 *
		 * @param pos
		 * @return
		 */
		virtual const R& get (const VecDi& pos) const = 0;


		/**
		 * Get index in data array of position vector.
		 *
		 * The grid is packed in a 1D array, so this method is required to
		 * get the index in that array of the D-dimensional position.
		 *
		 * @param pos
		 * @return
		 */
		UINT index (const VecDi& pos) const
		{
			return GridBase<T,D,R>::index(pos, this->dims(), this->offset());
		}

		/**
		 * Get index in data array of position vector.
		 *
		 * The grid is packed in a 1D array, so this method is required to
		 * get the index in that array of the D-dimensional position.
		 *
		 * @param pos
		 * @param dims
		 * @param offset
		 * @return
		 */
		static UINT index (
			const VecDi& pos, const VecDu& dims,
			const VecDi& offset = VecDi::Constant(0)
		)
		{
			UINT idx = 0;
			for (INT i = 0; i < dims.size(); i++)
			{
				INT u_pos = pos(i) - offset(i);
				for (INT j = i+1; j < dims.size(); j++)
				{
					u_pos *= dims(j);
				}
				idx += u_pos;
			}
			return idx;
		}

		/**
		 * Get position of index.
		 *
		 * Given an index in the 1D grid data array, calculate the position
		 * vector that it pertains to.
		 *
		 * @param idx
		 * @return
		 */
		VecDi index (const UINT& idx) const
		{
			return GridBase<T,D,R>::index(idx, this->dims(), this->offset());
		}

		/**
		 * Get position of index.
		 *
		 * Given an index and the dimensions and offset of a grid, calculate
		 * the position vector that the index pertains to in a representative
		 * 1D array.
		 *
		 * @param idx
		 * @param dims
		 * @param offset
		 * @return
		 */
		static VecDi index (
			UINT idx, const VecDu& dims, const VecDi& offset = VecDi::Zero()
		) {
/*
	Eg. 2D: row major order (3x4=12): (x,y)[idx] =>
		(0,0)[0], (0,1)[1], (0,2)[2],  (0,3)[3]
		(1,0)[4], (1,1)[5], (1,2)[6],  (1,3)[7]
		(2,0)[8], (2,1)[9], (2,2)[10], (2,3)[11]

	E.g. 3D:
	z = idx % Dz
	y = (idx/Dz) % Dy
	x = (idx/Dz)/Dy % Dx
*/
			VecDi pos(D);

			for (INT axis = dims.size()-1; axis >= 0; axis--)
			{
				pos(axis) = idx % dims(axis) + offset(axis);
				idx /= dims(axis);
			}

			return pos;
		}



		/**
		 * Retrieve a reference to the raw grid data array.
		 *
		 * @return
		 */
		ArrayData& data ()
		{
			return m_data;
		}

		/**
		 * Retrieve a reference to the raw grid data array (const version).
		 *
		 * @return
		 */
		const ArrayData& data () const
		{
			return m_data;
		}

		const iterator begin () const
		{
			return iterator(*this, 0);
		}

		const iterator end () const
		{
			return iterator(*this, this->data().size());
		}

		/**
		 * Set the dimensions of the grid and resize it.
		 *
		 * Values will be default initialised (or not at all).
		 *
		 * @param vec_NewDims
		 * @return
		 */
		virtual void dims (const VecDu& dims_new)
		{
			m_dims = dims_new;

			INT uGridBaseSize = m_dims(0);
			for (INT i = 1; i < m_dims.size(); i++)
			{
				uGridBaseSize *= m_dims(i);
			}
			m_data.resize(uGridBaseSize);
		}

		/**
		 * Get grid dimensions.
		 *
		 * @return
		 */
		const VecDu& dims () const
		{
			return m_dims;
		}

		/**
		 * Fill grid with a single value.
		 *
		 * @param val
		 */
		virtual void fill (const T& val)
		{
			this->data().fill(val);
		}

		/**
		 * Test if a position is inside the grid bounds.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		bool inside (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			const VecDu& dims = this->dims();
			const VecDi& offset = this->offset();

			for (INT i = 0; i < pos.size(); i++)
			{
				if (pos(i) >= (INT)dims(i) + offset(i))
					return false;
				if (pos(i) < offset(i))
					return false;
			}
			return true;
		}

		/**
		 * Get the neighbouring positions in the cardinal directions.
		 *
		 * Neighbour positions will be added to vout.
		 *
		 * If bcheck is true then duplicates will not be allowed in vout,
		 * using a linear search to ensure so.
		 *
		 * @param pos
		 * @param vout
		 * @param bcheck
		 */
		void neighs (
			const VecDi& pos, PosArray& vout, bool bcheck = false
		) const
		{
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Most likely all neighbours are valid.
			vout.reserve(2*D);
			// Position for look-around.
			VecDi vec_dir(pos);
			for (INT axis = 0; axis < dims.size(); axis++) {
				// Check if backward value is within GridBase.
				vec_dir(axis) -= 1;
				if (this->inside(vec_dir) && (!bcheck ||
					std::find(vout.begin(), vout.end(), vec_dir) == vout.end())
				)
				{
					vout.push_back(vec_dir);
				}
				// Check if forward value is within GridBase.
				vec_dir(axis) += 2;
				if (this->inside(vec_dir) && (!bcheck ||
					std::find(vout.begin(), vout.end(), vec_dir) == vout.end())
				)
				{
					vout.push_back(vec_dir);
				}
				vec_dir(axis) -= 1;
			}
		}

		/**
		 * Get the neighbouring positions in the cardinal directions.
		 *
		 * Neighbour positions will be added to vout.
		 *
		 * grid_check will be used as a lookup to ensure no duplicates are
		 * added to apos_out.
		 *
		 * @param pos
		 * @param apos_out
		 * @param grid_check
		 */
		void neighs (
			const VecDi& pos, const std::function<void(const VecDi&)>& fn
		) const {
			// Position for look-around.
			VecDi vec_dir(pos);
			for (INT axis = 0; axis < this->dims().size(); axis++)
			{
				// Check if backward value is within GridBase.
				vec_dir(axis) -= 1;
				if (this->inside(vec_dir))
					fn(vec_dir);

				// Check if forward value is within GridBase.
				vec_dir(axis) += 2;

				if (this->inside(vec_dir))
					fn(vec_dir);

				vec_dir(axis) -= 1;
			}
		}

		/**
		 * Forward difference gradient.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		VecDT gradF (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			// Value at this point.
			const T val_centre = (*this)(pos);
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-ahead.
			Eigen::Matrix<PosType, D, 1> vec_dir(pos);

			for (INT axis = 0; axis < dims.size(); axis++)
			{
				vec_dir(axis) += 1;
				vec_grad(axis) = (*this)(vec_dir) - val_centre;
				vec_dir(axis) -= 1;
			}

			return vec_grad / this->dx();
		}

		/**
		 * Backward difference gradient
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		VecDT gradB (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			// Value at this point.
			const T val_centre = (*this)(pos);
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-behind.
			Eigen::Matrix<PosType, D, 1> vec_dir(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				vec_dir(axis) -= 1;
				vec_grad(axis) = val_centre - (*this)(vec_dir);
				vec_dir(axis) += 1;
			}

			return vec_grad / this->dx();
		}

		/**
		 * Central difference gradient.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		VecDT gradC (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-around.
			Eigen::Matrix<PosType, D, 1> vec_dir(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				vec_dir(axis) -= 1;
				const T back = (*this)(vec_dir);
				vec_dir(axis) += 2;
				const T forward = (*this)(vec_dir);
				vec_dir(axis) -= 1;

				vec_grad(axis) = (forward - back) /  2;
			}

			return vec_grad / this->dx();
		}

		/**
		 * Safe gradient.
		 *
		 * Will calculate central, forward or backward difference along each
		 * axis, depending what grid values are available.
		 * That is, for grid points at the edge of the grid it will
		 * return forward/backward differences.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		VecDT grad (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			typedef Eigen::Matrix<PosType, D, 1>	VecDR;
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-around.
			VecDR pos_test(pos);

			// Central value.
			const T centre = (*this)(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				T back = centre;
				T forward = centre;
				UINT order = 0;
				// Check if backward value is within GridBase.
				pos_test(axis) -= 1;
				if (this->inside(pos_test)) {
					back = (*this)(pos_test);
					order++;
				}
				// Check if forward value is within GridBase.
				pos_test(axis) += 2;
				if (this->inside(pos_test)) {
					forward = (*this)(pos_test);
					order++;
				}
				pos_test(axis) -= 1;
				// Calculate central/forward/backward difference along this
				// axis.
				if (order != 0)
					vec_grad(axis) = (forward - back) / order;
				else
					vec_grad(axis) = 0;
			}

			return vec_grad / this->dx();
		}

		/**
		 * Entropy satisfying differencing.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		VecDT gradE (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			typedef Eigen::Matrix<PosType, D, 1> VecDp;
			// Value at this point.
			const T centre = (*this)(pos);
			// Reference to GridBase dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-around.
			VecDp pos_test(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				pos_test(axis) -= 1;
				T back = (*this)(pos_test);
				pos_test(axis) += 2;
				T forward = (*this)(pos_test);
				pos_test(axis) -= 1;

				back = std::min((centre - back), 0.0f);
				forward = std::max((forward - centre), 0.0f);

				vec_grad(axis) = (forward + back);
			}

			return vec_grad / this->dx();
		}

		/**
		 * Linear interpolation.
		 *
		 * @param vec_fpos
		 * @return
		 */
		T interp (const VecDf& vec_fpos) const
		{
			const VecDu& dims = this->dims();

			// Store all 2^d corners.
			std::vector< T > val_corners(1 << dims.size());

			// Get all corners of containing cell.
			for (UINT i = 0; i < val_corners.size(); i++)
			{
				// 0 = 00 => (x,y)
				// 1 = 01 => (x+1,y)
				// 2 = 10 => (x,y+1)
				// 3 = 11 => (x+1,y+1)

				VecDi pos_corner(dims.size());
				for (INT dim = 0; dim < pos_corner.size(); dim++)
				{
					INT pos = (INT)floor(vec_fpos(dim));
					const INT dir = (i >> dim) & 1;
					pos += dir;
					pos_corner(dim) = pos;
				}

				val_corners[i] = (*this)(pos_corner);
			}

			// Translate position vector into 'hypercube space',
			// so 0 <= v(x) <= 1.
			VecDf vec_dir(vec_fpos);
			for (INT dim = 0; dim < vec_dir.size(); dim++)
			{
				vec_dir(dim) = vec_dir(dim) - floor(vec_dir(dim));
			}

			// Repeatedly reduce along axes,
			// i.e. hypercube -> cube -> square -> line -> point
			while (val_corners.size() > 1)
			{
				val_corners = this->interp(val_corners, vec_dir);
			}

			return val_corners[0];
		}

		/**
		 * Curvature calculation based on difference of normals method.
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		T curv (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			typedef Eigen::Matrix<PosType, D, 1> VecDp;

			const T val_centre = (*this)(pos);
			const VecDu& dims = this->dims();
			VecDp vec_dir(pos);

			// Forward directed principal normal.
			VecDT n_forward(D);

			for (INT axis = 0; axis < dims.size(); axis++)
			{
				vec_dir(axis) += 1;

				const T val_axis = (*this)(vec_dir) - val_centre;
				T val_neighs_sq = 0;

				// Loop other dimensions to get central difference across them.
				for (
					INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++
				) {
					// Only getting differences across other axes.
					if (axis_neigh != axis) {
						// Central difference across this forward point.
						VecDp vec_dir_neigh(vec_dir);
						vec_dir_neigh(axis_neigh) -= 1;
						const T val_low = (*this)(vec_dir_neigh);
						vec_dir_neigh(axis_neigh) += 2;
						const T val_high = (*this)(vec_dir_neigh);

						const T val_neigh = (val_high - val_low) / 2;
						val_neighs_sq += val_neigh*val_neigh;
					}
				}

				n_forward(axis) =		val_axis /
								sqrt(val_axis*val_axis + val_neighs_sq);

				vec_dir(axis) -= 1;
			}


			// Backward directed principal normal.
			VecDT n_backward(D);

			for (INT axis = 0; axis < dims.size(); axis++)
			{
				vec_dir(axis) -= 1;
				VecDT vec_vals(D);

				const T val_axis = val_centre - (*this)(vec_dir);
				T val_neighs_sq = 0;

				// Loop other dimensions to get central difference across them.
				for (
					INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++
				) {
					// Only getting differences across other axes.
					if (axis_neigh != axis) {
						// Central difference across this backward point.
						VecDp vec_dir_neigh(vec_dir);
						vec_dir_neigh(axis_neigh) -= 1;
						const T val_low = (*this)(vec_dir_neigh);
						vec_dir_neigh(axis_neigh) += 2;
						const T val_high = (*this)(vec_dir_neigh);

						const T val_neigh = (val_high - val_low) / 2;
						val_neighs_sq += val_neigh*val_neigh;
					}
				}

				n_backward(axis) =		val_axis /
								sqrt(val_axis*val_axis + val_neighs_sq);

				vec_dir(axis) += 1;
			}

			const VecDT dn_by_dx = (n_forward - n_backward);

			T curvature = 0;
			for (INT axis = 0; axis < dims.size(); axis++)
			{
				curvature += dn_by_dx(axis);
			}
			curvature /= 2;

			return curvature;
		}

		/**
		 * Calculate divergence, i.e. d2f(x)/dx2
		 *
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		T divergence (const Eigen::Matrix<PosType, D, 1>& pos) const
		{
			const VecDT vec_grad_f = this->gradF(pos);
			const VecDT vec_grad_b = this->gradB(pos);
			const VecDT vec_grad_diff = vec_grad_b - vec_grad_f;

			// Component-wise sum.
			const T val = vec_grad_diff.sum();

			return val / (this->dx()*this->dx());
		}

		/**
		 * Get the value stored in the grid, circumventing subclass's mutation.
		 *
		 * @param pos
		 * @return
		 */
		T& get_internal (const VecDi& pos) {
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			assert_pos_bounds(pos, "get_internal: ");
			#endif
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}

		/**
		 * Get the value stored in the grid, circumventing subclass's mutation
		 * (const version).
		 *
		 * @param pos
		 * @return
		 */
		const T& get_internal (const VecDi& pos) const
		{
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			assert_pos_bounds(pos, "get_internal: ");
			#endif
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}

	protected:
		void assert_pos_bounds (
			const VecDi& pos, std::string title
		) const {
			if (!this->inside(pos))
			{
				using namespace Eigen;
				IOFormat fmt(
					StreamPrecision, DontAlignCols, ", ", ", ",
					"", "", "(", ")"
				);
				const VecDi& pos_min = this->offset();
				const VecDi& pos_max = (
					this->dims().template cast<INT>() + pos_min
					- VecDi::Constant(1)
				);
				std::stringstream err;
				err << title << pos.transpose().format(fmt)
					<< " is outside grid "
					<< pos_min.format(fmt) << "-" << pos_max.format(fmt)
					<< std::endl;
				std::string err_str = err.str();
				throw std::domain_error(err_str);
			}
		}

	public:

	#ifndef _TESTING
	protected:
	#endif


		/**
		 * Interpolate down one dimension.
		 *
		 * The values of val_corners_in are interpolated to one dimension
		 * smaller than they are currently (cube->square, square->line,
		 * line->point).
		 *
		 * @param val_corners_in
		 * @param vec_fpos
		 * @return
		 */
		std::vector<T> interp (
			const std::vector<T>& val_corners_in, const VecDf& vec_fpos
		) const
		{
			const size_t num_corners = val_corners_in.size();

			// Number of values returned.
			// This is a power of 2 less than input dimensions
			// (cube becomes square, square becomes line, line becomes point).
			const size_t num_out = num_corners >> 1;

			// The axis along which to interpolate.
			// This is computed from the dimensions of the original input and
			// the dimensions of the intended output.
			const size_t axis_idx = vec_fpos.size() - log2(num_corners);

			// The weighting to be used in interpolating each pair of points.
			// This is the position along the axis of interpolation.
			const FLOAT axis_pos = vec_fpos(axis_idx);

			std::vector<T> val_corners_out(num_out);

			for (size_t i = 0; i < num_out; i++)
			{
				const T low = val_corners_in[(i << 1)];
				const T high = val_corners_in[(i << 1) + 1];
				const T val = axis_pos*high + (1.0f-axis_pos)*low;
				val_corners_out[i] = val;
			}

			return val_corners_out;
		}

	};


	/**
	 * A standard n-dimensional grid for storing arbitrary values.
	 *
	 * @tparam T the type of data to store in the grid.
	 * @tparam D the number of dimensions of the grid.
	 */
	template <typename T, UINT D>
	class Grid : public GridBase<T, D, T>
	{
	protected:
		typedef GridBase<T, D, T>	Base;

	public:
		typedef typename Base::VecDu	VecDu;
		typedef typename Base::VecDi	VecDi;
		typedef typename Base::VecDf	VecDf;


		virtual ~Grid ()
		{}

		/**
		 * Initialise a zero-dimensional GridBase.
		 */
		Grid () : Base()
		{}

		/**
		 * Delegates to GridBase's constructor.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		Grid (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) : Base(dims, offset, delta)
		{}

		/**
		 * Override GridBase::get to simply return the value stored in the grid.
		 *
		 * @param pos
		 * @return
		 */
		virtual T& get (const VecDi& pos)
		{
			return this->get_internal(pos);
		}

		/**
		 * Override GridBase::get to simply return the value stored in the grid
		 * (const version).
		 *
		 * @param pos
		 * @return
		 */
		virtual const T& get (const VecDi& pos) const
		{
			return this->get_internal(pos);
		}
	};

}// End namespace felt.
#endif
