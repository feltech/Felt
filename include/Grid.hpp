#ifndef Lattice_hpp
#define Lattice_hpp

#include <inttypes.h>
#include <math.h>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace felt {

	typedef float FLOAT;
	typedef int INT;
	typedef size_t UINT;

	typedef Eigen::Matrix<FLOAT, 2, 1> Vec2f;
	typedef Eigen::Matrix<UINT, 2, 1> Vec2u;
	typedef Eigen::Matrix<INT, 2, 1> Vec2i;
	typedef Eigen::Matrix<FLOAT, 3, 1> Vec3f;
	typedef Eigen::Matrix<UINT, 3, 1> Vec3u;
	typedef Eigen::Matrix<INT, 3, 1> Vec3i;

	// Following signnum from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	static inline uint32_t log2(const uint32_t x) {
	  uint32_t y;
	  asm ( "\tbsr %1, %0\n"
		  : "=r"(y)
		  : "r" (x)
	  );
	  return y;
	}



	template <class T, UINT D>
	class Grid {

		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;
		typedef Eigen::Matrix<T, D, 1> VecDT;
		typedef Eigen::Array<T, 1, Eigen::Dynamic> ArrayData;

	protected:
		VecDi m_vec_offset;
		VecDu m_vec_dims;
		FLOAT m_dx;
		ArrayData m_vec_Data;

	public:

		/**
		 * Initialise a zero-dimensional grid.
		 */
		Grid () :
		m_vec_offset(VecDi::Zero()),
		m_vec_dims(VecDu::Zero()),
		m_dx(1)
		{
		}

		Grid (UINT x, UINT y) :
		m_dx(1)
		{
			this->init(Vec2u(x,y));
		}

		Grid (UINT x, UINT y, UINT z) :
		m_dx(1)
		{
			this->init(Vec3u(x,y,z));
		}

		/**
		 * Initialise a grid with given dimension, offset and delta x.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		Grid (const VecDu& dims, const VecDi& offset = VecDi::Zero(), const FLOAT& delta = 1) :
		m_dx(1)
		{
			this->init(dims, offset, delta);
		}

		void init (const VecDu& dims, const VecDi& offset = VecDi::Zero(), const FLOAT& delta = 1) {
			this->dx(delta);
			this->offset(offset);
			this->dims(dims);
		}

		/**
		 * Set grid offset.
		 *
		 * @return
		 */
		void offset (const VecDi& vec_offset) {
			m_vec_offset = vec_offset;
		}

		/**
		 * Get grid offset.
		 *
		 * @return
		 */
		const VecDi& offset () const {
			return m_vec_offset;
		}


		/**
		 * Get grid delta x.
		 * @return
		 */
		inline const FLOAT& dx () const {
			return m_dx;
		}

		/**
		 * Set grid delta x.
		 * @param delta
		 */
		void dx (const FLOAT& delta) {
			m_dx = delta;
		}

		/*
		 * Helpers for cases where D is 2 or 3.
		 */
		T& operator() (const INT& x, const INT& y) {
			return (*this)(VecDi(x,y));
		}
		const T& operator() (const INT& x, const INT& y) const {
			return (*this)(VecDi(x,y));
		}
		T operator() (const FLOAT& x, const FLOAT& y) {
			return (*this)(VecDi(x,y));
		}
		const T operator() (const FLOAT& x, const FLOAT& y) const {
			return (*this)(VecDi(x,y));
		}
		T& operator() (const INT& x, const INT& y, const INT& z) {
			return (*this)(VecDi(x,y,z));
		}
		const T& operator() (const INT& x, const INT& y, const INT& z) const {
			return (*this)(VecDi(x,y,z));
		}
		T operator() (const FLOAT& x, const FLOAT& y, const FLOAT& z) {
			return (*this)(VecDi(x,y,z));
		}
		const T operator() (const FLOAT& x, const FLOAT& y, const FLOAT& z) const {
			return (*this)(VecDi(x,y,z));
		}

		/**
		 * Get/set grid values.
		 * @param pos
		 * @return
		 */
		T& operator() (const VecDi& pos) {
			return this->data()(this->index(pos));
		}


		/**
		 * Get grid values.
		 * @param pos
		 * @return
		 */
		const T& operator() (const VecDi& pos) const
		{
			return this->data()(this->index(pos));
		}



		/**
		 * Get index of position.
		 * @param pos
		 * @return
		 */
		UINT index (const VecDi& pos) const
		{
			return Grid<T,D>::index(pos, this->dims(), this->offset());
		}


		/**
		 * Get index of position.
		 * @param pos
		 * @param dims
		 * @param offset
		 * @return
		 */
		static UINT index (const VecDi& pos, const VecDu& dims, const VecDi& offset = VecDi::Constant(0))
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
		 * @param idx
		 * @return
		 */
		VecDi index (UINT idx) const
		{
			return Grid<T,D>::index(idx, this->dims(), this->offset());
		}

		/**
		 * Get position of index.
		 * @param idx
		 * @param dims
		 * @param offset
		 * @return
		 */
		static VecDi index (UINT idx, const VecDu& dims, const VecDi& offset = VecDi::Zero())
		{
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

			for (INT axis = dims.size()-1; axis >= 0 ; axis--)
			{
				pos(axis) = idx % dims(axis) + offset(axis);
				idx /= dims(axis);
			}

			return pos;
		}


		/**
		 * Get interpolated grid value.
		 * @param pos
		 * @return
		 */
		const T operator() (const VecDf& pos) const {
			return this->interp(pos);
		}

		/**
		 * Get interpolated grid value.
		 * @param pos
		 * @return
		 */
		T operator() (const VecDf& pos) {
			return this->interp(pos);
		}

		/**
		 * Retrieve a reference to the data stored in grid.
		 * @return
		 */
		ArrayData& data () {
			return m_vec_Data;
		}
		/**
		 * Retrieve a reference to the data stored in grid.
		 * @return
		 */
		const ArrayData& data () const {
			return m_vec_Data;
		}


		/**
		 * Reshape grid.
		 *
		 * @param vec_NewDims
		 * @return
		 */
		const VecDu dims (const VecDu& dims_new) {
			VecDu dims_old = m_vec_dims;
			m_vec_dims = dims_new;

			INT uGridSize = m_vec_dims(0);
			for (INT i = 1; i < m_vec_dims.size(); i++) {
				uGridSize *= m_vec_dims(i);
			}
			m_vec_Data.resize(uGridSize);

			// Return old dimensions.
			return dims_old;
		}



		/**
		 * Get grid dimensions.
		 *
		 * @return
		 */
		const VecDu& dims () const {
			return m_vec_dims;
		}


		/**
		 * Fill with a single value.
		 *
		 * @param val
		 */
		void fill (const T& val) {
			this->data() = ArrayData::Constant(1, this->data().size(), val);
//			for (INT i = 0; i < this->data().size(); i++) {
//				this->data()(i) = val;
//			}
		}


		/**
		 * Inside/outside test.
		 *
		 * @param pos
		 * @return
		 */
		bool inside (const VecDi& pos) const {
			const VecDu& dims = this->dims();
			const VecDi& offset = this->offset();

			for (INT i = 0; i < pos.size(); i++) {
				if (pos(i) >= (INT)dims(i) + offset(i))
					return false;
				if (pos(i) < offset(i))
					return false;
			}
			return true;
		}

		void neighs (const VecDi& pos, std::vector<VecDi>& vout, bool bcheck = false) const
		{
			// Reference to grid dimensions.
			const VecDu& dims = this->dims();
			// Most likely all 6 neighbours are valid.
			vout.reserve(6);
			// Position for look-around.
			VecDi vec_dir(pos);
			for (INT axis = 0; axis < dims.size(); axis++) {
				// Check if backward value is within grid.
				vec_dir(axis) -= 1;
				if (this->inside(vec_dir) && (!bcheck || std::find(vout.begin(), vout.end(), vec_dir) == vout.end()))
				{
					vout.push_back(vec_dir);
				}
				// Check if forward value is within grid.
				vec_dir(axis) += 2;
				if (this->inside(vec_dir) && (!bcheck || std::find(vout.begin(), vout.end(), vec_dir) == vout.end()))
				{
					vout.push_back(vec_dir);
				}
				vec_dir(axis) -= 1;
			}
		}

//		void neighs (const VecDi& pos, std::unordered_set<VecDi, UINT (*) (const VecDi& a)>& vout) const
//		{
//			// Reference to grid dimensions.
//			const VecDu& dims = this->dims();
//			// Position for look-around.
//			VecDi vec_dir(pos);
//			for (INT axis = 0; axis < dims.size(); axis++) {
//				// Check if backward value is within grid.
//				vec_dir(axis) -= 1;
//				if (this->inside(vec_dir))
//				{
//					vout.insert(vec_dir);
//				}
//				// Check if forward value is within grid.
//				vec_dir(axis) += 2;
//				if (this->inside(vec_dir))
//				{
//					vout.insert(vec_dir);
//				}
//				vec_dir(axis) -= 1;
//			}
//		}


		void neighs (const VecDi& pos, std::vector<VecDi>& vout, Grid<bool,D>& grid_check) const
		{
			// Reference to grid dimensions.
			const VecDu& dims = this->dims();
			// Most likely all 6 neighbours are valid.
			vout.reserve(6);
			// Position for look-around.
			VecDi vec_dir(pos);
			for (INT axis = 0; axis < dims.size(); axis++) {
				// Check if backward value is within grid.
				vec_dir(axis) -= 1;
				if (this->inside(vec_dir) && !grid_check(vec_dir))
				{
					vout.push_back(vec_dir);
					grid_check(vec_dir) = true;
				}
				// Check if forward value is within grid.
				vec_dir(axis) += 2;
				if (this->inside(vec_dir) && !grid_check(vec_dir))
				{
					vout.push_back(vec_dir);
					grid_check(vec_dir) = true;
				}
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
			// Reference to grid dimensions.
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
			// Reference to grid dimensions.
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
			// Reference to grid dimensions.
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
		 * Will calculate central, forward or backward difference along each axis,
		 * depending what grid values are available.
		 * That is, for grid points at the edge of the grid it will return forward/backward
		 * differences.
		 * @param pos
		 * @return
		 */
		VecDT grad (const VecDi& pos) const
		{
			// Reference to grid dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-around.
			VecDi vec_dir(pos);

			// Central value.
			const T centre = (*this)(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				T back = centre;
				T forward = centre;
				UINT order = 0;
				// Check if backward value is within grid.
				vec_dir(axis) -= 1;
				if (this->inside(vec_dir)) {
					back = (*this)(vec_dir);
					order++;
				}
				// Check if forward value is within grid.
				vec_dir(axis) += 2;
				if (this->inside(vec_dir)) {
					forward = (*this)(vec_dir);
					order++;
				}
				vec_dir(axis) -= 1;
				// Calculate central/forward/backward difference along this axis.
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
			// Reference to grid dimensions.
			const VecDu& dims = this->dims();
			// Vector to store gradient calculation.
			VecDT vec_grad(dims.size());
			// Position for look-around.
			VecDp vec_dir(pos);

			for (INT axis = 0; axis < dims.size(); axis++) {
				vec_dir(axis) -= 1;
				T back = (*this)(vec_dir);
				vec_dir(axis) += 2;
				T forward = (*this)(vec_dir);
				vec_dir(axis) -= 1;

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
		T interp (const VecDf& vec_fpos) const {
			const VecDu& dims = this->dims();

			// Store all 2^d corners.
			std::vector< T > val_corners(1 << dims.size());

			// Get all corners of containing cell.
			for (UINT i = 0; i < val_corners.size(); i++) {
				// 0 = 00 => (x,y)
				// 1 = 01 => (x+1,y)
				// 2 = 10 => (x,y+1)
				// 3 = 11 => (x+1,y+1)

				VecDi pos_corner(dims.size());
				for (INT dim = 0; dim < pos_corner.size(); dim++) {
					INT pos = (INT)floor(vec_fpos(dim));
					const INT dir = (i >> dim) & 1;
					pos += dir;
					pos_corner(dim) = pos;
				}

				val_corners[i] = (*this)(pos_corner);
			}

			// Translate position vector into 'hypercube space', so 0 <= v(x) <= 1.
			VecDf vec_dir(vec_fpos);
			for (INT dim = 0; dim < vec_dir.size(); dim++) {
				vec_dir(dim) = vec_dir(dim) - floor(vec_dir(dim));
			}

			// Repeatedly reduce along axes, i.e. hypercube -> cube -> square -> line -> point
			while (val_corners.size() > 1) {
				val_corners = this->_interp(val_corners, vec_dir);
			}

			return val_corners[0];
		}

		/**
		 * Curvature calculation based on difference of normals method.
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		T curv (const Eigen::Matrix<PosType, D, 1>& pos) const {
			typedef Eigen::Matrix<PosType, D, 1> VecDp;

			const T val_centre = (*this)(pos);
			const VecDu& dims = this->dims();
			VecDp vec_dir(pos);

			// Forward directed principal normal.
			VecDT n_forward(D);

			for (INT axis = 0; axis < dims.size(); axis++) {
				vec_dir(axis) += 1;

				const T val_axis = (*this)(vec_dir) - val_centre;
				T val_neighs_sq = 0;

				// Loop other dimensions to get central difference across them.
				for (INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++) {
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

			for (INT axis = 0; axis < dims.size(); axis++) {
				vec_dir(axis) -= 1;
				VecDT vec_vals(D);

				const T val_axis = val_centre - (*this)(vec_dir);
				T val_neighs_sq = 0;

				// Loop other dimensions to get central difference across them.
				for (INT axis_neigh = 0; axis_neigh < dims.size(); axis_neigh++) {
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
			for (INT axis = 0; axis < dims.size(); axis++) {
				curvature += dn_by_dx(axis);
			}
			curvature /= 2;

			return curvature;
		}


		/**
		 * Calculate divergence, i.e. d2f(x)/dx2
		 * @param pos
		 * @return
		 */
		template <typename PosType>
		T divergence (const Eigen::Matrix<PosType, D, 1>& pos) const {
			const VecDT vec_grad_f = this->gradF(pos);
			const VecDT vec_grad_b = this->gradB(pos);
			const VecDT vec_grad_diff = vec_grad_b - vec_grad_f;

			// Component-wise sum.
			const T val = vec_grad_diff.sum();

			return val / (this->dx()*this->dx());
		}


	#ifndef _TESTING
	protected:
	#endif

		std::vector<T> _interp (const std::vector<T>& val_corners_in, const VecDf& vec_fpos) const {
			const size_t num_corners = val_corners_in.size();

			// Number of values returned.
			// This is a power of 2 less than input dimensions (cube becomes square,
			// square becomes line, line becomes point).
			const size_t num_out = num_corners >> 1;

			// The axis along which to interpolate.
			// This is computed from the dimensions of the original input and the dimensions
			// of the intended output.
			const size_t axis_idx = vec_fpos.size() - log2(num_corners);

			// The weighting to be used in interpolating each pair of points.
			// This is the position along the axis of interpolation.
			const FLOAT axis_pos = vec_fpos(axis_idx);

			std::vector<T> val_corners_out(num_out);

			for (size_t i = 0; i < num_out; i++) {
				const T low = val_corners_in[(i << 1)];
				const T high = val_corners_in[(i << 1) + 1];
				const T val = axis_pos*high + (1.0f-axis_pos)*low;
				val_corners_out[i] = val;
			}

			return val_corners_out;
		}

	};

}// End namespace felt.
#endif
