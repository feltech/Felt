#include <boost/numeric/ublas/io.hpp>
#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "Felt/Grid.hpp"
#include "Utils.hpp"

using namespace felt;

/**
 * Test the Grid class.
 */
BOOST_AUTO_TEST_SUITE(test_Grid)

	/**
	 * Initialsing grid size.
	 */
	BOOST_AUTO_TEST_CASE(dimensions)
	{
		//! [Initialsing grid size]
		// ==== Setup ====
		Grid <FLOAT, 3> grid(Vec3u(3, 7, 11));
		const Vec3u& size = grid.size();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL((UINT)size(0), 3);
		BOOST_CHECK_EQUAL((UINT)size(1), 7);
		BOOST_CHECK_EQUAL((UINT)size(2), 11);
		BOOST_CHECK_EQUAL(grid.data().size(), (3*7*11));

		// ==== Action ====
		grid.size(Vec3u(5, 11, 13));

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.data().size(), (5*11*13));
		//! [Initialsing grid size]
	}

	/**
	 * Getting/setting grid values.
	 */
	BOOST_AUTO_TEST_CASE(get_and_set)
	{
		//! [Get and set]
		// ==== Setup ====
		Grid<FLOAT, 3> grid(Vec3u(3, 7, 11));

		// ==== Action ====
		grid(Vec3i(0,0,0)) = 13.0f;
		grid(Vec3i(1,2,3)) = 17.0f;
		grid(Vec3i(2,6,10)) = 19.0f;

		// ==== Confirm ====
		const float check_get = grid.get(Vec3i(1,2,3));
		const float check_op = grid(Vec3i(1,2,3));
		const float check_data_start = grid.data()(0);
		const float check_data_end = grid.data()(grid.data().size() - 1);
		BOOST_CHECK_EQUAL(check_get, 17.0f);
		BOOST_CHECK_EQUAL(check_op, 17.0f);
		BOOST_CHECK_EQUAL(check_data_start, 13.0f);
		BOOST_CHECK_EQUAL(check_data_end, 19.0f);
		//! [Get and set]
	}

	/**
	 * Getting grid point indices.
	 */
	BOOST_AUTO_TEST_CASE(position_index)
	{
		//! [Position index]
		// ==== Setup ====

		using GridType = Grid<FLOAT, 2>;

		const Vec2u size(3, 4);
		const Vec2i offset(-1, -1);

		GridType grid(size, offset);

		/*
			Row major order: (x,y) => [
				(-1,-1),(-1,0),	(-1,1),	(-1,2)
				(0,-1),	(0,0),	(0,1),	(0,2)
				(1,-1),	(1,0),	(1,1),	(1,2)
			]
		*/

		// ==== Action ====

		const UINT index_from_pos_static = GridType::index(Vec2i(1, 0), size, offset);
		const Vec2i& pos_from_index_static = GridType::index(7, size, offset);

		const UINT index_from_pos = grid.index(Vec2i(1, 0));
		const Vec2i& pos_from_index = grid.index(7);

		// === Confirm ===

		BOOST_CHECK_EQUAL(index_from_pos_static, 9);
		BOOST_CHECK_EQUAL(pos_from_index_static, Vec2i(0, 2));

		BOOST_CHECK_EQUAL(index_from_pos, 9);
		BOOST_CHECK_EQUAL(pos_from_index, Vec2i(0, 2));
		//! [Position index]
	}


	/**
	 * Fill grid with a value.
	 */
	BOOST_AUTO_TEST_CASE(filling)
	{
		Grid <UINT,3> grid(Vec3u(3, 7, 11));

		grid.fill(7);

		INT sum = 0;
		for (INT i = 0; i < grid.data().size(); i++) {
			sum += grid.data()(i);
		}

		BOOST_CHECK_EQUAL(sum, 3*7*11 * 7);
	}

	/**
	 * Check position is within grid.
	 */
	BOOST_AUTO_TEST_CASE(inside_outside_check)
	{
		Grid <FLOAT,3> grid(Vec3u(3, 7, 11));

		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-1,0,0)), false);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(0,0,0)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(1,2,3)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(3,7,11)), false);
	}

	/**
	 * Offsetting the grid.
	 */
	BOOST_AUTO_TEST_CASE(offsetting)
	{
		//! [Offsetting the grid]
		// ==== Setup ====
		Grid <FLOAT,3> grid(Vec3u(7, 11, 13), Vec3i(-3,-3,-3));

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-2,0,0)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-4,0,0)), false);

		// ==== Action ====
		grid.offset(Vec3i(-1,-1,-1));

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-1,0,0)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-2,0,0)), false);

		// ==== Action ====
		grid(Vec3i(-1,-1,-1)) = 21.0f;
		grid(Vec3i(-1,0,-1)) = 23.0f;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL((FLOAT)grid.data()(0), 21.0f);
		BOOST_CHECK_EQUAL(grid(Vec3i(-1,0,-1)), 23.0f);
		//! [Offsetting the grid]
	}


	/**
	 * Forward differencing.
	 */
	BOOST_AUTO_TEST_CASE(grad_forward)
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradF(Vec2i(0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				const Vec2f vec_grad = grid.gradF(Vec2i(0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -1.0f);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), -1.0f);
			}
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));

			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradF(Vec3i(0,0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), -1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 1);
			}
		}
	}


	/**
	 * Backward differencing.
	 */
	BOOST_AUTO_TEST_CASE(grad_backward)
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradB(Vec2i(0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				{
					const Vec2f vec_grad = grid.gradB(Vec2i(0,0));

					BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 1.0f);
					BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 1.0f);
				}
			}
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));

			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradB(Vec3i(0,0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 1);
			}
		}
	}


	/**
	 * Central differencing.
	 */
	BOOST_AUTO_TEST_CASE(grad_central)
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradC(Vec2i(0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				{
					const Vec2f vec_grad = grid.gradC(Vec2i(0,0));

					BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 0.0f);
					BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0.0f);
				}
			}
		}


		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
/**
	Row major order: (x,y,z) =>
	(-1,-1,-1),	(-1,-1, 0),	(-1,-1, 1)
	(-1, 0,-1),	(-1, 0, 0),	(-1, 0, 1)
	(-1, 1,-1),	(-1, 1, 0),	(-1, 1, 1)

	(0,-1,-1),	(0,-1, 0),	(0,-1, 1)
	(0, 0,-1),	(0, 0, 0),	(0, 0, 1)
	(0, 1,-1),	(0, 1, 0),	(0, 1, 1)

	(1,-1,-1),	(1,-1, 0),	(1,-1, 1)
	(1, 0,-1),	(1, 0, 0),	(1, 0, 1)
	(1, 1,-1),	(1, 1, 0),	(1, 1, 1)
*/
			grid.data() <<
				0,	0,	0,
				0,	2,	0,
				0,	0,	0,

				0,  0,	0,
				0,	1,	2,
				0,	0,	0,

				0,	0,	0,
				0,	0,	0,
				0,	0,	0;
//			grid.fill(0);
//			grid(Vec3i(0,0,0)) = 1.0f;
//			grid(Vec3i(0,0,1)) = 2.0f;
//			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradC(Vec3i(0,0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 1);
			}

			// Safe differencing - checks boundaries and uses
			// central/forward/back as appropriate.
			{
				Vec3f vec_grad = grid.grad(Vec3i(0,0,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 1);

				vec_grad = grid.grad(Vec3i(0,-1,0));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 1);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 0);

				vec_grad = grid.grad(Vec3i(1,0,1));

				BOOST_CHECK_EQUAL((FLOAT)vec_grad(0), -2);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(1), 0);
				BOOST_CHECK_EQUAL((FLOAT)vec_grad(2), 0);

			}
		}
	}


	/**
	 * Divergence.
	 */
	BOOST_AUTO_TEST_CASE(divergence_d2f_by_dx2)
	{
		//! [Divergence] Divergence
		// ==== Setup ====
		Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
		grid.dx(2);
		grid.data() <<	1,	1,	1,
						1,	1,	1,
						1,	1,	1,

						1,	1,	1,
						1,	0,	1,
						1,	1,	1,

						1,	1,	1,
						1,	1,	1,
						1,	1,	1;

		// ==== Action ====
		const FLOAT d2f_by_dx2 = grid.divergence(Vec3i(0,0,0));

		// ==== Confirm ====
		BOOST_CHECK_CLOSE(d2f_by_dx2, -0.75f, 0.00001f);
		//! [Divergence]
	}

	/**
	 * Using delta x.
	 */
	BOOST_AUTO_TEST_CASE(delta_x)
	{
		using namespace boost::numeric::ublas;

		// Default and explicitly setting.
		{
			//! [Delta x setter]
			// ==== Setup ====
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(grid.dx(), 1.0f);

			// ==== Action ====
			grid.dx(2.0f);

			// ==== Confirm ====
			BOOST_CHECK_EQUAL(grid.dx(), 2.0f);
			//! [Delta x setter]
		}

		// Spatial derivatives with dx != 1.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
			grid.dx(2.0f);
			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			BOOST_CHECK_EQUAL((FLOAT)grid.gradF(Vec3i(0,0,0))(0), -0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradF(Vec3i(0,0,0))(1), -0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradF(Vec3i(0,0,0))(2), 0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradB(Vec3i(0,0,0))(0), -0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradB(Vec3i(0,0,0))(1), 0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradB(Vec3i(0,0,0))(2), 0.5);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradC(Vec3i(0,0,0))(0), -0.5f);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradC(Vec3i(0,0,0))(1), 0);
			BOOST_CHECK_EQUAL((FLOAT)grid.gradC(Vec3i(0,0,0))(2), 0.5f);

		}
	}

	/**
	 * Interpolation.
	 */
	BOOST_AUTO_TEST_CASE(interpolate1D)
	{
		Grid <FLOAT,1> grid = Grid <FLOAT,1>();
/**
		0----1
*/
		{
			std::vector<FLOAT> input = std::vector<FLOAT>(2);
			input[0] = 1.0f;
			input[1] = 0;

			Eigen::Matrix<FLOAT, 1, 1> pos;
			pos(0) = 0.3;

			std::vector< FLOAT > out = grid.interp(input, pos);

			BOOST_CHECK_EQUAL(out[0], 0.7f);
		}
	}
	BOOST_AUTO_TEST_CASE(interpolate2D)
	{
		Grid <FLOAT,2> grid = Grid <FLOAT,2>();
/**
		10----11
		|	   |
		|	   |
		00----01
*/
		{
			std::vector<FLOAT> input = std::vector<FLOAT>(4);
			input[0 /**00*/] = 2.0f;
			input[1 /**01*/] = 0;
			input[2 /**10*/] = 0.0f;
			input[3 /**11*/] = 1.0;

			Vec2f pos(0.8f, 0.5f);

			std::vector< FLOAT > output = grid.interp(input, pos);

			BOOST_CHECK_CLOSE(output[0], 0.4f, 0.00001f);
			BOOST_CHECK_CLOSE(output[1], 0.8f, 0.00001f);

			input = output;

			output = grid.interp(input, pos);
			BOOST_CHECK_CLOSE(output[0], 0.6f, 0.00001f);
		}
	}

	BOOST_AUTO_TEST_CASE(interpolate3D)
	{
/**
		  011----111
		 /|		  /|
		010----011 |
		| 100----|101
		|/		 |/
		000----001
*/
		Grid <FLOAT,3> grid = Grid <FLOAT,3>();

		{
			std::vector<FLOAT> input = std::vector<FLOAT>(8);
			input[0 /**000*/] = 0.0f;
			input[1 /**001*/] = 0.8f;
			input[2 /**010*/] = 1.0f;
			input[3 /**011*/] = 1.0f;
			input[4 /**100*/] = 0.0f;
			input[5 /**101*/] = 0.0f;
			input[6 /**110*/] = 1.0f;
			input[7 /**111*/] = 1.0f;

			Vec3f pos(0.5f, 0.75f, 0.5f);

			std::vector< FLOAT > output4 = grid.interp(input, pos);

			BOOST_CHECK_EQUAL(output4[0 /**00x*/], 0.4f);
			BOOST_CHECK_EQUAL(output4[1 /**01x*/], 1.0f);
			BOOST_CHECK_EQUAL(output4[2 /**10x*/], 0.0f);
			BOOST_CHECK_EQUAL(output4[3 /**11x*/], 1.0f);

			input = output4;
			std::vector< FLOAT > output2 = grid.interp(input, pos);
			BOOST_CHECK_CLOSE(output2[0 /**0yx*/], 0.85f, 0.00001f);
			BOOST_CHECK_CLOSE(output2[1 /**1yx*/], 0.75f, 0.00001f);

			input = output2;
			std::vector< FLOAT > output1 = grid.interp(input, pos);
			BOOST_CHECK_CLOSE(output1[0 /**zyx*/], 0.8f, 0.00001f);
		}
	}

	BOOST_AUTO_TEST_CASE(interpolate_grid_2D)
	{
		//! [Interpolation]
		// ==== Setup ====
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
		grid.fill(0);
		grid(Vec2i(-1,-1)) = 1.0f;
		grid(Vec2i(-1,0)) = 1.0f;
		grid(Vec2i(0,1)) = 2.0f;
		grid(Vec2i(1,1)) = 2.0f;

		// ==== Action ====
		// Via explicit function call.
		const FLOAT val1 = grid.interp(Vec2f(0.0f, 0.0f));
		const FLOAT val2 = grid.interp(Vec2f(-0.5f, -0.5f));
		const FLOAT val3 = grid.interp(Vec2f(0.5f, 0.5f));
		// Via getter function overload.
		const FLOAT val4 = grid.val(Vec2f(0.5f, 0.5f));
		// Via operator overload.
		const FLOAT val5 = grid(Vec2f(0.5f, 0.5f));

		// ==== Confirm ====
		BOOST_CHECK_CLOSE(val1, 0.0f, 0.00001f);
		BOOST_CHECK_CLOSE(val2, 0.5f, 0.00001f);
		BOOST_CHECK_CLOSE(val3, 1.0f, 0.00001f);
		BOOST_CHECK_CLOSE(val4, 1.0f, 0.00001f);
		BOOST_CHECK_CLOSE(val5, 1.0f, 0.00001f);
		//! [Interpolation]
	}


	/**
	 * Gradient interpolation.
	 */
	BOOST_AUTO_TEST_CASE(grad_forwardinterp)
	{
		Grid <FLOAT,2> grid(Vec2u(5, 5), Vec2i(-2,-2));
		grid.fill(0);
		grid(Vec2i(-1,-1)) = 1.0f;
		grid(Vec2i(-1,0)) = 1.0f;
		grid(Vec2i(0,1)) = 2.0f;
		grid(Vec2i(1,1)) = 2.0f;

		Vec2f pos = Vec2f(0.0f, 0.0f);
		Vec2f vec_grad = grid.gradF(pos);

		BOOST_CHECK_CLOSE((FLOAT)vec_grad(0), 0.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad(1), 2.0f, 0.00001f);
	}


	/**
	 * Entropy satisfying gradient.
	 */
	BOOST_AUTO_TEST_CASE(grad_entropy_satisfying)
	{
		// ==== Setup ====
		/**
			Row major order: (x,y) =>
			(-1,-1),(-1,0),	(-1,1)
			(0,-1),	(0,0),	(0,1)
			(1,-1),	(1,0),	(1,1)

		*/
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
		Vec2i pos(0, 0);
		Vec2f vec_grad;

		// grad +'ve
		grid.data() <<	0,	0,	0,
						0,	1,	3,
						0,	3,	0;

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		BOOST_CHECK_CLOSE((FLOAT)vec_grad(0), 1.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad(1), 1.0f, 0.00001f);


		// ==== Setup ====

		// grad -'ve
		grid.data() <<	0,	3,	0,
						3,	1,	0,
						0,	0,	0;

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		BOOST_CHECK_CLOSE((FLOAT)vec_grad(0), -1.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad(1), -1.0f, 0.00001f);


		// ==== Setup ====

		// div -'ve
		grid.data() <<	0,	2,	0,
						3,	1,	2,
						0,	3,	0;

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		BOOST_CHECK_CLOSE((FLOAT)vec_grad(0), 0.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad(1), 0.0f, 0.00001f);

		// ==== Setup ====

		// div +'ve
		grid.data() <<	0,	6,	0,
						6,	9,	1,
						0,	1,	0;

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		BOOST_CHECK_CLOSE((FLOAT)vec_grad(0), -5.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad(1), -5.0f, 0.00001f);
	}

	BOOST_AUTO_TEST_CASE(curvature)
	{
		// 2D.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
			grid.data() <<	1,	1,	1,
							1,	0,	1,
							1,	1,	1;

			const float singularity_curvature = grid.curv(Vec2i(0,0));

			BOOST_CHECK_CLOSE(singularity_curvature, 2, 0.00001f);

			grid.data() <<	 1,	 1,	 1,
							 0,	 0,  1,
							-1,	 0,	 1;

			const float corner_curvature = grid.curv(Vec2i(0,0));

			BOOST_CHECK_CLOSE(corner_curvature, 1, 0.00001f);
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1));
			grid.data() <<
				1,	1,	1,
				1,	1,	1,
				1,	1,	1,

				1,	1,	1,
				1,	0,	1,
				1,	1,	1,

				1,	1,	1,
				1,	1,	1,
				1,	1,	1;

			const float singularity_curvature_3D = grid.curv(Vec3i(0,0,0));

			BOOST_CHECK_CLOSE(singularity_curvature_3D, 3, 0.00001f);

			grid.data() <<
				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1;
			const float corner_curvature_3D = grid.curv(Vec3i(0,0,0));

			BOOST_CHECK_CLOSE(corner_curvature_3D, 1, 0.00001f);


			grid.data() <<
				 1,	 1,	 1,
				 1,	 1,	 1,
				 1,	 1,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				 0,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1;
			const float sharp_corner_curvature_3D = grid.curv(Vec3i(0,0,0));

			BOOST_CHECK_CLOSE(sharp_corner_curvature_3D, 1.5f, 0.00001f);
		}
	}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_LazyGrid)
	BOOST_AUTO_TEST_CASE(initialisation)
	{
		/// [LazyGrid initialisation]
		// ==== Setup ====
		LazyGrid <FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 7.0f);

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), 7.0f);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 7.0f);
		/// [LazyGrid initialisation]
	}

	BOOST_AUTO_TEST_CASE(activation)
	{
		/// [LazyGrid activation]
		// ==== Setup ====
		LazyGrid <FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 7.0f);

		// ==== Action ====
		grid.activate();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), true);
		BOOST_CHECK_EQUAL(grid.data().size(), 3*3*3);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 7.0f);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(0,1,1)), 7.0f);

		// ==== Action ====
		grid.get(Vec3i(1,1,1)) = 5.0f;

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 5.0f);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(0,1,1)), 7.0f);
		/// [LazyGrid activation]
	}

	BOOST_AUTO_TEST_CASE(deactivation)
	{
		/// [LazyGrid deactivation]
		// ==== Setup ====
		LazyGrid <FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 7.0f);

		// ==== Action ====
		grid.get(Vec3i(1,1,1)) = 5.0f;
		grid.background() = 3.0f;
		grid.deactivate();

		// ==== Confirm ====
		BOOST_CHECK_EQUAL(grid.is_active(), false);
		BOOST_CHECK_EQUAL(grid.data().size(), 0);
		BOOST_CHECK_EQUAL(grid.background(), 3.0f);
		BOOST_CHECK_EQUAL(grid.get(Vec3i(1,1,1)), 3.0f);
		/// [LazyGrid deactivation]
	}
BOOST_AUTO_TEST_SUITE_END()
