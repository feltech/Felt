#include <boost/numeric/ublas/io.hpp>
#include <boost/test/unit_test.hpp>
#include <omp.h>

#define _TESTING

#include "Felt.hpp"
#include "Grid.hpp"

using namespace felt;

/*
 * Test the Grid library.
 */
BOOST_AUTO_TEST_SUITE(test_Grid)

	/*
	 * Initialsing grid dimensions.
	 */
	BOOST_AUTO_TEST_CASE(dimensions)
	{
		Vec3u vec_Dims(3, 7, 11);

		Grid <FLOAT,3> grid(vec_Dims);

		Vec3u vec_DimsRet = grid.dims();

		BOOST_CHECK_EQUAL((size_t)vec_DimsRet(0), (size_t)3);
		BOOST_CHECK_EQUAL((size_t)vec_DimsRet(1), (size_t)7);
		BOOST_CHECK_EQUAL((size_t)vec_DimsRet(2), (size_t)11);

		BOOST_CHECK_EQUAL(grid.data().size(), (size_t)(3*7*11));

		grid.dims(Vec3u(5,11,13));
		BOOST_CHECK_EQUAL(grid.data().size(), (size_t)(5*11*13));
	}

	/*
	 * Getting/setting grid values.
	 */
	BOOST_AUTO_TEST_CASE(get_and_set)
	{
		Vec3u vec_Dims(3, 7, 11);

		Grid <FLOAT,3> grid(Vec3u(3, 7, 11));

		grid(Vec3i(0,0,0)) = 13.0f;
		grid(Vec3i(1,2,3)) = 17.0f;
		grid(Vec3i(2,6,10)) = 19.0f;

		BOOST_CHECK_EQUAL((FLOAT)grid.data()(0), 13.0f);
		BOOST_CHECK_EQUAL((FLOAT)grid(Vec3i(1,2,3)), 17.0f);
		BOOST_CHECK_EQUAL(grid.data()(grid.data().size()-1), (float)19.0f);

		const float test_const = grid(Vec3i(1,2,3));
		BOOST_CHECK_EQUAL(test_const, 17.0f);

		// Shortcuts
		{
			// 2D.
			Grid <FLOAT,2> grid(3, 7);
			grid.fill(2);
			BOOST_CHECK_EQUAL(grid(1,1), (float)2.0f);
		}
		{
			// 3D.
			Grid <FLOAT,3> grid(3, 7, 11);
			grid.fill(2);
			BOOST_CHECK_EQUAL(grid(1,1,1), (float)2.0f);
		}
	}

	/*
	 * Getting grid point indices.
	 */
	BOOST_AUTO_TEST_CASE(get_indices)
	{
		Grid <FLOAT,2> grid(Vec2u(3, 4), Vec2i(-1,-1));
/*
	Row major order: (x,y) =>
	(-1,-1),(-1,0),	(-1,1),	(-1,2)
	(0,-1),	(0,0),	(0,1),	(0,2)
	(1,-1),	(1,0),	(1,1),	(1,2)

*/
		BOOST_CHECK_EQUAL(grid.index(Vec2i(1,0)), 9);
		Vec2i pos = grid.index(7);
		BOOST_CHECK_EQUAL((INT)pos(0), 0);
		BOOST_CHECK_EQUAL((INT)pos(1), 2);
	}


	/*
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

	/*
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

	/*
	 * Offsetting the grid.
	 */
	BOOST_AUTO_TEST_CASE(offsetting)
	{
		Grid <FLOAT,3> grid(Vec3u(7, 11, 13), Vec3i(-3,-3,-3));

		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-2,0,0)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-4,0,0)), false);

		grid.offset(Vec3i(-1,-1,-1));
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-1,0,0)), true);
		BOOST_CHECK_EQUAL(grid.inside(Vec3i(-2,0,0)), false);

		grid(Vec3i(-1,-1,-1)) = 21.0f;
		grid(Vec3i(-1,0,-1)) = 23.0f;

		BOOST_CHECK_EQUAL((FLOAT)grid.data()(0), 21.0f);
		BOOST_CHECK_EQUAL(grid(Vec3i(-1,0,-1)), 23.0f);
	}


	/*
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


	/*
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


	/*
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
/*
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


	/*
	 * Divergence.
	 */
	BOOST_AUTO_TEST_CASE(divergence__d2f_by_dx2)
	{
		Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 2);

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

		const FLOAT d2f_by_dx2__negative = grid.divergence(Vec3i(0,0,0));
		BOOST_CHECK_CLOSE(d2f_by_dx2__negative, -0.75f, 0.00001f);
	}

	/*
	 * Using delta x.
	 */
	BOOST_AUTO_TEST_CASE(delta_x)
	{
		using namespace boost::numeric::ublas;

		// Default and explicitly setting.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
			BOOST_CHECK_EQUAL(grid.dx(), 1.0f);

			grid.dx(2.0f);
			BOOST_CHECK_EQUAL(grid.dx(), 2.0f);
		}

		// Setting at construction.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 3.0f);
			BOOST_CHECK_EQUAL(grid.dx(), 3.0f);
		}

		// Spatial derivatives with dx != 1.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 2.0f);
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

	/*
	 * Interpolation.
	 */
	BOOST_AUTO_TEST_CASE(interpolate1D)
	{
		Grid <FLOAT,1> grid = Grid <FLOAT,1>();
/*
		0----1
*/
		{
			std::vector<FLOAT> input = std::vector<FLOAT>(2);
			input[0] = 1.0f;
			input[1] = 0;

			Eigen::Matrix<FLOAT, 1, 1> pos;
			pos(0) = 0.3;

			std::vector< FLOAT > out = grid._interp(input, pos);

			BOOST_CHECK_EQUAL(out[0], 0.7f);
		}
	}
	BOOST_AUTO_TEST_CASE(interpolate2D)
	{
		Grid <FLOAT,2> grid = Grid <FLOAT,2>();
/*
		10----11
		|	   |
		|	   |
		00----01
*/
		{
			std::vector<FLOAT> input = std::vector<FLOAT>(4);
			input[0 /*00*/] = 2.0f;
			input[1 /*01*/] = 0;
			input[2 /*10*/] = 0.0f;
			input[3 /*11*/] = 1.0;

			Vec2f pos(0.8f, 0.5f);

			std::vector< FLOAT > output = grid._interp(input, pos);

			BOOST_CHECK_CLOSE(output[0], 0.4f, 0.00001f);
			BOOST_CHECK_CLOSE(output[1], 0.8f, 0.00001f);

			input = output;

			output = grid._interp(input, pos);
			BOOST_CHECK_CLOSE(output[0], 0.6f, 0.00001f);
		}
	}

	BOOST_AUTO_TEST_CASE(interpolate3D)
	{
/*
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
			input[0 /*000*/] = 0.0f;
			input[1 /*001*/] = 0.8f;
			input[2 /*010*/] = 1.0f;
			input[3 /*011*/] = 1.0f;
			input[4 /*100*/] = 0.0f;
			input[5 /*101*/] = 0.0f;
			input[6 /*110*/] = 1.0f;
			input[7 /*111*/] = 1.0f;

			Vec3f pos(0.5f, 0.75f, 0.5f);

			std::vector< FLOAT > output4 = grid._interp(input, pos);

			BOOST_CHECK_EQUAL(output4[0 /*00x*/], 0.4f);
			BOOST_CHECK_EQUAL(output4[1 /*01x*/], 1.0f);
			BOOST_CHECK_EQUAL(output4[2 /*10x*/], 0.0f);
			BOOST_CHECK_EQUAL(output4[3 /*11x*/], 1.0f);

			input = output4;
			std::vector< FLOAT > output2 = grid._interp(input, pos);
			BOOST_CHECK_CLOSE(output2[0 /*0yx*/], 0.85f, 0.00001f);
			BOOST_CHECK_CLOSE(output2[1 /*1yx*/], 0.75f, 0.00001f);

			input = output2;
			std::vector< FLOAT > output1 = grid._interp(input, pos);
			BOOST_CHECK_CLOSE(output1[0 /*zyx*/], 0.8f, 0.00001f);
		}
	}

	BOOST_AUTO_TEST_CASE(interpolate_grid_2D)
	{
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
		grid.fill(0);
		grid(Vec2i(-1,-1)) = 1.0f;
		grid(Vec2i(-1,0)) = 1.0f;
		grid(Vec2i(0,1)) = 2.0f;
		grid(Vec2i(1,1)) = 2.0f;

		Vec2f pos = Vec2f(0.0f, 0.0f);

		// Explicit call.
		FLOAT val = grid.interp(pos);

		BOOST_CHECK_CLOSE(val, 0.0f, 0.00001f);

		BOOST_CHECK_CLOSE(grid.interp(Vec2f(-0.5f, -0.5f)), 0.5f, 0.00001f);


		BOOST_CHECK_CLOSE(grid.interp(Vec2f(0.5f, 0.5f)), 1.0f, 0.00001f);

		pos = Vec2f(0.5f, 0.5f);
		// Via operator.
		val = grid(pos);
		BOOST_CHECK_CLOSE(val, 1.0f, 0.00001f);
	}


	/*
	 * Gradient interpolation.
	 */
	BOOST_AUTO_TEST_CASE(grad_forward_interp)
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


	/*
	 * Entropy satisfying gradient.
	 */
	BOOST_AUTO_TEST_CASE(grad_entropy_satisfying)
	{
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1));
		grid.fill(0);
		grid(Vec2i(0,0)) = 1.0f;

		grid(Vec2i(-1,0)) = 0.0f;
		grid(Vec2i(1,0)) = 2.0f;

		grid(Vec2i(0,-1)) = 0.5f;
		grid(Vec2i(0,1)) = 0.7f;


		Vec2i pos(0, 0);
		Vec2f vec_grad_y_entropy = grid.gradE(pos);

		BOOST_CHECK_CLOSE((FLOAT)vec_grad_y_entropy(0), 1.0f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad_y_entropy(1), 0.0f, 0.00001f);

		grid(Vec2i(-1,0)) = 0.9f;
		grid(Vec2i(1,0)) = 0.8f;

		grid(Vec2i(0,-1)) = 2.0f;
		grid(Vec2i(0,1)) = 0.0f;

		Vec2f vec_grad_x_entropy = grid.gradE(pos);

		BOOST_CHECK_CLOSE((FLOAT)vec_grad_x_entropy(0), 0, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad_x_entropy(1), -1.0f, 0.00001f);


/*
	Row major order: (x,y) =>
	(-1,-1),(-1,0),	(-1,1)
	(0,-1),	(0,0),	(0,1)
	(1,-1),	(1,0),	(1,1)

*/
		grid.data() <<	0,	1,		0,
						2,	1,		0,
						0,	1.5,	0;

		Vec2f vec_grad_no_entropy = grid.gradE(pos);

		BOOST_CHECK_CLOSE((FLOAT)vec_grad_no_entropy(0), 0.5f, 0.00001f);
		BOOST_CHECK_CLOSE((FLOAT)vec_grad_no_entropy(1), -1.0f, 0.00001f);
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
