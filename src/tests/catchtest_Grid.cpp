#include <boost/numeric/ublas/io.hpp>
#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include "Felt/Grid.hpp"
#include "Utils.hpp"

using namespace felt;

/**
 * Test the Grid class.
 */
SCENARIO("Grid")
{
	GIVEN("a 3x7x11 grid with no offset")
	{
		//! [Initialsing grid size]
		Grid<FLOAT, 3> grid(Vec3u(3, 7, 11), Vec3i::Zero(), 0);

		THEN("memory is allocated and the size is reported correctly")
		{
			const Vec3u &size = grid.size();

			CHECK((UINT) size(0) == 3);
			CHECK((UINT) size(1) == 7);
			CHECK((UINT) size(2) == 11);
			CHECK(grid.data().size() == (3 * 7 * 11));
		}

		WHEN("the size is updated")
		{
			grid.size(Vec3u(5, 11, 13));

			THEN("the memory is resized")
			{
				CHECK(grid.data().size() == (5 * 11 * 13));
			}
		}
		//! [Initialsing grid size]

		//! [Get and set]
		WHEN("some positions values are set")
		{
			// ==== Action ====
			grid(Vec3i(0, 0, 0))  = 13.0f;
			grid(Vec3i(1, 2, 3))  = 17.0f;
			grid(Vec3i(2, 6, 10)) = 19.0f;

			THEN("querying those positions returns the same values")
			{
				// ==== Confirm ====
				CHECK(grid.get(Vec3i(1, 2, 3)) == 17.0f);
				CHECK(grid(Vec3i(1, 2, 3)) == 17.0f);
			}

			THEN("expected elements of the underlying array contain those values")
			{
				CHECK(grid.data()[0] == 13.0f);
				CHECK(grid.data()[grid.data().size() - 1] == 19.0f);
			}
		}
		//! [Get and set]
	}

	//! [Position index]
	/**
	 * Getting grid point indices.
	 */
	GIVEN("a 3x4 grid with (-1,-1) offset")
	{
		using GridType = Grid<FLOAT, 2>;
		const Vec2u size(3, 4);
		const Vec2i offset(-1, -1);

		GridType grid(size, offset, 0);

/*
	Row major order: (x,y) => [
		(-1,-1),(-1,0),	(-1,1),	(-1,2)
		(0,-1),	(0,0),	(0,1),	(0,2)
		(1,-1),	(1,0),	(1,1),	(1,2)
	]
*/
		WHEN("we query for a point within the grid")
		{
			THEN("the index of a point in the data array is reported correctly")
			{
				CHECK(GridType::index(Vec2i(1, 0), size, offset) == 9);
				CHECK(grid.index(Vec2i(1, 0)) == 9);
			}

			THEN("the point represented by an index in the data array is reported correctly")
			{
				CHECK(grid.index(7) == Vec2i(0, 2));
				CHECK(GridType::index(7, size, offset) == Vec2i(0, 2));
			}
		}

		WHEN("we query for a point outside the grid (positive)")
		{
			THEN("the index of a point in the data array is reported correctly")
			{
				CHECK(GridType::index(Vec2i(2, -1), size, offset) == 0);
				CHECK(grid.index(Vec2i(2, -1)) == 0);
			}

			THEN("the point represented by an index in the data array is reported correctly")
			{
				CHECK(grid.index(0) == Vec2i(-1, -1));
				CHECK(GridType::index(0, size, offset) == Vec2i(-1, -1));
			}
		}
	}
	//! [Position index]


	GIVEN("a 7x7x7 grid with (-3,-3,-3) offset")
	{
		using GridType = Grid<FLOAT, 3>;
		GridType grid(Vec3u(7,7,7), Vec3i(-3,-3,-3), 0);

		WHEN("we modulo transform a position outside the grid in the positive direction")
		{
			Vec3f pos_mod = grid.mod(Vec3f(5,5,5));

			THEN("the position is transformed correctly")
			{
				CHECK(pos_mod == ApproxVec(Vec3f(-2,-2,-2)));
			}
		}

		WHEN("we modulo transform a position outside the grid in the negative direction")
		{
			Vec3f pos_mod = grid.mod(Vec3f(-5,-5,-5));

			THEN("the position is transformed correctly")
			{
				CHECK(pos_mod == ApproxVec(Vec3f(2,2,2)));
			}
		}
	}

	/**
	 * Fill grid with a value.
	 */
	WHEN("filling")
	{
		Grid <UINT,3> grid(Vec3u(3, 7, 11), Vec3i::Zero(), 0);

		grid.fill(7);

		INT sum = 0;
		for (INT i = 0; i < grid.data().size(); i++) {
			sum += grid.data()[i];
		}

		CHECK(sum == 3*7*11 * 7);
	}

	/**
	 * Check position is within grid.
	 */
	WHEN("inside_outside_check")
	{
		Grid <FLOAT,3> grid(Vec3u(3, 7, 11), Vec3i::Zero(), 0);

		CHECK(grid.inside(Vec3i(-1,0,0)) == false);
		CHECK(grid.inside(Vec3i(0,0,0)) == true);
		CHECK(grid.inside(Vec3i(1,2,3)) == true);
		CHECK(grid.inside(Vec3i(3,7,11)) == false);
	}

	/**
	 * Offsetting the grid.
	 */
	WHEN("offsetting")
	{
		//! [Offsetting the grid]
		// ==== Setup ====
		Grid <FLOAT,3> grid(Vec3u(7, 11, 13), Vec3i(-3,-3,-3), 0);

		// ==== Confirm ====
		CHECK(grid.inside(Vec3i(-2,0,0)) == true);
		CHECK(grid.inside(Vec3i(-4,0,0)) == false);

		// ==== Action ====
		grid.offset(Vec3i(-1,-1,-1));

		// ==== Confirm ====
		CHECK(grid.inside(Vec3i(-1,0,0)) == true);
		CHECK(grid.inside(Vec3i(-2,0,0)) == false);

		// ==== Action ====
		grid(Vec3i(-1,-1,-1)) = 21.0f;
		grid(Vec3i(-1,0,-1)) = 23.0f;

		// ==== Confirm ====
		CHECK((FLOAT)grid.data()[0] == 21.0f);
		CHECK(grid(Vec3i(-1,0,-1)) == 23.0f);
		//! [Offsetting the grid]
	}


	/**
	 * Forward differencing.
	 */
	WHEN("grad_forward")
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradF(Vec2i(0,0));

				CHECK((FLOAT)vec_grad(0) == 0);
				CHECK((FLOAT)vec_grad(1) == 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				const Vec2f vec_grad = grid.gradF(Vec2i(0,0));

				CHECK((FLOAT)vec_grad(0) == -1.0f);
				CHECK((FLOAT)vec_grad(1) == -1.0f);
			}
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);

			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradF(Vec3i(0,0,0));

				CHECK((FLOAT)vec_grad(0) == -1);
				CHECK((FLOAT)vec_grad(1) == -1);
				CHECK((FLOAT)vec_grad(2) == 1);
			}
		}
	}


	/**
	 * Backward differencing.
	 */
	WHEN("grad_backward")
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradB(Vec2i(0,0));

				CHECK((FLOAT)vec_grad(0) == 0);
				CHECK((FLOAT)vec_grad(1) == 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				{
					const Vec2f vec_grad = grid.gradB(Vec2i(0,0));

					CHECK((FLOAT)vec_grad(0) == 1.0f);
					CHECK((FLOAT)vec_grad(1) == 1.0f);
				}
			}
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);

			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradB(Vec3i(0,0,0));

				CHECK((FLOAT)vec_grad(0) == -1);
				CHECK((FLOAT)vec_grad(1) == 1);
				CHECK((FLOAT)vec_grad(2) == 1);
			}
		}
	}


	/**
	 * Central differencing.
	 */
	WHEN("grad_central")
	{
		// Basic 2D testing.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);

			grid.fill(0);

			{
				const Vec2f vec_grad = grid.gradC(Vec2i(0,0));

				CHECK((FLOAT)vec_grad(0) == 0);
				CHECK((FLOAT)vec_grad(1) == 0);
			}

			grid(Vec2i(0,0)) = 1.0f;

			{
				{
					const Vec2f vec_grad = grid.gradC(Vec2i(0,0));

					CHECK((FLOAT)vec_grad(0) == 0.0f);
					CHECK((FLOAT)vec_grad(1) == 0.0f);
				}
			}
		}


		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);
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
			grid.data() = {
				0,	0,	0,
				0,	2,	0,
				0,	0,	0,

				0,  0,	0,
				0,	1,	2,
				0,	0,	0,

				0,	0,	0,
				0,	0,	0,
				0,	0,	0
			};
//			grid.fill(0);
//			grid(Vec3i(0,0,0)) = 1.0f;
//			grid(Vec3i(0,0,1)) = 2.0f;
//			grid(Vec3i(-1,0,0)) = 2.0f;

			{
				const Vec3f vec_grad = grid.gradC(Vec3i(0,0,0));

				CHECK((FLOAT)vec_grad(0) == -1);
				CHECK((FLOAT)vec_grad(1) == 0);
				CHECK((FLOAT)vec_grad(2) == 1);
			}

			// Safe differencing - checks boundaries and uses
			// central/forward/back as appropriate.
			{
				Vec3f vec_grad = grid.grad(Vec3i(0,0,0));

				CHECK((FLOAT)vec_grad(0) == -1);
				CHECK((FLOAT)vec_grad(1) == 0);
				CHECK((FLOAT)vec_grad(2) == 1);

				vec_grad = grid.grad(Vec3i(0,-1,0));

				CHECK((FLOAT)vec_grad(0) == 0);
				CHECK((FLOAT)vec_grad(1) == 1);
				CHECK((FLOAT)vec_grad(2) == 0);

				vec_grad = grid.grad(Vec3i(1,0,1));

				CHECK((FLOAT)vec_grad(0) == -2);
				CHECK((FLOAT)vec_grad(1) == 0);
				CHECK((FLOAT)vec_grad(2) == 0);

			}
		}
	}


	/**
	 * Divergence.
	 */
	WHEN("divergence_d2f_by_dx2")
	{
		//! [Divergence] Divergence
		// ==== Setup ====
		Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);
		grid.dx(2);
		grid.data() = {	1,	1,	1,
						1,	1,	1,
						1,	1,	1,

						1,	1,	1,
						1,	0,	1,
						1,	1,	1,

						1,	1,	1,
						1,	1,	1,
						1,	1,	1 };

		// ==== Action ====
		const FLOAT d2f_by_dx2 = grid.divergence(Vec3i(0,0,0));

		// ==== Confirm ====
		CHECK(d2f_by_dx2 == Approx(-0.75f).epsilon(0.00001f));
		//! [Divergence]
	}

	/**
	 * Using delta x.
	 */
	WHEN("delta_x")
	{
		using namespace boost::numeric::ublas;

		// Default and explicitly setting.
		{
			//! [Delta x setter]
			// ==== Setup ====
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);

			// ==== Confirm ====
			CHECK(grid.dx() == 1.0f);

			// ==== Action ====
			grid.dx(2.0f);

			// ==== Confirm ====
			CHECK(grid.dx() == 2.0f);
			//! [Delta x setter]
		}

		// Spatial derivatives with dx != 1.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);
			grid.dx(2.0f);
			grid.fill(0);
			grid(Vec3i(0,0,0)) = 1.0f;
			grid(Vec3i(0,0,1)) = 2.0f;
			grid(Vec3i(-1,0,0)) = 2.0f;

			CHECK((FLOAT)grid.gradF(Vec3i(0,0,0))(0) == -0.5f);
			CHECK((FLOAT)grid.gradF(Vec3i(0,0,0))(1) == -0.5f);
			CHECK((FLOAT)grid.gradF(Vec3i(0,0,0))(2) == 0.5f);
			CHECK((FLOAT)grid.gradB(Vec3i(0,0,0))(0) == -0.5f);
			CHECK((FLOAT)grid.gradB(Vec3i(0,0,0))(1) == 0.5f);
			CHECK((FLOAT)grid.gradB(Vec3i(0,0,0))(2) == 0.5);
			CHECK((FLOAT)grid.gradC(Vec3i(0,0,0))(0) == -0.5f);
			CHECK((FLOAT)grid.gradC(Vec3i(0,0,0))(1) == 0);
			CHECK((FLOAT)grid.gradC(Vec3i(0,0,0))(2) == 0.5f);

		}
	}

	/**
	 * Interpolation.
	 */
	WHEN("interpolate1D")
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

			grid.interp(input, pos);

			CHECK(input[0] == 0.7f);
		}
	}
	WHEN("interpolate2D")
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

			grid.interp(input, pos);

			CHECK(input[0] == Approx(0.4f).epsilon(0.00001f));
			CHECK(input[1] == Approx(0.8f).epsilon(0.00001f));

			grid.interp(input, pos);
			CHECK(input[0] == Approx(0.6f).epsilon(0.00001f));
		}
	}

	WHEN("interpolate3D")
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

			grid.interp(input, pos);

			CHECK(input[0 /**00x*/] == 0.4f);
			CHECK(input[1 /**01x*/] == 1.0f);
			CHECK(input[2 /**10x*/] == 0.0f);
			CHECK(input[3 /**11x*/] == 1.0f);

			grid.interp(input, pos);
			CHECK(input[0 /**0yx*/] == Approx(0.85f).epsilon(0.00001f));
			CHECK(input[1 /**1yx*/] == Approx(0.75f).epsilon(0.00001f));

			grid.interp(input, pos);
			CHECK(input[0 /**zyx*/] == Approx(0.8f).epsilon(0.00001f));
			CHECK(input.size() == 1);
		}
	}

	WHEN("interpolate_grid_2D")
	{
		//! [Interpolation]
		// ==== Setup ====
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);
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
		CHECK(val1 == Approx(0.0f).epsilon(0.00001f));
		CHECK(val2 == Approx(0.5f).epsilon(0.00001f));
		CHECK(val3 == Approx(1.0f).epsilon(0.00001f));
		CHECK(val4 == Approx(1.0f).epsilon(0.00001f));
		CHECK(val5 == Approx(1.0f).epsilon(0.00001f));
		//! [Interpolation]
	}


	/**
	 * Gradient interpolation.
	 */
	WHEN("grad_forwardinterp")
	{
		Grid <FLOAT,2> grid(Vec2u(5, 5), Vec2i(-2,-2), 0);
		grid.fill(0);
		grid(Vec2i(-1,-1)) = 1.0f;
		grid(Vec2i(-1,0)) = 1.0f;
		grid(Vec2i(0,1)) = 2.0f;
		grid(Vec2i(1,1)) = 2.0f;

		Vec2f pos = Vec2f(0.0f, 0.0f);
		Vec2f vec_grad = grid.gradF(pos);

		CHECK((FLOAT)vec_grad(0) == Approx(0.0f).epsilon(0.00001f));
		CHECK((FLOAT)vec_grad(1) == Approx(2.0f).epsilon(0.00001f));
	}


	/**
	 * Entropy satisfying gradient.
	 */
	WHEN("grad_entropy_satisfying")
	{
		// ==== Setup ====
		/**
			Row major order: (x,y) =>
			(-1,-1),(-1,0),	(-1,1)
			(0,-1),	(0,0),	(0,1)
			(1,-1),	(1,0),	(1,1)

		*/
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);
		Vec2i pos(0, 0);
		Vec2f vec_grad;

		// grad +'ve
		grid.data() = {	0,	0,	0,
						0,	1,	3,
						0,	3,	0};

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		CHECK((FLOAT)vec_grad(0) == Approx(1.0f).epsilon(0.00001f));
		CHECK((FLOAT)vec_grad(1) == Approx(1.0f).epsilon(0.00001f));


		// ==== Setup ====

		// grad -'ve
		grid.data() = {	0,	3,	0,
						3,	1,	0,
						0,	0,	0};

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		CHECK((FLOAT)vec_grad(0) == Approx(-1.0f).epsilon(0.00001f));
		CHECK((FLOAT)vec_grad(1) == Approx(-1.0f).epsilon(0.00001f));


		// ==== Setup ====

		// div -'ve
		grid.data() = {	0,	2,	0,
						3,	1,	2,
						0,	3,	0};

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		CHECK((FLOAT)vec_grad(0) == Approx(0.0f).epsilon(0.00001f));
		CHECK((FLOAT)vec_grad(1) == Approx(0.0f).epsilon(0.00001f));

		// ==== Setup ====

		// div +'ve
		grid.data() = {	0,	6,	0,
						6,	9,	1,
						0,	1,	0};

		// ==== Action ====

		vec_grad = grid.gradE(pos);


		// ==== Confirm ====

		CHECK((FLOAT)vec_grad(0) == Approx(-5.0f).epsilon(0.00001f));
		CHECK((FLOAT)vec_grad(1) == Approx(-5.0f).epsilon(0.00001f));
	}

	WHEN("curvature")
	{
		// 2D.
		{
			Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);
			grid.data() = {	1,	1,	1,
							1,	0,	1,
							1,	1,	1};

			const float singularity_curvature = grid.curv(Vec2i(0,0));

			CHECK(singularity_curvature == Approx(2).epsilon(0.00001f));

			grid.data() = {	 1,	 1,	 1,
							 0,	 0,  1,
							-1,	 0,	 1};

			const float corner_curvature = grid.curv(Vec2i(0,0));

			CHECK(corner_curvature == Approx(1).epsilon(0.00001f));
		}

		// 3D.
		{
			Grid <FLOAT,3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);
			grid.data() = {
				1,	1,	1,
				1,	1,	1,
				1,	1,	1,

				1,	1,	1,
				1,	0,	1,
				1,	1,	1,

				1,	1,	1,
				1,	1,	1,
				1,	1,	1};

			const float singularity_curvature_3D = grid.curv(Vec3i(0,0,0));

			CHECK(singularity_curvature_3D == Approx(3).epsilon(0.00001f));

			grid.data() = {
				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1};
			const float corner_curvature_3D = grid.curv(Vec3i(0,0,0));

			CHECK(corner_curvature_3D == Approx(1).epsilon(0.00001f));


			grid.data() = {
				 1,	 1,	 1,
				 1,	 1,	 1,
				 1,	 1,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				 0,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1};
			const float sharp_corner_curvature_3D = grid.curv(Vec3i(0,0,0));

			CHECK(sharp_corner_curvature_3D == Approx(1.5f).epsilon(0.00001f));
		}
	}
}


SCENARIO("LazyGrid")
{
	GIVEN("a 3x3x3 grid")
	{
		LazyGrid<FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1, -1, -1), 7.0f);

		THEN("the grid is initialised as inactive with background value returned by queries")
		{
			/// [LazyGrid initialisation]
			// ==== Confirm ====
			CHECK(grid.is_active() == false);
			CHECK(grid.data().size() == 0);
			CHECK(grid.background() == 7.0f);
			CHECK(grid.get(Vec3i(1, 1, 1)) == 7.0f);
			/// [LazyGrid initialisation]
		}

		WHEN("the grid is activated")
		{
			/// [LazyGrid activation]
			// ==== Action ====
			grid.activate();

			THEN("memory is allocated and the grid filled with background value")
			{
				// ==== Confirm ====
				CHECK(grid.is_active() == true);
				CHECK(grid.data().size() == 3 * 3 * 3);

				CHECK(grid.get(Vec3i(1, 1, 1)) == 7.0f);
				CHECK(grid.get(Vec3i(0, 1, 1)) == 7.0f);
			}

			AND_WHEN("a value is modified")
			{
				// ==== Action ====
				grid.get(Vec3i(1, 1, 1)) = 5.0f;

				THEN("that position's value is updated")
				{
					// ==== Confirm ====
					CHECK(grid.get(Vec3i(1, 1, 1)) == 5.0f);
					CHECK(grid.get(Vec3i(0, 1, 1)) == 7.0f);
				}
			}
			/// [LazyGrid activation]
		}

		WHEN("the grid is deactivated after being modified")
		{
			/// [LazyGrid deactivation]
			// ==== Action ====
			grid.get(Vec3i(1, 1, 1)) = 5.0f;
			grid.background()        = 3.0f;
			grid.deactivate();

			THEN("memory is deallocated and the queries report the background value")
			{
				// ==== Confirm ====
				CHECK(grid.is_active() == false);
				CHECK(grid.data().size() == 0);
				CHECK(grid.background() == 3.0f);
				CHECK(grid.get(Vec3i(1, 1, 1)) == 3.0f);
				/// [LazyGrid deactivation]
			}
		}
	}
}
