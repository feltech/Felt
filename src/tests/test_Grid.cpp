#include "catch.hpp"
#include <omp.h>

#define _TESTING

#include "Felt/Grid.hpp"
#include "Utils.hpp"

using namespace felt;

/**
 * Test the Grid class.
 */
SCENARIO("Grid - basics")
{
	/// [Grid - basics: GIVEN 3x7x11]
	GIVEN("a 3x7x11 grid with no offset and background value of 0")
	{
		Grid<FLOAT, 3> grid(Vec3u(3, 7, 11), Vec3i::Zero(), 0);
		/// [Grid - basics: GIVEN 3x7x11]

		/// [Grid - basics: THEN memory is allocated]
		THEN("memory is allocated and the size is reported correctly")
		{
			CHECK(grid.size()(0) == 3);
			CHECK(grid.size()(1) == 7);
			CHECK(grid.size()(2) == 11);
			CHECK(grid.data().size() == (3 * 7 * 11));
		}
		/// [Grid - basics: THEN memory is allocated]

		THEN("we can test if locations lie within the grid")
		{
			CHECK(grid.inside(Vec3i(-1,0,0)) == false);
			CHECK(grid.inside(Vec3i(0,0,0)) == true);
			CHECK(grid.inside(Vec3i(1,2,3)) == true);
			CHECK(grid.inside(Vec3i(3,7,11)) == false);
			CHECK(grid.inside(Vec3f(0,-0.00001,0)) == false);
			CHECK(grid.inside(Vec3f(0,0,9.99999)) == true);
		}

		WHEN("some positions values are set")
		{
			grid(Vec3i(0, 0, 0))  = 13.0f;
			grid(Vec3i(1, 2, 3))  = 17.0f;
			grid(Vec3i(2, 6, 10)) = 19.0f;

			//! [Get and set]
			THEN("querying those positions returns the same values")
			{
				CHECK(grid.get(Vec3i(1, 2, 3)) == 17.0f);
				CHECK(grid(Vec3i(1, 2, 3)) == 17.0f);
			}

			THEN("expected elements of the underlying array contain those values")
			{
				CHECK(grid.data()[0] == 13.0f);
				CHECK(grid.data()[grid.data().size() - 1] == 19.0f);
			}
			//! [Get and set]

			AND_WHEN("the size is updated")
			{
				grid.size(Vec3u(5, 11, 13));

				THEN("the memory is resized")
				{
					CHECK(grid.data().size() == (5 * 11 * 13));
				}
				THEN("the memory layout has changed")
				{
					CHECK(grid.get(Vec3i(0, 0, 0)) == 13.0f); // Same.
					CHECK(grid.get(Vec3i(1, 2, 3)) == 0); // Now moved...
					CHECK(grid.get(Vec3i(0, 7, 11)) == 17.0f); // ...to here.
				}
			}
		}

		WHEN("the grid is filled with a value")
		{
			grid.fill(7);

			THEN("every point in the grid contains that value")
			{
				INT sum = 0;
				for (INT i = 0; i < grid.data().size(); i++) {
					sum += grid.data()[i];
				}

				CHECK(sum == 3*7*11 * 7);
			}
		}
	}

	GIVEN("a 7x11x13 grid with (-3,-3,-3) offset and background value of 0")
	{
		using GridType = Grid <FLOAT, 3>;
		Vec3u size(7, 11, 13);
		Vec3i offset(-3, -3, -3);
		GridType grid(size, offset, 0);

		//! [Position index]
		THEN("the index of a point in the data array is reported correctly")
		{
			CHECK(GridType::index(Vec3i(1, 0, -1), size, offset) == 613);
			CHECK(grid.index(Vec3i(1, 0, -1)) == 613);
		}

		THEN("the point represented by an index in the data array is reported correctly")
		{
			CHECK(grid.index(613) == Vec3i(1, 0, -1));
			CHECK(GridType::index(613, size, offset) == Vec3i(1, 0, -1));
		}
		//! [Position index]

		THEN("we can test if locations lie within the offset grid")
		{
			CHECK(grid.inside(Vec3i(-2,0,0)) == true);
			CHECK(grid.inside(Vec3i(-4,0,0)) == false);

			AND_WHEN("we change the offset")
			{
				// ==== Action ====
				grid.offset(Vec3i(-1,-1,-1));

				THEN("testing the same positions inside the grid yields a different result")
				{
						// ==== Confirm ====
					CHECK(grid.inside(Vec3i(-2,0,0)) == false);
					CHECK(grid.inside(Vec3i(-4,0,0)) == false);
				}
			}
		}

		WHEN("editing points in the offset grid")
		{
			grid(Vec3i(-3,-3,-3)) = 21.0f;
			grid(Vec3i(-1,0,-1)) = 23.0f;

			THEN("we can retrieve the values from the offset positions")
			{
				CHECK(grid.data()[0] == 21.0f);
				CHECK(grid(Vec3i(-1,0,-1)) == 23.0f);
			}
		}
	}

	/// [Grid - basics: GIVEN 2x2x2]
	GIVEN("a 2x2x2 grid with (-1,-1,-1) offset and background value of 7")
	{
		Grid<FLOAT, 3> grid(Vec3u(2, 2, 2), Vec3i(-1, -1, -1), 7);
		/// [Grid - basics: GIVEN 2x2x2]

		/// [Grid - basics: WHEN we cycle over all]
		WHEN("we cycle over all grid locations with a range based for loop")
		{
			FLOAT sum = 0;
			std::vector <Vec3i> visits;

			for (auto pos : grid)
			{
				visits.push_back(pos);
				sum += grid(pos);
			}
			/// [Grid - basics: WHEN we cycle over all]

			/// [Grid - basics: THEN we have visited all]
			THEN("we have visited all locations in the grid")
			{
				std::vector<Vec3i> expected_visits = {
					Vec3i(-1, -1, -1), Vec3i(-1, -1, 0), Vec3i(-1, 0, -1), Vec3i(-1, 0, 0),
					Vec3i(0, -1, -1), Vec3i(0, -1, 0), Vec3i(0, 0, -1), Vec3i(0, 0, 0)
				};

				CHECK(visits == expected_visits);
			}
			/// [Grid - basics: THEN we have visited all]

			/// [Grid - basics: THEN filled with the background value]
			THEN("the grid is found to be filled with the background value")
			{
				CHECK(sum == 7 * 2*2*2);
			}
			/// [Grid - basics: THEN filled with the background value]
		}
	}
}

SCENARIO("Grid - numerical methods")
{

	GIVEN("A 1D grid and input vector of values")
	{
		Grid <FLOAT,1> grid = Grid <FLOAT,1>();

		std::vector<FLOAT> input = { 1.0f, 0 };

		WHEN("we interpolate a distance of 0.3 between the vector values")
		{
			using Vec1f = Eigen::Matrix<FLOAT, 1, 1>;
			const Vec1f pos(0.3);

			grid.interp(input, pos);

			THEN("the input vector now contains the single interpolated value")
			{
				CHECK(input.size() == 1);
				CHECK(input[0] == 0.7f);
			}
		}
	}


	GIVEN("A 2D grid and input vector of values")
	{
		Grid <FLOAT,2> grid = Grid <FLOAT,2>();

		std::vector<FLOAT> input = std::vector<FLOAT>(4);
		input[0 /*00*/] = 2.0f;
		input[1 /*01*/] = 0;
		input[2 /*10*/] = 0.0f;
		input[3 /*11*/] = 1.0;

		WHEN("we bilinearly interpolate a distance of (0.8, 0.5) between the vector values")
		{
			const Vec2f pos(0.8f, 0.5f);

			grid.interp(input, pos);

			THEN("the input vector now contains the correct interpolated values")
			{
				CHECK(input.size() == 2);
				CHECK(input[0] == Approx(0.4f));
				CHECK(input[1] == Approx(0.8f));

				AND_WHEN("we interpolate along this line")
				{
					grid.interp(input, pos);

					THEN("the input vector now contains the final interpolated value")
					{
						CHECK(input.size() == 1);
						CHECK(input[0] == Approx(0.6f));
					}
				}
			}
		}
	}


	GIVEN("A 3D grid and input vector of values")
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

		std::vector<FLOAT> input = std::vector<FLOAT>(8);
		input[0 /**000*/] = 0.0f;
		input[1 /**001*/] = 0.8f;
		input[2 /**010*/] = 1.0f;
		input[3 /**011*/] = 1.0f;
		input[4 /**100*/] = 0.0f;
		input[5 /**101*/] = 0.0f;
		input[6 /**110*/] = 1.0f;
		input[7 /**111*/] = 1.0f;

		WHEN("we trilinearly interpolate a distance of (0.5, 0.75, 0.5) between the vector values")
		{
			Vec3f pos(0.5f, 0.75f, 0.5f);

			grid.interp(input, pos);

			THEN("the input vector now contains the correct interpolated values")
			{
				CHECK(input.size() == 4);
				CHECK(input[0 /**00x*/] == 0.4f);
				CHECK(input[1 /**01x*/] == 1.0f);
				CHECK(input[2 /**10x*/] == 0.0f);
				CHECK(input[3 /**11x*/] == 1.0f);

				AND_WHEN(
					"we bilinearly interpolate a distance of (0.8, 0.5) between these vector values"
				) {
					grid.interp(input, pos);

					THEN("the input vector now contains the correct interpolated values")
					{
						CHECK(input.size() == 2);
						CHECK(input[0 /**0yx*/] == Approx(0.85f));
						CHECK(input[1 /**1yx*/] == Approx(0.75f));

						AND_WHEN("we interpolate along this line")
						{
							grid.interp(input, pos);

							THEN("the input vector now contains the final interpolated value")
							{
								CHECK(input.size() == 1);
								CHECK(input[0 /**zyx*/] == Approx(0.8f));
							}
						}
					}
				}
			}
		}
	} // End GIVEN("A 3D grid and input vector of values").


	GIVEN("a 3x3 grid with (-1,-1) offset and background value of 0")
	{
		Grid <FLOAT,2> grid(Vec2u(3, 3), Vec2i(-1,-1), 0);

		THEN("the spatial resolution has a default of 1")
		{
			CHECK(grid.dx() == 1.0f);

			AND_WHEN("we change the spatial resolution to 2")
			{
				grid.dx(2.0f);

				THEN("the resolution is reported as 2")
				{
					CHECK(grid.dx() == 2.0f);
				}
			}
		}


		WHEN("we modify some values near the centre")
		{
			grid(Vec2i(-1,-1)) = 1.0f;
			grid(Vec2i(-1,0)) = 1.0f;
			grid(Vec2i(0,1)) = 2.0f;
			grid(Vec2i(1,1)) = 2.0f;

			AND_WHEN("we interpolate at some real locations using explicit function calls")
			{
				const FLOAT val1 = grid.interp(Vec2f(0.0f, 0.0f));
				const FLOAT val2 = grid.interp(Vec2f(-0.5f, -0.5f));
				const FLOAT val3 = grid.interp(Vec2f(0.5f, 0.5f));

				THEN("the interpolated values are correct")
				{
					CHECK(val1 == Approx(0.0f));
					CHECK(val2 == Approx(0.5f));
					CHECK(val3 == Approx(1.0f));
				}
			}

			AND_WHEN("we implicitly interpolate at a real location using value getter")
			{
				const FLOAT val = grid.val(Vec2f(0.5f, 0.5f));

				THEN("the interpolated value is correct")
				{
					CHECK(val == Approx(1.0f));
				}
			}

			AND_WHEN("we implicitly interpolate at a real location using bracket operator")
			{
				const FLOAT val = grid(Vec2f(0.5f, 0.5f));

				THEN("the interpolated value is correct")
				{
					CHECK(val == Approx(1.0f));
				}
			}
		}

		WHEN("we calculate the forward difference gradient at the centre")
		{
			const Vec2f grad = grid.gradF(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we calculate the backward difference gradient at the centre")
		{
			const Vec2f grad = grid.gradB(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we calculate the central difference gradient at the centre")
		{
			const Vec2f grad = grid.gradC(Vec2i(0,0));

			THEN("the gradient is zero")
			{
				CHECK(grad == Vec2f(0, 0));
			}
		}

		WHEN("we set the central grid location to 1")
		{
			grid(Vec2i(0,0)) = 1.0f;

			AND_WHEN("we calculate the forward difference gradient at the centre")
			{
				const Vec2f grad = grid.gradF(Vec2i(0,0));

				THEN("the gradient is negative")
				{
					CHECK(grad == Vec2f(-1, -1));
				}
			}

			AND_WHEN("we calculate the backward difference gradient at the centre")
			{
				const Vec2f grad = grid.gradB(Vec2i(0,0));

				THEN("the gradient is positive")
				{
					CHECK(grad == Vec2f(1, 1));
				}
			}

			AND_WHEN("we calculate the central difference gradient at the centre")
			{
				const Vec2f grad = grid.gradC(Vec2i(0,0));

				THEN("the gradient is zero")
				{
					CHECK(grad == Vec2f(0, 0));
				}
			}


			WHEN("we calculate the forward difference gradient at a real position")
			{
				Vec2f pos(-0.3f, 0.3f);

				Vec2f grad = grid.gradF(pos);

				THEN("the gradient is interpolated correctly")
				{
					CHECK(grad == ApproxVec(Vec2f(-0.28, -0.49)));
				}
			}
		}


		WHEN("we have a (positively directed) entropic flow")
		{
			grid.data() = {
			//<	  -	y  + >
				0,	0,	0,	// -
				0,	1,	3,	// x
				0,	3,	0	// +
			};

			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(1, 1));
			}
		}

		WHEN("we have a (negatively directed) entropic flow")
		{
			grid.data() = {
			//<	  -	y  + >
				0,	3,	0,	// -
				3,	1,	0,	// x
				0,	0,	0	// +
			};

			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(-1, -1));
			}
		}

		WHEN("we have a positively divergent flow")
		{
			grid.data() = {
			//<	  -	y  + >
				0,	2,	0,	// -
				2,	1,	3,	// x
				0,	3,	0	// +
			};
			THEN("the entropy satisfying gradient is zero")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(0, 0));
			}
		}

		WHEN("we have a negatively divergent flow")
		{
			grid.data() = {
			//<	  -	y  + >
				0,	6,	0,	// -
				6,	9,	1,	// x
				0,	1,	0	// +
			};
			THEN("the entropy satisfying gradient gives an upwind value")
			{
				CHECK(grid.gradE(Vec2i(0, 0)) == Vec2f(-5, -5));
			}
		}

		WHEN("we set up a positive divergence at the centre")
		{
			grid.data() = {
				1,	1,	1,
				1,	0,	1,
				1,	1,	1
			};

			THEN("the curvature is at its maximum")
			{
				CHECK(grid.curv(Vec2i(0,0)) == 2);
			}
		}

		WHEN("we set up a 'corner' at the centre")
		{
			grid.data() = {
				 1,	 1,	 1,
				 0,	 0,  1,
				-1,	 0,	 1
			};

			THEN("the curvature is 1")
			{
				CHECK(grid.curv(Vec2i(0,0)) == 1);
			}
		}

	} // End GIVEN("a 3x3 grid with (-1,-1) offset and background value of 0")


	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 0")
	{
		Grid <FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1,-1,-1), 0);

		WHEN("we set up a gradient about the centre")
		{
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

			THEN("the forward difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradF(Vec3i(0,0,0));
				CHECK(grad == Vec3f(-1, -1, 1));
			}

			THEN("the backward difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradB(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 1, 1));
			}

			THEN("the central difference at the centre is as expected")
			{
				const Vec3f grad = grid.gradC(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 0, 1));
			}

			THEN("the safe gradient at the centre uses the central difference")
			{
				const Vec3f grad = grid.grad(Vec3i(0,0,0));

				CHECK(grad == Vec3f(-1, 0, 1));
			}

			THEN("the safe gradient at the bottom face uses the forward difference")
			{
				const Vec3f grad = grid.grad(Vec3i(0,-1,0));

				CHECK(grad == Vec3f(0, 1, 0));
			}

			THEN("the safe gradient at the right-forward edge uses the backward difference")
			{
				const Vec3f grad = grid.grad(Vec3i(1,0,1));

				CHECK(grad == Vec3f(-2, 0, 0));
			}

			AND_WHEN("we decrease the spatial resolution")
			{
				grid.dx(2);

				THEN("the forward difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradF(Vec3i(0,0,0));
					CHECK(grad == Vec3f(-0.5, -0.5, 0.5));
				}

				THEN("the backward difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradB(Vec3i(0,0,0));

					CHECK(grad == Vec3f(-0.5, 0.5, 0.5));
				}

				THEN("the central difference at the centre is as expected")
				{
					const Vec3f grad = grid.gradC(Vec3i(0,0,0));

					CHECK(grad == Vec3f(-0.5, 0, 0.5));
				}
			}
		}

		WHEN("we set up a positive divergence at the centre")
		{
			grid.data() = {
				1,	1,	1,
				1,	1,	1,
				1,	1,	1,

				1,	1,	1,
				1,	0,	1,
				1,	1,	1,

				1,	1,	1,
				1,	1,	1,
				1,	1,	1
			};

			THEN("the divergence is calculated correctly")
			{
				CHECK(grid.divergence(Vec3i(0,0,0)) == Approx(6));
			}

			AND_WHEN("we decrease the spatial resolution")
			{
				grid.dx(2);

				THEN("the divergence is calculated correctly")
				{
					CHECK(grid.divergence(Vec3i(0,0,0)) == Approx(0.75));
				}
			}

			THEN("the mean curvature is at the maximum")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 3);
			}
		}

		WHEN("we set up a 'corner' along two dimensions")
		{
			grid.data() = {
				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1
			};

			THEN("the mean curvature is 1")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 1);
			}
		}

		WHEN("we set up a sharp 'corner' along all three dimensions")
		{
			grid.data() = {
				 1,	 1,	 1,
				 1,	 1,	 1,
				 1,	 1,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				 0,	 0,	 1,

				 1,	 1,	 1,
				 0,	 0,	 1,
				-1,	 0,	 1
			};

			THEN("the mean curvature is 1.5")
			{
				CHECK(grid.curv(Vec3i(0,0,0)) == 1.5);
			}
		}
	}
} // End SCENARIO("Grid - numerical methods")


SCENARIO("LazyGrid")
{
	GIVEN("a 3x3x3 grid with (-1,-1,-1) offset and background value of 7")
	{
		LazyGrid<FLOAT, 3> grid(Vec3u(3, 3, 3), Vec3i(-1, -1, -1), 7.0f);

		THEN("the grid is initialised as inactive with background value returned by queries")
		{
			/// [LazyGrid initialisation]
			CHECK(grid.is_active() == false);
			CHECK(grid.data().size() == 0);
			CHECK(grid.background() == 7.0f);
			CHECK(grid.get(Vec3i(1, 1, 1)) == 7.0f);
			/// [LazyGrid initialisation]
		}

		WHEN("the grid is activated")
		{
			/// [LazyGrid activation]
			grid.activate();

			THEN("memory is allocated and the grid filled with background value")
			{
				CHECK(grid.is_active() == true);
				CHECK(grid.data().size() == 3 * 3 * 3);

				CHECK(grid.get(Vec3i(1, 1, 1)) == 7.0f);
				CHECK(grid.get(Vec3i(0, 1, 1)) == 7.0f);
			}

			AND_WHEN("a value is modified")
			{
				grid.get(Vec3i(1, 1, 1)) = 5.0f;

				THEN("that position's value is updated")
				{
					CHECK(grid.get(Vec3i(1, 1, 1)) == 5.0f);
					CHECK(grid.get(Vec3i(0, 1, 1)) == 7.0f);
				}
			}
			/// [LazyGrid activation]
		}

		WHEN("the grid is deactivated after being modified")
		{
			/// [LazyGrid deactivation]
			grid.get(Vec3i(1, 1, 1)) = 5.0f;
			grid.background() = 3.0f;
			grid.deactivate();

			THEN("memory is deallocated and queries report the background value")
			{
				CHECK(grid.is_active() == false);
				CHECK(grid.data().size() == 0);
				CHECK(grid.background() == 3.0f);
				CHECK(grid.get(Vec3i(1, 1, 1)) == 3.0f);
				/// [LazyGrid deactivation]
			}
		}
	}
}
