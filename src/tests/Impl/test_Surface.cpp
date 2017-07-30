#include <unordered_set>

#include "../catch.hpp"
#include <Felt/Impl/Common.hpp>
#include "Utils.hpp"
#include <Felt/Impl/Surface.hpp>


namespace Felt
{
template <Dim D, LayerId L>
ListIdx layer_size(const Surface<D, L>& surface, const LayerId layer_id);

template <Dim D, LayerId L>
ListIdx layer_size(
	const Surface<D, L>& surface, const VecDi<D>& pos_child_, const LayerId layer_id
);
/// Alias std::to_string for compactness.
template <typename... Args>
auto s(Args&&... args) -> decltype(std::to_string(std::forward<Args>(args)...));

}

using namespace Felt;

using PosSet = std::unordered_set< Vec2i, matrix_hash<Vec2i> >;

SCENARIO("Surface - global updates")
{


GIVEN("a 2-layer 2D surface in a 7x7 isogrid with 3x3 spatial partitions")
{
	Surface<2, 2> surface(Vec2i(7, 7), Vec2i(3, 3));

	THEN("the isogrid is initialised correctly")
	{
		CHECK(surface.isogrid().size() == Vec2i(7, 7));
		CHECK(surface.isogrid().children().data().size() == 9);
		CHECK(surface.isogrid().children().get(Vec2i(0, 0)).size() == Vec2i(3, 3));
		CHECK(surface.isogrid().children().get(Vec2i(0, 0)).data().size() == 0);
		CHECK(surface.isogrid().size() == Vec2i(7, 7));
		// Grid is initialised to all points 'outside' the surface (since there is no surface yet).
		CHECK(surface.isogrid().get(Vec2i(0, 0)) == 3);
	}
}

GIVEN("a 2-layer 2D surface in a 9x9 isogrid with 3x3 partitions")
{
	// 2D surface with 2x narrow band layers, respectively.
	using SurfaceType = Surface<2, 2>;

	// Construct the surface.
	SurfaceType surface(Vec2i(9, 9), Vec2i(3, 3));


	WHEN("a singularity seed is created at the centre")
	{
		// Create seed point in the centre
		surface.seed(Vec2i(0, 0));

		INFO(stringify_grid_slice(surface.isogrid()));

		// Trivially check centre of seed is indeed a zero-level point (i.e. point
		// on the surface).
		THEN("the value at the centre of the grid is 0")
		{
			const FLOAT val_centre = surface.isogrid().get(Vec2i(0, 0));
			CHECK(val_centre == 0);
		}

		// Grid to use for checking data.
		Impl::Grid::Snapshot<FLOAT, 2> isogrid_check(Vec2i(5, 5), Vec2i::Zero(), 0);

		THEN("the surface data matches a singularity seed point")
		{
			// A 2D 2-layer singularity (seed) point should look like the following.
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.array() -= surface.isogrid().snapshot()->array();

			const Distance diff = isogrid_check.array().sum();

			CHECK(diff == 0);

			// Check appropriate points have been added to narrow band layers.
			CHECK(layer_size(surface, -2) == 0);
			CHECK(layer_size(surface, -1) == 0);
			CHECK(layer_size(surface, 0) == 1);
			CHECK(layer_size(surface, 1) == 4);
			CHECK(layer_size(surface, 2) == 8);
		}

		AND_WHEN("we expand the surface one unit outwards")
		{
			surface.update([](const auto& pos_, const auto& isogrid_) {
				(void)pos_; (void)isogrid_;
				return -1.0f;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data matches a surface of radius 1")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    2,    1,    0,   -1,    0,    1,    2,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.array() = isogrid_check.array() -
					surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();

				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 1);
				CHECK(layer_size(surface, 0) == 4);
				CHECK(layer_size(surface, 1) == 8);
				CHECK(layer_size(surface, 2) == 12);
			}

			AND_WHEN("iterating over layer 0 and recording each point hit")
			{
				PosSet pos_leafs;

				surface.isogrid().leafs(
					surface.layer_idx(0),
					[&pos_leafs](auto pos_) {
						pos_leafs.insert(pos_);
					}
				);

				THEN("the four zero layer points are recorded")
				{
					PosSet pos_leafs_expected{
						Vec2i(-1, 0), Vec2i(1, 0), Vec2i(0, 1), Vec2i(0, -1)
					};
					CHECK(pos_leafs == pos_leafs_expected);
				}
			}

			AND_WHEN("we expand by one unit again")
			{
				surface.update([](const auto& pos_, const auto& isogrid_) {
					(void)pos_; (void)isogrid_;
					return -1.0f;
				});


				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("the grid data matches a surface of radius 2")
				{
					isogrid_check.data() = {
						  3,    3,    3,    3,    2,    3,    3,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    2,    1,    0,    1,    2,    3,    3,
						  3,    2,    1,    0,   -1,    0,    1,    2,    3,
						  2,    1,    0,   -1,   -2,   -1,    0,    1,    2,
						  3,    2,    1,    0,   -1,    0,    1,    2,    3,
						  3,    3,    2,    1,    0,    1,    2,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    3,    3,    2,    3,    3,    3,    3
					};

					isogrid_check.array() = (
						isogrid_check.array() - surface.isogrid().snapshot()->array()
					);
					const Distance diff = isogrid_check.array().sum();

					CHECK(diff == Approx(0));
				}

				AND_WHEN("we expand by one unit 9 more times")
				{
					for (UINT i = 0; i < 9; i++)
					{
						surface.update([](const auto& pos_, const auto& isogrid_) {
							(void)pos_; (void)isogrid_;
							return -1.0f;
						});
					}

					INFO(stringify_grid_slice(surface.isogrid()));

					THEN("the surface data matches an area completely consumed by the surface")
					{
						isogrid_check.data() = {
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,
							 -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3,   -3
						};

						isogrid_check.array() = (
							isogrid_check.array() - surface.isogrid().snapshot()->array()
						);
						const Distance diff = isogrid_check.array().sum();

						CHECK(diff == Approx(0));

						CHECK(layer_size(surface, 0) == 0);
						CHECK(layer_size(surface, -1) == 0);
						CHECK(layer_size(surface, -2) == 0);
						CHECK(layer_size(surface, 1) == 0);
						CHECK(layer_size(surface, 2) == 0);
					}
				}
			}

			AND_WHEN("we contract the surface by 1 unit inwards")
			{
				surface.update([](const auto& pos_, const auto& isogrid_) {
					(void)pos_; (void)isogrid_;
					return 1.0f;
				});

				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("the surface has collapsed back to a singularity")
				{
					// A 2D 2-layer singularity (seed) point should look like the following.

					isogrid_check.data() = {
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    2,    3,    3,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    2,    1,    0,    1,    2,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    3,    3,    2,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3
					};

					isogrid_check.array() = isogrid_check.array() -
						surface.isogrid().snapshot()->array();
					INFO(stringify_grid_slice(*surface.isogrid().snapshot()));
					INFO(stringify_grid_slice(isogrid_check));

					const Distance diff = isogrid_check.array().sum();

					CHECK(diff == 0);

					CHECK(layer_size(surface, -2) == 0);
					CHECK(layer_size(surface, -1) == 0);
					CHECK(layer_size(surface, 0) == 1);
					CHECK(layer_size(surface, 1) == 4);
					CHECK(layer_size(surface, 2) == 8);
				}

				AND_WHEN("we contract the surface by 1 unit inwards again")
				{
					surface.update([](const auto& pos_, const auto& isogrid_) {
						(void)pos_; (void)isogrid_;
						return 1.0f;
					});

					INFO(stringify_grid_slice(surface.isogrid()));

					AND_WHEN("iterating over layer 0 and recording each point hit")
					{
						using PosSet = std::unordered_set< Vec2i, matrix_hash<Vec2i> >;
						PosSet pos_leafs;

						surface.isogrid().leafs(
							surface.layer_idx(0),
							[&pos_leafs](auto pos_) {
								pos_leafs.insert(pos_);
							}
						);

						THEN("there are no points recorded")
						{
							CHECK(pos_leafs.size() == 0);
						}
					}

					THEN("the surface data matches an area completely outside the surface")
					{
						isogrid_check.data() = {
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3
						};

						isogrid_check.array() = isogrid_check.array() -
							surface.isogrid().snapshot()->array();

						const Distance diff = isogrid_check.array().sum();

						CHECK(diff == 0);

						CHECK(layer_size(surface, -2) == 0);
						CHECK(layer_size(surface, -1) == 0);
						CHECK(layer_size(surface, 0) == 0);
						CHECK(layer_size(surface, 1) == 0);
						CHECK(layer_size(surface, 2) == 0);
					}
				}
			} // End when contract by 1.
		} // End when expand by 1.

		WHEN("we expand by 0.6")
		{
			surface.update([](const auto& pos_, const auto& isogrid_) {
				(void)pos_; (void)isogrid_;
				return -0.6f;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data and layers are as expected")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
					  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
					  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
					  3,  2.4,  1.4,  0.4, -0.6,  0.4,  1.4,  2.4,    3,
					  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
					  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
					  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				// snapshot() copies the spatially partitioned grid into a single
				// unpartitioned grid. array() returns an Eigen vector with the grid data, giving
				// access to arithmetic and other BLAS functions.
				isogrid_check.array() =
					isogrid_check.array() - surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 1);
				CHECK(layer_size(surface, 0) == 4);
				CHECK(layer_size(surface, 1) == 8);
				CHECK(layer_size(surface, 2) == 12);
			}

			AND_WHEN("we contract by 0.6")
			{
				surface.update([](const auto& pos_, const auto& isogrid_) {
					(void)pos_; (void)isogrid_;
					return 0.6f;
				});

				INFO(stringify_grid_slice(surface.isogrid()));

				THEN("the surface is once more a seed")
				{
					isogrid_check.data() = {
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    2,    3,    3,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    2,    1,    0,    1,    2,    3,    3,
						  3,    3,    3,    2,    1,    2,    3,    3,    3,
						  3,    3,    3,    3,    2,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3,
						  3,    3,    3,    3,    3,    3,    3,    3,    3
					};

					CHECK(layer_size(surface, -2) == 0);
					CHECK(layer_size(surface, -1) == 0);
					CHECK(layer_size(surface, 0) == 1);
					CHECK(layer_size(surface, 1) == 4);
					CHECK(layer_size(surface, 2) == 8);

					isogrid_check.array() =
						isogrid_check.array() - surface.isogrid().snapshot()->array();
					const Distance diff = isogrid_check.array().sum();
					CHECK(diff == Approx(0));
				}

				AND_WHEN("we contract by 0.6 again")
				{
					surface.update([](const auto& pos_, const auto& isogrid_) {
						(void)pos_; (void)isogrid_;
						return 0.6f;
					});

					INFO(stringify_grid_slice(surface.isogrid()));

					THEN("the surface has completely collapsed and all points are outside")
					{
						isogrid_check.data() = {
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3,
							  3,    3,    3,    3,    3,    3,    3,    3,    3
						};

						CHECK(layer_size(surface, -2) == 0);
						CHECK(layer_size(surface, -1) == 0);
						CHECK(layer_size(surface, 0) == 0);
						CHECK(layer_size(surface, 1) == 0);
						CHECK(layer_size(surface, 2) == 0);

						isogrid_check.array() =
							isogrid_check.array() - surface.isogrid().snapshot()->array();
						const Distance diff = isogrid_check.array().sum();
						CHECK(diff == Approx(0));
					}
				}
			}
		} // End WHEN we expand by 0.6
	}
}


GIVEN("a 2-layer 2D surface in a 21x21 isogrid with 2x2 partitions")
{
	using SurfaceType = Surface<2, 2>;
	SurfaceType surface(Vec2i{21, 21}, Vec2i{2, 2});
	// Grid to set values of manually, for checking against.
	Impl::Grid::Snapshot<Distance, 2> isogrid_check(Vec2i{21, 21}, Vec2i::Zero(), 0);

	WHEN("an initial seed is expanded such that the central partition is all inside")
	{
		// Create seed point and expand the narrow band.
		surface.seed(Vec2i(0, 0));

		for (UINT i = 0; i < 5; i++)
			surface.update([](const auto& pos_, const auto& isogrid_) {
				(void)pos_; (void)isogrid_;
				return -1.0f;
			});

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the isogrid values are as expected")
		{
			isogrid_check.data() = {
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,
				3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,
				3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -3,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -3,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,   -2,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    2,    1,    0,   -1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.array() = isogrid_check.array() - surface.isogrid().snapshot()->array();
			const Distance diff = isogrid_check.array().sum();
			CHECK(diff == Approx(0));
		}

		THEN("the central partition is deactivated")
		{
			// Surrounding partitions are still active.
			CHECK(surface.isogrid().children().get(Vec2i{0,1}).is_active() == true);
			CHECK(surface.isogrid().children().get(Vec2i{0,-1}).is_active() == true);
			CHECK(surface.isogrid().children().get(Vec2i{1,0}).is_active() == true);
			CHECK(surface.isogrid().children().get(Vec2i{-1,0}).is_active() == true);
			// But central partition has now been deactivated.
			CHECK(surface.isogrid().children().get(Vec2i{0,0}).is_active() == false);
		}
	}
}

} // End SCENARIO Surface - global up


SCENARIO("Surface - local updates")
{
	GIVEN("a 9x9 2-layer surface with 2x2 partitions initialised with a seed point in the centre")
	{
		using SurfaceType = Surface<2, 2>;
		SurfaceType surface(Vec2i{9, 9}, Vec2i{2, 2});

		// Grid to set values of manually, for checking against.
		Impl::Grid::Snapshot<Distance, 2> isogrid_check(Vec2i{9, 9}, Vec2i::Zero(), 0);

		// Create seed point and expand the narrow band.
		surface.seed(Vec2i(0, 0));

		WHEN("we contract the surface by 1 unit inwards")
		{
			surface.update_start();
			surface.delta(Vec2i(0, 0), 1.0f);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the surface data matches an area completely outside the surface")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.array() = isogrid_check.array() -
					surface.isogrid().snapshot()->array();

				const Distance diff = isogrid_check.array().sum();

				CHECK(diff == 0);

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 0);
				CHECK(layer_size(surface, 0) == 0);
				CHECK(layer_size(surface, 1) == 0);
				CHECK(layer_size(surface, 2) == 0);
			}
		}

		WHEN("we expand by 1 unit")
		{
			surface.update_start();
			surface.delta(Vec2i(0, 0), -1.0f);
			surface.update_end();

			THEN("the grid data matches a surface of radius 1")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    2,    1,    0,   -1,    0,    1,    2,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.array() = isogrid_check.array() -
					surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();

				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 1);
				CHECK(layer_size(surface, 0) == 4);
				CHECK(layer_size(surface, 1) == 8);
				CHECK(layer_size(surface, 2) == 12);
			}

			AND_WHEN("we modify a couple of points and calculate the affected narrow band points")
			{
				// Clean up from previous update.
				surface.update_start();
				// Add a couple of points that could affect the narrow band.
				surface.delta(Vec2i(0, 1), 0.3f);
				surface.delta(Vec2i(1, 0), 0.3f);
				surface.update_end();

				THEN("the affected narrow band points are as expected")
				{
					using PosArray = std::vector<Vec2i>;
					PosArray check_layers_pos[5];
					check_layers_pos[2 + -2] = PosArray();
					check_layers_pos[2 + -1] = PosArray({
						Vec2i(0, 0)
					});
					check_layers_pos[2 + 0] = PosArray({
					// We don't care for now about zero-layer points.
				//		Vec2i(0,1),
				//		Vec2i(1,0)
					});
					check_layers_pos[2 + 1] = PosArray({
						// For (0,1):
						Vec2i(-1, 1), Vec2i(1, 1), Vec2i(0, 2),
						// For (1,0):
						Vec2i(2, 0), Vec2i(1, -1)
					});

					check_layers_pos[2 + 2] = PosArray({
						// For (0,1):
						Vec2i(-2, 1), Vec2i(2, 1),
						Vec2i(-1, 2), Vec2i(1, 2),
						Vec2i(0, 3),
						// For (1,0):
						Vec2i(3, 0), Vec2i(1, -2), Vec2i(2, -1)
					});

					for (LayerId layer_id = -2; layer_id <= 2; layer_id++)
					{
						if (layer_id == 0)
							continue;

						const TupleIdx layer_idx = 2 + layer_id;
						ListIdx num_leafs = 0;
						surface.affected().leafs(layer_idx, [&num_leafs](auto) {
							num_leafs++;
						});


						INFO(
							"Layer " + s(layer_id) + " at index " + s(layer_idx) +
							" number of leafs " +
							s(num_leafs) + " == " +
							s(check_layers_pos[layer_idx].size())
						);
						CHECK(num_leafs == check_layers_pos[layer_idx].size());

						for (auto pos : check_layers_pos[layer_idx])
						{
							bool found = false;
							surface.affected().leafs(layer_idx, [&found, &pos](auto pos_) {
								found |= pos == pos_;
							});

							INFO(
								"Affected grid layer " + s(layer_id) + " at index " +
								s(layer_idx) + " should contain (" + s(pos(0)) + "," +
								s(pos(1)) + ")"
							);
							CHECK(found == true);
						}

						surface.affected().leafs(
							layer_idx,
							[&check_layers_pos, layer_id, layer_idx](const auto& pos_) {
								auto iter = std::find(
									check_layers_pos[layer_idx].begin(),
									check_layers_pos[layer_idx].end(), pos_
								);

								INFO(
									"Checking list layer " + s(layer_id) + " at index " +
									s(layer_idx) + " should contain (" + s(pos_(0)) + "," +
									s(pos_(1)) + ")"
								);
								CHECK(iter != check_layers_pos[layer_idx].end());
							}
						);
					}
				}
			}

			AND_WHEN("we cycle a square region partially containing the surface")
			{
				UINT num_visited = 0;
				PosSet pos_visits;
				surface.update(
					Vec2i(1, 0), Vec2i(3, 3),
					[&num_visited, &pos_visits](
						const Vec2i& pos, const auto& isogrid
					) -> Distance {
						(void)isogrid;
						num_visited++;
						pos_visits.insert(pos);
						return 0;
					}
				);

				THEN("we only visit the points in the region")
				{
					CHECK(num_visited == 1);
					CHECK(pos_visits == PosSet{Vec2i(1, 0)});
				}
			}

			AND_WHEN("we cycle a square region completely containing the surface")
			{
				UINT num_visited = 0;
				PosSet pos_visits;
				surface.update(
					Vec2i(-100, -100), Vec2i(100, 100),
					[&num_visited, &pos_visits](
						const Vec2i& pos, const auto& isogrid
					) -> Distance {
						(void)isogrid;
						num_visited++;
						pos_visits.insert(pos);
						return 0;
					}
				);

				THEN("we only visit the valid points")
				{
					CHECK(
						pos_visits ==
						(PosSet{Vec2i(1, 0), Vec2i(-1, 0), Vec2i(0, 1), Vec2i(0, -1)})
					);
					CHECK(num_visited == 4);
				}
			}
		} // WHEN we expand by 1 unit.

		AND_WHEN("we update square region containing the surface")
		{
			surface.update(
				Vec2i(-1, -1), Vec2i(1, 1),
				[](const Vec2i& pos, const auto& isogrid) -> Distance {
					(void)pos; (void)isogrid;
					return -0.6f;
				}
			);

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data and layers are as expected")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
					  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
					  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
					  3,  2.4,  1.4,  0.4, -0.6,  0.4,  1.4,  2.4,    3,
					  3,    3,  2.4,  1.4,  0.4,  1.4,  2.4,    3,    3,
					  3,    3,    3,  2.4,  1.4,  2.4,    3,    3,    3,
					  3,    3,    3,    3,  2.4,    3,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				// snapshot() copies the spatially partitioned grid into a single
				// unpartitioned grid. array() returns an Eigen vector with the grid data, giving
				// access to arithmetic and other BLAS functions.
				isogrid_check.array() =
					isogrid_check.array() - surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 1);
				CHECK(layer_size(surface, 0) == 4);
				CHECK(layer_size(surface, 1) == 8);
				CHECK(layer_size(surface, 2) == 12);
			}
		}
	} // End GIVEN 9x9 2-layer with 2x2 partitions with seed.
}// End SCENARIO Surface - local updates.


SCENARIO("Surface - layer interactions")
{

	GIVEN("a 16x9 2-layer surface with two small regions side-by-side")
	{
		Surface<2, 2> surface(Vec2i(16, 9), Vec2i::Constant(3));
		Impl::Grid::Snapshot<FLOAT, 2> isogrid_check(Vec2i(16, 9), Vec2i::Zero(), 0);

		// Create two seed points and expand the narrow band.
		surface.seed(Vec2i(-4, 0));
		surface.seed(Vec2i(4, 0));
		surface.update([](const auto& pos_, const auto& isogrid_) {
			(void)pos_; (void)isogrid_;
			return -1.0f;
		});

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("outermost layers in central partitions are as expected")
		{
			CHECK((layer_size<2,2>(surface, Vec2i(0,0), 2)) == 3);
			CHECK((layer_size<2,2>(surface, Vec2i(1,0), 2)) == 3);
		}

		THEN("the surface is in the expected state")
		{
			isogrid_check.data() = {
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
				  3,    2,    1,    0,   -1,    0,    1,    2,    3,    2,    1,    0,   -1,    0,    1,    2,
				  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
				  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
				  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
				  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
			};

			isogrid_check.array() = isogrid_check.array() - surface.isogrid().snapshot()->array();
			const Distance diff = isogrid_check.array().sum();
			CHECK(diff == Approx(0));

			CHECK(layer_size(surface, -2) == 0);
			CHECK(layer_size(surface, -1) == 2);
			CHECK(layer_size(surface, 0) == 8);
			CHECK(layer_size(surface, 1) == 16);
			CHECK(layer_size(surface, 2) == 24);
		}

		WHEN("we expand the subsurfaces towards one-another")
		{
			surface.update_start();
			surface.delta(Vec2i(-3, 0), -1.0f);
			surface.delta(Vec2i(3, 0), -1.0f);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the centremost partitions contain the expected number of outer layer points")
			{
				CHECK((layer_size<2,2>(surface, Vec2i(0,0), 2)) == 3);
				CHECK((layer_size<2,2>(surface, Vec2i(1,0), 2)) == 2);
			}

			THEN("the surface is in the expected state")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    2,    3,    3,    3,    3,    3,    2,    2,    3,    3,    3,
					  3,    3,    3,    2,    1,    1,    2,    3,    3,    3,    2,    1,    1,    2,    3,    3,
					  3,    3,    2,    1,    0,    0,    1,    2,    3,    2,    1,    0,    0,    1,    2,    3,
					  3,    2,    1,    0,   -1,   -1,    0,    1,    2,    1,    0,   -1,   -1,    0,    1,    2,
					  3,    3,    2,    1,    0,    0,    1,    2,    3,    2,    1,    0,    0,    1,    2,    3,
					  3,    3,    3,    2,    1,    1,    2,    3,    3,    3,    2,    1,    1,    2,    3,    3,
					  3,    3,    3,    3,    2,    2,    3,    3,    3,    3,    3,    2,    2,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.array() = isogrid_check.array() - surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 4);
				CHECK(layer_size(surface, 0) == 12);
				CHECK(layer_size(surface, 1) == 20);
				CHECK(layer_size(surface, 2) == 27);
			}
		}

		WHEN("we expand/contract the subsurfaces towards one-another")
		{
			for (UINT i = 0; i < 10; i++)
			{
				surface.update_start();
				surface.delta(Vec2i(-3, 0), -1.0f);
				surface.delta(Vec2i(3, 0), -1.0f);
				surface.update_end();

				surface.update_start();
				surface.delta(Vec2i(-3, 1), 1.0f);
				surface.delta(Vec2i(-2, 0), 1.0f);
				surface.delta(Vec2i(-3, -1), 1.0f);

				surface.delta(Vec2i(3, 1), 1.0f);
				surface.delta(Vec2i(2, 0), 1.0f);
				surface.delta(Vec2i(3, -1), 1.0f);
				surface.update_end();
			}

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("outermost layers in central partitions are as expected")
			{
				CHECK((layer_size<2,2>(surface, Vec2i(0,0), 2)) == 3);
				CHECK((layer_size<2,2>(surface, Vec2i(1,0), 2)) == 3);
			}

			THEN("the surface has returned to its inital state")
			{
				isogrid_check.data() = {
					  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
					  3,    2,    1,    0,   -1,    0,    1,    2,    3,    2,    1,    0,   -1,    0,    1,    2,
					  3,    3,    2,    1,    0,    1,    2,    3,    3,    3,    2,    1,    0,    1,    2,    3,
					  3,    3,    3,    2,    1,    2,    3,    3,    3,    3,    3,    2,    1,    2,    3,    3,
					  3,    3,    3,    3,    2,    3,    3,    3,    3,    3,    3,    3,    2,    3,    3,    3,
					  3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,    3
				};

				isogrid_check.array() =
					isogrid_check.array() - surface.isogrid().snapshot()->array();
				const FLOAT diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0));

				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 2);
				CHECK(layer_size(surface, 0) == 8);
				CHECK(layer_size(surface, 1) == 16);
				CHECK(layer_size(surface, 2) == 24);
			}
		}
	}

	GIVEN(
		"a 11x11x11 3-layer surface with 3x3x3 partitions initialised with a 1 unit radius surface"
	) {
		Surface<3,3> surface(Vec3i(11,11,11), Vec3i(3,3,3));
		surface.seed(Vec3i(0,0,0));
		surface.update([](const auto& pos, const auto& grid) {
			(void)pos; (void)grid;
			return -1.0f;
		});

		INFO(stringify_grid_slice(surface.isogrid()));

		THEN("the layer lists have the expected size")
		{
			CHECK(layer_size(surface, -3) == 0);
			CHECK(layer_size(surface, -2) == 0);
			CHECK(layer_size(surface, -1) == 1);
			CHECK(layer_size(surface, 0) == 6);
			CHECK(layer_size(surface, 1) == 18);
			CHECK(layer_size(surface, 2) == 38);
			CHECK(layer_size(surface, 3) == 66);
		}

		AND_WHEN("we expand and contract two nearby points using local update")
		{
			surface.update_start();
			surface.delta(Vec3i(0,1,0), -1);
			surface.update_end();

			surface.update_start();
			surface.delta(Vec3i(0,2,0), 1);
			surface.delta(Vec3i(1,1,0), 1);
			surface.delta(Vec3i(-1,1,0), 1);
			surface.delta(Vec3i(0,1,1), 1);
			surface.delta(Vec3i(0,1,-1), 1);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the layer lists have the same size as before")
			{
				CHECK(layer_size(surface, -3) == 0);
				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 1);
				CHECK(layer_size(surface, 0) == 6);
				CHECK(layer_size(surface, 1) == 18);
				CHECK(layer_size(surface, 2) == 38);
				CHECK(layer_size(surface, 3) == 66);
			}
		}
	}

	GIVEN(
		"a 12x12 3-layer grid with two seeds diagonally opposite, and the left seed expanded"
		" towards the right"
	) {
		using SurfaceType = Surface<2, 3>;
		Vec2i size(12, 12);

		SurfaceType surface(size, Vec2i(2, 2));

		// Grid to set values of manually, for checking against.
		Impl::Grid::Snapshot<Distance, 2> isogrid_check(size, Vec2i::Zero(), 0);

		// Create seed point and expand the narrow band.
		surface.seed(Vec2i(-2, -2));
		surface.seed(Vec2i(2, 2));

		INFO(stringify_grid_slice(surface.isogrid()));

//    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
//    4,    4,    4,    4,    3,    4,    4,    4,    4,    4,    4,    4,
//    4,    4,    4,    3,    2,    3,    4,    4,    4,    4,    4,    4,
//    4,    4,    3,    2,    1,    2,    3,    4,    4,    4,    4,    4,
//    4,    3,    2,    1,    0,    1,    2,    3,    4,    4,    4,    4,
//    4,    4,    3,    2,    1,    2,    3,    4,    3,    4,    4,    4,
//    4,    4,    4,    3,    2,    3,    4,    3,    2,    3,    4,    4,
//    4,    4,    4,    4,    3,    4,    3,    2,    1,    2,    3,    4,
//    4,    4,    4,    4,    4,    3,    2,    1,    0,    1,    2,    3,
//    4,    4,    4,    4,    4,    4,    3,    2,    1,    2,    3,    4,
//    4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//    4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4

		surface.update_start();
		surface.delta(Vec2i(-2, -2), -1);
		surface.update_end();
		surface.update_start();
		surface.delta(Vec2i(-1, -2), -1);
		surface.update_end();
		surface.update_start();
		surface.delta(Vec2i(0, -2), -1);
		surface.update_end();

		INFO(stringify_grid_slice(surface.isogrid()));

//	      4,    4,    4,    4,    3,    3,    3,    4,    4,    4,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    3,    4,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    2,    3,    4,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    1,    2,    3,    4,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,    0,    1,    2,    3,    4,
//	      4,    3,    2,    1,    0,    0,    0,    1,    2,    3,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    2,    2,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    3,    3,    2,    1,    0,    1,    2,    3,
//	      4,    4,    4,    4,    4,    4,    3,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4


//		Bad:
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    3,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    3,    3,    3,    2,    1,    2,    3,    4,
//	      4,    4,    4,    4,    4,    4,    4,    3,    2,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4

//		Also bad:
//	      4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
//	      4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
//	      4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
//	      4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
//	      4,    4,    4,    4,    3,    3,    3,    3,    3,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    3,    3,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
//	      4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4

		WHEN(
			"we simultaneously expand the left seed and contract the right, then expand the "
			"left again, using local updates"
		) {

			surface.update_start();
			surface.delta(Vec2i(1, -2), -1);
			surface.delta(Vec2i(2, 2), 1);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			surface.update_start();
			surface.delta(Vec2i(1, -1), -1);
			surface.update_end();

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data is as expected")
			{

				isogrid_check.data() = {
					  4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
					  4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
					  4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
					  4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
					  3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
					  4,    3,    2,    1,    0,    0,    0,   -1,    0,    1,    2,    3,
					  4,    4,    3,    2,    1,    1,    1,    0,    1,    2,    3,    4,
					  4,    4,    4,    3,    2,    2,    2,    1,    2,    3,    4,    4,
					  4,    4,    4,    4,    3,    3,    3,    2,    3,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4
				};

				isogrid_check.array() =
					isogrid_check.array() - surface.isogrid().snapshot()->array();
				const FLOAT diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0).epsilon(0.000001f));

				CHECK(layer_size(surface, -3) == 0);
				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 5);
				CHECK(layer_size(surface, 0) == 11);
				CHECK(layer_size(surface, 1) == 15);
				CHECK(layer_size(surface, 2) == 19);
				CHECK(layer_size(surface, 3) == 23);
			}
		}

		WHEN(
			"we simultaneously expand the left seed and contract the right, then expand the "
			"left again, using global updates"
		) {

			surface.update([](const auto& pos_, const auto& isogrid_) -> Distance {
				(void)isogrid_;
				if (pos_ == Vec2i(1, -2))
					return -1;
				if (pos_ == Vec2i(2, 2))
					return 1;
				return 0;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

			surface.update([](const auto& pos_, const auto& isogrid_) -> Distance {
				(void)isogrid_;
				if (pos_ == Vec2i(1, -1))
					return -1;
				return 0;
			});

			INFO(stringify_grid_slice(surface.isogrid()));

			THEN("the grid data is as expected")
			{
				isogrid_check.data() = {
					  4,    4,    4,    4,    3,    3,    3,    3,    4,    4,    4,    4,
					  4,    4,    4,    3,    2,    2,    2,    2,    3,    4,    4,    4,
					  4,    4,    3,    2,    1,    1,    1,    1,    2,    3,    4,    4,
					  4,    3,    2,    1,    0,    0,    0,    0,    1,    2,    3,    4,
					  3,    2,    1,    0,   -1,   -1,   -1,   -1,    0,    1,    2,    3,
					  4,    3,    2,    1,    0,    0,    0,   -1,    0,    1,    2,    3,
					  4,    4,    3,    2,    1,    1,    1,    0,    1,    2,    3,    4,
					  4,    4,    4,    3,    2,    2,    2,    1,    2,    3,    4,    4,
					  4,    4,    4,    4,    3,    3,    3,    2,    3,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    3,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,
					  4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4
				};

				isogrid_check.array() =
					isogrid_check.array() - surface.isogrid().snapshot()->array();
				const Distance diff = isogrid_check.array().sum();
				CHECK(diff == Approx(0).epsilon(0.000001f));

				CHECK(layer_size(surface, -3) == 0);
				CHECK(layer_size(surface, -2) == 0);
				CHECK(layer_size(surface, -1) == 5);
				CHECK(layer_size(surface, 0) == 11);
				CHECK(layer_size(surface, 1) == 15);
				CHECK(layer_size(surface, 2) == 19);
				CHECK(layer_size(surface, 3) == 23);
			}
		}
	}
}// End SCENARIO Surface - layer interactions


namespace Felt
{
template <Dim D, LayerId L>
ListIdx layer_size(const Surface<D, L>& surface, const LayerId layer_id)
{
	ListIdx size = 0;

	const TupleIdx layer_idx = surface.layer_idx(layer_id);
	const PosArray& pos_idxs_child = surface.isogrid().children().lookup().list(layer_idx);

	for (ListIdx list_idx = 0; list_idx < pos_idxs_child.size(); list_idx++)
	{
		const PosIdx pos_idx_child = pos_idxs_child[list_idx];

		size += surface.isogrid().children().get(pos_idx_child).lookup().list(layer_idx).size();
	}

	return size;
}

template <Dim D, LayerId L>
ListIdx layer_size(
	const Surface<D, L>& surface, const VecDi<D>& pos_child_, const LayerId layer_id
) {
	ListIdx size = 0;

	const TupleIdx layer_idx = surface.layer_idx(layer_id);
	const PosIdx pos_idx_child = surface.isogrid().children().index(pos_child_);
	return surface.isogrid().children().get(pos_idx_child).lookup().list(layer_idx).size();

}

/// Alias std::to_string for compactness.
template <typename... Args>
auto s(Args&&... args) -> decltype(std::to_string(std::forward<Args>(args)...)) {
	return std::to_string(std::forward<Args>(args)...);
}
} // Felt.

