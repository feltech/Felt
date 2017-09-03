Deformable surface simulation and visualisation library
-------------------------------------------------------

[![Build Status](https://travis-ci.org/feltech/Felt.svg?branch=master)](https://travis-ci.org/feltech/Felt)

**Documentation: http://feltech.github.io/Felt/**

Felt is an implementation of the
[sparse field level set](https://en.wikipedia.org/wiki/Level_set_(data_structures)#Sparse_field)
method for deformable surfaces, and the
[marching cubes](https://en.wikipedia.org/wiki/Marching_cubes) method of polygonising them.

The surface is embedded in a grid, where each grid node stores its distance to the surface. Thus
the surface can be extracted from the points where the distance is zero - it is an *implicit*
surface representation.

Felt makes heavy use of spatial partitioning to reduce memory footprint, provide localised
deformations and repolygonisation, and aid in parallelism.

Felt leverages C++14 and is header-only, with a dependencies on [boost](http://www.boost.org/) and
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) libraries.


##To do
The library is in a stable state, but the documentation needs work.

There are Doxygen annotations for almost all functions and classes, and quite extensive test cases.
I aim to finish tidying up the test cases, which now use the awesome
[Catch](https://github.com/philsquared/Catch) BDD style testing library, and use these as integral
parts of the documentation via snippets.

There is ongoing development to support the physics engine [Bullet](http://bulletphysics.org) via 
the amazing game engine [Urho3D](https://urho3d.github.io/) in my test-bed demo app
[UrFelt](https://github.com/feltech/UrFelt).  I aim to move the Bullet and Urho3D classes developed
there to form an optional part of the main Felt project.


##Example
A super simple example from the tests:

```cpp
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

		WHEN("we expand by 0.6")
		{
			surface.update([](const auto& pos_, const auto& isogrid_) {
				(void)pos_; (void)isogrid_;
				return -0.6f;
			});

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
```

##Legal
_Copyright David Feltell all rights reserved._
_Released under [GPL 3](http://www.gnu.org/licenses/gpl-3.0.en.html) (for now):_
