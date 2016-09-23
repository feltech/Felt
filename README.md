Deformable surface simulation and visualisation library
-------------------------------------------------------

Documentation: http://feltech.github.io/Felt/

Felt is an implementation of the
[sparse field level set](https://en.wikipedia.org/wiki/Level_set_(data_structures)#Sparse_field)
method for deformable surfaces, and the
[marching cubes](https://en.wikipedia.org/wiki/Marching_cubes) method of polygonising them.

The surface is embedded in a grid, where each grid node stores its distance to the surface. Thus
the surface can be extracted from the points where the distance is zero - it is an *implicit*
surface representation.

Felt makes heavy use of spatial partitioning to reduce memory footprint, provide localised
deformations and repolygonisation, and aid in parallelism.

Felt leverages C++14 and is header-only, with a dependency on the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) BLAS library.


##To do
The library is in a stable state, but the documentation needs work.

There are Doxygen annotations for almost all functions and classes, and quite extensive test cases.
I aim to finish tidying up the test cases, which now use the awesome
[Catch](https://github.com/philsquared/Catch) BDD style testing library, and use these as integral
parts of the documentation via snippets.

There is experimental support for the physics engine [Bullet](http://bulletphysics.org) via the
amazing game engine [Urho3D](https://urho3d.github.io/) in my test-bed demo app
[UrFelt](https://github.com/feltech/UrFelt).  I aim to move the Bullet and Urho3D classes developed
there to form an optional part of the main Felt project.


##Example
A super simple example from the tests:

```cpp
GIVEN("a 9x9 2-layer surface with a singularity seed at the centre")
{
	// 2D surface with 2x narrow band layers, respectively.
	using Surface_t = Surface<2, 2>;
	// Construct the surface.
	Surface_t surface(Vec2u(9, 9));

	// Grid to use for checking data.
	Grid<FLOAT, 2> isogrid_check(Vec2u(5, 5), Vec2i::Zero(), 0);

	// Create seed point in the centre
	surface.seed(Vec2i(0, 0));

	WHEN("we expand by 0.6")
	{
		// Expand the narrow band.
		surface.update([](auto& pos, auto& isogrid) {
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

			isogrid_check.vdata() = isogrid_check.vdata() - surface.isogrid().snapshot().vdata();
			const FLOAT diff = isogrid_check.vdata().sum();
			CHECK(diff == Approx(0));

			CHECK(surface.layer(-2).size() == 0);
			CHECK(surface.layer(-1).size() == 1);
			CHECK(surface.layer(0).size() == 4);
			CHECK(surface.layer(1).size() == 8);
			CHECK(surface.layer(2).size() == 12);
		}
	}
}
```

##Data structures

###Surfaces
The [Surface](@ref felt::Surface) is where most of the magic happens.  It leverages the tracking
grid structures to implement a spatially partitioned deformable surface.

###Rendering
Felt provides [Poly](@ref felt::Poly) objects held in a [PolyGrid](@ref felt::PolyGridBase),
to allow spatial partitions of the surface to be polygonised. Although surfaces
can have arbitrary dimension, polygonisation can only be done for 2D (lines) and 3D (triangles)
surfaces.

###Grids
The lowest level component in Felt is the [Grid](@ref felt::GridBase), with arbitrary
dimension and storing arbitrary types.

There are two types of simple grid, [Eager](@ref felt::EagerGridBase) and
[Lazy](@ref felt::LazyGrid). Lazy grids allow the internal memory to be released, returning a
background value for all queries.

Built on this are the [Partitioned](@ref felt::PartitionedGridBase) grid types, which spoof
other grid types, but internally use lazy subgrids to store the data, allowing sections of the
grid to be deactivated.

A key data structure in Felt is to couple a list of grid points with the grid itself, termed a
*tracking list*, which can keep track of active nodes.

[LookupGrid](@ref felt::LookupGridBase)s have one or more tracking lists, and have each grid node
store its index in the tracking list, so that tracking lists can be queryed (and hence updated) by
location. There can be multiple tracking lists for a grid, and grid nodes can hold a
[Single](@ref felt::EagerSingleLookupGrid) list index, or
[Multi](@ref felt::EagerMultiLookupGrid)ple list indices (one for each tracking
list).

Secondly, there are [TrackedGrid](@ref felt::TrackedGridBase)s, which compose a lookup grid to track
active grid nodes, as well as a grid to store arbitrary data.

There are then [Partitioned](@ref felt::TrackingPartitionedGridBase) variants of the lookup and
tracked grids.

##Legal
_Copyright David Feltell all rights reserved._
_Released under [GPL 3](http://www.gnu.org/licenses/gpl-3.0.en.html) (for now):_
