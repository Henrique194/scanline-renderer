# Scanline Renderer

Program to demonstrate scanline rendering using z-sorted spans.

In this implementation, polygon faces must not be interpenetrating. Also,
correct sorting is not guaranteed if two polygonal objects butt up against
each other. In other words, each polygonal object must be made of a continuous,
non-self-intersecting skin, and polygonal objects must not interpenetrate or
touch in order for proper sorting to result. More complex, slower sorting is
required to make those cases work reliably.

Polygon processing could be considerably more efficient if polygons shared
common edges and edges shared common vertices. Also, indirection to vertices
could be used to avoid having to copy all the vertices during every clip test.
Outcode-type testing could be used to determine completely clipped or unclipped
polygons ahead of time, avoiding the need to clip and copy entirely for such
polygons. Outcode-type tests work best in viewspace, with the frustum normalized
so that the field of view is 90 degrees, so simple compares, rather than dot
products, can be used to categorize points with respect to the frustum. See
*Computer Graphics*, by Foley & van Dam, or *Procedural Elements of Computer
Graphics*, by Rogers, for further information.

## Credits

This implementation is based on Chapter 67 of Michael Abrash's Graphics
Programming Black Book. The scanline renderer described here employs the same
span-sorting technique found in the test version of Quake (QTEST1.ZIP) that id
Software released in early March 1996, and the same data structures went on to
form the foundation of the final software renderer in Quake.
