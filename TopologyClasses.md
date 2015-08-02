# Topology Classes #



## Summary ##
`MGTopology` is an abstract base class of class `MGCellBase` and class `MGComplex`.

![http://mgcl.googlecode.com/svn/wiki/image/mg-topology.png](http://mgcl.googlecode.com/svn/wiki/image/mg-topology.png)

`MGTopology`:

> `MGTopology` is an abstract class that represents all of the topological structures.

`MGCellBase`:

> `MGCellBase` is the base class of whole Cell.
> `MGCellBase` holds only the binder Cell pointer (if a Cell belongs to a partnership) that represents the partner relationship.
> `MGCellBase` is a super class of classes `MGPVertex` and `MGCellNB`.

`MGPVertex`:

> `MGPVertex` is a special Cell that is a parameter Cell and also a boundary of an Edge simultaneously.
> `MGPVertex` is a subclass of `MGCellBase`.

> As member data, `MGPVertex` holds a parameter value of an `MGEdge` (of a curve that is a geometry of the Edge)
> and the Edge pointer that is the star Cell of the `MGPVertex`.
> The parameter value is an information of class parameter Cell,
> and the Edge pointer is an information of ordinal class Boundary's star Cell pointer.

> `MGPVertex` cannot be a binder Cell, can be only a parameter Cell.
> Binder Cell of an `MGPVertex` is `MGBVertex`.

`MGCellNB`:

> `MGCellNB` represents the theoretical Cell structures without boundary data out of whole Cell information,
> is a subclass of `MGCellBase`. Only an instance of class `MGCellNB` (or its subclasses) can be a member Cell of an `MGComplex`.
> Also, only an instance of class `MGCellNB` can be a binder Cell.
> `MGCellNB` is the super class of classes `MGBVertex`, `MGEdge`, and `MGCell`.

> As member data, `MGCellNB` has [geometry](GeometryClasses.md) (Cell's extent),
> the parent Complex pointer (if the `MGCellNB` is a member Cell of an `MGComplex`),
> and partner member Cells (if the `MGCellNB` is a binder Cell).

`MGBVertex`:

> `MGBVertex` is 0D manifold (point),
> has no boundaries,
> and is a binder Cell version of an `MGPVertex`.
> `MGBVertex` has no member data.
> The information of the super class `MGCellNB` suffices.

`MGEdge`:

> `MGEdge` is 1D manifold (line), does not have any general Boundary,
> but has instances of class `MGPVertex` as the boundary data.
> An array of two instances of class `MGPVertex` (start and end point) is the member data.
> An instance of `MGEdge` can be either a parameter Cell or a binder Cell.

`MGCell`:

> `MGCell` is an abstract class that represents manifolds whose dimensions are 2 or more than 2.
> Face, solid, or greater manifold dimension's objects are the examples.

> As a member data, Cell has array of Boundaries.
> The first one in the array is the outer boundary of the star Cell (if an outer boundary exists),
> and other boundaries are inner boundaries.
> We allow inactive boundaries are included in a Cell.

> Cell can be either a parameter Cell or a binder Cell.

`MGFace`:

> `MGFace` is a class of 2D manifold, whose boundary is loop.
> `MGFace` is a subclass of class `MGCell`.
> Class `MGFace` is a constituent of class `MGShell`.

> ![http://mgcl.googlecode.com/svn/wiki/image/trimmed-surface.png](http://mgcl.googlecode.com/svn/wiki/image/trimmed-surface.png)

`MGComplex`:

> `MGComplex` is two ordered sets of `MGCell`s,
> which are a set of parameter Cells and a set of binder Cells.
> Cells of different manifold dimensions may be members of a Complex.
> Thus a Complex can express so-called nonmanifold models.
> Geometries of all member Cells of a Complex, either binder Cells or parameter Cells,
> have the same world space coordinates,
> which may have different manifold dimensions.
> When a Complex is to express a body in our ordinal world, all Cells included in the body Complex have ordinal world coordinates of the body.
> When a Complex is to express a boundary of a Cell, all Cells included in the Complex have the star Cell's parameter space coordinates.

> `MGComplex` is a super class of `MGBoundary` class that is a boundary expression of `MGCell`.
> `MGComplex` does not have [geometry](GeometryClasses.md) data as a member data, but `MGCell` has.

`MGBoundarynD`:

> Boundary is an ordered set of same manifold dimension Cells, all of which must be orientable.
> An instance of class `MGBoundarynD` must not intersect with itself.
> Boundary class is a special Complex class, that is, is a subclass of Complex whose Cells all have the same manifold dimension.
> Boundary is a boundary of the star Cell.
> The manifold dimension is lower by one than that of the star Cell.

`MGLoop`:

> Loop is an ordered set of edges (parameter Cells), and is a boundary of a face.
> `MGLoop` is a subclass of class `MGBoundary`.
> An edge in a Loop is connected to at most one neighbor edge before and after the edge.
> When the start point of the first edge is connected to the end point of the last edge,
> the loop is closed, and is homeomorphic to a circle (that is a unit circle S^1).
> Coordinate data of edges in a loop that is a boundary of a face is the face's parameter `(u, v)`.

`MGShell`:

> Shell is an ordered set of Faces, which are connected each other through a binder Edge at their parameter Edge.
> Shell makes a boundary of a solid, and is a subclass of class Boundary.

## Terminologies of Topological Elements Relationship ##

![http://mgcl.googlecode.com/svn/wiki/image/topoelm-relation.png](http://mgcl.googlecode.com/svn/wiki/image/topoelm-relation.png)

`Star`:

> 自身が境界 (の一部) となっている 1 (以上) 高い多様体次元のセル。

`Boundary`:

> あるセルの境界を構成するセル群 (Complex)。

> 通常のモデルでは star/boudary は 1 だけ高い、または低い多様体次元のセルを求める。
> Boundary はあるセルの境界を構成するセル群 (Complex) である。

> 図では C11 の star Cell は face1 であり、
> face1 の boundary は C10, C11, C12, C13 である。

`Neighbour`:

> 隣接セル、すなわち同一セルの境界となっている隣の境界である。

> 隣接 (neighbour) とは、同じセル (一つ多様体次元が高い) の境界を構成するセルで、
> 自身のセルに隣接するものである。

> 図では C11 の neighbour は C10, C12 である。

`Partner`:

> パートナー (partner) とは自身のセルのとは異なるセルの境界を表現していて、
> 自身のセルと binder を通して接続されている。

> 図では C11 の partner は C2 である。
> binder B1 を通して接続されている。

## MGCell ##
最小位相要素 `MGGeometry` を構成要素に持ち、均一の多様体次元の Geometry を表現する。

  * `MGFace` - 2-dimensional manifold (trimmed surface).
  * `MGEdge` - 1-dimensional manifold.

![http://mgcl.googlecode.com/svn/wiki/image/mg-cell-fig.png](http://mgcl.googlecode.com/svn/wiki/image/mg-cell-fig.png)

## MGComplex ##
セルの有機的な結合。
`MGBoundary` の構成要素となる。

## MGBoundary ##
Class `MGBoundary` is a boundary of an instance of class `MGCell`.

  * Class `MGBoundarynD` represents an n-dimensional boundary.
  * Subclasses are as follows:
    * `MGPVertex`: 1-dimensional boundary.
    * `MGLoop`: 2-dimensional boundary.
    * `MGShell`: 3-dimensional boundary.

![http://mgcl.googlecode.com/svn/wiki/image/mg-boundary-fig.png](http://mgcl.googlecode.com/svn/wiki/image/mg-boundary-fig.png)