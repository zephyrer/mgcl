# Classes #



MGCL supports object-oriented programming.

  * [MGGeometry](GeometryClasses.md) - An abstract class for geometry classes.
  * [MGTopology](TopologyClasses.md) - An abstract class for topology classes.
  * [MGGel](Classes#MGGel.md) - An abstract class for class geometry, topology, attribute, and group.
  * [MGGroup](Classes#MGGroup.md) - Groups group elements (GoF's Composite design pattern)
  * [MGObject](Classes#MGObject.md) - An abstract class for `MGGeometry` and `MGTopology`.
  * Base Classes - Form geometry classes. `MGEReal` (extended real number), `MGKnot` (knot), etc.
  * Transformation Classes - Represents transformations. `MGVector`, `MGMatrix`, `MGTransform`, etc.
  * B-Representation Classes - `MGKnotVector` (knot vector), `MGBPointSeq` (coefficient of a B-Rep curve), etc.
  * Serialization Classes - Classes to serialize MGCL objects.

# Class Hierarchy Diagram #

![http://mgcl.googlecode.com/svn/wiki/image/object-list.png](http://mgcl.googlecode.com/svn/wiki/image/object-list.png)

## MGGel ##
An abstract class `MGGel` represents a group element.
Gel is the abbreviation of Group ELement.

`MGGel` is designed to be stored in `MGGroup` as an element.
The subclasses of `MGGel` are as follows:

  * `MGGroup`
  * `MGObject`
  * `MGAttrib`

### MGGroup ###
`MGGroup` is a class which contains instances of `MGGel`.

![http://mgcl.googlecode.com/svn/wiki/image/mg-group.png](http://mgcl.googlecode.com/svn/wiki/image/mg-group.png)

Because `MGGroup` is also a subclass of `MGGroup`, its structure is recursive (GoF's Composite design pattern).
The ownership of instances contained in a group is transferred to the group.

Main functions defined in this class as follows:

  * List manipulation methods (`pop`, `push`, `begin`, `end`, `front`, `back`)

### MGObject ###

Class `MGObject` is the superclass of `MGGeometry` and `MGTopology`.

Main functions defined in this class as follows:

  * Serialization of instances via class `MGOfstream` and `MGIfstream`.
  * Dump the instance.
  * Clone (deep copy) the instance
  * Transformation.

## Base Classes ##
Base classes form objects of subclasses of `MGObject`.

  * `MGEReal` - Extended real number; positive infinity, negative infinity and real numbers.
  * `MGInterval`, `MGBox` - 1 or N-dimensional ranges; bounds are defined by `MGEReal` values.
  * `MGNDDArray` - A non-decreasing sequence that consists of double (real) numbers.
  * `MGPosition` - Coordinates in N-dimensional space.
  * `MGKnot` - A knot; defined by a real number and multiplicity.
  * `MGKnotArray` - A sequence of knots.
  * `MGVector` - Vector in N-dimensional space; used transformations as well.

## Transformation Classes ##

Transformations (rigid motion, rotation, scaling, etc.) are performed by
using one of the classes `MGVector`, `MGMatrix`, and `MGTransf` or their combination.
Scaling can be performed by simply multiplying a double value to an instance of `MGObject`.

All geometry classes (except for `MGPPRep`: Piecewise Polynomial Representation) have operators
for transformations:
`operator+`, `operator-`, `operator+=`, `operator-=`,
`operator*`, `operator/`, `operator*=`, `operator/=`.

MGCL の幾何オブジェクトはすべて、座標変換によりオブジェクトを変換してもパラメータ値は変わりません。
幾何オブジェクト `g` と座標変換 `T` に関して次の式が成立します：

![http://mgcl.googlecode.com/svn/wiki/image/expr1.png](http://mgcl.googlecode.com/svn/wiki/image/expr1.png)

`eval()` は位置座標の評価であり、幾何オブジェクト `g` のあるパラメータ値 `t` の座標値を求めてから
変換 `T` を加えた点と、 `g` に座標変換 `T` を加えた後の幾何オブジェクトのパラメータ値 `t` の座標値は同じ点です。

### MGVector ###

`MGVector` is the superclass of `MGUnit_vector`.
Space dimension can be any non-negative number.

![http://mgcl.googlecode.com/svn/wiki/image/mg-vector.png](http://mgcl.googlecode.com/svn/wiki/image/mg-vector.png)

MGVector features:

  * Provides so many constructors.
  * Calculates length of a vector.
  * Computes scalar product of two vectors (`operator%`).
  * Computes vector product of two vectors (`operator*`).
> > Generally, vector product is defined in 3 dimensional space,
> > MGCL extends vector product to N dimensional space.
  * Computes triple product.
  * Tests if two vectors are perpendicular or parallel.
  * Calculates angles between two vectors, by any form of angle, sine, or cosine.
  * Computes interpolation of two vectors, by linear or rotation.
  * Tests if a vector is zero, has unit length.
  * Translate an instance of `MGObject`.

### MGMatrix, MGTransform ###
`MGMatrix` is a matix of N by N, where N is the space dimension.

`MGMatrix` provides transformation around the origin,
while `MGTransf` general transformation.

![http://mgcl.googlecode.com/svn/wiki/image/mg-matrix.png](http://mgcl.googlecode.com/svn/wiki/image/mg-matrix.png)

![http://mgcl.googlecode.com/svn/wiki/image/mg-transf.png](http://mgcl.googlecode.com/svn/wiki/image/mg-transf.png)

`MGMatrix`, `MGTransform` features:

  * Provides a number of constructors.
  * Computes the transposed matrix.
  * Computes the determinant of a matrix.
  * Scales objects.
  * Rotates objects by given angle around a given vector.
  * Transforms a given vector to X or Y axis (and vise versa)
  * Transforms two orthogonal vectors to X and Y axes (and vise versa)
  * Mirrors objects.

# B-Representatioal Curves/Surfaces Component Classes #

  * `MGKnotVector` - A knot vector; order `k` (`== degree - 1`) and B-Rep dimension `n`.


> ![http://mgcl.googlecode.com/svn/wiki/image/mg-knot-vector.png](http://mgcl.googlecode.com/svn/wiki/image/mg-knot-vector.png)

  * `MGBPointSeq` - A sequence of control points for a curve.

> ![http://mgcl.googlecode.com/svn/wiki/image/mg-b-point-seq.png](http://mgcl.googlecode.com/svn/wiki/image/mg-b-point-seq.png)

  * `MGSPointSeq` - A sequence of control points for a surface.

> ![http://mgcl.googlecode.com/svn/wiki/image/mg-s-point-seq.png](http://mgcl.googlecode.com/svn/wiki/image/mg-s-point-seq.png)

# Serialization Classes #

There are two serialization classes in MGCL:

  * `MGIfstream` - Input file stream.
  * `MGOfstream` - Output file stream.

These classes support the following functions:

  * Input or output for objects of all classes of geometry and topology.
  * Big endian style binary data.
  * Prevent the same objects from writing or reading twice to the stream.