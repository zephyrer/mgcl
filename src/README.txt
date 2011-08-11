======================================================================
MGCL V8.1 README
======================================================================

Configuration
=============

All libary files in this web page are built with Microsoft
Visual Studio .NET 2003 (msvc71).

MGCL has the following configurations:

 * Release (): Multi-threaded DLL, release mode
 * Debug (d): Debug Multi-threaded DLL, debug mode
 * ReleaseWithoutMFC (nomfc): Multi-threaded DLL, release mode
 * DebugWithoutMFC (nomfc-d): Multi-threaded DLL, debug mode

'nomfc' configurations are built with pre-processor symbol
'MGCL_NO_MFC' defined, which means the subpackages do not use
mgGL and mgIges. This is intended to save build time and
target module size by omiting enormous MFC object linking.

Settings
========

All the compiled libraries are placed in 'lib' directory.
To link MGCL, you may need to:

 * add lib/ to LIB search path
 * add include/ to INCLUDE search path
 * and link mgclXXX.lib (opengl32.lib and gdiplus.lib)

Note that if you use 'nomfc' libraries, you must define
pre-processor symbol 'MGCL_NO_MFC' in your project.

