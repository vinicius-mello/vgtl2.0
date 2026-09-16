/*! 
\mainpage 

The VGTL library is an experimental c++ template library for representation
of triangulations and multi-triangulations in arbitrary dimensions. 


\section Installation

VGTL is a template library, so is not necessary to build a lib file,
just include the most suitable header files for the task at hand. 
A good c++ compiler is required (gcc >= 3.4 and Visual C++ >= 2005 
must work). The Boost Library is also required, but just Boost.Array,
Boost.MultiArray, Boost.Utility, Boost.Random and Boost.LexicalCast.
GSL and GLUT are used in sample applications. To build the applications,
change the *DIR variables in the file configure.mak and run make in
top directory.

\section Overview

Our main goal is to represent triangulations and multi-triangulations. 
A multi-triangulation is a kind of multi-resolution triangulation. The 
design goals of VGTL are two: decouple algorithms from data structures
and be generic with respect to dimension, i.e., triangle and tetrahedra 
meshes, and simplicial meshes of higher dimension, have the same interface.

We implement four data structures for meshes: a list-based structure 
for single resolution meshes (\c lc<d,Data>), a vector-based structure
for single resolution meshes (\c vc<d,Data>), a list-based structure
for multi-resolution meshes (\c lmt<d,sd,Data>) and a vector-based structure
for multi-resolution meshes (\c nmt<d,sd,Data>). As usual, list-based 
structures are usefull when deletions are needed, otherwise vector-
based structures are more efficient. 

All data structures above are parameterized by dimension \c d and 
by an extra storage class \c Data. Suppose that in your application 
we need to represent a triangle mesh with a point and normal for each
vertex, a spring coefficient for each edge and a color for each triangle.
You could do like this:

\code 
	using namespace vgtl;

	//mydata is parameterized by dimension
	template <dim_t n>
	struct mydata {
	};

	//We use template specialization to set the attributes
	//for each dimension
	template <>
	struct mydata<0> {
		vec<3,double> normal;
		point<3,double> point;
	};

	template <>
	struct mydata<1> {
		double k;
	};

	template <>
	struct mydata<2> {
		vec<3,double> color;
	};
 
	typedef lc<2,mydata> T;

	T t;
\endcode

When \c t is instantiated, some template magic happens such that 
\c T contains 3 STL lists, one for each dimension. Se usassemos 
o container \c vc, teriamos 3 vetores STL.

Os elementos das listas podem ser percorridos usando iteradores:

\code 

	\\The iterator type for n-dimensional simplices is Simplex_it(T,n)
	Simplex_it(T,1) ei, eend;
	for(simplices(t,ei,eend);ei!=eend;++ei) {
		\\The type of n-dimensional simplices is Simplex(T,n)
		Simplex(T,1) e=*ei;
		\\do something with e
	}

	\\Some macros are helpfull:
	\\Vertex(T)=Simplex(T,1)
	\\Edge(T)=Simplex(T,1)
	\\Facet(T)=Simplex(T,n-1)
	\\Cell(T)=Simplex(T,n)

\endcode
There are a few low-level functions to retrive incidence relations:

\code 
	
	Cell(T) c1,c2;
	\\Suppose that t=c1+c2 where c1=<a,b,c,d> and c2=<a',b,c,d>	
	\\the face of c1 obtained by removing element in position:
	Facet(T) f=face_op(t,c1,1); 
	\\now f=<a,c,d>
	Edge(T) e=face_op(t,f,0);
	\\now e=<c,d>

	

	


\endcode
*/

