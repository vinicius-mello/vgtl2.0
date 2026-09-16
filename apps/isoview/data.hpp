#ifndef APPS_MESH_HPP
#define APPS_MESH_HPP

#include <vector>
//#include <vgtl/top/model/lmt.hpp>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/top/model/lc.hpp>
#include <vgtl/top/euclidean.hpp>
#include <vgtl/alg/point.hpp>
#include <vgtl/alg/packed_point.hpp>
#include <vgtl/cg/GL.hpp>

typedef vgtl::point<3,double> point3;
typedef vgtl::packed_point<3,double> ppoint3;
typedef vgtl::vec<3,double> vec3;
typedef vgtl::vec<4,double> vec4;

namespace vgtl {
				
	template <vgtl::dim_t k>
	struct mesh_data;

}

typedef vgtl::nmt<3,1,vgtl::mesh_data> T;
//typedef vgtl::lmt<3,1,vgtl::mesh_data> T;

namespace vgtl {

	using std::vector;
	
	template <dim_t k>
	struct mesh_data {
	};
	
	template <>
  struct mesh_data<0> {
    mesh_data() : value(1.0) {}
		point3 p;
    double value;
	};

	template <>
	struct mesh_data<2> {
    mesh_data() : x(0), err(0) {par[0]=0.3333;par[1]=0.3333;}
    point<2,double> par;
    double * x;
    double err;
		vec3 normal;
	};

  template <>
  struct mesh_data<1> {
    mesh_data() : x(0) {par[0]=0.5;}
    point<1,double> par;
    double * x;
    //vec<3,double> n;
  };

	template <>
	struct mesh_data<3> {
    mesh_data() : x(0), vtx(0), nml(0) {par[0]=0.25;par[1]=0.25;par[2]=0.25;}
    point<3,double> par;
    double * x;
    double err;
    point<3,double> * vtx;
    vec<3,double> * nml;
	};

	template <>
	struct mesh_data<4> {
		vector<Cell(T)> ini_cv;
	};

  inline
  double scalar_value(const T& t, Vertex(T) v) {
    return attr(t,v)->value;
  }

  inline
  void scalar_value_set(T& t, Vertex(T) v, double s) {
    attr(t,v)->value=s;
  }

  inline
  int signal(const T& t, Vertex(T) v) {
    return (attr(t,v)->value<0) ? -1 : 1;
  }

  inline
  void signal_set(T& t, Vertex(T) v, int s) {
    attr(t,v)->value=abs(attr(t,v)->value)*((s<0) ? -1.0 : 1.0);
  }

  inline
  double * signal_ptr(T& t, Vertex(T) v) {
    return &attr(t,v)->value;
  }

  inline
	point3 euclidean_point(const T& t, Vertex(T) v) {
		return attr(t,v)->p;
	}
	
	inline
	void euclidean_point_set(T& t, Vertex(T) v, const point3& p) {
		attr(t,v)->p=p;
	}

  inline
  void euclidean_point_set(T& t, Vertex(T) v, const vgtl::array<int,3>& p) {
    attr(t,v)->p[0]=p[0];
    attr(t,v)->p[1]=p[1];
    attr(t,v)->p[2]=p[2];
  }

  inline
  double cell_error(const T& t, Cell(T) c) {
    return attr(t,c)->err;
  }

  inline
  void cell_error_set(T& t, Cell(T) c, double e) {
    attr(t,c)->err=e;
  }

  inline
  double facet_error(const T& t, Facet(T) f) {
    return attr(t,f)->err;
  }

  inline
  void facet_error_set(T& t, Facet(T) f, double e) {
    attr(t,f)->err=e;
  }

  inline
  point<3,double> * cell_vertices(const T& t, Cell(T) c) {
    return attr(t,c)->vtx;
  }

  inline
  void cell_vertices_set(T& t, Cell(T) c, point<3,double> * vtx) {
    attr(t,c)->vtx=vtx;
  }

  inline
  vec<3,double> * cell_normals(const T& t, Cell(T) c) {
    return attr(t,c)->nml;
  }

  inline
  void cell_normals_set(T& t, Cell(T) c, vec<3,double> * nml) {
    attr(t,c)->nml=nml;
  }

  inline
	vec3 normal(const T& t, Facet(T) f) {
	/*	vgtl::array<point3,3> p;
		euclidean_points(t,f,p);
		vec3 v=cross(p[1]-p[0],p[2]-p[0]);
		normalize(v);
		return v;*/
		return attr(t,f)->normal;
	}

  template <class T, dim_t k>
  inline
  double * diffpars_x(const T& t, Simplex(T,k) s) {
    return attr(t,s)->x;
  }

  template <class T, dim_t k>
  inline
  void diffpars_x_set(T& t, Simplex(T,k) s, double * x) {
    if(x==0) {
      double sum=0;
      for(dim_t i=0;i<=k;++i) sum+=diffpars_x(t,s)[i];
      for(dim_t i=0;i<k;++i) attr(t,s)->par[i]=diffpars_x(t,s)[i]/sum;
    } else {
      x[k]=1;
      for(dim_t i=0;i<k;++i) {
        x[i]=attr(t,s)->par[i];
        x[k]-=x[i];
      }
    }
    attr(t,s)->x=x;
  }

  template <class T, dim_t k>
  inline
  point<k+1,double> diffpars(const T& t, Simplex(T,k) s) {
    point<k+1,double> p;
    p[k]=1;
    if(diffpars_x(t,s)==0) {
      for(dim_t i=0;i<k;++i) {
        p[i]=attr(t,s)->par[i];
        p[k]-=p[i];
      }
    } else {
      double sum=0;
      for(dim_t i=0;i<=k;++i) sum+=diffpars_x(t,s)[i];
      for(dim_t i=0;i<k;++i) {
        p[i]=diffpars_x(t,s)[i]/sum;
        p[k]-=p[i];
      }
    }
    return p;
  }

  template <class T, dim_t k>
  inline
  void diffpars_set(T& t, Simplex(T,k) s, const point<k+1,double>& p) {
    for(dim_t i=0;i<k;++i) attr(t,s)->par[i]=p[i];
  }

  template <class T, dim_t k>
  inline
  double * diffpars_ptr(T& t, Simplex(T,k) s) {
    return &attr(t,s)->par[0];
  }

  inline
	vector<Cell(T)>& initial_cells(T& t) {
		return t.ini_cv;
	}
	
  inline
	const vector<Cell(T)>& initial_cells(const T& t) {
		return t.ini_cv;
	}
	
	inline
	void normal_set(T& t, Facet(T) f, const vec3& v) {
		attr(t,f)->normal=v;
	}

	struct new_simplices {
    new_simplices()  {}
		void apply(T& t, Simplex(T,0) s) {
      signal_set(t,s,1.0);

			/*Edge(T) e=split_simplex(t,s);
			if(empty(e)) return;
			euclidean_ppoint_set(t,s,mediator(euclidean_ppoint(t,face_op(t,e,1)),
											euclidean_ppoint(t,face_op(t,e,0))));*/
		}
		void apply(T& t, Simplex(T,2) s) { 
			vgtl::array<point3,3> p;
			euclidean_points(t,s,p);
			vec3 v=cross(p[1]-p[0],p[2]-p[0]);
			normalize(v);
			normal_set(t,s,v);
		}
		void apply(T& t, Simplex(T,3) s) { 
		}
		void apply(T& t, Simplex(T,1) s) { 
		}
		Vertex(T) add_edge_vertex(T& t, Edge(T) e) {
			Vertex(T) v=add(t);
			vgtl::array<point3,2> ps;
			euclidean_points(t,e,ps);
			euclidean_point_set(t,v,barycenter(ps));
			return v;
		}
		Vertex(T) barycentric_vertex(T& t, Edge(T) e) {
			Vertex(T) v=add(t);
			vgtl::array<point3,2> ps;
			euclidean_points(t,e,ps);
			euclidean_point_set(t,v,barycenter(ps));
			return v;
		}
		Vertex(T) barycentric_vertex(T& t, Facet(T) f) {
			Vertex(T) v=add(t);
			vgtl::array<point3,3> ps;
			euclidean_points(t,f,ps);
			euclidean_point_set(t,v,barycenter(ps));
			return v;
		}
		Vertex(T) barycentric_vertex(T& t, Cell(T) c) {
			Vertex(T) v=add(t);
			vgtl::array<point3,4> ps;
			euclidean_points(t,c,ps);
			euclidean_point_set(t,v,barycenter(ps));
			return v;
		}
	};
	
}

#endif //APPS_MESH_HPP
