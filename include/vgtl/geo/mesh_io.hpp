#ifndef VGTL_MESH_IO_HPP
#define VGTL_MESH_IO_HPP

/*! \file
 * \brief Mesh IO
 */

#include <fstream>
#include <vector>
#include <map>
#include <vgtl/alg/point.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/euclidean.hpp>
#include <vgtl/top/add_simplex.hpp>

namespace vgtl {
	
	using std::vector;
	using std::map;
	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::endl;

	/*! \addtogroup top 
	 * @{
	 */

	template <class T, class Apply>
	void 
	read_mesh(ifstream& in, T& t, Apply& app, vector<Vertex(T)>& vv);

	template <class T, class Apply>
	void 
	read_mesh(ifstream& in, T& t, Apply& app) {
		vector<Vertex(T)> vv;
		read_mesh(in,t,app,vv);
	}

	template <class T, class Apply>
	void 
	read_mesh(ifstream& in, T& t, Apply& app, vector<Vertex(T)>& vv) {
		complex_builder<T> cb(t);
		int nv;
		in>>nv;
		for(int i=0;i<nv;++i) {
			vv.push_back(add(t));
			point<Dim(T)> p;
			in>>p;
			euclidean_point_set(t,vv[i],p);
			app.apply(t,vv[i]);
		}
		int nc;
		in>>nc;
		for(int i=0;i<nc;++i) {
			array<int,Dim(T)+1> a;
			in>>a;
			array<Vertex(T),Dim(T)+1> vs;
			for(int j=0;j<=Dim(T);++j) vs[j]=vv[a[j]];
			Cell(T) c=add(cb,vs);
			array<point<Dim(T)>,Dim(T)+1> ps;
			euclidean_points(t,c,ps);
			orientation_set(t,c,((bracket(ps)<0)?-1:1));
		}
		apply(cb,app);
	}
	
	template <class T>
	void 
	write_mesh(ofstream& out, const T& t, map<Vertex(T),int>& mv);
					
	template <class T>
	void 
	write_mesh(ofstream& out, const T& t) {
		map<Vertex(T),int> mv;
		write_mesh(out,t,mv);
	}

	template <class T>
	void 
	write_mesh(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		int nv=0;
		Vertex_it(T) vi, vend;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			if(!is_current(t,v)) continue;
			mv[v]=nv;++nv;
		}
		out<<nv<<endl;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			out<<euclidean_point(t,v)<<endl;
		}
		int nc=0;
		Cell_it(T) ci, cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			++nc;
		}
		out<<nc<<endl;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			array<Vertex(T),Dim(T)+1> vs;
			vertices(t,c,vs);
			for(int j=0;j<=Dim(T);++j) {
				out<<mv[vs[j]];
				if(j!=Dim(T)) out<<" ";
				else out<<endl;
			}
		}
	}

	template <class T, class Apply>
	void 
	read_mesh_iso(ifstream& in, T& t, Apply& app,vector<Vertex(T)>& vv);

	template <class T, class Apply>
	void 
	read_mesh_iso(ifstream& in, T& t, Apply& app) {
		vector<Vertex(T)> vv;
		read_mesh_iso(in,t,app,vv);
	}
	
	template <class T, class Apply>
	void 
	read_mesh_iso(ifstream& in, T& t, Apply& app,vector<Vertex(T)>& vv) {
		complex_builder<T> cb(t);
		int nv;
		in>>nv;
		for(int i=0;i<nv;++i) {
			vv.push_back(add(t));
			point<Dim(T)> p;
			in>>p;
			euclidean_point_set(t,vv[i],p);
			app.apply(t,vv[i]);
			double s;
			in>>s;
			scalar_value_set(t,vv[i],s);
		}
		int nc;
		in>>nc;
		for(int i=0;i<nc;++i) {
			array<int,Dim(T)+1> a;
			in>>a;
			sort(a.begin(),a.end());
			array<Vertex(T),Dim(T)+1> vs;
			for(int j=0;j<=Dim(T);++j) vs[j]=vv[a[j]];
			Cell(T) c=add(cb,vs);
			array<point<Dim(T)>,Dim(T)+1> ps;
			euclidean_points(t,c,ps);
			orientation_set(t,c,((bracket(ps)<0)?-1:1));
		}
		apply(cb,app);
	}
	
	template <class T>
	void 
	write_mesh_iso(ofstream& out, const T& t, map<Vertex(T),int>& mv);
	
	template <class T>
	void 
	write_mesh_iso(ofstream& out, const T& t) {
		map<Vertex(T),int> mv;
		write_mesh_iso(out,t,mv);
	}
	
	template <class T>
	void 
	write_mesh_iso(ofstream& out, const T& t, map<Vertex(T),int>& mv) {
		int nv=0;
		Vertex_it(T) vi, vend;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			if(!is_current(t,v)) continue;
			mv[v]=nv;++nv;
		}
		out<<nv<<endl;
		for(simplices(t,vi,vend);vi!=vend;++vi) {
			Vertex(T) v=*vi;
			out<<euclidean_point(t,v)<<" ";
			out<<scalar_value(t,v)<<endl;
		}
		int nc=0;
		Cell_it(T) ci, cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			++nc;
		}
		out<<nc<<endl;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			array<Vertex(T),Dim(T)+1> vs;
			vertices(t,c,vs);
			for(int j=0;j<=Dim(T);++j) {
				out<<mv[vs[j]];
				if(j!=Dim(T)) out<<" ";
				else out<<endl;
			}
		}
	}

}

#endif // VGTL_MESH_IO_HPP
