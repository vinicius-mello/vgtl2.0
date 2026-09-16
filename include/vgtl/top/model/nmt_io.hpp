#ifndef VGTL_MODEL_NMT_IO_HPP
#define VGTL_MODEL_NMT_IO_HPP

/*! \file
 * \brief Models a Retangular Triangulation
 */

#include <vgtl/top/model/nmt.hpp>
#include <vgtl/top/vertices_inv.hpp>
#include <vgtl/top/maubach.hpp>
#include <fstream>

namespace vgtl {

	using std::ifstream;
	using std::ofstream;
	using std::ios;
	using std::endl;
	using std::vector;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename enable_if<equal_c<k,0>,void>::type
	write_simplex(ofstream& out, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k> v) {
		out<<v.desc<<" ";
	}
						
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename disable_if<equal_c<k,0>,void>::type
	write_simplex(ofstream& out, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k> s) {
		array<simplex_descriptor<nmt<d,sd,Data>,0>,k+1> vs;
		vertices(t,s,vs);
		for(int i=0;i<=k;++i) write_simplex(out,t,vs[i]);
	}
						
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename enable_if<equal_c<k,0>,void>::type
	read_simplex(ifstream& in, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k>& v);
					
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename disable_if<equal_c<k,0>,void>::type
	read_simplex(ifstream& in, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k>& s);

	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename enable_if<equal_c<k,0>,void>::type
	read_simplex(ifstream& in, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k>& v) {
		in>>v.desc;
	}
						
	template <dim_t d, dim_t sd, dim_t k,
	          template <dim_t> class Data>
	typename disable_if<equal_c<k,0>,void>::type
	read_simplex(ifstream& in, const nmt<d,sd,Data>& t,
								simplex_descriptor<nmt<d,sd,Data>,k>& s) {
		array<simplex_descriptor<nmt<d,sd,Data>,0>,k+1> vs;
		for(int i=0;i<=k;++i) read_simplex(in,t,vs[i]);
		s=vertices_inv(t,vs);
	}
						
	template <dim_t d, dim_t sd,
	          template <dim_t> class Data>
	void write_splits(ofstream& out, const nmt<d,sd,Data>& t) {
		out<<d<<endl;
		int nsplits=0;
		simplex_iterator<nmt<d,sd,Data>,0> i,end;
		for(simplices(t,i,end);i!=end;++i) {
			if(!is_initial(t,*i)) nsplits++;
		}
		out<<nsplits<<endl;
		for(simplices(t,i,end);i!=end;++i) {
			simplex_descriptor<nmt<d,sd,Data>,0> v=*i;
			simplex_descriptor<nmt<d,sd,Data>,sd> s=split_simplex(t,v);
			if(empty(s)) continue;
			write_simplex(out,t,s);
			out<<endl;
		}
	}

	template <dim_t d, 
	          template <dim_t> class Data,
						class Apply>
	void read_splits(ifstream& in, nmt<d,1,Data>& t, Apply& app) {
		int dim;
		in>>dim;
		if(dim!=d) return; //throw exception
		int nsplits;
		in>>nsplits;
		do_nothing dn;
		for(int i=0;i<nsplits;++i) {
			simplex_descriptor<nmt<d,1,Data>,1> e;
			read_simplex(in,t,e);
			maubach_subdivide(t,top_cell(t,e),app);
		}
	}

	/*! @} */

}

#endif // VGTL_MODEL_NMT_IO_HPP
