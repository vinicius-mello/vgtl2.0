#ifndef VGTL_TUBULAR_HPP
#define VGTL_TUBULAR_HPP

#include <map>
#include <vgtl/top/fill.hpp>
#include <vgtl/top/iso/iso.hpp>

namespace vgtl {

	namespace iso {

		template <class Iso, dim_t k>
		bool cross(const Iso& t, Simplex(Iso,k) s) {
			int m,p;
			pair_tie(p,m)=dim(t,s);
			if(m==-1 || p==-1) return false;
			else return true;
		}

		template <class Iso, dim_t kk>
		bool cross(const Iso& t, const array<Vertex(Iso),kk>& vs) {
			int m,p;
			pair_tie(p,m)=dim(t,vs);
			if(m==-1 || p==-1) return false;
			else return true;
		}

		template <class Iso, class CellIt>
		void tubular(complex_buffer<Iso>& cb, CellIt begin, CellIt end) {
			for(CellIt ci=begin;ci!=end;++ci) {
				Cell(Iso) ce=*ci;
				if(!is_current(cb.t,ce)) continue;
				if(cross(cb.t,ce)) fill(cb,ce);
			}
		}

		template <class Iso>
		void tubular(complex_buffer<Iso>& cb) {
			Cell_it(Iso) ci, cend;
			for(simplices(cb.t,ci,cend);ci!=cend;++ci) {
				Cell(Iso) ce=*ci;
				if(!is_current(cb.t,ce)) continue;
				if(cross(cb.t,ce)) fill(cb,ce);
			}
		}
	
		template <class Iso>
		void tubular_v(complex_buffer<Iso>& cb) {
			tubular(cb);
			Vertex_it(Iso) vi, vend;
			for(simplices(cb.t,vi,vend);vi!=vend;++vi) {
				Vertex(Iso) v=*vi;
				if(!is_current(cb.t,v)) continue;
				put(cb,v);
			}
		}
	
		template <class Iso>
		void tubular_vc(complex_buffer<Iso>& cb) {
			Cell_it(Iso) ci, cend;
			for(simplices(cb.t,ci,cend);ci!=cend;++ci) {
				Cell(Iso) ce=*ci;
				if(!is_current(cb.t,ce)) continue;
				if(cross(cb.t,ce)) fill_vc(cb,ce);
			}
		}
	
	}

}

#endif // VGTL_TUBULAR_HPP
