// Build (from this directory): g++ -O2 -std=c++17 -I../../../include -I$HOME/code/lpt mesh_mem.cpp $HOME/code/lpt/lpt.cpp -o mesh_mem && ./mesh_mem nmt 12 && ./mesh_mem glpt 12
// Mesh-storage microbenchmark: uniform Maubach refinement of Gaifullin's
// 108-cell CP^2 seed to depth D, pure topology (no coordinates, no curve),
// with (a) vgtl's explicit nmt<4> and (b) the pointerless glpt_tree.
// Usage: mesh_mem nmt|glpt D
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <vgtl/top/model/nmt.hpp>
#include <vgtl/utl/array_cons.hpp>
#include <vgtl/top/add_simplex.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/top/do_nothing.hpp>
#include "glpt_tree.hpp"

using namespace vgtl;
typedef vgtl::nmt<4> T;

long rss_kb(const char* key) {
	FILE* f=fopen("/proc/self/status","r"); char line[256]; long v=-1;
	while(fgets(line,sizeof line,f)) if(!strncmp(line,key,strlen(key))) { v=atol(line+strlen(key)+1); break; }
	fclose(f); return v;
}

struct app : do_nothing {
	using do_nothing::apply;
	Vertex(T) add_edge_vertex(T& t, Edge(T) e) { return add(t); }
};

int main(int argc,char**argv) {
	std::string mode=argv[1]; int D=atoi(argv[2]);
	long base=rss_kb("VmRSS:");
	if(mode=="nmt") {
		static T t;
		std::vector<Vertex(T)> V(15);
		for(int i=0;i<15;++i) V[i]=add(t);
		{
			complex_builder<T> cb(t);
			for(int c=0;c<GLPT_NCELLS;++c) {
				vgtl::array<Vertex(T),5> vs;
				for(int k=0;k<5;++k) vs[k]=V[glpt_gaifullin_cells[c][k]];
				Cell(T) cv=add(cb,vs); orientation_set(t,cv,1); level_set(t,cv,0);
			}
		}
		for(int pass=0;;++pass) {
			std::vector<Cell(T)> todo;
			Cell_it(T) i,e; for(simplices(t,i,e);i!=e;++i) if(is_current(t,*i) && level(t,*i)<D) todo.push_back(*i);
			if(todo.empty()) break;
			for(size_t k=0;k<todo.size();++k) if(is_current(t,todo[k]) && level(t,todo[k])<D) { app a; maubach_subdivide(t,todo[k],a); }
		}
		long cur=0; { Cell_it(T) i,e; for(simplices(t,i,e);i!=e;++i) if(is_current(t,*i)) ++cur; }
		size_t n[5]={simplex_container<0>(t).size(),simplex_container<1>(t).size(),simplex_container<2>(t).size(),simplex_container<3>(t).size(),simplex_container<4>(t).size()};
		size_t sz[5]={sizeof(ms<T,0>),sizeof(ms<T,1>),sizeof(ms<T,2>),sizeof(ms<T,3>),sizeof(ms<T,4>)};
		double bytes=0; for(int k=0;k<5;++k) bytes+=double(n[k])*sz[k];
		printf("nmt D=%d current_cells=%ld stored: v=%zu e=%zu t=%zu f=%zu c=%zu  recsize=%zu/%zu/%zu/%zu/%zu B\n",D,cur,n[0],n[1],n[2],n[3],n[4],sz[0],sz[1],sz[2],sz[3],sz[4]);
		printf("  records=%.1f MB (%.1f B/current cell)  RSS=%ld KB  peak=%ld KB  (%.1f B/cell over baseline)\n",bytes/1e6,bytes/cur,rss_kb("VmRSS:"),rss_kb("VmHWM:"),1024.0*(rss_kb("VmHWM:")-base)/cur);
	} else {
		glpt_tree tr; tr.seed_all_roots();
		for(;;) {
			std::vector<glpt> todo;
			tr.for_each_leaf([&](const glpt& c){ if(c.simplex_level()<D) todo.push_back(c); });
			if(todo.empty()) break;
			for(auto& c:todo) tr.bisect(c);
			tr.clear_recent();
		}
		printf("glpt D=%d leaves=%zu buckets=%zu load=%.2f table=%.1f MB (%.1f B/leaf)  RSS=%ld KB peak=%ld KB (%.1f B/leaf over baseline)\n",D,tr.leaf_count(),tr.bucket_count(),tr.load_factor(),tr.bucket_count()*8/1e6,tr.bucket_count()*8.0/tr.leaf_count(),rss_kb("VmRSS:"),rss_kb("VmHWM:"),1024.0*(rss_kb("VmHWM:")-base)/tr.leaf_count());
	}
}
