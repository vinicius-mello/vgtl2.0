#ifndef VGTL_MT_HPP
#define VGTL_MT_HPP

#include <queue>
#include <vector>

/*! \file
 * \brief Generic functions related to MT Concept
 */

namespace vgtl {
	
	using std::queue;
	using std::vector;

	/*! \addtogroup top 
	 * @{
	 */

	//! Tests if a simplex belongs to the initial triangulation
	/*! Extends the \c is_initial requirement to all simplices, not
	 * only vertices.
	 */
	template <class T, dim_t k>
	typename disable_if<vertex_c<T,k>,bool>::type
	is_initial(const T& t,
		       	 Simplex(T,k) s) {
		array<Vertex(T),k+1> vs;
		vertices(t,s,vs);
		for(dim_t i=0; i<=k; ++i) if(!is_initial(t,vs[i])) return false;
		return true;
	}
		
	template <class T>
	inline
	bool was_subdivided(const T& t, SplitSimplex(T) s) {
		return !empty(split_vertex(t,s));
	}

	template <class T>
	inline
	bool is_initial(const T& t, Vertex(T) v) {
		return empty(split_simplex(t,v));
	}

	template <class T>
	inline
	bool has_children(const T& t, Cell(T) c) {
		return !empty(child(t,c,0));
	}

	template <class T>
	inline
	bool has_parent(const T& t, Cell(T) c) {
		return !empty(parent(t,c));
	}

	// A la T. Lewiner
	template <class T>
  class mt_cell_iterator
  {
  public :
    mt_cell_iterator( const T &_t ) : t(_t) { reset() ; }
    ~mt_cell_iterator() {}
    void reset() { 
			empty_set(_curr);
			for(typename vector<Cell(T)>::const_iterator i=initial_cells(t).begin();
				i!=initial_cells(t).end();++i ) _que.push(*i);
		}
		bool next() {
			while( !_que.empty() ) {
				_curr = _que.front();
				_que.pop();
				if(is_current(t,_curr)) {
					return true;
				} else if(has_children(t,_curr)) {
          for(int i=0;i<=SplitDim(T);++i) _que.push(child(t,_curr,i));
        }
			}
			return false;
		}
    Cell(T) curr() { return _curr; }
    void push(Cell(T) f) { _que.push(f); }
    Cell(T) operator *() { return curr(); }
    bool  operator ++() { return next() ; }
  private :
    const T &t;
    queue<Cell(T)> _que;
    Cell(T) _curr;
  };

	/*! @}*/
}

#endif // VGTL_MT_HPP
