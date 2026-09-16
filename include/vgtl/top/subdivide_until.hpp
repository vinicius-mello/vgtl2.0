#ifndef VGTL_SUBDIVIDE_UNTIL_HPP
#define VGTL_SUBDIVIDE_UNTIL_HPP

/*! \file
 * \brief Subdivide until 
 */

#include <vgtl/top/maubach.hpp>

namespace vgtl {

	using std::queue;
	
	/*! \addtogroup top 
	 * @{
	 */
	
	template <class T>
	struct put_children_on_queue {
		queue<Cell(T)> qu;
		put_children_on_queue() : qu() {}
		template <dim_t k>
		void
		apply(T& t, Simplex(T,k) s, 
			typename disable_if<cell_c<T,k>,void *>::type =0) {
		}
		template <dim_t k>
		void
		apply(T& t, Simplex(T,k) s, 
			typename enable_if<cell_c<T,k>,void *>::type =0) {
			qu.push(child(t,s,0));
			qu.push(child(t,s,1));
		}
	};
	
	//! 
	template <class T, class ApplyNew, class Condition>
	void
	maubach_subdivide_until(T& t, Condition& cond, ApplyNew& app)  {
		put_children_on_queue<T> pq;
		Cell_it(T) ci, cend;
		for(simplices(t,ci,cend);ci!=cend;++ci) {
			Cell(T) c=*ci;
			if(!is_current(t,c)) continue;
			pq.qu.push(c);
		}
		while(!pq.qu.empty()) {
			Cell(T) c=pq.qu.front();
			pq.qu.pop();
			if(!is_current(t,c)) continue;
			if(cond(t,c)) continue;
			maubach_subdivide(t,c,app,pq);
		}
	}

	template <class T>
	struct level_greater_or_equal {
		int n;
		level_greater_or_equal(int _n) : n(_n) {}
		bool operator()(const T& t, Cell(T) c) {
			return level(t,c)>=n;
		}
	};

	template <class T, class ApplyNew>
	void
	maubach_subdivide_to_level(T& t, int n, ApplyNew& app)  {
		level_greater_or_equal<T> lv(n);
		maubach_subdivide_until(t,lv,app);
	}

	/*! @} */

}

#endif // VGTL_SUBDIVIDE_UNTIL_HPP
