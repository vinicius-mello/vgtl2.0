#ifndef VGTL_SIMULATED_ANNEALING_HPP
#define VGTL_SIMULATED_ANNEALING_HPP

#include <cmath>

namespace vgtl {

  using std::exp;
  using std::pow;
  using std::min;

	template <class RandomGenerator>
	class boltzmann_criteria {
    RandomGenerator r;
		double K;
		public:
    boltzmann_criteria(double _K, RandomGenerator& _r) : K(_K), r(_r) {}
		bool operator()(double chg, double T) {
			return r() < exp(-chg/(K*T));
		}
	};

	template <class RandomGenerator>
	class tsallis_criteria {
    RandomGenerator r;
		double K;
		double _q;
		public:
    tsallis_criteria(double q, double _K, RandomGenerator& _r) : _q(1.0-q),K(_K),r(_r) {}
		bool operator()(double chg, double T) {
			double ap=min(1.0,pow(1.0-(1.0/(K*T))*_q*chg,1.0/_q));
			return r() < ap;
		}
	};

	class geometric_cooling_schedule {
		double Ti, Tf;
		double cooling_factor;
		double _T;
		public:
		geometric_cooling_schedule(double _Ti, double _Tf, double cf) : 
						Ti(_Ti), Tf(_Tf), cooling_factor(cf) {}
		double T() const {
			return _T;
		}
		void initialize_temperature() {
			_T=Ti;
		}
		void decrease_temperature() {
			_T*=cooling_factor;
		}
		bool freezed() const {
			return _T<=Tf;
		}
	};
	
	class simple_tries_handler {
		int i;
		int max_tries;
		public:
		simple_tries_handler(int m) : max_tries(m) {}
		void start() {
			i=0;
		}
		bool enough() {
			++i;
			return i>=max_tries;
		}
		void local_hit() {
		}
		void global_hit() {
		}
		void escape() {
		}
	};

	template <class ConfigurationGenerator,
		class CoolingSchedule, class TriesHandler, class AcceptanceCriteria>
	double simulated_annealing(ConfigurationGenerator& cg,
									CoolingSchedule& cs,
									TriesHandler& tr,
									AcceptanceCriteria& ac) {
		cs.initialize_temperature();
  	double min_energy=cg.energy();
  	double cur_energy=min_energy;
  	do {
			tr.start();
			do {
				double chg_energy=cg.change_configuration();
    		if(cur_energy+chg_energy<min_energy) {
					cg.save_configuration();
					min_energy=cur_energy+chg_energy;
					tr.global_hit();
				}
				if((chg_energy<0)) { /*new-old<0*/
					cg.accept_configuration();
					cur_energy=cur_energy+chg_energy;
					tr.local_hit();
				} else if(ac(chg_energy,cs.T())) {
					cg.accept_configuration();
					cur_energy=cur_energy+chg_energy;
					tr.escape();
				}
			} while(!tr.enough());
			cs.decrease_temperature();
		} while(!cs.freezed());
		cg.restore_configuration();
		return min_energy;
	}

}

#endif // VGTL_SIMULATED_ANNEALING_HPP
