#ifndef CG_ARCBALL_HPP
#define CG_ARCBALL_HPP

/*! \file
 * \brief Simple Arcball
 */

#include <vgtl/alg/vec.hpp>
#include <vgtl/alg/quat.hpp>

namespace vgtl {

	/*! \addtogroup cg 
	 * @{
	 */
	
	//! Simple Arcball
	/*!
	 * \b Example:
	 * \include arcball.cpp 
	 * \b Output:
	 * \include arcball.out
	 */
  class arcball {
      vec<2> now, down;
      quat<double> qnow, qdrag;
      int w_, h_;
      vec<3> project_to_sphere(const vec<2>& a) {
          double r=dot(a,a);
          vec<3> b;
          b[0]=a[0];b[1]=a[1];b[2]=0;
          if(r>1)
              normalize(b);
          else
              b[2]=sqrt(1-r);
          return b;
      }
  public:
      arcball() : qnow(1.0, 0.0, 0.0, 0.0), qdrag(1.0, 0.0, 0.0, 0.0) {
      }
      void init() {
          qnow.assign(0.0);
          qdrag.assign(0.0);
          qnow[0]=qdrag[0]=1.0;
      }
      void resize(int w, int h) {
          w_=w;
          h_=h;
      }

      void start_motion(int x, int y, int w, int h) {
          down[0]=2.0*x/w-1.0;
          down[1]=1.0-2.0*y/h;
      }
      void start_motion(int x, int y) {
          start_motion(x,y,w_,h_);
      }

      void stop_motion(int x,int y) {
      }
      void motion(int x, int y, int w, int h) {
          vec<3> from, to;
          now[0]=2.0*x/w-1.0;
          now[1]=1.0-2.0*y/h;
          from=project_to_sphere(down);
          to=project_to_sphere(now);
          qdrag=vgtl::rotation(from,to);
          qnow=qdrag*qnow;
          down=now;
      }
      void motion(int x, int y) {
          arcball::motion(x,y,w_,h_);
      }

      quat<double> rotation() const {
          return qnow;
      }
      void set_rotation(const quat<double>& q) {
          qnow=q;
      }
  };

	/* @} */

}

#endif //CG_ARCBALL_HPP
