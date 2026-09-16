#ifndef VGTL_CG_HANDLES_HPP
#define VGTL_CG_HANDLES_HPP

#include <string>
#include <map>
#include <vgtl/cg/GL.hpp>
#include <iostream>

using std::string;
using std::map;
using std::cout;
using std::endl;

namespace vgtl {

  class handle_controller;
  class handle {
  protected:
    int id_;
    handle_controller * hc;
  public:
    handle() {}
    virtual int id() {return id_;}
    virtual void move(int x, int y)=0;
    virtual void select()=0;
    virtual void release()=0;
    virtual void draw_proxy()=0;
    virtual void draw()=0;
    virtual ~handle() {};
    friend class handle_controller;
  };

  class handle_controller {
  protected:
    map<string,handle *> handle_by_name;
    map<int,handle *> handle_by_id;
    typedef map<int,handle *> hmap;
    int next_id;
    bool perspective_;
    array<vec<4,double>,4> mv_;
    array<vec<4,double>,4> pr_;
    point<2,int> pos_;
    vec<2,int> size_;
  public:
    handle_controller() {next_id=1;}
    handle * add_handle(const string& name, handle * h) {
      h->id_=next_id;
      h->hc=this;
      next_id++;
      handle_by_name[name]=h;
      handle_by_id[h->id()]=h;
      return h;
    }
    void delete_handle(const string& name) {
      if(handle_by_name.count(name)==1) {
        handle * h=handle_by_name[name];
        handle_by_name.erase(name);
        handle_by_id.erase(h->id());
        delete h;
      }
    }
    void delete_range(const string& first, const string& last) {
      map<string,handle *>::iterator itf, itl;
      itf=handle_by_name.lower_bound(first);
      itl=handle_by_name.upper_bound(last);
      for(map<string,handle *>::iterator i=itf;i!=itl;++i) {
        handle * h=i->second;
        handle_by_id.erase(h->id());
        delete h;
      }
      handle_by_name.erase(itf,itl);
    }
    handle * get(const string& name) {
      if(handle_by_name.count(name)==1) {
        return handle_by_name[name];
      }
      return 0;
    }
    void draw_handles() {
      glGetModelviewMatrix(mv_);
      glGetProjectionMatrix(pr_);
      glGetViewport(pos_,size_);
      for(hmap::iterator i=handle_by_id.begin(); i!=handle_by_id.end();++i) {
        i->second->draw();
      }
    }
    point<3,double> project(const point<3,double>& p) {
      return gluProject(p,mv_,pr_,pos_,size_);
    }
    point<3,double> unproject(const point<3,double>& p) {
      return gluUnProject(p,mv_,pr_,pos_,size_);
    }
    const array<vec<4,double>,4>& modelview() {
      return mv_;
    }
    point<3,double> plane_unproject(const point<2,double>& p,
                                    const array<point<3,double>,3>& ps) {
      array<point<2,double>,3> qs;
      vgtl::project(2,project(ps[0]),qs[0]);
      vgtl::project(2,project(ps[1]),qs[1]);
      vgtl::project(2,project(ps[2]),qs[2]);
      point<3,double> w=barycentric_coordinates(p,qs);
      if(perspective()) {
        point<3,double> z;
        z[0]=mv_[0][2]*ps[0][0]+mv_[1][2]*ps[0][1]+mv_[2][2]*ps[0][2]+mv_[3][2];
        z[1]=mv_[0][2]*ps[1][0]+mv_[1][2]*ps[1][1]+mv_[2][2]*ps[1][2]+mv_[3][2];
        z[2]=mv_[0][2]*ps[2][0]+mv_[1][2]*ps[2][1]+mv_[2][2]*ps[2][2]+mv_[3][2];
        w=_(w[0]/z[0],w[1]/z[1],w[2]/z[2]);
        w=barycentric_project(w);
      }
      return w;
    }
    double line_unproject(const point<2,double>& p,
                          const array<point<3,double>,2>& ps) {
      array<point<2,double>,2> qs;
      vgtl::project(2,project(ps[0]),qs[0]);
      vgtl::project(2,project(ps[1]),qs[1]);
      vec<2,double> v=qs[1]-qs[0];
      double l=sqrt(dot(v,v));
      normalize(v);
      double t=dot(p-qs[0],v)/l;
      if(perspective()) {
        point<2,double> z;
        z[0]=mv_[0][2]*ps[0][0]+mv_[1][2]*ps[0][1]+mv_[2][2]*ps[0][2]+mv_[3][2];
        z[1]=mv_[0][2]*ps[1][0]+mv_[1][2]*ps[1][1]+mv_[2][2]*ps[1][2]+mv_[3][2];
        t=(t*z[0])/(t*z[0]+(1.0-t)*z[1]);
      }
      return t;
    }
    bool tetra_intersect(const point<2,double>& p,
                                    const array<point<3,double>,4>& ps,
                                    point<3,double>& r, double& t) {
      bool inter=false;
      array<point<3,double>,3> fs;
      array<point<2,double>,3> qs;
      for(dim_t j=0;j<4;++j) {
        for(dim_t i=0; i<j; ++i) fs[i]=ps[i];
        for(dim_t i=j+1; i<4; ++i) fs[i-1]=ps[i];
        if(j%2==1) swap(fs[1],fs[2]);
        vgtl::project(2,project(fs[0]),qs[0]);
        vgtl::project(2,project(fs[1]),qs[1]);
        vgtl::project(2,project(fs[2]),qs[2]);
        if(bracket(qs)<0) continue;
        point<3,double> w=barycentric_coordinates(p,qs);
        if((w[0]>=0.0)&&(w[1]>=0.0)&&(w[2]>=0.0)) {
          inter=true;
          if(perspective()) {
            point<3,double> z;
            z[0]=mv_[0][2]*ps[0][0]+mv_[1][2]*ps[0][1]+mv_[2][2]*ps[0][2]+mv_[3][2];
            z[1]=mv_[0][2]*ps[1][0]+mv_[1][2]*ps[1][1]+mv_[2][2]*ps[1][2]+mv_[3][2];
            z[2]=mv_[0][2]*ps[2][0]+mv_[1][2]*ps[2][1]+mv_[2][2]*ps[2][2]+mv_[3][2];
            w=_(w[0]/z[0],w[1]/z[1],w[2]/z[2]);
            w=barycentric_project(w);
          }
          r=barycentric_combination(w,fs);
          point<3,double> a=project(r);
          t=a[2];
        }
      }
      return inter;
    }
    void perspective_on() {
      perspective_=true;
    }
    void perspective_off() {
      perspective_=false;
    }
    bool perspective() {
      return perspective_;
    }
    handle * pick_handle(int x, int y) {
      GLuint select_buffer[256];
      GLint hits;
      GLint viewport[4];
      handle * selected=0;
      glGetIntegerv(GL_VIEWPORT, viewport);
      //array<vec<4,double>,4> pr;
      //glGetProjectionMatrix(pr);
      glSelectBuffer(256, select_buffer);
      glRenderMode(GL_SELECT);
      glInitNames();
      glPushName(0);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluPickMatrix((GLdouble) x,(GLdouble) y,
                    3.0, 3.0, viewport);
      glMultMatrix(pr_);
      
      glMatrixMode(GL_MODELVIEW);
      for(hmap::iterator i=handle_by_id.begin(); i!=handle_by_id.end();++i) {
        glLoadName(i->first);
        i->second->draw_proxy();
      }
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();

      hits=glRenderMode(GL_RENDER);
      if(hits>0) {
        if(handle_by_id.count(select_buffer[3])==1)
          selected=handle_by_id[select_buffer[3]];
      }
      return selected;
    }

  };

  template <typename Scalar>
  class signal_handle : public handle {
  protected:
    point<3,double> p_;
    double size_;
    Scalar * data_;
    GLUquadric * sph;
  public:
    signal_handle() {
      sph=gluNewQuadric();
    }
    signal_handle(const point<3,double>& _p, double _size, Scalar * _data)
      : p_(_p), size_(_size), data_(_data) {
      sph=gluNewQuadric();
    }
    virtual void move(int x, int y) {};
    virtual void select() {
      *data_=-(*data_);
    }

    virtual void release() {
    }
    virtual void draw_proxy() {
      draw_sphere();
    }

    void draw_sphere() {
      glPushMatrix();
      glTranslate(p_);
      gluSphere(sph,size_,8,8);
      glPopMatrix();
    }
    virtual void draw() {
      if((*data_)<0)
        glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Red[0]);
      else
        glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Blue[0]);
      draw_sphere();
    }
    virtual ~signal_handle() {
      gluDeleteQuadric(sph);
    }

  };

  template <typename Scalar>
  class segment_handle : public handle {
  protected:
    point<3,double> p0_;
    point<3,double> p1_;
    double size_;
    Scalar * data_;
    GLUquadric * sph;
    bool restricted_;
  public:
    segment_handle()
      : restricted_(true) {
      sph=gluNewQuadric();
    }
    segment_handle(const point<3,double>& _p0,
                   const point<3,double>& _p1,
                   double _size, Scalar * _data,
                   bool _restricted=true)
      : p0_(_p0), p1_(_p1), size_(_size), data_(_data), restricted_(_restricted) {
      sph=gluNewQuadric();
    }
    virtual void move(int x, int y) {
      point<2,double> p=_((double)x,(double)y);
      double t=hc->line_unproject(p,_(p0_,p1_));
      if(restricted_)
        t= (t<0.0) ? 0.0 : ((t>1.0) ? 1.0 : t);
      data_[0]=(Scalar)t;
    }
    virtual void select() {
    }
    virtual void release() {
    }
    virtual void draw_proxy() {
      draw_sphere();
    }
    void draw_sphere() {
      glPushMatrix();
      glTranslate(p0_+(*data_)*(p1_-p0_));
      gluSphere(sph,size_,8,8);
      glPopMatrix();
    }
    virtual void draw() {
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Green[0]);
      draw_sphere();
      glBegin(GL_LINES);
      glVertex(p0_);
      glVertex(p1_);
      glEnd();
    }
    virtual ~segment_handle() {
      gluDeleteQuadric(sph);
    }

  };

  class triangle_handle : public handle {
  protected:
    point<3,double> p0_;
    point<3,double> p1_;
    point<3,double> p2_;
    double size_;
    double * data_;
    GLUquadric * sph;
    bool restricted_;
  public:
    triangle_handle() : restricted_(true) {
      sph=gluNewQuadric();
    }
    triangle_handle(const point<3,double>& _p0,
                   const point<3,double>& _p1,
                   const point<3,double>& _p2,
                   double _size, double * _data,
                   bool _restricted=true)
      : p0_(_p0), p1_(_p1), p2_(_p2), size_(_size), data_(_data), restricted_(_restricted){
      sph=gluNewQuadric();
    }
    virtual void move(int x, int y) {
      point<2,double> p=_((double)x,(double)y);
      point<3,double> w=hc->plane_unproject(p,_(p0_,p1_,p2_));
      if(restricted_)
        w=barycentric_snap(w);
      data_[0]=w[0];
      data_[1]=w[1];
    }
    virtual void select() {
    }
    virtual void release() {
    }
    virtual void draw_proxy() {
      draw_sphere();
    }
    void draw_sphere() {
      glPushMatrix();
      point<3,double> p=_(data_[0],data_[1],1.0-data_[0]-data_[1]);
      p=barycentric_combination(p,_(p0_,p1_,p2_));
      glTranslate(p);
      gluSphere(sph,size_,8,8);
      glPopMatrix();
    }
    virtual void draw() {
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Green[0]);
      draw_sphere();
      glBegin(GL_LINE_LOOP);
      glVertex(p0_);
      glVertex(p1_);
      glVertex(p2_);
      glEnd();
    }
    virtual ~triangle_handle() {
      gluDeleteQuadric(sph);
    }

  };

  class circle_handle : public triangle_handle {
    point<3,double> center_;
    vec<3,double> dir_;
    double radius_;
    vec<3,double> a,b;
  public:
    circle_handle(const point<3,double>& _center,
                   const vec<3,double>& _dir,
                   double _radius,
                   double _size, double * _data,
                   bool _restricted=true)
      : triangle_handle(), center_(_center), dir_(_dir), radius_(_radius)

    {
      size_=_size;
      data_=_data;
      restricted_=_restricted;
      normalize(dir_);
      vec<3,double> v=_(0.23423,0.56456,-0.7123);
      a=cross(dir_,v);
      normalize(a);
      b=cross(dir_,a);
      p0_=center_+radius_*a;
      p1_=center_+radius_*b;
      p2_=center_+(-radius_)*(a+b);
    }
    virtual void move(int x, int y) {
      point<2,double> p=_((double)x,(double)y);
      point<3,double> w=hc->plane_unproject(p,_(p0_,p1_,p2_));
      point<3,double> q=barycentric_combination(w,_(p0_,p1_,p2_));
      vec<3,double> v=q-center_;
      double l=sqrt(dot(v,v));
      normalize(v);
      if(restricted_)
        l=radius_;
      else
        l=(l<=radius_) ? l : radius_;
      q=center_+l*v;
      data_[0]=q[0];
      data_[1]=q[1];
      data_[2]=q[2];
    }
    virtual void select() {
    }
    virtual void release() {
    }
    virtual void draw_proxy() {
      draw_sphere();
    }
    void draw_sphere() {
      glPushMatrix();
      point<3,double> p=_(data_[0],data_[1],data_[2]);
      glTranslate(p);
      gluSphere(sph,size_,8,8);
      glPopMatrix();
    }
    virtual void draw() {
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Green[0]);
      draw_sphere();
      glBegin(GL_LINE_LOOP);
      for(double th=0;th<(2.0*M_PI);th+=M_PI/90)
        glVertex(center_+radius_*cos(th)*a+radius_*sin(th)*b);
      glEnd();
    }
    virtual ~circle_handle() {
    }

  };
}
#endif // VGTL_CG_HANDLES_HPP
