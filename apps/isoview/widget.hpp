#ifndef GLWIDGET_H
#define GLWIDGET_H
#include <QGLWidget>
#include <QFutureWatcher>
#include "app.hpp"
#include <vgtl/cg/trackball.hpp>
#include <vgtl/cg/handles.hpp>
#include <vgtl/alg/bb.hpp>

typedef Vertex(T) vtx;
typedef Edge(T) edg;
typedef Facet(T) fct;

class GLWidget : public QGLWidget {
  Q_OBJECT
public:
  GLWidget(QWidget* wid) :
      QGLWidget(QGLFormat(QGL::DepthBuffer|QGL::Rgba|QGL::DoubleBuffer),wid,0), dragging(false),
      light_dragging(false), smoothing(false)
  {}
protected:
  void paintGL();
  void initializeGL();
  void resizeGL(int w, int h);
  void mousePressEvent(QMouseEvent * ev);
  void mouseReleaseEvent(QMouseEvent * ev);
  void mouseMoveEvent(QMouseEvent * ev);
  void mouseDoubleClickEvent(QMouseEvent *ev);
  void keyPressEvent(QKeyEvent * ev);
  Cell(T) pick_cell(int x, int y);
  void draw_cell(Cell(T) c);
  void setControls(Cell(T) s);
  void removeControls();
  void setVertices(Cell(T) c);
public slots:
  void signalChanged(vtx v);
  void edgeChanged(edg e);
  void faceChanged(fct f);
  void contextMenu(const QPoint &);
  void splitCell();
  void setIsoView();
  void setCrossView();
  void setTetraView();
  void setCurvedView();
  void optimize();
  void optimizeThread();
  void optimizeFinished();
private:
  vgtl::trackball model_track;
  vgtl::trackball light_track;
  bool dragging;
  bool light_dragging;
  bool perspective;
  vgtl::handle_controller hc;
  vgtl::handle * selected;
  Cell(T) sel_cell;
  enum { tetraMode, crossMode, isoMode, curvedMode } viewMode;
  std::vector<vgtl::array<vgtl::index_t,3> > tri_idx;
  std::vector<vgtl::point<2,double> > tri_bcoord;
  std::vector<vgtl::array<vgtl::index_t,4> > quad_idx;
  std::vector<vgtl::point<2,double> > quad_bcoord;
  unsigned int nsub;
  vgtl::bbform<4,double> bb2;
  vgtl::bbform<4,double> bb3;
  vgtl::bbform<4,double> bb4;
  QFutureWatcher<void> watcher;
  QFuture<void> future;
  bool smoothing;
  //int kc3,kc2,kc1;
};


class SignalHandle : public QObject, public vgtl::signal_handle<double>  {
  Q_OBJECT
protected:
  Vertex(T) v_;
public:
  SignalHandle(T& t, Vertex(T) v, double _size) {
    v_=v;
    size_=_size;
    p_=vgtl::euclidean_point(t,v);
    data_=vgtl::signal_ptr(t,v);
  }
  virtual void move(int x, int y) {
    vgtl::signal_handle<double>::move(x,y);
  };
  virtual void select() {
    vgtl::signal_handle<double>::select();
    emit valueChanged(v_);
  }
  virtual void release() {
  }
  virtual void draw_proxy() {
    vgtl::signal_handle<double>::draw_proxy();
  }

  void draw_sphere() {
    vgtl::signal_handle<double>::draw_sphere();
  }
  virtual void draw() {
    vgtl::signal_handle<double>::draw();
  }
  virtual ~SignalHandle() {
  }
signals:
  void valueChanged(vtx v);
};

class FaceHandle : public QObject, public vgtl::triangle_handle  {
  Q_OBJECT
protected:
  Facet(T) f_;
public:
  FaceHandle(T& t, Facet(T) f, double _size) {
    f_=f;
    size_=_size;
    vgtl::array<vgtl::point<3,double>,3> vs;
    vgtl::euclidean_points(t,f,vs);
    p0_=vs[0];
    p1_=vs[1];
    p2_=vs[2];
    data_=vgtl::diffpars_ptr(t,f);
  }
  virtual void move(int x, int y) {
    vgtl::triangle_handle::move(x,y);
  };
  virtual void select() {
    vgtl::triangle_handle::select();
  }
  virtual void release() {
    emit valueChanged(f_);
  }
  virtual void draw_proxy() {
    vgtl::triangle_handle::draw_proxy();
  }
  virtual void draw() {
    vgtl::triangle_handle::draw();
  }
  virtual ~FaceHandle() {
  }
signals:
  void valueChanged(fct v);
};

class EdgeHandle : public QObject, public vgtl::segment_handle<double>  {
  Q_OBJECT
protected:
  Edge(T) e_;
public:
  EdgeHandle(T& t, Edge(T) e, double _size) {
    e_=e;
    size_=_size;
    vgtl::array<vgtl::point<3,double>,2> vs;
    vgtl::euclidean_points(t,e,vs);
    p0_=vs[1];
    p1_=vs[0];
    data_=vgtl::diffpars_ptr(t,e);
  }
  virtual void move(int x, int y) {
    vgtl::segment_handle<double>::move(x,y);
  };
  virtual void select() {
    vgtl::segment_handle<double>::select();
  }
  virtual void release() {
    emit valueChanged(e_);
  }
  virtual void draw_proxy() {
    vgtl::segment_handle<double>::draw_proxy();
  }
  virtual void draw() {
    vgtl::segment_handle<double>::draw();
  }
  virtual ~EdgeHandle() {
  }
signals:
  void valueChanged(edg e);
};

#endif // WIDGET_H

