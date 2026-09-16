#include "widget.hpp"
#include <QMouseEvent>
#include <QMenu>
#include <QtConcurrentRun>
#include <QFutureWatcher>
#include <vgtl/cg/GL.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/top/subdivide_until.hpp>
#include <vgtl/top/iso/tubular.hpp>
#include <vgtl/geo/iso.hpp>
#include <vgtl/top/complex_buffer.hpp>
#include <vgtl/comb/triangular.hpp>
#include <vgtl/geo/diffpars.hpp>
#include <vgtl/alg/linalg.hpp>
#include <stack>
#include <iostream>
#include <sstream>

void smooth(vgtl::complex_buffer<T>& cb, int it);


template <class T>
inline std::string to_string (const T& t)
{
   std::stringstream ss;
   ss << t;
   return ss.str();
}

using namespace std;
using namespace vgtl;

void smooth(complex_buffer<T>& cb, int it);

double clipxu=1.0, clipxb=0.0;
double clipyu=1.0, clipyb=0.0;
double clipzu=1.0, clipzb=0.0;

void GLWidget::setControls(Cell(T) s) {
  vgtl::array<Vertex(T),4> vs;
  vertices(*app::t,s,vs);
  int l=level(*app::t,s)/6;
  for(int i=0;i<4;++i) {
    string face("ctrl.face");
    FaceHandle * fh=new FaceHandle(*app::t,face_op(*app::t,s,i),0.01/pow(1.25,l));
    hc.add_handle(face+to_string<int>(i),fh);
    connect(fh,SIGNAL(valueChanged(fct)),this,SLOT(faceChanged(fct)));
  }
  vgtl::array<vgtl::array<dim_t,2>,6> ae;
  ae[0]=_<dim_t>(0,1);
  ae[1]=_<dim_t>(0,2);
  ae[2]=_<dim_t>(0,3);
  ae[3]=_<dim_t>(1,2);
  ae[4]=_<dim_t>(1,3);
  ae[5]=_<dim_t>(2,3);
  for(int i=0;i<6;++i) {
    string edge("ctrl.edge");
    EdgeHandle * eh=new EdgeHandle(*app::t,face_ind(*app::t,s,ae[i]),0.01/pow(1.25,l));
    hc.add_handle(edge+to_string<int>(i),eh);
    connect(eh,SIGNAL(valueChanged(edg)),this,SLOT(edgeChanged(edg)));
  }
  for(int i=0;i<4;++i) {
    string vert("ctrl.vert");
    SignalHandle * sh=new SignalHandle(*app::t,vs[i],0.02/pow(1.25,l));
    hc.add_handle(vert+to_string<int>(i),sh);
    connect(sh,SIGNAL(valueChanged(vtx)),this,SLOT(signalChanged(vtx)));
  }
}

void GLWidget::removeControls() {
  empty_set(sel_cell);
  hc.delete_range("ctrl.edge0","ctrl.edge5");
  hc.delete_range("ctrl.vert0","ctrl.vert3");
  hc.delete_range("ctrl.face0","ctrl.face3");
}


void init() {
  GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };

  glLightfv (GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv (GL_LIGHT0, GL_SPECULAR, light_specular);

  glShadeModel(GL_SMOOTH);
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE,&vgtl::White[0]);
  glEnable (GL_CULL_FACE);
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_NORMALIZE);
}

void draw_face(Facet(T) f, double o) {
  vec<3,double> n=o*normal(*app::t,f);
  vgtl::array<Vertex(T),3> vs;
  vertices(*app::t,f,vs);
  vgtl::array<point<3,double>,3> ps;
  double s[3];
  for(int i=0;i<3;++i) {
    s[i]=signal(*app::t,vs[i]);
    ps[i]=euclidean_point(*app::t,vs[i]);
  }
  if(o<0) {
    swap(ps[1],ps[2]);
    swap(s[1],s[2]);
  }
  bool crs=true;
  if(s[0]*s[1]<0) {
    if(s[1]*s[2]<0) {
      swap(ps[0],ps[1]);
      swap(s[0],s[1]);
      swap(ps[1],ps[2]);
      swap(s[1],s[2]);
    }
  } else  if(s[1]*s[2]<0) {
    swap(ps[0],ps[2]);
    swap(s[0],s[2]);
    swap(ps[1],ps[2]);
    swap(s[1],s[2]);
  } else {
    crs=false;
  }
  glPolygonMode(GL_FRONT,GL_LINE);
  glDisable(GL_LIGHTING);
  glColor(vgtl::White);
  glBegin(GL_TRIANGLES);
  glVertex(ps[0]);
  glVertex(ps[1]);
  glVertex(ps[2]);
  glEnd();
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0,1.0);
  glEnable(GL_LIGHTING);
  glPolygonMode(GL_FRONT,GL_FILL);
  vec<4,float> red=_(0.9f,0.2f,0.3f,1.0f);
  vec<4,float> blu=_(0.0f,0.4f,0.9f,1.0f);
  if(!crs) {
    if(s[0]<0)
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&red[0]);
    else
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&blu[0]);
    glBegin(GL_TRIANGLES);
    glNormal(n);
    glVertex(ps[0]);
    glVertex(ps[1]);
    glVertex(ps[2]);
    glEnd();
  } else {
    if(s[0]<0)
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&red[0]);
    else
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&blu[0]);
    point<3,double> a,b;
    a=ps[0]+0.5*(ps[1]-ps[0]);
    b=ps[0]+0.5*(ps[2]-ps[0]);
    glBegin(GL_TRIANGLES);
    glNormal(n);
    glVertex(ps[0]);
    glVertex(a);
    glVertex(b);
    glEnd();
    if(s[0]<0)
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&blu[0]);
    else
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&red[0]);
    glBegin(GL_QUADS);
    glNormal(n);
    glVertex(a);
    glVertex(ps[1]);
    glVertex(ps[2]);
    glVertex(b);
    glEnd();
    glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Yellow[0]);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex(a);
    glVertex(b);
    glEnd();
    glLineWidth(1.0);
  }
  glDisable(GL_POLYGON_OFFSET_FILL);
}

void GLWidget::draw_cell(Cell(T) s) {
  double o=orientation(*app::t,s);
  if(viewMode==tetraMode||viewMode==crossMode) {
    for(int i=0;i<=3;++i,o*=-1) {
         draw_face(face_op(*app::t,s,i),o);
    }
  } else if(viewMode==curvedMode) {
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Yellow[0]);
    glMaterialfv(GL_BACK,GL_DIFFUSE,&vgtl::Red[0]);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3,GL_DOUBLE,0,&cell_vertices(*app::t,s)[0][0]);
    glNormalPointer(GL_DOUBLE,0,&cell_normals(*app::t,s)[0][0]);
    pair<int,int> di=iso::dim(*app::t,s);
    if(di.first==1) {
      glDrawElements(GL_QUADS,4*quad_idx.size(),GL_UNSIGNED_INT,&quad_idx[0][0]);
    } else if(di.first==0||di.first==2) {
      glDrawElements(GL_TRIANGLES,3*tri_idx.size(),GL_UNSIGNED_INT,&tri_idx[0][0]);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glEnable(GL_CULL_FACE);
  } else if(viewMode==isoMode) {
    list<vgtl::array<Edge(T),3> > ls;
    list<int> ori;
    iso::triangulate(*app::t,s,ls,ori);
    list<vgtl::array<Edge(T),3> >::iterator li;
    list<int>::iterator lo=ori.begin();
    for(li=ls.begin();li!=ls.end();++li,++lo) {
      vgtl::array<Edge(T),3> ae=*li;
      vgtl::array<point<3,double>, 3> ap;
      for(int i=0;i<3;++i) {
        vgtl::array<point<3,double>,2> vs;
        euclidean_points(*app::t,ae[i],vs);
        ap[i]=vs[0]+0.5*(vs[1]-vs[0]);
      }
      if(*lo<0) swap(ap[1],ap[2]);
      vec<3,double> n=cross(ap[1]-ap[0],ap[2]-ap[0]);
      normalize(n);
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Yellow[0]);
      glBegin(GL_TRIANGLES);
      glNormal(n);
      glVertex(ap[0]);
      glVertex(ap[1]);
      glVertex(ap[2]);
      glEnd();
      glMaterialfv(GL_FRONT,GL_DIFFUSE,&vgtl::Red[0]);
      glBegin(GL_TRIANGLES);
      glNormal(-n);
      glVertex(ap[0]);
      glVertex(ap[2]);
      glVertex(ap[1]);
      glEnd();
    }
  }
}

void GLWidget::setVertices(Cell(T) c) {
  point<3,double> * vtx_ptr=cell_vertices(*app::t,c);
  vec<3,double> * nml_ptr=cell_normals(*app::t,c);
  if(vtx_ptr!=0) {
    delete [] vtx_ptr;
    delete [] nml_ptr;
    cell_vertices_set(*app::t,c,0);
    cell_normals_set(*app::t,c,0);
  }
  vgtl::array<Vertex(T),4> vs;
  vertices(*app::t,c,vs);
  pair<int,int> di=iso::dim(*app::t,vs);
  vector<point<2,double> > * vtx_vec_ptr;
  if(di.first==1) {
    vtx_ptr=new point<3,double>[(nsub+1)*(nsub+1)];
    nml_ptr=new vec<3,double>[(nsub+1)*(nsub+1)];
    vtx_vec_ptr=&quad_bcoord;
  } else if((di.first==0)||(di.first==2)) {
    vtx_ptr=new point<3,double>[(nsub+2)*(nsub+1)/2];
    nml_ptr=new vec<3,double>[(nsub+2)*(nsub+1)/2];
    vtx_vec_ptr=&tri_bcoord;
  } else return;
  vector<point<2,double> >& vtx_vec=*vtx_vec_ptr;
  cell_vertices_set(*app::t,c,vtx_ptr);
  cell_normals_set(*app::t,c,nml_ptr);

  diffpars(*app::t,c,bb2,Edge(T)());
  diffpars(*app::t,c,bb3,Facet(T)());
  diffpars(*app::t,c,bb4,Cell(T)());
  vgtl::array<point<3,double>,4> ps;
  euclidean_points(*app::t,c,ps);
  vgtl::array<point<4,double>,3> bps;
  iso::basis(*app::t,c,vs,bps);

  for(index_t i=0;i<vtx_vec.size();++i) {
    vgtl::array<vec<4,double>,2> bvs;
    bvs[0]=bps[1]-bps[0];
    bvs[1]=bps[2]-bps[0];
    vgtl::array<vec<4,double>,4> J;
    point<2,double> be=vtx_vec[i];
    point<4,double> bp=bps[0]+be[0]*bvs[0]+be[1]*bvs[1];
    for(int k=0;k<app::kc3;++k) {
      bp=jacobian(bb4,bp,J);
      bvs[0]=apply(J,bvs[0]);
      bvs[1]=apply(J,bvs[1]);
    }
    for(int k=0;k<app::kc2;++k) {
      bp=jacobian(bb3,bp,J);
      bvs[0]=apply(J,bvs[0]);
      bvs[1]=apply(J,bvs[1]);
    }
    for(int k=0;k<app::kc1;++k) {
      bp=jacobian(bb2,bp,J);
      bvs[0]=apply(J,bvs[0]);
      bvs[1]=apply(J,bvs[1]);
    }
    point<3,double> bc;
    bc=barycentric_combination(bp,ps);
    vtx_ptr[i]=bc;
    vec<3,double> du;
    du=barycentric_combination(bvs[0],ps);
    vec<3,double> dv;
    dv=barycentric_combination(bvs[1],ps);
    vec<3,double> n=cross(du,dv);
    vgtl::normalize(n);
    nml_ptr[i]=n;
  }
}


void GLWidget::signalChanged(vtx v) {
  list<Cell(T)> st;
  star(*app::t,v,back_inserter(st));
  list<Cell(T)>::iterator li;
  for(li=st.begin();li!=st.end();++li) {
    Cell(T) c=*li;
    setVertices(c);
  }
}

void GLWidget::edgeChanged(edg e) {
  list<Cell(T)> st;
  star(*app::t,e,back_inserter(st));
  list<Cell(T)>::iterator li;
  for(li=st.begin();li!=st.end();++li) {
    Cell(T) c=*li;
    setVertices(c);
  }
}

void GLWidget::faceChanged(fct f) {
  list<Cell(T)> st;
  star(*app::t,f,back_inserter(st));
  list<Cell(T)>::iterator li;
  for(li=st.begin();li!=st.end();++li) {
    Cell(T) c=*li;
    setVertices(c);
  }
}

void GLWidget::splitCell() {
  if(!empty(sel_cell)) {
    maubach_subdivide(*app::t,sel_cell,app::ns);
    removeControls();
    repaint();
  }
}

void GLWidget::setIsoView() {
  viewMode=isoMode;
  repaint();
}

void GLWidget::setTetraView() {
  viewMode=tetraMode;
  repaint();
}

void GLWidget::setCrossView() {
  viewMode=crossMode;
  repaint();
}

void GLWidget::setCurvedView() {
  viewMode=curvedMode;
  repaint();
}

void GLWidget::optimizeFinished() {
  smoothing=false;
  repaint();
}

void GLWidget::optimizeThread() {
  future=QtConcurrent::run(this,&GLWidget::optimize);
  smoothing=true;
  removeControls();
  watcher.setFuture(future);
  repaint();
}

void GLWidget::optimize() {
  complex_buffer<T> cb(*app::t);
  iso::tubular(cb);
  //buffer_set<3>(cb).clear();
  smooth(cb,10);
  for(set<Cell(T)>::iterator ci=buffer_set<3,T>(cb).begin();
      ci!=buffer_set<3,T>(cb).end(); ++ci) {
    Cell(T) c=*ci;
    setVertices(c);
  }
  //repaint();
}

void GLWidget::contextMenu(const QPoint &p) {
  QMenu menu(this);
  if(!empty(sel_cell)) menu.addAction("Split", this, SLOT(splitCell()));
  menu.addSeparator();
  if(viewMode!=isoMode) menu.addAction("Iso View", this, SLOT(setIsoView()));
  if(viewMode!=tetraMode) menu.addAction("Tetra View", this, SLOT(setTetraView()));
  if(viewMode!=crossMode) menu.addAction("Cross View", this, SLOT(setCrossView()));
  if(viewMode!=curvedMode) menu.addAction("Curved View", this, SLOT(setCurvedView()));
  menu.addSeparator();
  if(!smoothing) menu.addAction("Smooth", this, SLOT(optimizeThread()));
  menu.exec(mapToGlobal(p));
}

void GLWidget::initializeGL() {
  init();
  hc.perspective_on();
  this->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(this, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(contextMenu(const QPoint &)));
  hc.add_handle("clip.xb",new segment_handle<double>(_(0.0,0.0,-0.1),_(1.0,0.0,-0.1),0.02,&clipxb));
  hc.add_handle("clip.xu",new segment_handle<double>(_(0.0,0.0,-0.1),_(1.0,0.0,-0.1),0.02,&clipxu));
  hc.add_handle("clip.yb",new segment_handle<double>(_(-0.1,0.0,-0.1),_(-0.1,1.0,-0.1),0.02,&clipyb));
  hc.add_handle("clip.yu",new segment_handle<double>(_(-0.1,0.0,-0.1),_(-0.1,1.0,-0.1),0.02,&clipyu));
  hc.add_handle("clip.zb",new segment_handle<double>(_(-0.1,0.0,0.0),_(-0.1,0.0,1.0),0.02,&clipzb));
  hc.add_handle("clip.zu",new segment_handle<double>(_(-0.1,0.0,0.0),_(-0.1,0.0,1.0),0.02,&clipzu));
  maubach_subdivide_to_level(*app::t,3,app::ns);
  viewMode=tetraMode;
  nsub=4;
  map<vgtl::array<index_t,3>,index_t> trv;
  vgtl::array<index_t,3> c=_<index_t>(nsub,0,0);
  index_t i=0;
  do {
    trv[c]=i;++i;
    point<2,double> bcoord=_<double>(((double)c[0])/nsub,((double)c[1])/nsub);
    tri_bcoord.push_back(bcoord);
  } while(next_triangular(nsub,&c[0],&c[3]));
  c=_<index_t>(nsub-1,0,0);
  i=0;
  do {
    vgtl::array<index_t,3> ti;
    for(dim_t i=0;i<3; ++i) {
      c[i]+=1;
      ti[i]=trv[c];
      c[i]-=1;
    }
    tri_idx.push_back(ti);
  } while(next_triangular(nsub-1,&c[0],&c[3]));
  c=_<index_t>(nsub-2,0,0);
  do {
    vgtl::array<index_t,3> ti;
    for(dim_t i=0;i<3; ++i) {
      c[0]++;c[1]++;c[2]++;
      c[i]-=1;
      ti[i]=trv[c];
      c[i]+=1;
      c[0]--;c[1]--;c[2]--;
    }
    tri_idx.push_back(ti);
  } while(next_triangular(nsub-2,&c[0],&c[3]));
  i=0;
  map<vgtl::array<index_t,2>,index_t> qdv;
  for(index_t k=0;k<=nsub;++k) {
    for(index_t j=0;j<=nsub;++j) {
      qdv[_<index_t>(j,k)]=i;++i;
      point<2,double> bcoord=_<double>(((double)j)/nsub,((double)k)/nsub);
      quad_bcoord.push_back(bcoord);
    }
  }
  for(index_t k=0;k<nsub;++k) {
    for(index_t j=0;j<nsub;++j) {
      vgtl::array<index_t,4> qi;
      qi[0]=qdv[_<index_t>(j,k)];
      qi[1]=qdv[_<index_t>(j+1,k)];
      qi[2]=qdv[_<index_t>(j+1,k+1)];
      qi[3]=qdv[_<index_t>(j,k+1)];
      quad_idx.push_back(qi);
    }
  }
  app::kc1=2; app::kc2=4; app::kc3=6;
  connect(&watcher, SIGNAL(finished()), this, SLOT(optimizeFinished()));
}

void draw_background() {
  glDisable(GL_LIGHTING);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glPolygonMode(GL_FRONT,GL_FILL);
  glBegin(GL_QUADS);
  glColor(0.05*vgtl::White+0.95*vgtl::Black);
  glVertex(_(-1.0,-1.0,1.0));
  glVertex(_(1.0,-1.0,1.0));
  glColor(0.2*vgtl::White+0.8*vgtl::Black);
  glVertex(_(1.0,1.0,1.0));
  glVertex(_(-1.0,1.0,1.0));
  glEnd();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_LIGHTING);
}

void draw_rays() {
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  for(double x=-1.0;x<=1.0;x+=0.5)
    for(double y=-1.0;y<=1.0;y+=0.5) {
      glVertex(_(x,y,-2.0));
      glVertex(_(x,y,2.0));
    }
  glEnd();
  glEnable(GL_LIGHTING);
}

void draw_ball() {
  glDisable(GL_LIGHTING);
  glBegin(GL_LINE_LOOP);
  for(double theta=0;theta<2*M_PI;theta+=M_PI/30) {
    glVertex(_(cos(theta),0.0,sin(theta)));
  }
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(double theta=0;theta<2*M_PI;theta+=M_PI/30) {
    glVertex(_(0.0,cos(theta),sin(theta)));
  }
  glEnd();
  glBegin(GL_LINE_LOOP);
  for(double theta=0;theta<2*M_PI;theta+=M_PI/30) {
    glVertex(_(cos(theta),sin(theta),0.0));
  }
  glEnd();
  glEnable(GL_LIGHTING);
}



Cell(T) GLWidget::pick_cell(int x, int y) {
  makeCurrent();
  Cell(T) selected;
  if(smoothing) return selected;
  vgtl::array<point<3,double>,4> vert;
  point<3,double> r;
  double t;
  double min_t=1.0;
  stack<Cell(T)> sc;
  Cell_it(T) ci,cend;
  for(simplices(*app::t,ci,cend);ci!=cend;++ci) {
    Cell(T) c=*ci;
    sc.push(c);
  }
  while(!sc.empty()) {
    Cell(T) c=sc.top();
    sc.pop();
    vgtl::array<point<3,double>,4> ps;
    euclidean_points(*app::t,c,ps);
    point<3,double> bc=barycenter(ps);
    if(((viewMode==crossMode||viewMode==isoMode||viewMode==curvedMode)&&(!iso::cross(*app::t,c)))||
       (!((bc[0]>=clipxb) && (bc[0]<=clipxu) && (bc[1]>=clipyb) && (bc[1]<=clipyu) && (bc[2]>=clipzb) && (bc[2]<=clipzu)))) {
      continue;
    }
    if(!is_current(*app::t,c)) {
      if(has_children(*app::t,c)) {
        sc.push(child(*app::t,c,0));
        sc.push(child(*app::t,c,1));
      }
    } else {
      euclidean_points(*app::t,c,vert);
      double o=orientation(*app::t,c);
      if(o<0) swap(vert[2],vert[3]);
      bool inter=hc.tetra_intersect(_((double)x,(double)y),vert,r,t);
      if(inter && (t<min_t)) {
        min_t=t;
        selected=c;
      }
    }
  }

  return selected;
}


void GLWidget::paintGL() {
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(_(0.0,0.0,3.0),_(0.0,0.0,0.0),_(0.0,1.0,0.0));
  glClearColor(vgtl::Black);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  draw_background();
  glPushMatrix();
  light_track.transform();
  GLfloat light_position[] = { 0.0, 0.0, 1.0, 0.0 };
  glLightfv (GL_LIGHT0, GL_POSITION, light_position);
  glPopMatrix();
  glPushMatrix();
  glRotate(model_track.rotation());
  if(dragging) glColor(vgtl::White);
  else glColor(0.6*vgtl::White);
  draw_ball();
  glPopMatrix();
  glPushMatrix();
  glRotate(light_track.rotation());
  if(light_dragging) {
    glColor(vgtl::Yellow);
    draw_rays();
  }
  glPopMatrix();

  model_track.transform();
  glTranslate(_(-0.5,-0.5,-0.5));
  stack<Cell(T)> sc;
  Cell_it(T) ci,cend;
  for(simplices(*app::t,ci,cend);ci!=cend;++ci) {
    Cell(T) c=*ci;
    sc.push(c);
  }
  while(!sc.empty()) {
    Cell(T) c=sc.top();
    sc.pop();
    vgtl::array<point<3,double>,4> ps;
    euclidean_points(*app::t,c,ps);
    point<3,double> bc=barycenter(ps);
    if(((viewMode==crossMode||viewMode==isoMode||viewMode==curvedMode)&&(!iso::cross(*app::t,c)))||
       (!((bc[0]>=clipxb) && (bc[0]<=clipxu) && (bc[1]>=clipyb) && (bc[1]<=clipyu) && (bc[2]>=clipzb) && (bc[2]<=clipzu)))) {
      if(sel_cell==c) {
        removeControls();
      }
      continue;
    }
    if(!is_current(*app::t,c)) {
      if(has_children(*app::t,c)) {
        sc.push(child(*app::t,c,0));
        sc.push(child(*app::t,c,1));
      }
    } else {
      draw_cell(c);
    }
  }
  hc.draw_handles();
}

void GLWidget::resizeGL(int w, int h) {
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glLoadIdentity();
  double ar=static_cast<double>(w)/static_cast<double>(h);
  if(hc.perspective()) {
    gluPerspective(60,ar,0.5,10);
  }
  else {
    if(w<h)
      glOrtho(-1.0,1.0,-1.0/ar,1.0/ar,0.5,10);
    else
      glOrtho(-ar,ar,-1.0,1.0,0.5,10);
  }
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(_(0.0,0.0,3.0),_(0.0,0.0,0.0),_(0.0,1.0,0.0));
  model_track.resize();
  light_track.resize();
}

void GLWidget::mousePressEvent(QMouseEvent * ev) {
  makeCurrent();
  selected=hc.pick_handle(ev->x(),this->height()-ev->y());
  if((ev->modifiers() & Qt::ShiftModifier) && (ev->modifiers() & Qt::ControlModifier))
    light_track.start_motion(ev->x(),ev->y());
  else
    model_track.start_motion(ev->x(),ev->y());
  dragging=false;
}

void GLWidget::mouseReleaseEvent(QMouseEvent * ev) {
  makeCurrent();
  if(!dragging) {
    if(selected) selected->select();
    else {
      Cell(T) sel=pick_cell(ev->x(),this->height()-ev->y());
      if(!empty(sel)) {
        removeControls();
        sel_cell=sel;
        setControls(sel_cell);
      }
    }
  } else {
    if(selected) selected->release();
  }
  dragging=false;
  light_dragging=false;
  repaint();
}

void GLWidget::mouseMoveEvent(QMouseEvent * ev) {
  makeCurrent();
  if(selected)
    selected->move(ev->x(),this->height()-ev->y());
  else if((ev->modifiers() & Qt::ShiftModifier) && (ev->modifiers() & Qt::ControlModifier)) {
    light_track.move_rotation(ev->x(),ev->y());
    light_dragging=true;
  }
  else if(ev->modifiers() & Qt::ShiftModifier)
      model_track.move_scaling(ev->x(),ev->y());
  else if(ev->modifiers() & Qt::ControlModifier)
      model_track.move_pan(ev->x(),ev->y());
  else if(ev->modifiers() & Qt::AltModifier)
      model_track.move_zoom(ev->x(),ev->y());
  else model_track.move_rotation(ev->x(),ev->y());
  dragging=true;
  repaint();
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent * ev) {
  makeCurrent();
}

void GLWidget::keyPressEvent(QKeyEvent * ev) {
  makeCurrent();
  if(ev->key()==Qt::Key_O) {
    hc.perspective_off();
    resizeGL(this->width(),this->height());
  }
  else if(ev->key()==Qt::Key_P) {
    hc.perspective_on();
    resizeGL(this->width(),this->height());
  }
  repaint();
}
