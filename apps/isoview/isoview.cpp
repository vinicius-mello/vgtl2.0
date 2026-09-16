#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <list>
#include <vgtl/utl/pair_tie.hpp>
#include "data.hpp"
#include <vgtl/top/initial.hpp>
#include <vgtl/top/star.hpp>
#include <vgtl/top/maubach.hpp>
#include <vgtl/geo/mesh_io.hpp>
#include <vgtl/top/barycentric.hpp>
#include <vgtl/top/reset.hpp>
#include <vgtl/geo/cube.hpp>
#include <vgtl/cg/arcball.hpp>
#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QStatusBar>
#include <QtGui>
#include <QGLWidget>
#include "widget.hpp"
 
using namespace std;
using namespace vgtl;


namespace app {
				
	new_simplices ns;
  T * t;
  bool shrink;
  int kc3,kc2,kc1;


  //Cell(T) selected;
//  void select(Cell(T) s) {
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Red);
//    }
//    selected=s;
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Blue);
//      print_info(this->t,selected);
//    }
//    post_redisplay();
//  }
//  void refine() {
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Red);
//      maubach_subdivide(this->t,selected,app::ns);
//      Cell(T) c0,c1;
//      c0=child(this->t,selected,0);
//      c1=child(this->t,selected,1);
//      selected=c0;
//      color_set(this->t, selected, vgtl::Blue);
//      post_redisplay();
//    }
//  }
//  void weld() {
//    if(!empty(selected)) {
//      if(level(this->t,selected)==0) return;
//      color_set(this->t, selected, vgtl::Red);
//      do_nothing dn;
//      maubach_weld(this->t,selected,dn);
//      empty_set(selected);
//      post_redisplay();
//    }
//  }
//  void hide() {
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Black);
//      empty_set(selected);
//      post_redisplay();
//    }
//  }
//  void hide_neighboors() {
//    if(!empty(selected)) {
//      vgtl::array<Vertex(T),4> vs;
//      vertices(this->t,selected,vs);
//      for(int i=0; i<=3; ++i) {
//        list<Cell(T)> st;
//        star(this->t,vs[i],back_inserter(st));
//        for(list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j)
//          color_set(this->t,*j,vgtl::Black);
//      }
//      color_set(this->t, selected, vgtl::Blue);
//      post_redisplay();
//    }
//  }
//  void unhide_neighboors() {
//    if(!empty(selected)) {
//      vgtl::array<Vertex(T),4> vs;
//      vertices(this->t,selected,vs);
//      for(int i=0; i<=3; ++i) {
//        list<Cell(T)> st;
//        star(this->t,vs[i],back_inserter(st));
//        for(list<Cell(T)>::iterator j=st.begin(); j!=st.end(); ++j)
//          color_set(this->t,*j,vgtl::Red);
//      }
//      color_set(this->t, selected, vgtl::Blue);
//      post_redisplay();
//    }
//  }
//  void unhide_all() {
//    Cell_it(T) si, send;
//    for(simplices(this->t, si,send);si!=send;++si) {
//      color_set(this->t, *si, vgtl::Red);
//    }
//    empty_set(selected);
//    post_redisplay();
//  }
//  void initial() {
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Red);
//    }
//    vgtl::initial(this->t);
//    empty_set(selected);
//    post_redisplay();
//  }
//  void final() {
//    if(!empty(selected)) {
//      color_set(this->t, selected, vgtl::Red);
//    }
//    vgtl::final(this->t);
//    empty_set(selected);
//    post_redisplay();
//  }
//  void barycentric() {
//    //if(refined) return;
//    reset_subdivision(this->t);
//    barycentric_subdivision(this->t,app::ns);
//    empty_set(selected);
//    post_redisplay();
//  }
//  void tshrink() {
//    shrink=!shrink;
//    post_redisplay();
//  }

}
	


class MainWindow : public QMainWindow
{
public:
  MainWindow() {
    menuBar()->addMenu("&File");
    menuBar()->addMenu("&Edit");
    menuBar()->addMenu("&View");
    menuBar()->addMenu("&Help");

    // add status bar message
    statusBar()->showMessage("Welcome to meshview!");

    GLWidget * gl = new GLWidget(this);
    //gl->setAlignment( Qt::AlignLeft | Qt::AlignTop );
    //gl->setFrameStyle( 0 );
    setCentralWidget( gl );
  }

};

int main(int argc, char * argv[]) {
  QApplication a(argc,argv);
  app::t=new T();
  set_cube(*app::t,app::ns);
  initial_cells_set(*app::t);
  //MainWindow window;
  //window.show();
  GLWidget * window=new GLWidget(0);
  window->show();
  //ifstream in(argv[1],ios::binary|ios::in);
  //read_mesh(in,app::t,app::ns);

  return a.exec();
}




