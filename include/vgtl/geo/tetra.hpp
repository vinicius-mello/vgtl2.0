#ifndef VGTL_GEO_TETRA_HPP
#define VGTL_GEO_TETRA_HPP

/*! \file
 * \brief Geometric tetrahedron functions
 */


#include <vgtl/alg/plucker.hpp>

namespace vgtl {

  /*! \addtogroup geo
   * @{
   */

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Source code for the paper
  //
  //  "Fast Ray--Tetrahedron Intersection"
  //
  //  (c) Nikos Platis, Theoharis Theoharis 2003
  //
  //  Department of Informatics and Telecommunications,
  //  University of Athens, Greece
  //  {nplatis|theotheo}@di.uoa.gr
  //
  ////////////////////////////////////////////////////////////////////////////////


  // Computes the parametric distance tEnter and tLeave
  // of enterPoint and leavePoint from orig
  // (To enhance performance of the intersection algorithm, this code
  // should in practice be incorporated to the function below.)

  template <class Scalar>
  void compute_parametric_dist(
      const point<3,Scalar>& orig, const vec<3,Scalar>& dir,
      const point<3,Scalar>& enterPoint, const point<3,Scalar>& leavePoint,
      Scalar& tEnter, Scalar& tLeave)
  {
      Scalar one=1;
      if (dir[0])  {
          Scalar invDirx = one / dir[0];
          tEnter = (enterPoint[0] - orig[0]) * invDirx;
          tLeave = (leavePoint[0] - orig[0]) * invDirx;
      }
      else if (dir[1])  {
          Scalar invDiry = one / dir[1];
          tEnter = (enterPoint[1] - orig[1]) * invDiry;
          tLeave = (leavePoint[1] - orig[1]) * invDiry;
      }
      else  {
          double invDirz = one / dir[2];
          tEnter = (enterPoint[2] - orig[2]) * invDirz;
          tLeave = (leavePoint[2] - orig[2]) * invDirz;
      }
  }


  template <class Scalar>
  int Sign(Scalar x, Scalar Zero=0.000001) {
      return (x>Zero ? 1 : (x<-Zero ? -1 : 0));
  }

  // Ray--tetrahedron intersection algorithm using Pluecker coordinates

  template <class Scalar>
  bool ray_tetra_plucker(
      const point<3,Scalar>& orig, const vec<3,Scalar>& dir,
      const array<point<3,Scalar>,4>& vert,
      int& enterFace, int& leaveFace,
      point<3,Scalar>& enterPoint, point<3,Scalar>& leavePoint,
      Scalar& uEnter1, Scalar& uEnter2, Scalar& uLeave1, Scalar& uLeave2,
      Scalar& tEnter, Scalar& tLeave)
  {
      enterFace = -1;
      leaveFace = -1;

      Scalar uAB = 0, uAC = 0, uDB = 0, uDC = 0, uBC = 0, uAD = 0;
                      // Keep the compiler happy about uninitialized variables
      int signAB = -2, signAC = -2, signDB = -2,
          signDC = -2, signBC = -2, signAD = -2;

      // In the following: A,B,C,D=vert[i], i=0,1,2,3.
      plucker<Scalar> plRay(orig, dir);

      int nextSign = 0;

      Scalar one=1;

      // Examine face ABC
      uAB = plRay * plucker<Scalar>(vert[0], vert[1]);
      signAB = Sign(uAB);

      uAC = plRay * plucker<Scalar>(vert[0], vert[2]);
      signAC = Sign(uAC);

      if ((signAC == -signAB)  ||  (signAC == 0)  ||  (signAB == 0))  {
          // Face ABC may intersect with the ray
          uBC = plRay * plucker<Scalar>(vert[1], vert[2]);
          signBC = Sign(uBC);

          int signABC = signAB;
          if (signABC == 0)  {
              signABC = -signAC;
              if (signABC == 0)  {
                  signABC = signBC;
              }
          }

          if ((signABC != 0)  &&
              ((signBC == signABC)  ||  (signBC == 0)))  {
              // Face ABC intersects with the ray
              Scalar invVolABC = one / (uAB + uBC - uAC);
              if (signABC == 1)  {
                  enterFace = 3;
                  uEnter1 = -uAC * invVolABC;
                  uEnter2 =  uAB * invVolABC;
                  enterPoint =
                      barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                              _(vert[0],vert[1],vert[2]));

                  nextSign = -1;
              }
              else  {
                  leaveFace = 3;
                  uLeave1 = -uAC * invVolABC;
                  uLeave2 =  uAB * invVolABC;
                  leavePoint =
                      barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                              _(vert[0],vert[1],vert[2]));

                  nextSign = 1;
              }

              // Determine the other intersecting face between BAD, CDA, DCB
              // Examine face BAD
              uAD = plRay * plucker<Scalar>(vert[0], vert[3]);
              signAD = Sign(uAD);

              if ((signAD == nextSign)  ||  (signAD == 0))  {
                  // Face BAD may intersect with the ray
                  uDB = plRay * plucker<Scalar>(vert[3], vert[1]);
                  signDB = Sign(uDB);

                  if ((signDB == nextSign)  ||
                      ((signDB == 0)  &&
                       ((signAD != 0)  ||  (signAB != 0))))  {
                      // Face BAD intersects with the ray
                      Scalar invVolBAD = one / (uAD + uDB - uAB);
                      if (nextSign == 1)  {
                          enterFace = 2;
                          uEnter1 =  uDB * invVolBAD;
                          uEnter2 = -uAB * invVolBAD;
                          enterPoint =
                            barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                                    _(vert[1],vert[0],vert[3]));

                          compute_parametric_dist(orig, dir,
                                                enterPoint, leavePoint,
                                                tEnter, tLeave);
                          return true;
                      }
                      else  {
                          leaveFace = 2;
                          uLeave1 =  uDB * invVolBAD;
                          uLeave2 = -uAB * invVolBAD;
                          leavePoint =
                              barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                                      _(vert[1],vert[0],vert[3]));

                          compute_parametric_dist(orig, dir,
                                                enterPoint, leavePoint,
                                                tEnter, tLeave);
                          return true;
                      }
                  }
              }

              // Face BAD does not intersect with the ray.
              // Determine the other intersecting face between CDA, DCB
              uDC = plRay * plucker<Scalar>(vert[3], vert[2]);
              signDC = Sign(uDC);

              if ((signDC == -nextSign)  ||
                  ((signDC == 0)  &&  ((signAD != 0)  ||  (signAC != 0))))  {
                  // Face CDA intersects with the ray
                  Scalar invVolCDA = one / (uAC - uDC - uAD);
                  if (nextSign == 1)  {
                      enterFace = 1;
                      uEnter1 =  uAC * invVolCDA;
                      uEnter2 = -uDC * invVolCDA;

                      enterPoint =
                          barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                                  _(vert[2],vert[3],vert[0]));

                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
                  else  {
                      leaveFace = 1;
                      uLeave1 =  uAC * invVolCDA;
                      uLeave2 = -uDC * invVolCDA;
                      leavePoint =
                        barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                              _(vert[2],vert[3],vert[0]));

                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
              }
              else  {
                  // Face DCB intersects with the ray
                  if (signDB == -2)  {
                      uDB = plRay * plucker<Scalar>(vert[3], vert[1]);
                  }

                  Scalar invVolDCB = one / (uDC - uBC - uDB);
                  if (nextSign == 1)  {
                      enterFace = 0;
                      uEnter1 = -uDB * invVolDCB;
                      uEnter2 =  uDC * invVolDCB;
                      enterPoint =
                        barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                              _(vert[3],vert[2],vert[1]));

                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
                  else  {
                      leaveFace = 0;
                      uLeave1 = -uDB * invVolDCB;
                      uLeave2 =  uDC * invVolDCB;
                      leavePoint =
                        barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                            _(vert[3],vert[2],vert[1]));
                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
              }
          }
      }

      // Examine face BAD
      uAD = plRay * plucker<Scalar>(vert[0], vert[3]);
      signAD = Sign(uAD);

      if ((signAD == -signAB)  ||  (signAB == 0)  ||  (signAD == 0))  {
          // Face BAD may intersect with the ray
          uDB = plRay * plucker<Scalar>(vert[3], vert[1]);
          signDB = Sign(uDB);

          int signBAD = -signAB;
          if (signBAD == 0)  {
              signBAD = signAD;
              if (signBAD == 0)  {
                  signBAD = signDB;
              }
          }

          if ((signBAD != 0)  &&
              ((signDB == signBAD)  ||  (signDB == 0)))  {
              // Face BAD intersects with the ray
              Scalar invVolBAD = one / (uAD + uDB - uAB);
              if (signBAD == 1)  {
                  enterFace = 2;
                  uEnter1 =  uDB * invVolBAD;
                  uEnter2 = -uAB * invVolBAD;
                  enterPoint =
                    barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                        _(vert[1],vert[0],vert[3]));
                  nextSign = -1;
              }
              else  {
                  leaveFace = 2;
                  uLeave1 =  uDB * invVolBAD;
                  uLeave2 = -uAB * invVolBAD;
                  leavePoint =
                    barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                      _(vert[1],vert[0],vert[3]));

                  nextSign = 1;
              }

              // Determine the other intersecting face between CDA, DCB
              uDC = plRay * plucker<Scalar>(vert[3], vert[2]);
              signDC = Sign(uDC);

              if ((signDC == -nextSign)  ||
                  ((signDC == 0)  &&  ((signAD != 0)  ||  (signAC != 0))))  {
                  // Face CDA intersects with the ray
                  Scalar invVolCDA = one / (uAC - uDC - uAD);
                  if (nextSign == 1)  {
                      enterFace = 1;
                      uEnter1 =  uAC * invVolCDA;
                      uEnter2 = -uDC * invVolCDA;
                      enterPoint =
                        barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                          _(vert[2],vert[3],vert[0]));

                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
                  else  {
                      leaveFace = 1;
                      uLeave1 =  uAC * invVolCDA;
                      uLeave2 = -uDC * invVolCDA;
                      leavePoint =
                        barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                        _(vert[2],vert[3],vert[0]));
                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
              }
              else  {
                  // Face DCB intersects with the ray
                  if (signBC == -2)  {
                      uBC = plRay * plucker<Scalar>(vert[1], vert[2]);
                  }

                  Scalar invVolDCB = one / (uDC - uBC - uDB);
                  if (nextSign == 1)  {
                      enterFace = 0;
                      uEnter1 = -uDB * invVolDCB;
                      uEnter2 =  uDC * invVolDCB;
                      enterPoint =
                        barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                        _(vert[3],vert[2],vert[1]));
                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
                  else  {
                      leaveFace = 0;
                      uLeave1 = -uDB * invVolDCB;
                      uLeave2 =  uDC * invVolDCB;
                      leavePoint =
                        barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                      _(vert[3],vert[2],vert[1]));
                      compute_parametric_dist(orig, dir,
                                            enterPoint, leavePoint,
                                            tEnter, tLeave);
                      return true;
                  }
              }
          }
      }

      // Examine face CDA
      if ((-signAD == signAC)  ||  (signAC == 0)  ||  (signAD == 0))  {
          // Face CDA may intersect with the ray
          uDC = plRay * plucker<Scalar>(vert[3], vert[2]);
          signDC = Sign(uDC);

          int signCDA = signAC;
          if (signCDA == 0)  {
              signCDA = -signAD;
              if (signCDA == 0)  {
                  signCDA = -signDC;
              }
          }

          if ((signCDA != 0)  &&
              ((-signDC == signCDA)  ||  (signDC == 0)))  {
              // Face CDA intersects with the ray
              // Face DCB also intersects with the ray
              Scalar invVolCDA = one / (uAC - uDC - uAD);

              if (signBC == -2)  {
                  uBC = plRay * plucker<Scalar>(vert[1], vert[2]);
              }
              if (signDB == -2)  {
                  uDB = plRay * plucker<Scalar>(vert[3], vert[1]);
              }
              Scalar invVolDCB = one / (uDC - uBC - uDB);

              if (signCDA == 1)  {
                  enterFace = 1;
                  uEnter1 =  uAC * invVolCDA;
                  uEnter2 = -uDC * invVolCDA;
                  enterPoint =
                    barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                  _(vert[2],vert[3],vert[0]));

                  leaveFace = 0;
                  uLeave1 = -uDB * invVolDCB;
                  uLeave2 =  uDC * invVolDCB;
                  leavePoint =
                    barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                                _(vert[3],vert[2],vert[1]));
                  compute_parametric_dist(orig, dir,
                                        enterPoint, leavePoint,
                                        tEnter, tLeave);
                  return true;
              }
              else  {
                  leaveFace = 1;
                  uLeave1 =  uAC * invVolCDA;
                  uLeave2 = -uDC * invVolCDA;
                  leavePoint =
                    barycentric_combination(_((one-uLeave1-uLeave2),uLeave1,uLeave2),
                              _(vert[2],vert[3],vert[0]));

                  enterFace = 0;
                  uEnter1 = -uDB * invVolDCB;
                  uEnter2 =  uDC * invVolDCB;
                  enterPoint =
                    barycentric_combination(_((one-uEnter1-uEnter2),uEnter1,uEnter2),
                                _(vert[3],vert[2],vert[1]));
                  compute_parametric_dist(orig, dir,
                                        enterPoint, leavePoint,
                                        tEnter, tLeave);
                  return true;
              }
          }
      }

      // Three faces do not intersect with the ray, the fourth will not.
      return false;
  }



}

#endif // VGTL_GEO_TETRA_HPP
