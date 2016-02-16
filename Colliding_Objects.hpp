
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include "CG_VECTOR.h"
#include "CG_SYSTEM.h"
#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#ifndef __COLLIDING_OBJECTS__
#define __COLLIDING_OBJECTS__
namespace PhysBAM{
  template<class TV> class Colliding_Objects{
      
      typedef typename TV::SCALAR T;
      typedef typename TV::template REBIND<int>::TYPE TV_INT;
      enum workaround1{d=TV::m};
      //I am using box for now
      //TODO change to Collision Geometry Collection
      ARRAY<BOX<TV>,int> colliding_Boxes;
      static const T mu = 0.1;
      
  public:
      Colliding_Objects(){};
      ~Colliding_Objects(){};
      void Collide(const TV x_in,const TV v_in,TV& x_out,TV& v_out)const{
          T min_distance = 0.0;
          x_out = x_in;
          v_out = v_in;
          for (int i = 1; i <= colliding_Boxes.Size(); ++i) {
              T tmp_distance = colliding_Boxes(i).Signed_Distance(x_in);
              if (tmp_distance < min_distance) {
                  
                  min_distance = tmp_distance;
                  
                  TV v_co = TV();
                  
                  TV v_rel = v_in - v_co;
                  
                  TV n = colliding_Boxes(i).Normal(x_in);
                  //std::cout << n << std::endl;
                  T v_n = TV::Dot_Product(n,v_rel);
                  if(v_n < 0){
                      TV v_t = v_rel - v_n * n;
                      if(v_t.Magnitude() <= - mu * v_n){
                          v_rel = TV();
                      }else{
                          v_rel = v_t + v_t * mu * v_n / v_t.Magnitude();
                      }
                      //std::cout << "v_rel: " << v_rel << std::endl;
                      x_out = x_in - min_distance * n;
                      v_out = v_rel + v_co;
                  }
              }
          }
      }
      void AddBox(const TV& minimum_corner,const TV& maximum_corner){
          colliding_Boxes.Resize(colliding_Boxes.Size() + 1);
          colliding_Boxes(colliding_Boxes.Size()) = BOX<TV>(minimum_corner,maximum_corner);
      }
  };
}
#endif
