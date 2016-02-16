
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
#include "CG_VECTOR.h"
#include "CG_SYSTEM.h"
#include "Colliding_Objects.hpp"
#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#ifndef __SNOW_SCENE__
#define __SNOW_SCENE__
namespace PhysBAM{
  template<class TV> class Snow_Scene{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV::m};

  public:
      //the following are simulation parameters:
      static const T massThreshold = 1e-9;
      static const T youngs_modulus = 1.4e5;
      static const T poissons_ratio = 0.2;
      static const T densities = 4e2;
      static const T theta_c = 2.5e-2;
      static const T theta_s = 7.5e-3;
      static const T sigma = 10;//hardening coeffcient
      static const T alpha = 0.95;
      static const T gravity = 9.81;
      const T maximum_dt;
      const int number_of_frames;  // Total number of frames
      const T frame_time;          // Frame (snapshot) interval
      T mu0,lambda0;
      
      //end parameters
      STREAM_TYPE stream_type;
      GRID<TV> uniform_grid;
      ARRAY<TV,TV_INT> node_velocities;
      ARRAY<TV,TV_INT> node_velocities_old;
      ARRAY<TV,TV_INT> node_force;
      ARRAY<T,TV_INT> mass;
      ARRAY<TV,TV_INT> node_position;
      ARRAY<int,TV_INT> node_index_table_3d_to_1d;
      ARRAY<TV_INT,int> node_index_table_1d_to_3d;
      ARRAY<T,int> particle_mass;
      ARRAY<T,int> particle_volume;
      ARRAY<T,int> particle_density;

      ARRAY<MATRIX<T,3>,int> particle_velocity_gradient;
      int first_frame,last_frame;
      std::string output_directory;

      GEOMETRY_PARTICLES<TV> snow_particles;
      FREE_PARTICLES<TV>* snow_particles_view;

      ARRAY<ARRAY<T,TV_INT>,int> weight_table;
      ARRAY<ARRAY<TV,TV_INT>,int> weight_gradient_table;
      ARRAY<ARRAY<TV_INT,TV_INT>,int> node_index_table;//stores the indices of the nodes that are contributes to the particle
      
      ARRAY<MATRIX<T,3>,int> Fe;
      ARRAY<MATRIX<T,3>,int> Fe_hat;
      ARRAY<MATRIX<T,3>,int> Fp;
      ARRAY<MATRIX<T,3>,int> dF;
      //ARRAY<MATRIX<T,3>,int> P;//the stress tensor
      
      DEFORMABLE_GEOMETRY_COLLECTION<TV> collection;

      Colliding_Objects<TV> colliding_Objects;
      //std::vector<GRID<TV>::INDEX> adj_nodes;//the adjacent nodes' index for each particle;

      Snow_Scene(STREAM_TYPE st_in);
      void Create_Ball(TV position,T radius,TV v = TV());
      void Rasterize_To_Grid();
      void Compute_Volumes();

      void Initialize_Grid(TV_INT counts,RANGE<TV> domain);
      void Initialize_Fields();
      void Reset_Fields();
      void Clear_Node_Force();
      void Initilize_Deformation_Gradients();
      void Compute_Fe_Hat(const ARRAY<TV,TV_INT>& x_hat);
      void Add_Node_Force();
      void Add_Node_Body_Force();
      void Explicit_Time_Integral_For_Grid(const T dt);
      void Semi_Implicit_Time_Integral_For_Grid(const T time, const T dt);
      void Update_Deformation_Gradient(const T dt);
      void Update_Particle_Position(const T dt);
      void Update_Particle_Velocity(const T dt);
      void Collide_Particles();
      void Write_Output_Files(const int frame);
      void Read_Output_Files(const int frame);
      void Ap(const ARRAY<TV>& v,ARRAY<TV>& result);
      T get_Mass(const int index)const{
          return mass(node_index_table_1d_to_3d(index));
      };
      static MATRIX<T,3> dR(const MATRIX<T,3>& F,const MATRIX<T,3>& dF){
          MATRIX<T,3> U,V;
          DIAGONAL_MATRIX<T,3> E;
          F.Fast_Singular_Value_Decomposition(U,E,V);
          MATRIX<T,3> R,S,W;
          R = U * V.Transposed();
          S = V * E * V.Transposed();
          W = R.Transposed() * dF;
          
          //std::cout << F - R * S << std::endl;
          
          MATRIX<T,3> A;
          TV b;
          A(1,1) = S(1,3); A(1,2) = S(2,3); A(1,3) = -(S(1,1) + S(2,2)); b(1) = (W(1,2) - W(2,1));
          A(2,1) = S(1,2); A(2,2) = -(S(1,1) + S(3,3)); A(2,3) = S(2,3); b(2) = (W(3,1) - W(1,3));
          A(3,1) = -(S(2,2) + S(3,3)); A(3,2) = S(1,2); A(3,3) = S(1,3); b(3) = (W(2,3) - W(3,2));
          
          TV r = A.Inverse() * b;
          MATRIX<T,3> rx;
          rx(1,1) = 0; rx(1,2) = - r(3); rx(1,3) = r(2);
          rx(2,1) = r(3); rx(2,2) = 0; rx(2,3) = -r(1);
          rx(3,1) = -r(2); rx(3,2) = r(1); rx(3,3) = 0;
          
          MATRIX<T,3> dR = R * rx;
          
          //this following is just sanity check
          MATRIX<T,3> dS = W - rx * S;
          MATRIX<T,3> tmp = dF - (dR * S + R * dS);
          
          //std::cout << tmp << std::endl;
          bool error = false;
          for (int i = 0; i < 9; ++i) {
              //std::cerr << "dR conputation is wrong, check the algorithm please!: " << tmp.x[i] << std::endl;
              if(tmp.x[i] > 1e-4){
                  error = true;
              }
          }
          if(error){
              std::cerr << "dR conputation is wrong, check the algorithm please!: " << tmp << std::endl;
              std::cerr << "F: " << F << std::endl;
              std::cerr << "dF: " << dF << std::endl;
          }
          return dR;
          
      };
      static T N_function(T x){
          if (x >= -1 && x <= 1) {
              return T(0.5) * sign(x) * x * x * x - x * x + T(2.0)/T(3.0);
          }else if((x >= 1 && x <= 2) || (x <= -1 && x >= -2)){
              return T(-1)/T(6.0) * sign(x) * x * x * x + x * x - 2 * x * sign(x) + T(4)/T(3);
          }else{
              return 0;
          }
      }
      static T dN_function(T x){
          if (x >= -1 && x <= 1) {
              return T(3) / T(2) * sign(x) * x * x - 2 * x;
          }else if((x >= 1 && x <= 2) || (x <= -1 && x >= -2)){
              return T(-1)/T(2.0) * sign(x) * x * x + 2 * x - 2 * sign(x);
          }else{
              return 0;
          }
      }
      T Weight_function(TV p){
          STATIC_ASSERT(d == 3);
          T h = uniform_grid.Minimum_Edge_Length();
          return N_function(T(1) / h * p(1)) * N_function(T(1) / h * p(2)) * N_function(T(1) / h * p(3));
      }
      TV grad_Weight_function(TV p){
          STATIC_ASSERT(d == 3);
          T h = uniform_grid.Minimum_Edge_Length();
          return TV(dN_function(T(1) / h * p(1)) / h * N_function(T(1) / h * p(2)) * N_function(T(1) / h * p(3)),
                    dN_function(T(1) / h * p(2)) / h * N_function(T(1) / h * p(1)) * N_function(T(1) / h * p(3)),
                    dN_function(T(1) / h * p(3)) / h * N_function(T(1) / h * p(1)) * N_function(T(1) / h * p(2)));
      }
      static MATRIX<T,3> Card_Product(const TV& x1,const TV& x2){
          //this function performs a card product from a 3 x 1 vector to a 1 x 3 vector, generate a 3 x 3 matrix;
          MATRIX<T,3> result;
          for (int i = 1; i <= 3; ++i) {
              for (int j = 1;j <= 3;++j){
                  result(i,j) = x1(i) * x2(j);
              }
          }
          return result;
      };
      MATRIX<T,3> P(MATRIX<T,3> Fe,MATRIX<T,3> Fp){
          MATRIX<T,3> Ue,Ve;
          DIAGONAL_MATRIX<T,3> Ee;
          Fe.Fast_Singular_Value_Decomposition(Ue,Ee,Ve);
          MATRIX<T,3> Re = Ue * Ve.Transposed();
          MATRIX<T,3> Se = Ve * Ee * Ve.Transposed();
          T mu = mu0 * std::exp(sigma * (1 - Fp.Determinant()));
          T lambda = lambda0 * std::exp(sigma * (1 - Fp.Determinant()));
          return Re * 2 * mu * (Se - MATRIX<T,3>::Identity_Matrix()) + lambda * (Fe.Determinant() - 1) * Fe.Determinant() * Fe.Inverse_Transposed();
      };
      ~Snow_Scene()
      {
      };
  };
}
#endif
