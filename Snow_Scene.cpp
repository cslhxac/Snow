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
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "Snow_Scene.hpp"
using namespace PhysBAM;
template<class TV>
Snow_Scene<TV>::Snow_Scene(STREAM_TYPE st_in)
:stream_type(st_in),output_directory("out"),collection(snow_particles),maximum_dt(1e-3),number_of_frames(60),frame_time(1e-2)
{
    lambda0 = youngs_modulus * poissons_ratio / ((1 + poissons_ratio) * (1 - 2 * poissons_ratio));
    mu0 = youngs_modulus / (2 * (1 + poissons_ratio));
    
    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
    
    snow_particles_view = FREE_PARTICLES<TV>::Create(snow_particles);
    
    snow_particles.Store_Velocity();
    
    int scale = 100;
    RANGE<TV> range(TV(),TV::All_Ones_Vector());TV_INT counts = TV_INT::All_Ones_Vector() * scale + 1;;
    std::cout << "init s" << std::endl;
    
    MATRIX<T,3> test;
    test(1,2) = 1;
    std::cout << test << std::endl;
    //std::cout << MATRIX<T,3>() << std::endl;
    //std::cout << TV() << std::endl;
    
    Initialize_Grid(counts,range);
    //std::cout << uniform_grid.Minimum_Edge_Length() << std::endl;
    //std::cout << uniform_grid.Maximum_Edge_Length() << std::endl;
    Initialize_Fields();
    snow_particles.Store_Velocity();
    /*snow_particles.array_collection -> Add_Elements(10);
     for(int i =0;i < 10;++i){
     snow_particles.X(i) = TV(0,T(1-i)/T(1-10),0.5);
     snow_particles_view -> nodes.Append(i);
     }*/
    
    collection.Add_Structure(snow_particles_view);
    std::cout << "init e" << std::endl;
    
    /*MATRIX<T,3> tmpF;
    tmpF (1,1) = 1;tmpF (1,2) = 2;tmpF (1,3) = 3;
    tmpF (2,1) = 4;tmpF (2,2) = 5;tmpF (2,3) = 6;
    tmpF (3,1) = 7;tmpF (3,2) = 8;tmpF (3,3) = 9;
    
    MATRIX<T,3> tmpdF;
    tmpdF (1,1) = 1;tmpdF (1,2) = 2;tmpdF (1,3) = 3;
    tmpdF (2,1) = 3;tmpdF (2,2) = 2;tmpdF (2,3) = 1;
    tmpdF (3,1) = 1;tmpdF (3,2) = 2;tmpdF (3,3) = 3;
    
    dR(tmpF,tmpdF);
    
    int k;
    std::cin >> k ;*/
    /*MATRIX<T,3> tmp;
    tmp (2,1) = 10;
    TV tmp_tv(1,2,3);
    
    std::cout << tmp * tmp_tv << std::endl;*/
};

template<class TV>
void Snow_Scene<TV>::Create_Ball(TV position,T radius,TV v){
    T h = uniform_grid.Minimum_Edge_Length();
    //This function should iterate through the center of the grids that intercects the ball and add a particle to it.
    for(typename GRID<TV>::CELL_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        RANGE<TV> bound = iterator.Bounding_Box();
        T weight = 0; //counting how many coners in the ball radius.
        //the following part works for only 3d
        STATIC_ASSERT(d==3);
        ARRAY<TV,VECTOR<int,3> > corners;
        bound.Corners(corners);
        for(int i=0;i<=1;i++){
            for(int j=0;j<=1;j++){
                for(int k=0;k<=1;k++){
                    if((corners(i,j,k) - position).Magnitude() <= radius){
                        ++weight;
                    }
                }
            }
        }
        if(weight > 0){
            snow_particles.array_collection -> Add_Elements(1);
            //std::cout << snow_particles.X.Size() << ":" << snow_particles.V.Size() << std::endl;
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities );
            snow_particles.X(snow_particles.X.Size()) = bound.Center();
            
            snow_particles.V(snow_particles.X.Size()) = v;
            //if(std::isnan(snow_particles.V(snow_particles.X.Size())(1))) std::cout << "here!!!" << std::endl;
            snow_particles_view -> nodes.Append(snow_particles.X.Size());
            /*
            snow_particles.array_collection -> Add_Elements(6);
            //std::cout << snow_particles.X.Size() << ":" << snow_particles.V.Size() << std::endl;
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            particle_mass.Append(T(weight) / T(8.0) * h * h * h * densities / 6.0);
            
            snow_particles.X(snow_particles.X.Size()) = bound.Center() - TV(-h,0,0);
            snow_particles.X(snow_particles.X.Size() - 1) = bound.Center() - TV(h,0,0);
            snow_particles.X(snow_particles.X.Size() - 2) = bound.Center() - TV(0,h,0);
            snow_particles.X(snow_particles.X.Size() - 3) = bound.Center() - TV(0,-h,0);
            snow_particles.X(snow_particles.X.Size() - 4) = bound.Center() - TV(0,0,h);
            snow_particles.X(snow_particles.X.Size() - 5) = bound.Center() - TV(0,0,-h);
            
            snow_particles.V(snow_particles.X.Size()) = v;
            snow_particles.V(snow_particles.X.Size() - 1) = v;
            snow_particles.V(snow_particles.X.Size() - 2) = v;
            snow_particles.V(snow_particles.X.Size() - 3) = v;
            snow_particles.V(snow_particles.X.Size() - 4) = v;
            snow_particles.V(snow_particles.X.Size() - 5) = v;
            //if(std::isnan(snow_particles.V(snow_particles.X.Size())(1))) std::cout << "here!!!" << std::endl;
            snow_particles_view -> nodes.Append(snow_particles.X.Size() - 5);
            snow_particles_view -> nodes.Append(snow_particles.X.Size() - 4);
            snow_particles_view -> nodes.Append(snow_particles.X.Size() - 3);
            snow_particles_view -> nodes.Append(snow_particles.X.Size() - 2);
            snow_particles_view -> nodes.Append(snow_particles.X.Size() - 1);
            snow_particles_view -> nodes.Append(snow_particles.X.Size());*/
        }
    }
}
template<class TV>
void Snow_Scene<TV>::Rasterize_To_Grid(){
    //std::cout << snow_particles.X.Size() << std::endl;
    //This function Rasterize Velocity and Mass onto the grid from the particles.
    //const T h = uniform_grid.Minimum_Edge_Length();
    //when rasterize to grid, the grid position need to be reset
    //reset grid position
    
    //Here we can regenerate grid based on the boundary of the particles!
    Reset_Fields();
    //end reset
    weight_gradient_table.Resize(snow_particles.X.Size());
    weight_table.Resize(snow_particles.X.Size());
    node_index_table.Resize(snow_particles.X.Size());
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        TV p_position = snow_particles.X(i);
        
        weight_table(i).Resize(1,5,1,5,1,5,true,false,0);
        weight_gradient_table(i).Resize(1,5,1,5,1,5,true,false,TV());
        node_index_table(i).Resize(1,5,1,5,1,5,true,false,TV_INT());
        TV_INT closest_node = uniform_grid.Closest_Node(p_position);
        
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    node_index_table(i)(TV_INT(x + 3,y + 3,z + 3)) = closest_node + TV_INT(x,y,z);
                    if(uniform_grid.Node_Indices().Inside(closest_node + TV_INT(x,y,z),0)){
                        
                        TV grid_position = uniform_grid.Node(closest_node + TV_INT(x,y,z));
                        //std::cout << uniform_grid.Node(closest_node + TV_INT(x,y,z)) << ":" << uniform_grid.X(closest_node + TV_INT(x,y,z)) << std::endl;
                        (weight_table(i))(TV_INT(x + 3,y + 3,z + 3)) = Weight_function(p_position - grid_position);
                        (weight_gradient_table(i))(TV_INT(x + 3,y + 3,z + 3)) = grad_Weight_function(p_position - grid_position);
                        
                        mass(closest_node + TV_INT(x,y,z)) += (weight_table(i))(TV_INT(x + 3,y + 3,z + 3)) * particle_mass(i);
                        
                        //std::cout << mass(closest_node + TV_INT(x,y,z)) << ":" << particle_mass(i) << ":" << p_position - grid_position << ":" << Weight_function(p_position - grid_position) << std:: endl;
                        node_velocities(closest_node + TV_INT(x,y,z)) += (weight_table(i))(TV_INT(x + 3,y + 3,z + 3)) * particle_mass(i) * snow_particles.V(i);
                        //if(std::isnan(snow_particles.V(i)(1))) std::cout << "here!!!" << std::endl;
                        //if(std::isnan(node_velocities(closest_node + TV_INT(x,y,z))(1))) std::cout << "here!!!" << std::endl;
                    }
                }
            }
        }
    }
    //now normalize the velocity
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        if(mass(iterator.Node_Index()) >= massThreshold){
            node_velocities(iterator.Node_Index()) /= mass(iterator.Node_Index());
            node_velocities_old(iterator.Node_Index()) = node_velocities(iterator.Node_Index());
        }
    }
}
template<class TV>
void Snow_Scene<TV>::Compute_Volumes(){
    //Compute the Volumes and densities for each particle
    particle_volume.Resize(snow_particles.X.Size());
    particle_density.Resize(snow_particles.X.Size());
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        particle_volume(i) = 0;
        particle_density(i) = 0;
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        T h = uniform_grid.Minimum_Edge_Length();
                        particle_density(i) += mass(node_index) * weight_table(i)(TV_INT(x + 3,y + 3,z + 3))/h/h/h;
                    }
                }
            }
        }
        particle_volume(i) = particle_mass(i) / particle_density(i);
    }
}

template<class TV>
void Snow_Scene<TV>::Initialize_Grid(TV_INT counts,RANGE<TV> domain){
    uniform_grid.Initialize(counts,domain,false);
    node_velocities.Resize(uniform_grid.Node_Indices());
    node_velocities_old.Resize(uniform_grid.Node_Indices());
    mass.Resize(uniform_grid.Node_Indices());
    node_force.Resize(uniform_grid.Node_Indices());
    node_position.Resize(uniform_grid.Node_Indices());
    node_index_table_3d_to_1d.Resize(uniform_grid.Node_Indices());
    //std::cout << "node count: " << uniform_grid.Node_Indices().Size() << std::endl;
    node_index_table_1d_to_3d.Resize((uniform_grid.Node_Indices().Edge_Lengths() + TV_INT(1,1,1)).Product());
};
template<class TV>
void Snow_Scene<TV>::Initialize_Fields(){
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        mass(iterator.Node_Index()) = 0;
        node_velocities(iterator.Node_Index()) = TV();
        node_force(iterator.Node_Index()) = TV();
        node_position(iterator.Node_Index()) = uniform_grid.Node(iterator.Node_Index());
    }
};
template<class TV>
void Snow_Scene<TV>::Reset_Fields(){
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        mass(iterator.Node_Index()) = 0;
        node_velocities(iterator.Node_Index()) = TV();
        node_velocities_old(iterator.Node_Index()) = TV();
        node_force(iterator.Node_Index()) = TV();
        
        node_position(iterator.Node_Index()) = uniform_grid.Node(iterator.Node_Index());
    }
};
template<class TV>
void Snow_Scene<TV>::Clear_Node_Force(){
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        node_force(iterator.Node_Index()) = TV();
    }
};
template<class TV>
void Snow_Scene<TV>::Initilize_Deformation_Gradients(){
    //make sure all the paticle mass are generated;
    Fe.Remove_All();
    Fe_hat.Remove_All();
    Fp.Remove_All();
    
    Fe.Resize(particle_mass.Size(),true,false,MATRIX<T,3>::Identity_Matrix());
    Fe_hat.Resize(particle_mass.Size(),true,false,MATRIX<T,3>::Identity_Matrix());
    Fp.Resize(particle_mass.Size(),true,false,MATRIX<T,3>::Identity_Matrix());
};
template<class TV>
void Snow_Scene<TV>::Compute_Fe_Hat(const ARRAY<TV,TV_INT>& x_hat){
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        Fe_hat(i) = MATRIX<T,3>();
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        Fe_hat(i) += Card_Product(x_hat(node_index) - /*uniform_grid.Node*/node_position(node_index),(weight_gradient_table(i))(TV_INT(x + 3,y + 3,z + 3)));
                    }
                }
            }
        }
        Fe_hat(i) += MATRIX<T,3>::Identity_Matrix();
        Fe_hat(i) *= Fe(i);
    }
}
template<class TV>
void Snow_Scene<TV>::Add_Node_Force(){
    Compute_Fe_Hat(node_position);
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        node_force(node_index) -= particle_volume(i) / Fp(i).Determinant() / Fe(i).Determinant() * P(Fe_hat(i),Fp(i)) * (Fe(i).Transposed()) * (weight_gradient_table(i))(TV_INT(x + 3,y + 3,z + 3));
                    }
                }
            }
        }
    }
    TV tmp;
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        tmp += node_force(iterator.Node_Index());
    }
    std::cout << "NetForce " << tmp << std::endl;
}
template<class TV>
void Snow_Scene<TV>::Add_Node_Body_Force(){
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        node_force(iterator.Node_Index()) -= TV::Axis_Vector(2) * mass(iterator.Node_Index()) * gravity;
    }
}
template<class TV>
void Snow_Scene<TV>::Explicit_Time_Integral_For_Grid(const T dt){
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        //We have NAN here!! causing particles to disapear!
        //std::cout << node_velocities(iterator.Node_Index()) << std::endl;
        if(mass(iterator.Node_Index()) >= massThreshold){
            node_position(iterator.Node_Index()) += T(0.5) * dt * node_velocities(iterator.Node_Index());
            node_velocities(iterator.Node_Index()) += dt * node_force(iterator.Node_Index()) / mass(iterator.Node_Index());
            node_position(iterator.Node_Index()) += T(0.5) * dt * node_velocities(iterator.Node_Index());
        }
    }
};
template<class TV>
void Snow_Scene<TV>::Update_Deformation_Gradient(const T dt){
    particle_velocity_gradient.Resize(snow_particles.X.Size());
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        particle_velocity_gradient(i) = MATRIX<T,3>();
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        particle_velocity_gradient(i) += Card_Product(node_velocities(node_index), (weight_gradient_table(i))(TV_INT(x + 3,y + 3,z + 3)));
                    }
                }
            }
        }
        MATRIX<T,3> Fe_tmp = (MATRIX<T,3>::Identity_Matrix() + dt * particle_velocity_gradient(i)) * Fe(i);
        MATRIX<T,3> F_tmp = Fe_tmp * Fp(i);
        MATRIX<T,3> U,V;
        DIAGONAL_MATRIX<T,3> E;
        Fe_tmp.Fast_Singular_Value_Decomposition(U,E,V);
        DIAGONAL_MATRIX<T,3> E_clamped;
        E_clamped(1,1) = clamp(E(1,1),1 - theta_c,1 + theta_s);
        E_clamped(2,2) = clamp(E(2,2),1 - theta_c,1 + theta_s);
        E_clamped(3,3) = clamp(E(3,3),1 - theta_c,1 + theta_s);
        Fe(i) = U * E_clamped * V.Transposed();
        Fp(i) = V * E_clamped.Inverse() * U.Transposed() * F_tmp;
    }
};
template<class TV>
void Snow_Scene<TV>::Update_Particle_Velocity(const T dt){
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        TV V_flip = snow_particles.V(i);
        snow_particles.V(i) = TV();
        TV V_pic = TV();
        for (int x = - 2; x <= 2; ++x) {
            for (int y = - 2; y <= 2; ++y) {
                for (int z = - 2; z <= 2; ++z) {
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        V_pic += node_velocities(node_index) * weight_table(i)(TV_INT(x + 3,y + 3,z + 3));
                        V_flip += (node_velocities(node_index) - node_velocities_old(node_index)) * weight_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    }
                }
            }
        }
        /*std::cout << (1 - alpha) * V_pic + alpha * V_flip << std::endl;
         std::cout << V_pic << std::endl;
         std::cout << V_flip << std::endl;
         std::cout << dt * gravity << std::endl;
         std::cout << std::endl;*/
        snow_particles.V(i) = (1 - alpha) * V_pic + alpha * V_flip;
    }
};
template<class TV>
void Snow_Scene<TV>::Collide_Particles(){
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        TV x_tmp = snow_particles.X(i);
        TV v_tmp = snow_particles.V(i);
         colliding_Objects.Collide(x_tmp,v_tmp,snow_particles.X(i),snow_particles.V(i));
    }
}
template<class TV>
void Snow_Scene<TV>::Update_Particle_Position(const T dt){
    for (int i = 1; i <= snow_particles.X.Size(); ++i) {
        snow_particles.X(i) += snow_particles.V(i) * dt;
    }
};
template<class TV>
void Snow_Scene<TV>::Write_Output_Files(const int frame){
    std::string f = STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocities",node_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",uniform_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",uniform_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",mass);
    
    collection.Write(stream_type,output_directory,frame,0,true);
};
template<class TV>
void Snow_Scene<TV>::Read_Output_Files(const int frame){};

template<class TV>
void Snow_Scene<TV>::Semi_Implicit_Time_Integral_For_Grid(const T time, const T dt){
    int nNodes = (uniform_grid.Node_Indices().Edge_Lengths() + TV_INT(1,1,1)).Product();
    //std::cout << "Checking nNodes: "<< nNodes << std::endl;
    ARRAY<TV> rhs(nNodes);
    ARRAY<TV> V_next(nNodes);
    int i = 1;
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next(),++i){
        if(mass(iterator.Node_Index()) >= massThreshold){
            TV rhs_tmp = node_velocities(iterator.Node_Index()) + dt * node_force(iterator.Node_Index()) / mass(iterator.Node_Index());
            TV x_tmp;
            colliding_Objects.Collide(node_position(iterator.Node_Index()),rhs_tmp,x_tmp,rhs(i));
            rhs(i) *= mass(iterator.Node_Index());
            //std::cout << "Node Force: "<< node_force(iterator.Node_Index()) << std::endl;
        }
        V_next(i) = node_velocities(iterator.Node_Index());
        node_index_table_3d_to_1d(iterator.Node_Index()) = i;
        node_index_table_1d_to_3d(i) = iterator.Node_Index();
    }
    //std::cout << "Checking nNodes: "<< i << std::endl;
    //Temporary vectors, required by Conjugate Gradients
    ARRAY<TV> temp_q(nNodes),temp_s(nNodes),temp_r(nNodes),
    temp_k(nNodes),temp_z(nNodes);
    
    // Encapsulate all vectors in CG-mandated format
    CG_VECTOR<T> cg_x(V_next),cg_b(rhs),
    cg_q(temp_q),cg_s(temp_s),cg_r(temp_r),cg_k(temp_k),cg_z(temp_z);
    
    // Generate CG-formatted system object
    CG_SYSTEM<T> cg_system(*this,time,dt);
    
    // Generate Conjugate Gradients solver object
    CONJUGATE_GRADIENT<T> cg;
    cg.print_residuals = false;
    cg.print_diagnostics = true;
    cg.restart_iterations = 20;
    
    // Solve linear system using CG
    cg.Solve(cg_system,
             cg_x,cg_b,cg_q,cg_s,cg_r,cg_k,cg_z,
             1e-7,0,50);
    
    for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
        if(mass(iterator.Node_Index()) >= massThreshold){
            //Hack
            /*if (V_next(node_index_table_3d_to_1d(iterator.Node_Index())).Magnitude() > 1e2) {
                V_next(node_index_table_3d_to_1d(iterator.Node_Index())).Normalize();
                V_next(node_index_table_3d_to_1d(iterator.Node_Index())) *= 1e2;
            }*/
            //end Hack
            node_position(iterator.Node_Index()) += T(0.5) * (node_velocities(iterator.Node_Index()) + V_next(node_index_table_3d_to_1d(iterator.Node_Index())));
            node_velocities(iterator.Node_Index()) = V_next(node_index_table_3d_to_1d(iterator.Node_Index()));
        }
    }
};

template<class TV>
void Snow_Scene<TV>::Ap(const ARRAY<TV>& v,ARRAY<TV>& result){
    dF.Remove_All();
    dF.Resize(snow_particles.X.Size(),true,false,MATRIX<T,3>());
    MATRIX<T,3> Ap_tmp;
    
    result.Fill(TV());
    
    for (int i = 1; i <= snow_particles.X.Size(); ++i){
        dF(i) = MATRIX<T,3>();
        for (int x = - 2; x <= 2; ++x){
            for (int y = - 2; y <= 2; ++y){
                for (int z = - 2; z <= 2; ++z){
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        dF(i) += Card_Product(v(node_index_table_3d_to_1d(node_index)), weight_gradient_table(i)(TV_INT(x + 3,y + 3,z + 3))) * Fe(i);
                    }
                }
            }
        }
    }
    for (int i = 1; i <= snow_particles.X.Size(); ++i){
        MATRIX<T,3> Ap_tmp = MATRIX<T,3>();
        T mu = mu0 * std::exp(sigma * (1 - Fp(i).Determinant()));
        T lambda = lambda0 * std::exp(sigma * (1 - Fp(i).Determinant()));
        T J = Fe(i).Determinant();
        T dJ = J * MATRIX<T,3>::Inner_Product(Fe(i).Inverse_Transposed(),dF(i));
        MATRIX<T,3> dF_Inverse_Transposed = -Fe(i).Inverse_Transposed() * dF(i).Transposed() * Fe(i).Inverse_Transposed();
        MATRIX<T,3> dJF_Inverse_Transposed = dJ * Fe(i).Inverse_Transposed() + J * dF_Inverse_Transposed;
        Ap_tmp = 2 * mu * (dF(i) - dR(Fe(i),dF(i))) + lambda * J * Fe(i).Inverse_Transposed() * dJ + lambda * (J - 1) * dJF_Inverse_Transposed;
    
        for (int x = - 2; x <= 2; ++x){
            for (int y = - 2; y <= 2; ++y){
                for (int z = - 2; z <= 2; ++z){
                    TV_INT node_index = node_index_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    if(uniform_grid.Node_Indices().Inside(node_index,0)){
                        result(node_index_table_3d_to_1d(node_index)) += particle_volume(i) * Ap_tmp * Fe(i).Transposed() * weight_gradient_table(i)(TV_INT(x + 3,y + 3,z + 3));
                    }
                }
            }
        }
    }
};
template class Snow_Scene<VECTOR<double,3> >;

