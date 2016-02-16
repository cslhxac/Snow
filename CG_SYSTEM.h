// Copyright (c) 2011-2014, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include "CG_VECTOR.h"
#include "Snow_Scene.hpp"
#include<iostream>
#include <stdlib.h>
#include <time.h>

namespace PhysBAM{
    
    template<class T>
    class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
    {
        typedef KRYLOV_SYSTEM_BASE<T> BASE;
        typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
        
        Snow_Scene<VECTOR<T,3> >& layout;
        const T time;
        const T dt;
        
    public:
        CG_SYSTEM(Snow_Scene<VECTOR<T,3> >& layout_input,const T time_input,const T dt_input)
        :BASE(false,false),layout(layout_input),time(time_input),dt(dt_input) {}
        
        void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
        {
            const ARRAY<VECTOR<T,3>,int>& v_array = CG_VECTOR<T>::Array(v);
            ARRAY<VECTOR<T,3>,int>& result_array = CG_VECTOR<T>::Array(result);
            T beta = 0.5;//use trapezoidal rule
            
            result_array.Fill(VECTOR<T,3>());
            layout.Ap(v_array,result_array);
            //std::cout << "v_array.Size()" << v_array.Size() << std::endl;
            for(int p = 1; p <= v_array.Size(); p++){
                //if(layout.get_Mass(p) >= Snow_Scene<VECTOR<T,3> >::massThreshold){
                result_array(p) = layout.get_Mass(p) * v_array(p) + beta * dt * dt * result_array(p);
                //}else{
                //result_array(p) = TV();
                //}
            }
            //testing Ap is linear.
            /*std::srand(666);
            ARRAY<VECTOR<T,3>,int> v_tmp1(v_array.Size()),v_tmp2(v_array.Size()),v_tmp3(v_array.Size());
            ARRAY<VECTOR<T,3>,int> r_tmp1(v_array.Size()),r_tmp2(v_array.Size()),r_tmp3(v_array.Size());
            r_tmp1.Fill(VECTOR<T,3>());
            r_tmp2.Fill(VECTOR<T,3>());
            r_tmp3.Fill(VECTOR<T,3>());
            for(int p = 1; p <= v_array.Size(); p++){
                v_tmp1(p) = VECTOR<T,3>(std::rand() % 100 / 100.0,std::rand() % 100 / 100.0,std::rand() % 100 / 100.0);
                v_tmp2(p) = VECTOR<T,3>(std::rand() % 100 / 100.0,std::rand() % 100 / 100.0,std::rand() % 100 / 100.0);
                v_tmp3(p) = v_tmp1(p) + v_tmp2(p);
            }
            layout.Ap(v_tmp1,r_tmp1);
            layout.Ap(v_tmp2,r_tmp2);
            layout.Ap(v_tmp3,r_tmp3);
            for(int p = 1; p <= v_array.Size(); p++){
                VECTOR<T,3> tmp_V = r_tmp1(p) + r_tmp2(p) - r_tmp3(p);
                std::cout << tmp_V << std::endl;
                for (int i = 1; i <= 3; ++i) {
                    if(tmp_V(i) > 1e-6){
                        std::cout << "error: Ap is not linear!" << std::endl;
                    }
                }
            }*/
        }
        
        double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const
        {
            const ARRAY<VECTOR<T,3>,int>& x_array=CG_VECTOR<T>::Array(x);
            const ARRAY<VECTOR<T,3>,int>& y_array=CG_VECTOR<T>::Array(y);
            
            double result=0.;
            for(int i=1;i<=x_array.m;i++)
                result+=VECTOR<T,3>::Dot_Product(x_array(i),y_array(i));
            return result;
        }
        
        T Convergence_Norm(const VECTOR_BASE& x) const
        {
            const ARRAY<VECTOR<T,3>,int>& x_array=CG_VECTOR<T>::Array(x);
            
            T result=0.;
            for(int i=1;i<=x_array.m;i++)
                result=std::max(result,x_array(i).Magnitude());
            return result;
        }
        
        void Project(VECTOR_BASE& x) const
        {
            //std::cout << "Project" << std::endl;
            ARRAY<VECTOR<T,3>,int>& x_array=CG_VECTOR<T>::Array(x);
            
            //layout.Clear_Values_Of_Kinematic_Particles(x_array);
        }
        
        void Set_Boundary_Conditions(VECTOR_BASE& v) const
        {
            ARRAY<VECTOR<T,3>,int>& v_array=CG_VECTOR<T>::Array(v);
            //std::cout << "Set_Boundary_Conditions" << std::endl;
            //layout.Set_Kinematic_Velocities(time+dt,v_array);
        }
        
        void Project_Nullspace(VECTOR_BASE& x) const {} // Just a stub (for solids)
        
        void Apply_Preconditioner(const VECTOR_BASE& r,VECTOR_BASE& z) const{};
    };
    
}

