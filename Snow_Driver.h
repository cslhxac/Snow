#ifndef __SNOW_DRIVER__
#define __SNOW_DRIVER__

#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include"Snow_Scene.hpp"
#include <PhysBAM_Tools/Log/LOG.h>
using namespace PhysBAM;
namespace PhysBAM{
    template<class TV> class Snow_Scene;
    
    template<class T>
    class Snow_Driver{
    public:
        typedef VECTOR<T,3> TV;
        Snow_Scene<TV> scene;
        T time;
        Snow_Driver(STREAM_TYPE st_in):scene(st_in){
            
        }
        ~Snow_Driver(){
            
        }
        void Run()
        {
            FILE_UTILITIES::Create_Directory(scene.output_directory);
            FILE_UTILITIES::Create_Directory(scene.output_directory + "/common");
            //std::cout << "1" << std::endl;
            scene.Create_Ball(TV(0.5,0.2,0.5),0.1,TV(0,-10,0));
            //scene.Create_Ball(TV(0.5,0.25,0.5),0.2,TV(0,0,0));
            //scene.Create_Ball(TV(0.65,0.8,0.2),0.1,TV(-10,0,0));
            scene.colliding_Objects.AddBox(TV(0,-0.1,0),TV(1,0.1,1));
            //std::cout << "2" << std::endl;
            scene.Reset_Fields();
            //std::cout << "3" << std::endl;
            //std::cout << "4" << std::endl;
            scene.Rasterize_To_Grid();
            Write_Frame(0);
            /*layout.Initialize();layout.Write_Output(0);
            */
            time = 0;
            for(int frame = 1; frame <= scene.number_of_frames; frame++){
                std::cout << "Frame: " << frame << std::endl;
                Simulate_Frame(frame);
                Write_Frame(frame);
            }
        }
        void Write_Frame(const int frame){
            FILE_UTILITIES::Create_Directory(scene.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
            scene.Write_Output_Files(frame);
            FILE_UTILITIES::Write_To_Text_File(scene.output_directory+"/common/last_frame",frame,"\n");
        }
        void Simulate_Frame(const int frame)
        {
            T frame_end_time = scene.frame_time * (T)frame, dt_max, dt;
            
            for(;time < frame_end_time; time += dt){
                dt_max = scene.maximum_dt;
                dt = std::min(dt_max,(T)1.001*(frame_end_time-time));
                Simulate_Time_Step(time,dt, frame == 1);
            }
        };
        void Simulate_Time_Step(const T time,const T dt, bool first_frame)
        {
            //std::cout << dt << std::endl;
            scene.Reset_Fields();
            scene.Rasterize_To_Grid();
            if(first_frame){
                scene.Compute_Volumes();
                scene.Initilize_Deformation_Gradients();
            }
            scene.Clear_Node_Force();
            scene.Add_Node_Force();
            scene.Add_Node_Body_Force();
            //scene.Explicit_Time_Integral_For_Grid(dt);
            scene.Semi_Implicit_Time_Integral_For_Grid(time,dt);
            scene.Update_Deformation_Gradient(dt);
            scene.Update_Particle_Velocity(dt);
            scene.Update_Particle_Position(dt);
            //T tmp;
            //std::cin >> tmp;
        }

    };
}
#endif

