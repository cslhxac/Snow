#ifndef __SNOW_SCENE__
#define __SNOW_SCENE__
#include <PhysBAM_Tools/Grids_Uniform_Array/ARRAYS_ND.h>
#include <PhysBAM_Tools/PhysBAM_Tools/Grids_Uniform/UNIFORM_NODE_ITERATOR.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <vector>
namespace PhysBAM{
  template<class TV> class Snow_Scene<TV>{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV:m};

  public:
    STREAM_TYPE stream_type;
    GRID<TV> uniform_grid;
    ARRAY<T,TV_INT> node_velocity;
    ARRAY<T,TV_INT> mass;
    std::vector<GRID<TV>::INDEX> adj_nodes;//the adjacent nodes' index for each particle;

    Snow_Scene();
    ~Snow_Scene();
    void Write_Output_Files(const int frame){
          std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocities",node_velocities);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",grid);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mass",mass);
    }
    virtual void Read_Output_Files(const int frame);
  };
}
