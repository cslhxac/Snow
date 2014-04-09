#ifndef __SNOW_SCENE__
#define __SNOW_SCENE__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <iostream> 
#include <vector>
#include <string>
namespace PhysBAM{
  template<class TV> class Snow_Scene{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    enum workaround1{d=TV::m};

  public:
    STREAM_TYPE stream_type;
    GRID<TV> uniform_grid;
    ARRAY<TV,TV_INT> node_velocities;
    ARRAY<T,TV_INT> mass;
    int first_frame,last_frame;
    std::string output_directory;

    GEOMETRY_PARTICLES<TV> snow_particles;
    FREE_PARTICLES<TV>* snow_particles_view;

    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection;

    //std::vector<GRID<TV>::INDEX> adj_nodes;//the adjacent nodes' index for each particle;

    Snow_Scene(STREAM_TYPE st_in)
      :stream_type(st_in),output_directory("out"),collection(snow_particles)
    {
      Initialize_Geometry_Particle();
      Initialize_Read_Write_Structures();

      snow_particles_view = FREE_PARTICLES<TV>::Create(snow_particles);

      snow_particles.Store_Velocity();

      int scale = 100;
      RANGE<TV> range(TV(),TV::All_Ones_Vector()*0.5);range.max_corner(2)=1;TV_INT counts=TV_INT::All_Ones_Vector()*scale/2;counts(2)=scale;
      std::cout << "init s" << std::endl;
      Initialize_Grid(counts,range);
      node_velocities.Resize(uniform_grid.Node_Indices());
      mass.Resize(uniform_grid.Node_Indices());
      Initialize_Fields();

      snow_particles.array_collection -> Add_Elements(10);
      for(int i =0;i < 10;++i){
	snow_particles.X(i) = TV(0,T(1-i)/T(1-10),0.5);
	snow_particles_view -> nodes.Append(i);
      }

      collection.Add_Structure(snow_particles_view);
      std::cout << "init e" << std::endl;
    };
    ~Snow_Scene()
    {
    };
    void Write_Output_Files(const int frame){
          std::string f = STRING_UTILITIES::string_sprintf("%d",frame);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/velocities",node_velocities);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",uniform_grid);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",uniform_grid);
	  FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",mass);

	  collection.Write(stream_type,output_directory,frame,0,true);
    }

    void Initialize_Grid(TV_INT counts,RANGE<TV> domain){
      uniform_grid.Initialize(counts,domain,false);
    };
    void Initialize_Fields() 
    {
      for(typename GRID<TV>::NODE_ITERATOR iterator(uniform_grid);iterator.Valid();iterator.Next()){
	mass(iterator.Node_Index()) = 1;
	node_velocities(iterator.Node_Index()) = TV();
      }
    }

    void Read_Output_Files(const int frame){};
  };
}
#endif
