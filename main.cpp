#include "Snow_Scene.hpp"

using namepace PhysBAM;

int main(int argc,char *argv[]){
  typedef float T;
  typedef float RW;
  STREAM_TYPE stream_type((RW()));
  Snow_Scene<VECTOR<T,3>> scene();
    
  FILE_UTILITIES::Create_Directory(example.output_directory);
  FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",0));
  FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
  FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),"");
  FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",0,"\n");
  FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",0,"\n");

  scene.Write_Output_Files(0);
}

