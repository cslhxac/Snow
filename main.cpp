#include "Snow_Scene.hpp"
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>

using namespace PhysBAM;

int main(int argc,char *argv[]){
  typedef float T;
  typedef float RW;
  STREAM_TYPE stream_type((RW()));
  
  Snow_Scene<VECTOR<T,3> >* scene = new Snow_Scene<VECTOR<T,3> >(stream_type);
  
  FILE_UTILITIES::Create_Directory(scene -> output_directory);
  FILE_UTILITIES::Create_Directory(scene -> output_directory+STRING_UTILITIES::string_sprintf("/%d",0));
  FILE_UTILITIES::Create_Directory(scene -> output_directory+STRING_UTILITIES::string_sprintf("/%d",1));
  FILE_UTILITIES::Create_Directory(scene -> output_directory+"/common");
  FILE_UTILITIES::Write_To_Text_File(scene -> output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",0),"");
  FILE_UTILITIES::Write_To_Text_File(scene -> output_directory+"/common/first_frame",0,"\n");
  FILE_UTILITIES::Write_To_Text_File(scene -> output_directory+"/common/last_frame",1,"\n");

  scene -> Write_Output_Files(0);
  scene -> Write_Output_Files(1);
}

