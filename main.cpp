#include "Snow_Driver.h"
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>

using namespace PhysBAM;

int main(int argc,char *argv[]){
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
  
    Snow_Driver<T> driver(stream_type);
    
    driver.Run();
  
}

