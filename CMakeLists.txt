#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.20)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
#set (CMAKE_CXX_COMPILER_FORCED ON)
#using intel compiler
set (CMAKE_CXX_COMPILER "mpicxx")

#set as inlcude dir the current directory
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(project_sources
		src/TBG_LMC.cpp
		src/StandardFunctions.cpp
		src/Functions.cpp
        )

set(header_files
		include/Constants.h
		include/StandardFunctions.h
		include/Functions.h
		)

include_directories("./include")

#add compiler flags
set(CXX_flags "-O2 -std=c++11 -D 'MKL_ILP64'")
set(CMAKE_CXX_FLAGS	"${CMAKE_CXX_FLAGS} ${CXX_flags}")
#add MKL libraries
set(MKLROOT "/share/apps/spack/envs/fgci-centos7-generic/software/intel-mkl/2019.3.199/li7bk52/compilers_and_libraries_2019.3.199/linux/mkl")
set(MKL_LIBRARY_FLAGS "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ") 
set(MKL_LIBRARY_FLAGS "${MKL_LIBRARY_FLAGS} ${MKLROOT}/lib/intel64/libmkl_intel_thread.a")
set(MKL_LIBRARY_FLAGS "${MKL_LIBRARY_FLAGS} ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl")


#set project name
project(TBG_LMC)
#add an executable - using variables
add_executable(${PROJECT_NAME}.out ${project_sources} ${header_files})

target_link_libraries(TBG_LMC.out PUBLIC ${MKL_LIBRARY_FLAGS})
                                   

