cmake_minimum_required(VERSION 2.8.5)
project(FortranExercise)
enable_language(Fortran)
# Release会采用-O3优化
set(CMAKE_BUILD_TYPE Debug)
 
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
endif()

#这当前项目里lib目录下的所有文件(包含子目录下的文件)添加到SRCS这个变量里
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} lib/*)
# 将SRCS编译为mylib这个动态链接库
add_library(mylib SHARED ${SRCS})
 
file(GLOB_RECURSE SRCS RELATIVE ${PROJECT_SOURCE_DIR} src/*)
add_executable(myprogram ${SRCS})
# 让myprogram依赖mylib
add_dependencies(myprogram mylib)
# 链接mylib
target_link_libraries(myprogram mylib)

# 参考链接 https://www.lainme.com/doku.php/blog/2015/12/%E7%94%A8cmake%E6%9D%A5%E7%BC%96%E8%AF%91fortran
