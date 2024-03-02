cmake_minimum_required(VERSION 3.15)
project(Project)
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-fexec-charset=GBK)
set (CMAKE_EXE_LINKER_FLAGS \"-static\")
set(BinPath ${PROJECT_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${BinPath})

aux_source_directory(${PROJECT_SOURCE_DIR}/src cpp)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/link_prediciton lp)

set(Include ${PROJECT_SOURCE_DIR}/include)
include_directories(${Include})

add_executable(main main.cpp ${lp} ${cpp})
