cmake_minimum_required(VERSION 2.8)
message(STATUS "Start running cmake...")

enable_language(CXX)
if (UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif ()

if (LINUX)
    message(STATUS "Linux Env.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
elseif (APPLE)
    message(STATUS "Apple Env.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -O3")
endif ()

find_package(Boost REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Boost_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Boost_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Boost_EXE_LINKER_FLAGS}")
endif ()



add_subdirectory(source)
project(SimplForman3D CXX)
add_executable(${PROJECT_NAME} "source/main.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} libFORMAN)
