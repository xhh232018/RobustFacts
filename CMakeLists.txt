cmake_minimum_required(VERSION 3.0)
project(Perturbation)

SET(BOOST_ROOT ~/boost_1_63_0)
find_package(Boost 1.63.0 REQUIRED COMPONENTS serialization)
find_library(boost_serialization BOOST_ROOT)
set(CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS} )

#output dir
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

set(SRC_List GraphManager.h GraphManager.cpp QueryInterface.h QueryInterface.cpp QueryPerturbation.h QueryPerturbation.cpp DataPerturbation.h DataPerturbation.cpp)

add_executable(ptbClean main.cpp ${SRC_List})
target_link_libraries(ptbClean ${Boost_LIBRARIES})
