cmake_minimum_required(VERSION 3.27)
project(TDSPP)

set(CMAKE_CXX_STANDARD 17)

# Link to boost library folder
link_directories("/home/delus/Documents/cpp_library/boost_1_73_0/libs")

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra")

# Set working directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(TDSPP
        TDSPP.cpp)
#        TDSPP_8_19.cpp)
