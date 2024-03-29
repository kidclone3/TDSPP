cmake_minimum_required(VERSION 3.27)
project(TDSPP)

set(CMAKE_CXX_STANDARD 17)
set(BOOST_DIR "/home/delus/Downloads/boost_1_73_0/libs")

## Link to boost library folder
link_directories(${BOOST_DIR})

# Add flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra")


# Disable flags when debugging
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Set working directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(TDSPP
        TDSPP.cpp
        Output.cpp
        Output.h
        SummaryOutput.cpp
        SummaryOutput.h
        Graph.cpp
        Graph.h)
#        TDSPP_8_19.cpp)
