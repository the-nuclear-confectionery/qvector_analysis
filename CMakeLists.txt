cmake_minimum_required(VERSION 3.10)
project(v2calculator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core Hist MathCore Gpad RIO Tree)
include(${ROOT_USE_FILE})

set(SOURCES
    src/main.cpp
    src/qvector.cpp
    src/event.cpp
)

add_executable(v2calculator ${SOURCES})

target_include_directories(v2calculator PRIVATE
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(v2calculator ${ROOT_LIBRARIES})
