cmake_minimum_required(VERSION 3.10)
project(ana)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- ROOT ---
find_package(ROOT REQUIRED COMPONENTS Core Hist MathCore Gpad RIO Tree)
include(${ROOT_USE_FILE})

# --- YAML-CPP ---
find_package(yaml-cpp REQUIRED)

set(SOURCES
    src/main.cpp
    src/event.cpp
    src/observables.cpp
)

add_executable(ana ${SOURCES})

target_include_directories(ana PRIVATE
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(ana
    ${ROOT_LIBRARIES}
    yaml-cpp
)
