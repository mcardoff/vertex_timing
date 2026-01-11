# Define the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(ClusteringProject CXX)

# Find the ROOT and Boost libraries
# Assumes ROOT and Boost are installed and discoverable by find_package
find_package(ROOT REQUIRED COMPONENTS RIO Core Gpad Graf Hist Tree)
find_package(Boost COMPONENTS filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")

# Add your executable, specifying the path to the source file
add_executable(clustering_dt src/clustering_dt.cxx)

# Specify the include directories for your header files
target_include_directories(clustering_dt
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

# Link the executable with the required libraries
# ROOT libraries are linked using the ${ROOT_LIBRARIES} variable
# Boost libraries are linked using Boost::filesystem
target_link_libraries(clustering_dt
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)
