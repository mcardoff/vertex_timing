# Define the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(ClusteringProject CXX)

# Set the C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find the ROOT and Boost libraries
# Assumes ROOT and Boost are installed and discoverable by find_package
find_package(ROOT REQUIRED COMPONENTS RIO Core Gpad Graf Hist Tree)
find_package(Boost COMPONENTS filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")

# Add your executables, specifying the path to the source files
add_executable(clustering_dt src/clustering_dt.cxx)
add_executable(hgtd_matching src/hgtd_matching.cxx)
add_executable(extract_error_metrics src/extract_error_metrics.cxx)
add_executable(test_ml_model src/test_ml_model.cxx)

# Specify the include directories for your header files
target_include_directories(clustering_dt
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(hgtd_matching
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(extract_error_metrics
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(test_ml_model
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

# Link the executables with the required libraries
# ROOT libraries are linked using the ${ROOT_LIBRARIES} variable
# Boost libraries are linked using Boost::filesystem
target_link_libraries(clustering_dt
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(hgtd_matching
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(extract_error_metrics
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(test_ml_model
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)
