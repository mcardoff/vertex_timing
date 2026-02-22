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
add_executable(generate_rpt src/generate_rpt.cxx)
add_executable(test_vbf_selection test_vbf_selection.cxx)
add_executable(evaluate_ml_lowtrack src/evaluate_ml_lowtrack.cxx)
add_executable(export_training_data src/export_training_data.cxx)
add_executable(distcut_sweep src/distcut_sweep.cxx)

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

target_include_directories(generate_rpt
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(test_vbf_selection
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(evaluate_ml_lowtrack
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
	${ROOT_INCLUDE_DIRS}
)

target_include_directories(export_training_data
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${ROOT_INCLUDE_DIRS}
)

target_include_directories(distcut_sweep
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

target_link_libraries(generate_rpt
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(test_vbf_selection
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(evaluate_ml_lowtrack
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(export_training_data
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)

target_link_libraries(distcut_sweep
    PRIVATE
        ${ROOT_LIBRARIES}
        Boost::filesystem
)
