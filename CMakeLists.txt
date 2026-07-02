# Define the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name
project(ClusteringProject CXX)

# Set the C++ standard to 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- LCG / cvmfs dependency discovery (e.g. UChicago AF under `lsetup root`) ---
# `lsetup` puts each LCG package's lib dir on LD_LIBRARY_PATH but does NOT set
# CMAKE_PREFIX_PATH, so find_package can't locate ROOT's transitive deps (Vdt)
# or the matching Boost library.  Derive each package prefix from the cvmfs lib
# dirs on LD_LIBRARY_PATH and prepend them.  No-op on non-cvmfs machines (e.g.
# the macOS dev box, where LD_LIBRARY_PATH has no /cvmfs/.../lib entries).
if(DEFINED ENV{LD_LIBRARY_PATH})
  string(REPLACE ":" ";" _ld_dirs "$ENV{LD_LIBRARY_PATH}")
  foreach(_d ${_ld_dirs})
    if(_d MATCHES "^/cvmfs/.*/lib$")
      get_filename_component(_prefix "${_d}" DIRECTORY)
      list(APPEND CMAKE_PREFIX_PATH "${_prefix}")
    endif()
  endforeach()
  if(CMAKE_PREFIX_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
  endif()
endif()

# Find the ROOT and Boost libraries
# Assumes ROOT and Boost are installed and discoverable by find_package
# Imt: needed for TTreeProcessorMT (util/ttreemt_prototype.cxx). EnableImplicitMT()
# already works without listing it explicitly (pulled in transitively), but list it
# so the dependency is correct against other ROOT configs.
find_package(ROOT REQUIRED COMPONENTS RIO Core Gpad Graf Hist Tree Imt)
# TMVA excluded: HGTD_SORT/CONE_BDT are disabled. Re-add when re-enabling those scores.
find_package(Boost COMPONENTS filesystem REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -g")

# Export compile commands so clangd (lsp-mode, clangd, etc.) can resolve
# ROOT and Boost include paths without manual configuration.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add your executables, specifying the path to the source files
add_executable(clustering_dt src/clustering_dt.cxx)
add_executable(generate_rpt util/generate_rpt.cxx)
add_executable(rpt_v2 util/rpt_v2.cxx)
add_executable(rpt_v3 util/rpt_v3.cxx)
add_executable(rpt_v4 util/rpt_v4.cxx)
add_executable(rpt_v5 util/rpt_v5.cxx)
add_executable(export_training_data util/export_training_data.cxx)
add_executable(failure_decomposition util/failure_decomposition.cxx)
add_executable(maxpt_jet_test util/maxpt_jet_test.cxx)
add_executable(track_dt util/track_dt.cxx)
add_executable(ttreemt_prototype util/ttreemt_prototype.cxx)

# Specify the include directories for your header files
foreach(_tgt
    clustering_dt
    generate_rpt
    rpt_v2
    rpt_v3
    rpt_v4
    rpt_v5
    export_training_data
    failure_decomposition
    maxpt_jet_test
    track_dt
    ttreemt_prototype
)
    target_include_directories(${_tgt}
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${ROOT_INCLUDE_DIRS}
    )
    target_link_libraries(${_tgt}
        PRIVATE
            ${ROOT_LIBRARIES}
            Boost::filesystem
    )
endforeach()
