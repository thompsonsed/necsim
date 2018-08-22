cmake_minimum_required(VERSION 3.8)
project(NECSimCMD)

set(CMAKE_CXX_STANDARD 14)
set(SOURCE_DIR_NECSIM ${CMAKE_CURRENT_SOURCE_DIR})
include_directories( ${SOURCE_DIR_NECSIM})
set(SOURCE_FILES
        ${SOURCE_DIR_NECSIM}/ConfigFileParser.cpp
        ${SOURCE_DIR_NECSIM}/DispersalCoordinator.cpp
        ${SOURCE_DIR_NECSIM}/Logger.cpp
        ${SOURCE_DIR_NECSIM}/Logging.cpp
        ${SOURCE_DIR_NECSIM}/LogFile.cpp
        ${SOURCE_DIR_NECSIM}/LogFile.h
        ${SOURCE_DIR_NECSIM}/main.cpp
        ${SOURCE_DIR_NECSIM}/Landscape.cpp
        ${SOURCE_DIR_NECSIM}/ProtractedTree.cpp
        ${SOURCE_DIR_NECSIM}/Setup.cpp
        ${SOURCE_DIR_NECSIM}/Community.cpp
        ${SOURCE_DIR_NECSIM}/ReproductionMap.cpp
        ${SOURCE_DIR_NECSIM}/Metacommunity.cpp
        ${SOURCE_DIR_NECSIM}/Filesystem.cpp
        ${SOURCE_DIR_NECSIM}/DataMask.cpp
        ${SOURCE_DIR_NECSIM}/DataPoint.cpp
        ${SOURCE_DIR_NECSIM}/SpeciesList.cpp
        ${SOURCE_DIR_NECSIM}/TreeNode.cpp
        ${SOURCE_DIR_NECSIM}/SpatialTree.cpp
        ${SOURCE_DIR_NECSIM}/Tree.cpp
        ${SOURCE_DIR_NECSIM}/CPLCustomHandlerNecsim.cpp
        ${SOURCE_DIR_NECSIM}/DoubleComparison.cpp
        ${SOURCE_DIR_NECSIM}/NRrand.h
        ${SOURCE_DIR_NECSIM}/Matrix.h
        ${SOURCE_DIR_NECSIM}/SimulationTemplates.h)

find_package(Boost 1.5.7 COMPONENTS system filesystem REQUIRED)
find_library(gdal 2.1.0 REQUIRED)
find_library(SQL_DIR sqlite3)
include_directories(${Boost_INCLUDE_DIR})
include_directories(/usr/local/include)
include_directories(../)
link_directories(/usr/local/lib)
add_executable(NECSimCMD ${SOURCE_FILES})
add_definitions(-Dwith_gdal -Dverbose)
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG=1)
endif()
target_link_libraries(NECSimCMD "${SQL_DIR}")
target_link_libraries(NECSimCMD gdal)

target_link_libraries(NECSimCMD ${Boost_LIBRARIES})