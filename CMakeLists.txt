cmake_minimum_required(VERSION 3.8)
project(necsimCMD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(SOURCE_DIR_NECSIM ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SOURCE_DIR_NECSIM})
set(SOURCE_FILES
        ${SOURCE_DIR_NECSIM}/ConfigParser.cpp
        ${SOURCE_DIR_NECSIM}/DispersalCoordinator.cpp
        ${SOURCE_DIR_NECSIM}/Logger.cpp
        ${SOURCE_DIR_NECSIM}/Logging.cpp
        ${SOURCE_DIR_NECSIM}/LogFile.cpp
        ${SOURCE_DIR_NECSIM}/LogFile.h
        ${SOURCE_DIR_NECSIM}/main.cpp
        ${SOURCE_DIR_NECSIM}/Landscape.cpp
        ${SOURCE_DIR_NECSIM}/ProtractedTree.cpp
        ${SOURCE_DIR_NECSIM}/setup.cpp
        ${SOURCE_DIR_NECSIM}/Community.cpp
        ${SOURCE_DIR_NECSIM}/ActivityMap.cpp
        ${SOURCE_DIR_NECSIM}/Metacommunity.cpp
        ${SOURCE_DIR_NECSIM}/AnalyticalSpeciesAbundancesHandler.cpp
        ${SOURCE_DIR_NECSIM}/SimulatedSpeciesAbundancesHandler.cpp
        ${SOURCE_DIR_NECSIM}/SpeciesAbundancesHandler.cpp
        ${SOURCE_DIR_NECSIM}/file_system.cpp
        ${SOURCE_DIR_NECSIM}/DataMask.cpp
        ${SOURCE_DIR_NECSIM}/DataPoint.cpp
        ${SOURCE_DIR_NECSIM}/SpeciesList.cpp
        ${SOURCE_DIR_NECSIM}/TreeNode.cpp
        ${SOURCE_DIR_NECSIM}/SpatialTree.cpp
        ${SOURCE_DIR_NECSIM}/SQLiteHandler.cpp
        ${SOURCE_DIR_NECSIM}/Tree.cpp
        ${SOURCE_DIR_NECSIM}/cpl_custom_handler.cpp
        ${SOURCE_DIR_NECSIM}/custom_exceptions.h
        ${SOURCE_DIR_NECSIM}/double_comparison.cpp
        ${SOURCE_DIR_NECSIM}/neutral_analytical.cpp
        ${SOURCE_DIR_NECSIM}/parameters.cpp
        ${SOURCE_DIR_NECSIM}/Xoroshiro256plus.h
        ${SOURCE_DIR_NECSIM}/RNGController.h
        ${SOURCE_DIR_NECSIM}/Matrix.h
        ${SOURCE_DIR_NECSIM}/Map.h
        ${SOURCE_DIR_NECSIM}/SimulationTemplates.h
        ${SOURCE_DIR_NECSIM}/SimParameters.h
        ${SOURCE_DIR_NECSIM}/GillespieCalculator.cpp
        ${SOURCE_DIR_NECSIM}/GenericTree.h
        ${SOURCE_DIR_NECSIM}/MapLocation.h
        ${SOURCE_DIR_NECSIM}/MapLocation.cpp
        )

find_package(Boost 1.5.7 COMPONENTS system filesystem REQUIRED)
if (UNIX)
    find_package(GDAL)
else ()
    find_package(GDAL 2.1.0 REQUIRED)
endif ()
find_library(SQL_DIR sqlite3)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GDAL_INCLUDE_DIR})
include_directories(/usr/local/include)
include_directories(../)
link_directories(/usr/local/lib)
add_executable(necsimCMD ${SOURCE_FILES})
#add_definitions(-Dwith_gdal -Dverbose)
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG=1)
endif ()
target_link_libraries(necsimCMD "${SQL_DIR}")
target_link_libraries(necsimCMD gdal)

target_link_libraries(necsimCMD ${Boost_LIBRARIES})
