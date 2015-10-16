cmake_minimum_required(VERSION 2.8)
message( "\nConfiguring tools:" )

# dependencies
find_package( Cplex REQUIRED )

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include/
)

# autodiscover tool sources and add tools
file(GLOB TOOL_SRCS *.cpp)
foreach(tool_src ${TOOL_SRCS})
  
    get_filename_component(tool_name ${tool_src} NAME_WE)
    add_executable( ${tool_name} ${tool_src} )
    target_link_libraries( ${tool_name} pgmlink ${Boost_LIBRARIES} ${CPLEX_LIBRARIES} ${LEMON_LIBRARIES} ${HDF5_LIBRARIES})
  
endforeach(tool_src)
