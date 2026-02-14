#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(B1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
find_package(Geant4 REQUIRED ui_all vis_all analysis)
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)



#----------------------------------------------------------------------------
# Add the executable, use our local headers, and link it to the Geant4 libraries
#
add_executable(exampleB1 exampleB1.cc ${sources} ${headers})
target_include_directories(exampleB1 PRIVATE include)
target_link_libraries(exampleB1 PRIVATE ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  exampleB1.in
  exampleB1.out
  geometry.mac
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# Copy macros from SOURCE to executable dir so edits to .mac take effect without re-running cmake
add_custom_command(TARGET exampleB1 POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${PROJECT_SOURCE_DIR}/geometry.mac
    ${PROJECT_SOURCE_DIR}/init_vis.mac
    ${PROJECT_SOURCE_DIR}/vis.mac
    ${PROJECT_SOURCE_DIR}/run1.mac
    ${PROJECT_SOURCE_DIR}/run2.mac
    $<TARGET_FILE_DIR:exampleB1>
  COMMENT "Copying macro files from project source to executable directory"
)

# Copy Geant4 and Qt DLLs so exampleB1.exe can run directly (double-click)
get_filename_component(GEANT4_BIN "${Geant4_DIR}/../../../bin" ABSOLUTE)
find_package(Qt5 COMPONENTS Core QUIET)
set(QT_BIN_FOR_DEPLOY "")
if(Qt5_FOUND)
  get_filename_component(QT_BIN_FOR_DEPLOY "${Qt5Core_DIR}/../../../bin" ABSOLUTE)
endif()
add_custom_command(TARGET exampleB1 POST_BUILD
  COMMAND ${CMAKE_COMMAND}
    "-DDEST_DIR=$<TARGET_FILE_DIR:exampleB1>"
    "-DG4_BIN=${GEANT4_BIN}"
    "-DQT_BIN=${QT_BIN_FOR_DEPLOY}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/CopyDlls.cmake"
  COMMENT "Copying Geant4 and Qt DLLs to executable directory"
)
