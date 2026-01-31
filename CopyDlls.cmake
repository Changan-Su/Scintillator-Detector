# Copy Geant4 and Qt DLLs to exe dir so it can run standalone
# Usage: cmake -DDEST_DIR=<path> -DG4_BIN=<path> -DQT_BIN=<path> -P CopyDlls.cmake

if(NOT DEST_DIR)
  message(FATAL_ERROR "DEST_DIR not set")
endif()

# Copy Geant4 DLLs
if(G4_BIN AND EXISTS "${G4_BIN}")
  file(GLOB G4_DLLS "${G4_BIN}/G4*.dll"
       "${G4_BIN}/msvcp*.dll" "${G4_BIN}/vcruntime*.dll" "${G4_BIN}/concrt*.dll")
  foreach(_dll ${G4_DLLS})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_dll}" "${DEST_DIR}"
      RESULT_VARIABLE _r)
    if(_r)
      message(WARNING "Failed to copy ${_dll}")
    endif()
  endforeach()
  list(LENGTH G4_DLLS _n)
  message(STATUS "Copied ${_n} Geant4 DLLs to ${DEST_DIR}")
endif()

# Copy Qt DLLs (exampleB1 uses G4visQt3D which depends on Qt)
if(QT_BIN AND EXISTS "${QT_BIN}")
  file(GLOB QT_DLLS "${QT_BIN}/Qt5*.dll" "${QT_BIN}/D3Dcompiler_*.dll"
       "${QT_BIN}/opengl32sw.dll" "${QT_BIN}/libEGL*.dll" "${QT_BIN}/libGLES*.dll")
  foreach(_dll ${QT_DLLS})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_dll}" "${DEST_DIR}"
      RESULT_VARIABLE _r OUTPUT_QUIET ERROR_QUIET)
  endforeach()
  # Copy platforms plugin (required for Qt GUI)
  set(PLATFORMS_SRC "${QT_BIN}/../plugins/platforms")
  set(PLATFORMS_DST "${DEST_DIR}/platforms")
  if(EXISTS "${PLATFORMS_SRC}/qwindows.dll")
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${PLATFORMS_DST}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${PLATFORMS_SRC}/qwindows.dll" "${PLATFORMS_DST}/")
    message(STATUS "Qt DLLs and platforms plugin copied to ${DEST_DIR}")
  endif()
endif()
