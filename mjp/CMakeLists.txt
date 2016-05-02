
cmake_minimum_required(VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
set(PROJECT_NAME mjp)
project(${PROJECT_NAME})

# Determine the platform.
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  set(OS_LINUX 1)
  set(OS_POSIX 1)
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
  set(OS_MACOSX 1)
  set(OS_POSIX 1)  
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
#    'cmake -G "Visual Studio 12 Win64" ..
  set(OS_WINDOWS 1)
endif()

if (OS_LINUX)
  if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11 -fPIC -ggdb -O0 -Wall -Wno-write-strings -Wno-unused-local-typedefs)
    # add_definitions(-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
  endif()

  add_definitions(-D_LINUX -D_x64)
elseif(OS_WINDOWS)
  macro(get_WIN32_WINNT version)
    if (WIN32 AND CMAKE_SYSTEM_VERSION)
    set(ver ${CMAKE_SYSTEM_VERSION})
    string(REPLACE "." "" ver ${ver})
    string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})  
    set(${version} "0x${ver}")
    endif()
  endmacro()

  get_WIN32_WINNT(ver)
  add_definitions(-D_WIN32_WINNT=${ver} -D_x64 -D_CRT_SECURE_NO_WARNINGS)
  
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
endif()


set(DRAGONFLY_ROOT ..)
set(MJP_DIR ${DRAGONFLY_ROOT}/mjp)
set(UFI_DIR ${DRAGONFLY_ROOT}/ufi)
set(BASE_DIR ${DRAGONFLY_ROOT}/base)

# Access *custom* cmake find macros.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

FIND_PACKAGE(PythonInterp 2.7)
FIND_PACKAGE(PythonLibs 2.7)

include_directories(${PYTHON_INCLUDE_DIR})

FIND_PACKAGE(NumPy)
INCLUDE_DIRECTORIES(${PYTHON_NUMPY_INCLUDE_DIR})
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

include_directories(BEFORE 
  ${CMAKE_SOURCE_DIR}
  ${DRAGONFLY_ROOT}/base
  ${DRAGONFLY_ROOT}/ufi
  ${DRAGONFLY_ROOT}/mjp
  ${DRAGONFLY_ROOT}/third_party/include
)

link_directories(BEFORE
#  /opt/googletest/lib
)

set (SOURCES ${SOURCES} 
  ${MJP_DIR}/ND_matrix.h
  ${MJP_DIR}/ND_matrix_inline.h
  ${MJP_DIR}/ND_matrix_template.h
  ${MJP_DIR}/ND_matrix.cpp
  ${MJP_DIR}/particle_packing.h
  ${MJP_DIR}/particle_packing_inline.h
  ${MJP_DIR}/particle_packing_template.h
  ${MJP_DIR}/particle_packing.cpp

  ${UFI_DIR}/python_util.h
  ${UFI_DIR}/python_util_inline.h
  ${UFI_DIR}/python_util_template.h
  ${UFI_DIR}/python_util.cpp
  ${UFI_DIR}/python_util_module.h
  ${UFI_DIR}/python_util_module_template.h
  ${UFI_DIR}/python_util_module.cpp
  ${UFI_DIR}/scanner_util.h
  ${UFI_DIR}/scanner_util.cpp
  ${UFI_DIR}/scanner_util_module.cpp
  ${UFI_DIR}/scanner_util_module.h

# *Temporarily* remove "commUtil" dependency.
# (Will be upgraded to use "std::iostream" instead of current handle-based system.)
#  ${BASE_DIR}/commUtil.h
#  ${BASE_DIR}/commUtil_inline.h
#  ${BASE_DIR}/commUtil_template_forward.h
#  ${BASE_DIR}/commUtil_template.h
#  ${BASE_DIR}/commUtil.cpp
  ${BASE_DIR}/commUtil_stub.h
  ${BASE_DIR}/commUtil_stub_template.h
  ${BASE_DIR}/commUtil_stub.cpp

  ${BASE_DIR}/ghostIterator.h
  ${BASE_DIR}/ghostIterator_template.h
  ${BASE_DIR}/gmm_ext.h
  ${BASE_DIR}/linalgUtil.h
  ${BASE_DIR}/linalgUtil_inline.h
  ${BASE_DIR}/linalgUtil_template.h
  ${BASE_DIR}/linalgUtil.cpp

  ${BASE_DIR}/numberTraits.h

  ${BASE_DIR}/factorial.h
  ${BASE_DIR}/factorial.cpp
  ${BASE_DIR}/ntuple.h
  ${BASE_DIR}/ntuple_template.h
  ${BASE_DIR}/parallelUtil.h
  ${BASE_DIR}/parallelUtil_inline.h
  ${BASE_DIR}/parallelUtil_template.h
  ${BASE_DIR}/parallelUtil.cpp
)
add_definitions(-DEMBEDDED_PYTHON)

if (OS_LINUX)
  set (LIBS ${LIBS}
    pthread
    gomp
    rt
#    gtest
  )
elseif(OS_WINDOWS)
  set (LIBS ${LIBS}
#    gtest
  )
endif()

# MESSAGE( ${LIBS} )
add_executable(${PROJECT_NAME} ${SOURCES})

if(OS_WINDOWS)
  target_link_libraries(${PROJECT_NAME}
    ${LIBS}
    ${PYTHON_LIBRARY}
  )

else()
  target_link_libraries(${PROJECT_NAME}
    ${LIBS}
    ${PYTHON_LIBRARY}    
  )
endif()

