cmake_minimum_required(VERSION 3.15)
project(SIEVE_simulator VERSION 1.3.0)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_XOPEN_SOURCE=500")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g -D_XOPEN_SOURCE=500")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

# require the presence of some libraries
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# for test purpose, use "test" after TARGET_GROUP
# for running purpose, use "production" after TARGET_GROUP
set(TARGET_GROUP production CACHE STRING "Group to build")

set(EXE_NAME SIEVE_simulator)

# set directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/external)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# the main entrance
if (TARGET_GROUP STREQUAL production)
    add_executable(${EXE_NAME}
            ${SRC_DIR}/cellcoal.c
            ${SRC_DIR}/eigen.c
            ${SRC_DIR}/hashMap.c
            ${SRC_DIR}/signatures.c
            )
    target_link_libraries(${EXE_NAME} PRIVATE Threads::Threads OpenMP::OpenMP_C)
    if (UNIX AND NOT APPLE)
        target_link_libraries(${EXE_NAME} PRIVATE -lm -ldl)
    endif ()
# elseif (TARGET_GROUP STREQUAL test)
else()
    message(FATAL_ERROR "Given TARGET_GROUP unknown")
endif()
