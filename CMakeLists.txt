cmake_minimum_required(VERSION 3.15)
project(mySimulator VERSION 1.3.0)

#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_XOPEN_SOURCE=500")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g -D_XOPEN_SOURCE=500")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

# require the presence of some libraries
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(dlib REQUIRED)
#find_package(Boost REQUIRED)

#include_directories(${DLIB_INCLUDE_DIRECTORY})
#include_directories(${Boost_INCLUDE_DIRS})

#message(STATUS "dlib included from ${DLIB_INCLUDE_DIRECTORY}")
#message(STATUS "boost included from ${Boost_INCLUDE_DIRS}")

# for test purpose, use "test" after TARGET_GROUP
# for running purpose, use "production" after TARGET_GROUP
set(TARGET_GROUP production CACHE STRING "Group to build")

set(EXE_NAME SIEVE_simulator)

# set directories
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
#set(TEST_DIR ${CMAKE_SOURCE_DIR}/test)
set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/external)
#set(TEST_DIR ${CMAKE_SOURCE_DIR}/test)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# the main entrance
if (TARGET_GROUP STREQUAL production)
    #    add_library(statisticalTools SHARED ${SRC_DIR}/statisticalTools.cpp)
    #    target_compile_options(statisticalTools PRIVATE -fPIC)
    #    target_link_libraries(statisticalTools PRIVATE dlib::dlib Boost::boost)

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
elseif (TARGET_GROUP STREQUAL test)
    #    include(CTest)
    #
    #    add_subdirectory(${EXTERNAL_DIR})
    #    add_subdirectory(${TEST_DIR})
    #
    #    add_library(Unity STATIC
    #            ${EXTERNAL_DIR}/Unity-2.4.3/src/unity.c
    #            )
    #
    #    target_include_directories(Unity PUBLIC
    #            ${EXTERNAL_DIR}/Unity-2.4.3/src
    #            )
    #
    #    # test 1
    #    add_executable(testNumOverlapMDI
    #            ${TEST_DIR}/testNumOverlapMDI.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testNumOverlapMDI PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testNumOverlapMDI PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test1, testNumOverlapMDI)
    #    # end test 1
    #
    #    # test 2
    #    #    add_executable(testNumOverlapMD
    #    #            ${TEST_DIR}/testNumOverlapMD.c
    #    #            ${SRC_DIR}/hashMap.c
    #    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #    #            )
    #    #
    #    #    if (UNIX AND NOT APPLE)
    #    #        target_link_libraries(testNumOverlapMD PRIVATE Unity -lm -ldl)
    #    #    else()
    #    #        target_link_libraries(testNumOverlapMD PRIVATE Unity)
    #    #    endif ()
    #    #
    #    #    add_test(test2, testNumOverlapMD)
    #    # end test 2
    #
    #    # test 3
    #    add_executable(testGetMuDelInsNr_1
    #            ${TEST_DIR}/testGetMuDelInsNr_1.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetMuDelInsNr_1 PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetMuDelInsNr_1 PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test3, testGetMuDelInsNr_1)
    #    # end test 3
    #
    #    # test 4
    #    add_executable(testGetMuDelInsNr_2
    #            ${TEST_DIR}/testGetMuDelInsNr_2.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetMuDelInsNr_2 PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetMuDelInsNr_2 PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test4, testGetMuDelInsNr_2)
    #    # end test 4
    #
    #    # test 5
    #    add_executable(testGetMuDelInsNr_3
    #            ${TEST_DIR}/testGetMuDelInsNr_3.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetMuDelInsNr_3 PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetMuDelInsNr_3 PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test5, testGetMuDelInsNr_3)
    #    # end test 5
    #
    #    # test 6
    #    add_executable(testGetMuDelInsNr_4
    #            ${TEST_DIR}/testGetMuDelInsNr_4.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetMuDelInsNr_4 PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetMuDelInsNr_4 PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test6, testGetMuDelInsNr_4)
    #    # end test 6
    #
    #    # test 7
    #    #    add_executable(testGetMuDelNr
    #    #            ${TEST_DIR}/testGetMuDelNr.c
    #    #            ${SRC_DIR}/hashMap.c
    #    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #    #            )
    #    #
    #    #if (UNIX AND NOT APPLE)
    #    #        target_link_libraries(testGetMuDelNr PRIVATE Unity -lm -ldl)
    #    #    else()
    #    #        target_link_libraries(testGetMuDelNr PRIVATE Unity)
    #    #    endif ()
    #
    #    #    add_test(test7, testGetMuDelNr)
    #    # end test 7
    #
    #    # test 8
    #    add_executable(testGetBinaryKeyMDI
    #            ${TEST_DIR}/testGetBinaryKeyMDI.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetBinaryKeyMDI PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetBinaryKeyMDI PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test8, testGetBinaryKeyMDI)
    #    # end test 8
    #
    #    # test 9
    #    #    add_executable(testGetBinaryKeyMD
    #    #            ${TEST_DIR}/testGetBinaryKeyMD.c
    #    #            ${SRC_DIR}/hashMap.c
    #    #            ${SRC_DIR}/finiteMuInDelSimulator.c
    #    #            )
    #    #
    #    #if (UNIX AND NOT APPLE)
    #    #        target_link_libraries(testGetBinaryKeyMD PRIVATE Unity -lm -ldl)
    #    #    else()
    #    #        target_link_libraries(testGetBinaryKeyMD PRIVATE Unity)
    #    #    endif ()
    #
    #    #    add_test(9, testGetBinaryKeyMD)
    #    # end test 9
    #
    #    # test 10
    #    #    add_executable(testGetValIntArrToIntMD
    #    #            ${TEST_DIR}/testGetValIntArrToIntMD.c
    #    #            ${SRC_DIR}/hashMap.c
    #    #            ${SRC_DIR}/finiteMuInDelSimulator.c
    #    #            )
    #    #
    #    #   if (UNIX AND NOT APPLE)
    #    #        target_link_libraries(testGetValIntArrToIntMD PRIVATE Unity -lm -ldl)
    #    #    else()
    #    #        target_link_libraries(testGetValIntArrToIntMD PRIVATE Unity)
    #    #    endif ()
    #
    #    #    add_test(test10, testGetValIntArrToIntMD)
    #    # end test 10
    #
    #    # test 11
    #    add_executable(testGetValIntArrToIntMDI
    #            ${TEST_DIR}/testGetValIntArrToIntMDI.c
    #            ${SRC_DIR}/hashMap.c
    #            ${SRC_DIR}/finiteMuInDelSimulator.h
    #            )
    #
    #    if (UNIX AND NOT APPLE)
    #        target_link_libraries(testGetValIntArrToIntMDI PRIVATE Unity -lm -ldl)
    #    else ()
    #        target_link_libraries(testGetValIntArrToIntMDI PRIVATE Unity)
    #    endif ()
    #
    #    add_test(test11, testGetValIntArrToIntMDI)
    # end test 11
else()
    message(FATAL_ERROR "Given TARGET_GROUP unknown")
endif()
