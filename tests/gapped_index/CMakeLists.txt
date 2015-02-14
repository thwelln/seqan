# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /tests/gappedIndex/CMakeLists.txt
#
# CMakeLists.txt file for the gappedIndex module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_tests_gappedIndex)
message (STATUS "Configuring tests/gapped_index")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_gapped_index
                test_gapped_index.cpp
                test_gapped_index.h
                test_gapped_index_find.h
		test_gapped_index_stree.h
		test_gapped_index_construction.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_gappedIndex ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_gappedIndex COMMAND $<TARGET_FILE:test_gappedIndex>)
