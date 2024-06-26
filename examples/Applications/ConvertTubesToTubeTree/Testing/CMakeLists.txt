##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

include_regular_expression( "^.*$" )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1 - I-Junction
itk_add_test(
  NAME ${MODULE_NAME}Test1
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Input-IJunc.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Tree.tre )

# Test1 - Compare - I-Junction
itk_add_test(
  NAME ${MODULE_NAME}Test1-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Tree.tre
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Tree.tre}
    -d 0.01 )

set_tests_properties( ${MODULE_NAME}Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}Test1 )

# Test2 - Y-Junction
itk_add_test(
  NAME ${MODULE_NAME}Test2
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2-Input-YJunc.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2-Tree.tre )

# Test2 - Compare - Y-Junction
itk_add_test(
  NAME ${MODULE_NAME}Test2-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2-Tree.tre
    -B 2 DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2-Tree.tre} DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2-Tree2.tre}
    -d 0.01 )

set_tests_properties( ${MODULE_NAME}Test2-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}Test2 )
