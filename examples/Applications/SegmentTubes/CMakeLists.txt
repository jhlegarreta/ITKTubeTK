##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
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

set( MODULE_NAME SegmentTubes )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${CMAKE_CURRENT_LIST_DIR}/../../../docs/TubeTKLogo.h
  TARGET_LIBRARIES
    TubeTK ${ITK_LIBRARIES}
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
