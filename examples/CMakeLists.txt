# The header files
SET(HEADERS
  ${CMAKE_SOURCE_DIR}/include/itkImageToPathFilter.h
  ${CMAKE_SOURCE_DIR}/include/itkSingleImageCostFunction.h
  ${CMAKE_SOURCE_DIR}/include/itkArrivalFunctionToPathFilter.h
  ${CMAKE_SOURCE_DIR}/include/itkSpeedFunctionToPathFilter.h
  ${CMAKE_SOURCE_DIR}/include/itkImageToPathFilter.hxx
  ${CMAKE_SOURCE_DIR}/include/itkSingleImageCostFunction.hxx
  ${CMAKE_SOURCE_DIR}/include/itkArrivalFunctionToPathFilter.hxx
  ${CMAKE_SOURCE_DIR}/include/itkSpeedFunctionToPathFilter.hxx
)

# Add this as include directory
INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
)

# Main library
ADD_EXECUTABLE(MinimalPathMain main.cxx ${HEADERS} )
TARGET_LINK_LIBRARIES(MinimalPathMain ${ITK_LIBRARIES})

ADD_EXECUTABLE(MinimalPathExamples example.cxx ${HEADERS} )
TARGET_LINK_LIBRARIES(MinimalPathExamples ${ITK_LIBRARIES})
