cmake_minimum_required(VERSION 3.10.2)
project(MinimalPathExtractionExamples)

find_package(ITK REQUIRED
  COMPONENTS MinimalPathExtraction
    ITKIOImageBase
    ITKImageFunction
    ITKOptimizers
    ITKPath
    ITKIOMeta
    ITKIOPNG
    ITKIOJPEG
    )

include(${ITK_USE_FILE})

add_executable(MinimalPathExamples example.cxx)
target_link_libraries(MinimalPathExamples ${ITK_LIBRARIES})

add_executable(interptest interpolator.cxx)
target_link_libraries(interptest ${ITK_LIBRARIES})
