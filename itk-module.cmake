# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.rst" DOCUMENTATION)

itk_module(MinimalPathExtraction
  ENABLE_SHARED
  DEPENDS
    ITKCommon
    ITKOptimizers
    ITKPath
    ITKFastMarching
  TEST_DEPENDS
    ITKTestKernel
    ITKIOSpatialObjects
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
