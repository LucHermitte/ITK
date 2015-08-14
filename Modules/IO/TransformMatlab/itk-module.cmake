set(DOCUMENTATION "This module contains the classes for the input and output
of itkTransform object in  Matlab format.")

itk_module(ITKIOTransformMatlab
  PRIVATE_DEPENDS
    ITKIOTransformBase
  TEST_DEPENDS
    ITKTestKernel
    ITKIOTransformBase
  DESCRIPTION
    "${DOCUMENTATION}"
)
