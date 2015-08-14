set(DOCUMENTATION "This module contains the classes for the input and output
of itkTransform object in MINC format (.xfm).")

itk_module(ITKIOTransformMINC
  DEPENDS
    ITKIOTransformBase
    ITKMINC
  PRIVATE_DEPENDS
    ITKIOMINC
  TEST_DEPENDS
    ITKTestKernel
    ITKDisplacementField
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
)
