set(ITKIOXML_SRC
  itkXMLFile.cxx
  itkDOMNode.cxx
  itkDOMNodeXMLReader.cxx
  itkDOMNodeXMLWriter.cxx
  itkStringTools.cxx
  itkFancyString.cxx
)

add_library(ITKIOXML ${ITKIOXML_SRC})
itk_module_link_dependencies()
itk_module_target(ITKIOXML)
