cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
set(MODULE_NAME "LITTPlanV2") # Do not use 'project()'
set(MODULE_TITLE "${MODULE_NAME}")

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/${MODULE_NAME}")
  set(MODULE_CATEGORY "Examples")
  set(MODULE_STATUS "Beta")
  set(MODULE_CONTRIBUTORS "Erol Yeniaras (MD Anderson Cancer Center), Jean-Christophe Fillion-Robin (Kitware)")
  set(MODULE_DESCRIPTION "This is ...")
  set(MODULE_ICONURL "http://http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Base/QTGUI/Resources/Icons/ExtensionDefaultIcon.png?revision=19432&view=co")
  set(MODULE_SCREENSHOTURLS "")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
set(MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT")

set(MODULE_INCLUDE_DIRECTORIES
  )

set(MODULE_SRCS
  qSlicerLITTPlanV2IO.cxx
  qSlicerLITTPlanV2IO.h
  qSlicerLITTPlanV2Module.cxx
  qSlicerLITTPlanV2Module.h
  qSlicerLITTPlanV2ModuleWidget.cxx
  qSlicerLITTPlanV2ModuleWidget.h
  )

set(MODULE_MOC_SRCS
  qSlicerLITTPlanV2IO.h
  qSlicerLITTPlanV2Module.h
  qSlicerLITTPlanV2ModuleWidget.h
  )

set(MODULE_UI_SRCS
  Resources/UI/qSlicerLITTPlanV2Module.ui
  )

set(MODULE_TARGET_LIBRARIES
  vtkSlicerAnnotationsModuleMRML
  )

set(MODULE_RESOURCES
  Resources/qSlicerLITTPlanV2Module.qrc
  )

#-----------------------------------------------------------------------------
slicerMacroBuildQtModule(
  NAME ${MODULE_NAME}
  TITLE ${MODULE_TITLE}
  EXPORT_DIRECTIVE ${MODULE_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  SRCS ${MODULE_SRCS}
  MOC_SRCS ${MODULE_MOC_SRCS}
  UI_SRCS ${MODULE_UI_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RESOURCES ${MODULE_RESOURCES}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
