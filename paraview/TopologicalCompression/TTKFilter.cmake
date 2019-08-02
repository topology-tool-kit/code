# Allows to disable each filter
option(TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER "Build the TopologicalCompression filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER)

if(${TTK_BUILD_TOPOLOGICALCOMPRESSION_FILTER})
  ttk_register_pv_filter(ttkTopologicalCompression TopologicalCompression.xml)
endif()