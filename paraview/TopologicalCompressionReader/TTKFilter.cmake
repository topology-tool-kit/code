# Allows to disable each filter
option(TTK_BUILD_TOPOLOGICALCOMPRESSIONREADER_FILTER "Build the TopologicalCompressionReader filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TOPOLOGICALCOMPRESSIONREADER_FILTER)

if(${TTK_BUILD_TOPOLOGICALCOMPRESSIONREADER_FILTER})
  ttk_register_pv_filter(ttkTopologicalCompressionReader TopologicalCompressionReader.xml)
endif()