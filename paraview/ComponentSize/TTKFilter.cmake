option(TTK_BUILD_COMPONENTSIZE_FILTER "Build the ComponentSize filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_COMPONENTSIZE_FILTER)

if(${TTK_BUILD_COMPONENTSIZE_FILTER})
  ttk_register_pv_filter(ttkComponentSize ComponentSize.xml)
endif()