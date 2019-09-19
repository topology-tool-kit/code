option(TTK_BUILD_POINTDATASELECTOR_FILTER "Build the PointDataSelector filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_POINTDATASELECTOR_FILTER)

if(${TTK_BUILD_POINTDATASELECTOR_FILTER})
  ttk_register_pv_filter(ttkPointDataSelector PointDataSelector.xml)
endif()