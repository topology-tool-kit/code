option(TTK_BUILD_PLANARGRAPHLAYOUT_FILTER "Build the PlanarGraphLayout filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_PLANARGRAPHLAYOUT_FILTER)

if(${TTK_BUILD_PLANARGRAPHLAYOUT_FILTER})
  ttk_register_pv_filter(ttkPlanarGraphLayout PlanarGraphLayout.xml)
endif()