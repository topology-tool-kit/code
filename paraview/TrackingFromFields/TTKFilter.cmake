option(TTK_BUILD_TRACKINGFROMFIELDS_FILTER "Build the TrackingFromFields filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_TRACKINGFROMFIELDS_FILTER)

if(${TTK_BUILD_TRACKINGFROMFIELDS_FILTER})
  ttk_register_pv_filter(ttkTrackingFromFields TrackingFromFields.xml)
endif()