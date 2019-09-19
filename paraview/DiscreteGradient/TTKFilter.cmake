option(TTK_BUILD_DISCRETEGRADIENT_FILTER "Build the DiscreteGradient filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_DISCRETEGRADIENT_FILTER)

if(${TTK_BUILD_DISCRETEGRADIENT_FILTER})
  ttk_register_pv_filter(ttkDiscreteGradient DiscreteGradient.xml)
endif()