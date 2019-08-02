# Allows to disable each filter
option(TTK_BUILD_CELLDATACONVERTER_FILTER "Build the CellDataConverter filter" ${TTK_ENABLE_FILTER_DEFAULT})
mark_as_advanced(TTK_BUILD_CELLDATACONVERTER_FILTER)

if(${TTK_BUILD_CELLDATACONVERTER_FILTER})
  ttk_register_pv_filter(ttkCellDataConverter CellDataConverter.xml)
endif()