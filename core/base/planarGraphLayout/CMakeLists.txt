ttk_add_base_library(planarGraphLayout
  SOURCES
    PlanarGraphLayout.cpp
  HEADERS
    PlanarGraphLayout.h
    MergeTreeVisualization.h
  DEPENDS
    common
    ftmTree
)

if (TTK_ENABLE_GRAPHVIZ AND Graphviz_FOUND)
  target_compile_definitions(planarGraphLayout PRIVATE TTK_ENABLE_GRAPHVIZ)
  target_link_libraries(planarGraphLayout PRIVATE Graphviz::cgraph)
  target_link_libraries(planarGraphLayout PRIVATE Graphviz::gvc)
endif()
