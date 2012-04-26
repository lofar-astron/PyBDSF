find_package(Sphinx REQUIRED)

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

set(SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

add_custom_target(bdsm_doc ALL
    ${SPHINX_EXECUTABLE} -q -b html -d "${SPHINX_CACHE_DIR}" "${SPHINX_SOURCE_DIR}" "${SPHINX_HTML_DIR}")