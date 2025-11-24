cmake_minimum_required(VERSION 3.10)

project(shielding_qt LANGUAGES CXX)

# C++ standardı
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Qt için otomatik MOC/UIC/RCC
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Geant4: tüm UI ve görselleştirme modülleri
find_package(Geant4 REQUIRED ui_all vis_all)

# Sadece Qt5 Widgets yeterli
find_package(Qt5 REQUIRED COMPONENTS Widgets)

# Geant4 ayarları
include(${Geant4_USE_FILE})

# Include dizinleri
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR}
)

# Tüm .cc dosyalarını topla
file(GLOB SHIELDING_QT_SOURCES
    ${PROJECT_SOURCE_DIR}/src/*.cc
)

# Artık kullanmadığımız / sorun çıkaran DoseResultsDialog.cc'yi listeden çıkar
list(REMOVE_ITEM SHIELDING_QT_SOURCES
    ${PROJECT_SOURCE_DIR}/src/DoseResultsDialog.cc
)

# Uygulama
add_executable(shielding_qt
    ${SHIELDING_QT_SOURCES}
)

# Linkleme
target_link_libraries(shielding_qt
    PRIVATE
        ${Geant4_LIBRARIES}
        Qt5::Widgets
)

# Çıktı dizini (opsiyonel)
set_target_properties(shielding_qt PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin"
)

install(TARGETS shielding_qt DESTINATION bin)
