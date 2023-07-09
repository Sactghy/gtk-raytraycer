cmake_minimum_required(VERSION 3.5)

project(gr3d LANGUAGES C)

find_package(Threads REQUIRED)

add_executable(gr3d main.c)

include_directories(/usr/include/gtk-4.0 /usr/include/glib-2.0 /usr/lib64/glib-2.0/include)
include_directories(/usr/lib64/glib-2.0/include /usr/include/cairo /usr/include/pango-1.0)
include_directories(/usr/include/harfbuzz /usr/include/gdk-pixbuf-2.0 /usr/lib64/graphene-1.0/include)
include_directories(/usr/include/graphene-1.0)

target_link_libraries(gr3d m Threads::Threads)
target_link_libraries(gr3d gtk-4 pangocairo-1.0 pango-1.0 harfbuzz cairo-gobject cairo)
target_link_libraries(gr3d graphene-1.0 glib-2.0 gobject-2.0 gio-2.0 gdk_pixbuf-2.0)
