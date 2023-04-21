cmake_minimum_required(VERSION 3.5)

project(gtk-some3d LANGUAGES C)

find_package(Threads REQUIRED)

add_executable(gtk-some3d main.c)

include_directories(/usr/include/gtk-4.0 /usr/include/glib-2.0 /usr/lib64/glib-2.0/include)
include_directories(/usr/lib64/glib-2.0/include /usr/include/cairo /usr/include/pango-1.0)
include_directories(/usr/include/harfbuzz /usr/include/gdk-pixbuf-2.0 /usr/lib64/graphene-1.0/include)
include_directories(/usr/include/graphene-1.0)

target_link_libraries(gtk-some3d m Threads::Threads)
target_link_libraries(gtk-some3d gtk-4 pangocairo-1.0 pango-1.0 harfbuzz cairo-gobject cairo)
target_link_libraries(gtk-some3d graphene-1.0 glib-2.0 gobject-2.0 gio-2.0 gdk_pixbuf-2.0)
