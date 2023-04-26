#include <time.h>
#include <threads.h>
#include <unistd.h>
#include <gtk/gtk.h>
#include "geometry.h"

int xx = 600, yy = 600, xs = 0, num_threads = 0, curX, curY, kkey = 0;
double cc1 = 0.0, cc2 = 0.0, tcmp = 0.0;
double cx = 0.998629534755, sy = 0.052335956243;
const double pi180 = M_PI / 180;
bool nomouse = false;

GtkWidget *window, *dwRect;
GtkEventController *keycon;
GtkGesture *clickcon;
GdkPixbuf *pixels; guchar *pix;
GtkWidget *sld_Sa, *sld_S0, *sld_Sb;
GdkDisplay* dd;
GdkSurface* ss;
GdkSeat *seat;
GdkDevice *mouse_device;
double xmm, ymm;

Matrix44 cam = (Matrix44){ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
Vec3d from = (Vec3d){ 100, 100, 100 }, to = (Vec3d){ 0, 0, 0 }, fvec, orig;
Vec2d cursor = (Vec2d){ 0, 0 };
Sphere wSph, wSpc, bsp0, bsp1;
Point *p0;
Line l0;
Sphere *bndVols[2]; int bcnt = 2;

struct { double d, o1, o2; Vec3d c; } typedef DCOO;

struct { int y0; } typedef thPrm;

int cmpDCOO( const void *a, const void *b )
{
    DCOO aa = *(const DCOO *)a;
    DCOO bb = *(const DCOO *)b;
    if ( aa.d > bb.d ) return 1; else return 0;
}

int outPut( void *prm_ ) {

    thPrm *prm = prm_;

    guchar *p = pix + prm->y0 * xs;

    double pX, pY; Vec3d dir;

    for ( int x = 0; x <= xx; x++ ) { p += 4;

        pX = ( 2 * ( ( (double)x + 0.5 ) / (double)xx ) - 1 );
        pY = ( 1 - 2 * ( ( (double)prm->y0 + 0.5 ) / (double)yy ) );

        Vec3d tmp = (Vec3d){ pX, pY, -1 }; multDirMatrix( &cam, &tmp , &dir ); normalize( &dir );

        DCOO ray_res[64]; int rcnt = 0; double dist1 = 0, dist2 = 0, o2 = 1.0;
        Vec3d t_col, res_col; t_col = res_col = (Vec3d){0,0,0};

        for ( int i = 0; i < bcnt; i++ ) {

         if ( intersectbSph( bndVols[i], &orig, &dir, &dist1, &dist2, &res_col, &o2 ) == 1 )
         { //ray_res[rcnt].o1 = bndVols[i]->opq;  ray_res[rcnt].d = dist1;
          // ray_res[rcnt].o2 = 0.3; ray_res[rcnt].c = res_col; rcnt++;
            for ( int o = 0; o < bndVols[i]->owns.cnt; o++ ) {

                if ( bndVols[i]->owns.otc[o] == line )
                if ( intersectLn( (Line*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &o2 )  == 1 )
                      { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist1;
                        ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }

                if ( bndVols[i]->owns.otc[o] == point )
                if ( intersectPnt( (Point*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &o2 )  == 1 )
                      { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist1;
                        ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }

         } } } qsort( &ray_res, rcnt, sizeof(*ray_res), cmpDCOO );

        intersectWSph( &wSph, &from, &dir, &pX, &pY, &res_col, &o2 );

        ray_res[rcnt].c = res_col; ray_res[rcnt].d = pX;
        ray_res[rcnt].o1 = wSph.opq; ray_res[rcnt].o2 = o2; rcnt++;

        if ( wSph.opq < 1.0 ) { intersectWSpc( &wSpc, &from, &dir, &pX, &pY, &res_col, &o2 );

            ray_res[rcnt].c = res_col; ray_res[rcnt].d = pX;
            ray_res[rcnt].o1 = wSph.opq; ray_res[rcnt].o2 = o2; rcnt++; }

        if ( rcnt > 1 ) { double prevOpq = 1.0, tempOpq = 0.0; res_col = (Vec3d){0,0,0};

              for ( int r = 0; r < rcnt; r++ ) { prevOpq = 1 - tempOpq;

                  if ( ray_res[r].o1 == 1.0 ) { tempOpq += ray_res[r].o2;

                  multscl( &ray_res[r].c, ray_res[r].o2 * prevOpq, &t_col );

                  } else { tempOpq += prevOpq * ray_res[r].o1 * ray_res[r].o2;

                  multscl( &ray_res[r].c, prevOpq * ray_res[r].o1 * ray_res[r].o2, &t_col ); }

                  res_col.x += t_col.x; res_col.y += t_col.y; res_col.z += t_col.z;

                  if ( tempOpq > 0.99 ) break; }

              res_col.x = ( res_col.x > 255 ) ? 255 : res_col.x;
              res_col.y = ( res_col.y > 255 ) ? 255 : res_col.y;
              res_col.z = ( res_col.z > 255 ) ? 255 : res_col.z;
        }

        p[0] = (guchar)res_col.x; p[1] = (guchar)res_col.y; p[2] = (guchar)res_col.z;


    }

    thrd_exit(EXIT_SUCCESS);
}

static gboolean drawFrame( GtkWidget *widget, GdkFrameClock *fclock, gpointer udata )
{
    cc1 = 100.0 * clock() / CLOCKS_PER_SEC;

    gdk_surface_get_device_position( ss, mouse_device, &xmm, &ymm, NULL );

    int dx = curX - xmm, dy = curY - ymm; curX = xmm; curY = ymm;

    dy = ( dy < 3 && dy > -3 ) ? 0 : (double)(dy) / 3.0;
    dx = ( dx < 3 && dx > -3 ) ? 0 : (double)(dx) / 3.0;

    Vec3d atpos; subvec( &to, &from, &atpos ); normalize( &atpos );

    if ( dy != 0 && !nomouse ) { atpos.y += dy * pi180; addvec( &from, &atpos, &to ); }

    if ( dx != 0 && !nomouse ) { double xsin = sin( dx * pi180 ), xcos = cos( dx * pi180 ),
                                 tx = ( atpos.x * xcos ) + ( atpos.z * xsin ),
                                 tz = ( atpos.z * xcos ) - ( atpos.x * xsin );
                                 atpos.x = tx; atpos.z = tz; addvec( &from, &atpos, &to ); }

    Vec3d forward, right, up, tmp;
    subvec( &from, &to, &forward ); normalize( &forward );
    tmp = (Vec3d){ 0, 1, 0 }; normalize( &tmp ); cross( &tmp, &forward, &right );
    cross( &forward, &right, &up );
    orig = from; subvec( &to, &from, &fvec );

    switch ( kkey )
    {
       case 111 : from.x = from.x - forward.x * 7; to.x = to.x - forward.x * 7;
                  from.y = from.y - forward.y * 7; to.y = to.y - forward.y * 7;
                  from.z = from.z - forward.z * 7; to.z = to.z - forward.z * 7; break;

       case 116 : from.x = from.x + forward.x * 7; to.x = to.x + forward.x * 7;
                  from.y = from.y + forward.y * 7; to.y = to.y + forward.y * 7;
                  from.z = from.z + forward.z * 7; to.z = to.z + forward.z * 7; break;

       case 114 : from.x = from.x + right.x * 7; to.x = to.x + right.x * 7;
                  from.y = from.y + right.y * 7; to.y = to.y + right.y * 7;
                  from.z = from.z + right.z * 7; to.z = to.z + right.z * 7; break;

       case 113 : from.x = from.x - right.x * 7; to.x = to.x - right.x * 7;
                  from.y = from.y - right.y * 7; to.y = to.y - right.y * 7;
                  from.z = from.z - right.z * 7; to.z = to.z - right.z * 7; break;
    }

    cam.m[0][0] = right.x;   cam.m[0][1] = right.y;   cam.m[0][2] = right.z;
    cam.m[1][0] = up.x;      cam.m[1][1] = up.y;      cam.m[1][2] = up.z;
    cam.m[2][0] = forward.x; cam.m[2][1] = forward.y; cam.m[2][2] = forward.z;
    cam.m[3][0] = from.x;    cam.m[3][1] = from.y;    cam.m[3][2] = from.z;

    int yl = 0; while ( yl < yy ) { thrd_t t[num_threads]; thPrm p[num_threads];

        int ty = ( yy - yl - num_threads > 0 ) ? 0 : yy - yl - num_threads;

        for ( int i = 0; i < num_threads + ty; i++ ) {

                p[i].y0 = yl + i;

                thrd_create( &t[i], (thrd_start_t) &outPut, &p[i] ); }

        for ( int i = 0; i < num_threads + ty; i++ ) thrd_join( t[i], NULL );

        yl += num_threads; }

    //Vec3d fvec; subvec( &p0.pos, &orig, &fvec ); printf( " :::: %f\n", norm( &fvec ) ); fflush(stdout);

    for ( int i = 0; i < 169; i++ ) { Vec3d ps = (Vec3d){0,0,0};
    rotate( &p0[i].pos, &ps, cx, sy, 0 ); rotate( &p0[i].pos, &ps, cx, sy, 2 ); }

    gtk_picture_set_pixbuf( GTK_PICTURE (widget), pixels );

    cc2 = 100.0 * clock() / CLOCKS_PER_SEC;

    //printf("%f - %f = %f ::: %f ::: %f \n", cc2, cc1, cc2 - cc1, gdk_frame_clock_get_fps( fclock ), fade ); fflush(stdout);

    return true;
}

static void reserved( )
{
    //
}

static gboolean key_pressed( GtkEventControllerKey* self, guint val, guint code, GdkModifierType mod )
{
    // if ( !kkey ) printf("code : %i | val : %i | mod : %i", code, val, mod ); fflush( stdout );
    kkey = code;
    return TRUE;
}

static gboolean key_released( GtkEventControllerKey* self, guint val, guint code, GdkModifierType mod )
{
    // printf(" - RELEASED\n"); fflush( stdout );
    if ( kkey == 65 ) nomouse = !nomouse;
    kkey = 0;
    return FALSE;
}

static void m_pressed( GtkGestureClick* self, gint n_press, gdouble x, gdouble y, gpointer user_data)
{
    printf("%f : %f | npress : %i \n", x, y, n_press); fflush(stdout);
}

static gboolean sldVChng( GtkRange* self, gdouble value )
{
    if ( (long int) self == (long int) sld_Sa ) scwvec = value;
    if ( (long int) self == (long int) sld_Sb ) sccvec = value / 10;
    if ( (long int) self == (long int) sld_S0 ) wSph.opq = (double)value / 100.0;

    return false;
}

static void activate( GtkApplication *app, gpointer udata )
{
    srand(4124);
    xs = xx * 4; num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    tcmp = num_threads * 0.0015;
    cc2 = 1.0 * clock() / CLOCKS_PER_SEC + 1; cc1 = cc2 - tcmp - 0.01;
    wSph.col = (Vec3d){50,60,70}; wSph.pos = (Vec3d){0,0,0}; wSph.rad2 = 13699999; wSph.opq = 1.0;
    wSpc.col = (Vec3d){250,160,170}; wSpc.pos = (Vec3d){0,0,0}; wSpc.rad2 = 26999999; wSpc.opq = 1.0;
    bsp0.col = (Vec3d){130,130,130}; bsp0.pos = to; bsp0.rad2 = 56*56; bsp0.opq = 1.0;

    bndVols[0] = &bsp0; bndVols[1] = &bsp1;
    bsp0.owns.obj = malloc( sizeof( unsigned long int[169] ) );
    bsp0.owns.otc = malloc( sizeof( enum oType[169] ) );
    bsp0.owns.cnt = 169;

    p0 = malloc( sizeof( Point[169] ) );
    for ( int i = 0; i < 169; i++ ) {

    bsp0.owns.obj[i] = (unsigned long int)&p0[i]; bsp0.owns.otc[i] = point;
    guchar x = 0, y = 0, z = 0; int res = 9999999;
    while ( res > 54 * 54 || res < 36 * 36 ) { x = rand(); y = rand(); z = rand();
    res = (int)( ( (int)(0b00111111 & x)-28 ) * ( (int)(0b00111111 & x)-28 ) +
                 ( (int)(0b00111111 & y)-28 ) * ( (int)(0b00111111 & y)-28 ) +
                 ( (int)(0b00111111 & z)-28 ) * ( (int)(0b00111111 & z)-28 ) );
    }
    p0[i].col = (Vec3d){x^(0b01111000 & z),y ^ (0b00111110 & x),(0b00111111 & z) ^ y};
    p0[i].pos = (Vec3d){(int)(0b00111111 & x)-28,(int)(0b00111111 & y)-28,(int)(0b00111111 & z)-28};
    p0[i].rad = 5.0; p0[i].opq =  1.0 / p0[i].rad; p0[i].otp = point;

    }

    l0.p0 = (Vec3d){290,20,130}; l0.p1 = (Vec3d){260,50,90}; l0.col = (Vec3d){60,200,110}; l0.opq = 1.0;
    subvec( &l0.p1, &l0.p0, &l0.n0 ); normalize( &l0.n0 ); l0.otp = line;

    bsp1.pos = (Vec3d){250,0,100}; bsp1.col = (Vec3d){130,130,130}; bsp1.rad2 = 60*60; bsp1.opq = 1.0;
    bsp1.owns.obj = malloc( sizeof( unsigned long int[128] ) );
    bsp1.owns.otc = malloc( sizeof( enum oType[128] ) );
    bsp1.owns.obj[0] = (unsigned long int)&l0; bsp1.owns.otc[0] = line;
    bsp1.owns.cnt = 1;

    GtkWidget *button, *grid, *spr, *cgrid, *tbar;

    window = gtk_application_window_new( app );
    gtk_window_set_resizable( GTK_WINDOW (window), FALSE );
    gtk_window_set_title( GTK_WINDOW (window), "[ - : - ]" );
    gtk_window_set_default_size( GTK_WINDOW (window), xx+50, yy+170 );
    keycon = gtk_event_controller_key_new();
    g_signal_connect_object( keycon, "key-pressed", G_CALLBACK (key_pressed), window, 0 );
    g_signal_connect_object( keycon, "key-released", G_CALLBACK (key_released), window, 0 );
    gtk_widget_add_controller( window, keycon );

    tbar = gtk_header_bar_new();
    gtk_header_bar_set_show_title_buttons( GTK_HEADER_BAR (tbar), FALSE );
    gtk_window_set_titlebar( GTK_WINDOW (window), tbar );

    grid = gtk_grid_new(); gtk_window_set_child( GTK_WINDOW (window), grid );
    gtk_widget_set_halign( grid, GTK_ALIGN_CENTER );
    gtk_widget_set_valign( grid, GTK_ALIGN_CENTER );

    gtk_grid_set_baseline_row( GTK_GRID (grid), 1 );
    gtk_grid_set_row_baseline_position( GTK_GRID (grid), 3, 15 );

    gtk_grid_set_row_spacing( GTK_GRID (grid), 10 );
    gtk_grid_set_column_spacing( GTK_GRID (grid), 10 );

    cgrid = gtk_grid_new(); gtk_grid_attach ( GTK_GRID (grid), cgrid, 0, 3, 1, 1 );
    gtk_widget_set_halign ( cgrid, GTK_ALIGN_END );
    gtk_widget_set_valign ( cgrid, GTK_ALIGN_END );

    button = gtk_button_new_with_label( "x" );
    g_signal_connect_swapped( button, "clicked", G_CALLBACK (gtk_window_destroy), window );
    gtk_widget_set_can_focus( button, FALSE );
    gtk_header_bar_pack_end( GTK_HEADER_BAR (tbar), button );

    button = gtk_button_new_with_label( "+" );
    g_signal_connect( button, "clicked", G_CALLBACK (reserved), NULL );
    gtk_widget_set_can_focus( button, FALSE );
    gtk_header_bar_pack_end( GTK_HEADER_BAR (tbar), button );

   // gtk_grid_attach( GTK_GRID (cgrid), button, 0, 0, 1, 11 );

    pixels = gdk_pixbuf_new( GDK_COLORSPACE_RGB, TRUE, 8, xx, yy );
    pix = gdk_pixbuf_get_pixels( pixels );
    for ( int x = 0; x < xx; x++ ) { guchar rr = rand() & 0b00001111;
        for ( int y = 0; y < yy; y++ ) { guchar *p = pix + y * xs + x * 4;
            p[0] = 60 + rr; rr += 110;
            p[1] = 30 + rr; rr += 50;
            p[2] = 10 + rr; p[3] = 252; } }

    sld_Sa = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 1, 200, 1 );
    g_object_set( sld_Sa, "width-request", 165, NULL );
    gtk_grid_attach ( GTK_GRID (cgrid), sld_Sa, 1, 0, 1, 1 );
    g_signal_connect( GTK_WIDGET (sld_Sa), "change-value", G_CALLBACK (sldVChng), NULL );
    gtk_widget_set_focusable( GTK_WIDGET (sld_Sa), false );

    sld_S0 = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 0, 100, 1 );
    g_object_set( sld_S0, "width-request", 165, NULL );
    gtk_grid_attach ( GTK_GRID (cgrid), sld_S0, 1, 1, 1, 1 );
    gtk_range_set_value ( GTK_RANGE (sld_S0), 100 );
    g_signal_connect( GTK_WIDGET (sld_S0), "change-value", G_CALLBACK (sldVChng), NULL );
    gtk_widget_set_focusable( GTK_WIDGET (sld_S0), false );

    sld_Sb = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 1, 256, 1 );
    g_object_set( sld_Sb, "width-request", 165, NULL );
    gtk_grid_attach ( GTK_GRID (cgrid), sld_Sb, 0, 0, 1, 1 );
    g_signal_connect( GTK_WIDGET (sld_Sb), "change-value", G_CALLBACK (sldVChng), NULL );
    gtk_widget_set_focusable( GTK_WIDGET (sld_Sb), false );

    dwRect = gtk_picture_new_for_pixbuf( pixels );
    //gtk_picture_set_can_shrink( GTK_PICTURE (dwRect), FALSE );
    gtk_grid_attach( GTK_GRID (grid), dwRect, 0, 7, 1, 1 );
    gtk_widget_add_tick_callback( GTK_WIDGET (dwRect), drawFrame, NULL, NULL );
    clickcon = gtk_gesture_click_new();
    g_signal_connect_object( clickcon, "pressed", G_CALLBACK (m_pressed), window, 0 );
    gtk_widget_add_controller( dwRect, GTK_EVENT_CONTROLLER (clickcon) );

    seat = gdk_display_get_default_seat( gdk_display_get_default() );
    mouse_device = gdk_seat_get_pointer( seat );
    dd = gdk_display_get_default();
    ss = gdk_surface_new_toplevel( dd );
    gtk_widget_set_visible( window, TRUE );

    gdk_surface_get_device_position( ss, mouse_device, &xmm, &ymm, NULL );
    curX = xmm; curY = ymm;

}

int main()
{
    GtkApplication *app; int status;
    app = gtk_application_new( "org.gtk.example", G_APPLICATION_FLAGS_NONE );
    g_signal_connect( app, "activate", G_CALLBACK (activate), NULL );
    status = g_application_run( G_APPLICATION (app), 0, 0 );

    free( bsp0.owns.obj ); free ( p0 );
    printf("Exit\n"); fflush(stdout);

    g_object_unref( app ); g_object_unref( pixels );
    return status;
}
