#include <time.h>
#include <threads.h>
#include <unistd.h>
#include <gtk/gtk.h>
#include "geometry.h"

int xx = 600, yy = 600, xs = 0, num_threads = 0, curX, curY, kkey = 0;
double cc1 = 0.0, cc2 = 0.0, tcmp = 0.0, aa0 = 0;
double cos3 = 0.998629534755, sin3 = 0.052335956243, cos1 = 0.999847695156, sin1 = 0.017452406437;
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
Vec3d from = (Vec3d){ 100, 200, 400 }, to = (Vec3d){ 100,-75,200 }, fvec, orig;
Vec2d cursor = (Vec2d){ 0, 0 };
Sphere wSph, wSpc, bsp0, bsp1, bsp2;
Point *p0;
Line *l0;
Curve sp;
Ellipse el0, el1, el2, el3;
Sphere *bndVols[3]; int bcnt = 3;

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

        DCOO ray_res[64]; int rcnt = 0; double dist1, dist2, o1 = 1.0, o2 = 1.0;
        Vec3d t_col, res_col, res_col1; t_col = res_col = res_col1 = (Vec3d){0,0,0};

        for ( int i = 0; i < bcnt; i++ ) { dist1 = 0; dist2 = 0;

         if ( bndVols[i]->isvis && intersectbSph( bndVols[i], &orig, &dir, &dist1, &dist2, &res_col, &o2 ) == 1 )
         { //ray_res[rcnt].o1 = bndVols[i]->opq;  ray_res[rcnt].d = dist1;
         // ray_res[rcnt].o2 = 0.3; ray_res[rcnt].c = res_col; rcnt++;
            for ( int o = 0; o < bndVols[i]->owns.cnt; o++ ) {

                if ( bndVols[i]->owns.otc[o] == line )
                if ( intersectLn( (Line*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &o2 ) == 1 )
                      { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist1;
                        ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }

                if ( bndVols[i]->owns.otc[o] == point )
                if ( intersectPnt( (Point*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &o2 ) == 1 )
                      { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist1;
                        ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }

                if ( bndVols[i]->owns.otc[o] == curve )
                if ( intersectSp( (Curve*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &o2 ) == 1 )
                      { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = 1;
                        ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }

                if ( bndVols[i]->owns.otc[o] == ellipse )
                if ( intersectEl( (Ellipse*)bndVols[i]->owns.obj[o], &orig, &dir, &dist1, &dist2, &res_col, &res_col1, &o1, &o2 ) == 1 )
                      { if( dist1 != INFINITY ) { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist1;
                          ray_res[rcnt].o2 = o2; ray_res[rcnt].c = res_col; rcnt++; }
                    if( dist2 != INFINITY ) { ray_res[rcnt].o1 = bndVols[i]->opq; ray_res[rcnt].d = dist2;
                          ray_res[rcnt].o2 = o1; ray_res[rcnt].c = res_col1; rcnt++; }}

         } } } qsort( &ray_res, rcnt, sizeof(*ray_res), cmpDCOO );

        intersectWSph( &wSph, &from, &dir, &pX, &pY, &res_col, &o2 );

        ray_res[rcnt].c = res_col; ray_res[rcnt].d = pX;
        ray_res[rcnt].o1 = 1; ray_res[rcnt].o2 = 1; rcnt++;

        if ( wSph.opq < 1.0 ) { intersectWSpc( &wSpc, &from, &dir, &pX, &pY, &res_col, &o2 );

            ray_res[rcnt-1].c.x = ray_res[rcnt-1].c.x * wSph.opq + res_col.x * (1- wSph.opq);
            ray_res[rcnt-1].c.y = ray_res[rcnt-1].c.y * wSph.opq +res_col.y * (1 -wSph.opq);
            ray_res[rcnt-1].c.z = ray_res[rcnt-1].c.z * wSph.opq + res_col.z * (1- wSph.opq);
            res_col.x = ray_res[rcnt-1].c.x;
            res_col.y = ray_res[rcnt-1].c.y;
            res_col.z = ray_res[rcnt-1].c.z;


            //ray_res[rcnt].d = pX;
            //ray_res[rcnt].o1 = wSph.opq; ray_res[rcnt].o2 = o2;
        }

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

    //Vec3d aa22,to22; multscl(&atpos,100,&aa22);addvec( &from, &aa22, &to22 );
    Vec3d forward, right, up, tmp;//bsp2.pos = to22; el0.p0 = to22; el1.p0 = to22; el2.p0 = to22; el3.p0 = to22;
    subvec( &from, &to, &forward ); normalize( &forward );
    tmp = (Vec3d){ 0, 1, 0 }; normalize( &tmp ); cross( &tmp, &forward, &right );
    cross( &forward, &right, &up );
    orig = from; subvec( &to, &from, &fvec );

    switch ( kkey ) {

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

    Vec3d ps = (Vec3d){0,0,0};

    for ( int i = 0; i < 169; i++ ) {
    rotate( &p0[i].pos, &ps, cos3, sin3, 0 ); rotate( &p0[i].pos, &ps, cos3, sin3, 2 ); }

    el0.ang+=1; aa0-=1; an0 = aa0*M_PI / 180;
    double qDegRadX, qDegRadY, qDegRadZ;
    qDegRadX = 0;//( el0.ang * M_PI ) / 180;
    qDegRadY = 0;//( el0.ang * M_PI ) / 180;
    qDegRadZ = ( el0.ang * M_PI ) / 180;

    el0.rmx.m[0][0] = cos(qDegRadX) * cos(qDegRadY); el0.rmx.m[0][1] = ( cos(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) - ( sin(qDegRadX) * cos(qDegRadZ) ); el0.rmx.m[0][2] = ( cos(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) + ( sin(qDegRadX) * sin(qDegRadZ) ); el0.rmx.m[0][3] = 0;
    el0.rmx.m[1][0] = sin(qDegRadX) * cos(qDegRadY); el0.rmx.m[1][1] = ( sin(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) + ( cos(qDegRadX) * cos(qDegRadZ) ); el0.rmx.m[1][2] = ( sin(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) - ( cos(qDegRadX) * sin(qDegRadZ) ); el0.rmx.m[1][3] = 0;
    el0.rmx.m[2][0] = -1.0 * sin(qDegRadY); el0.rmx.m[2][1] = cos(qDegRadY) * sin(qDegRadZ); el0.rmx.m[2][2] = cos(qDegRadY) * cos(qDegRadZ); el0.rmx.m[2][3] = 0;
    el0.rmx.m[3][0] = 0; el0.rmx.m[3][1] = 0; el0.rmx.m[3][2] = 0; el0.rmx.m[3][3] = 1;

    el1.rmx.m[0][0] = cos(qDegRadX) * cos(qDegRadY); el1.rmx.m[0][1] = ( cos(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) - ( sin(qDegRadX) * cos(qDegRadZ) ); el1.rmx.m[0][2] = ( cos(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) + ( sin(qDegRadX) * sin(qDegRadZ) ); el1.rmx.m[0][3] = 0;
    el1.rmx.m[1][0] = sin(qDegRadX) * cos(qDegRadY); el1.rmx.m[1][1] = ( sin(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) + ( cos(qDegRadX) * cos(qDegRadZ) ); el1.rmx.m[1][2] = ( sin(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) - ( cos(qDegRadX) * sin(qDegRadZ) ); el1.rmx.m[1][3] = 0;
    el1.rmx.m[2][0] = -1.0 * sin(qDegRadY); el1.rmx.m[2][1] = cos(qDegRadY) * sin(qDegRadZ); el1.rmx.m[2][2] = cos(qDegRadY) * cos(qDegRadZ); el1.rmx.m[2][3] = 0;
    el1.rmx.m[3][0] = 0; el1.rmx.m[3][1] = 0; el1.rmx.m[3][2] = 0; el1.rmx.m[3][3] = 1;

    el2.rmx.m[0][0] = cos(qDegRadX) * cos(qDegRadY); el2.rmx.m[0][1] = ( cos(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) - ( sin(qDegRadX) * cos(qDegRadZ) ); el2.rmx.m[0][2] = ( cos(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) + ( sin(qDegRadX) * sin(qDegRadZ) ); el2.rmx.m[0][3] = 0;
    el2.rmx.m[1][0] = sin(qDegRadX) * cos(qDegRadY); el2.rmx.m[1][1] = ( sin(qDegRadX) * sin(qDegRadY) * sin(qDegRadZ) ) + ( cos(qDegRadX) * cos(qDegRadZ) ); el2.rmx.m[1][2] = ( sin(qDegRadX) * sin(qDegRadY) * cos(qDegRadZ) ) - ( cos(qDegRadX) * sin(qDegRadZ) ); el2.rmx.m[1][3] = 0;
    el2.rmx.m[2][0] = -1.0 * sin(qDegRadY); el2.rmx.m[2][1] = cos(qDegRadY) * sin(qDegRadZ); el2.rmx.m[2][2] = cos(qDegRadY) * cos(qDegRadZ); el2.rmx.m[2][3] = 0;
    el2.rmx.m[3][0] = 0; el2.rmx.m[3][1] = 0; el2.rmx.m[3][2] = 0; el2.rmx.m[3][3] = 1;

    el3.rmx.m[0][0] = el2.rmx.m[0][0]; el3.rmx.m[0][1] = el2.rmx.m[0][1]; el3.rmx.m[0][2] = el2.rmx.m[0][2]; el3.rmx.m[0][3] = 0;
    el3.rmx.m[1][0] = el2.rmx.m[1][0]; el3.rmx.m[1][1] = el2.rmx.m[1][1]; el3.rmx.m[1][2] = el2.rmx.m[1][2]; el3.rmx.m[1][3] = 0;
    el3.rmx.m[2][0] = el2.rmx.m[2][0]; el3.rmx.m[2][1] = el2.rmx.m[2][1]; el3.rmx.m[2][2] = el2.rmx.m[2][2]; el3.rmx.m[2][3] = 0;
    el3.rmx.m[3][0] = 0; el3.rmx.m[3][1] = 0; el3.rmx.m[3][2] = 0; el3.rmx.m[3][3] = 1;

    rrx.m[0][0] = cos(an0) * cos(an0); rrx.m[0][1] = ( cos(an0) * -sin(an0) * sin(an0) ) - ( sin(an0) * cos(an0) ); rrx.m[0][2] = ( cos(an0) * sin(an0) * cos(an0) ) + ( sin(an0) * sin(an0) ); rrx.m[0][3] = 0;
    rrx.m[1][0] = sin(an0) * -cos(an0); rrx.m[1][1] = ( sin(an0) * sin(an0) * -sin(an0*2) ) + ( cos(an0) * cos(an0) ); rrx.m[1][2] = ( sin(an0) * sin(an0) * cos(an0) ) - ( cos(an0) * sin(an0) ); rrx.m[1][3] = 0;
    rrx.m[2][0] = -1.0 * sin(an0); rrx.m[2][1] = cos(an0) * sin(an0); rrx.m[2][2] = -cos(an0) * cos(an0); rrx.m[2][3] = 0;
    rrx.m[3][0] = 0; rrx.m[3][1] = 0; rrx.m[3][2] = 0; rrx.m[3][3] = 1;

    for ( int i = 0; i < 90; i++ ) {

    rotate( &l0[i].p0, &bsp1.pos, cos1, sin1, 1 ); rotate( &l0[i].p0, &bsp1.pos, cos1, sin1, 0 ); //rotate( &l0[i].p0, &bsp1.pos, cos1, sin1, 2 );
    rotate( &l0[i].p1, &bsp1.pos, cos1, sin1, 1 ); rotate( &l0[i].p1, &bsp1.pos, cos1, sin1, 0 ); //rotate( &l0[i].p1, &bsp1.pos, cos1, sin1, 2 );
    rotate( &l0[i].n0, &ps, cos1, sin1, 1 ); rotate( &l0[i].n0, &ps, cos1, sin1, 0 ); //rotate( &l0[i].n0, &ps, cos1, sin1, 2 );
    }


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


    bndVols[0] = &bsp0; bndVols[1] = &bsp1; bndVols[2] = &bsp2;
    bsp0.isvis = 1; bsp1.isvis = 1; bsp2.isvis = 1;

    bsp0.col = (Vec3d){130,130,130}; bsp0.pos = (Vec3d){0,0,0}; bsp0.rad2 = 56*56; bsp0.opq = 1.0;
    bsp0.owns.obj = malloc( sizeof( unsigned long int[169] ) );
    bsp0.owns.otc = malloc( sizeof( enum oType[169] ) );
    bsp0.owns.cnt = 169; bsp0.isvis = 1;

    p0 = malloc( sizeof( Point[169] ) );
    for ( int i = 0; i < 169; i++ ) {

    bsp0.owns.obj[i] = (unsigned long int)&p0[i]; bsp0.owns.otc[i] = point;
    guchar x = 0, y = 0, z = 0; int res = 9999999;
    while ( res > 54 * 54 || res < 36 * 36 ) { x = rand(); y = rand(); z = rand();
    res = (int)( ( (int)( 0b00111111 & x ) - 28 ) * ( (int)( 0b00111111 & x ) - 28 ) +
                 ( (int)( 0b00111111 & y ) - 28 ) * ( (int)( 0b00111111 & y ) - 28 ) +
                 ( (int)( 0b00111111 & z ) - 28 ) * ( (int)( 0b00111111 & z ) - 28 ) ); }

    p0[i].col = (Vec3d){ x ^ ( 0b01111000 & z ), y ^ ( 0b00111110 & x ), ( 0b00111111 & z ) ^ y };
    p0[i].pos = (Vec3d){(int)( 0b00111111 & x ) - 28,(int)( 0b00111111 & y ) - 28,(int)( 0b00111111 & z ) - 28};
    p0[i].rad = 5.0; p0[i].opq =  1.0 / p0[i].rad; p0[i].otp = point;

    }

    bsp1.pos = (Vec3d){450,0,-200}; bsp1.col = (Vec3d){130,130,130}; bsp1.rad2 = 51*51; bsp1.opq = 1.0;
    bsp1.owns.obj = malloc( sizeof( unsigned long int[90] ) );
    bsp1.owns.otc = malloc( sizeof( enum oType[90] ) );
    bsp1.owns.cnt = 90;

    double de = 0.0, ispm = 1.0, im, ck = 1.0;
    l0 = malloc( sizeof( Line[90] ) );
    for ( int i = 0; i < 90; i++ ) {

    bsp1.owns.obj[i] = (unsigned long int)&l0[i]; bsp1.owns.otc[i] = line;
    l0[i].p0 = (Vec3d){ 0, 25 - de, 0 }; l0[i].p1 = (Vec3d){ 0, 35 + de, 0 }; l0[i].n0 = (Vec3d){ 0, 1, 0 };
    Vec3d def = (Vec3d){ 0, 0, 0 }, cen = (Vec3d){ 0, 30, 0 };

    im = 3.0 * ispm; de += im;
    if ( de >= 15.0 || de <= 0.0 ) { ispm *= -1; }

    if ( i < 45 ) ck += 2; else ck -= 2;
    double cosF = cos( ( i * 4 * M_PI ) / 180 ), sinF = sin( ( i * 4 * M_PI ) / 180 );

    rotate( &l0[i].p0, &cen, cosF , sinF , 0 );
    rotate( &l0[i].p1, &cen, cosF , sinF , 0 );
    rotate( &l0[i].n0, &def, cosF , sinF , 0 );

    rotate( &l0[i].p0, &def, cosF , sinF , 2 );
    rotate( &l0[i].p1, &def, cosF , sinF , 2 );
    rotate( &l0[i].n0, &def, cosF , sinF , 2 );

    l0[i].p0.x += bsp1.pos.x; l0[i].p0.y += bsp1.pos.y; l0[i].p0.z += bsp1.pos.z;
    l0[i].p1.x += bsp1.pos.x; l0[i].p1.y += bsp1.pos.y; l0[i].p1.z += bsp1.pos.z;
    l0[i].col = (Vec3d){ 196 + ck /2, 160 - ck, 150 }; l0[i].vis = 1; l0[i].opq =  1.0; p0[i].otp = line;

    }

    bsp2.col = (Vec3d){130,130,130}; bsp2.pos = (Vec3d){300,-75,200}; bsp2.rad2 = 78*78; bsp2.opq = 1.0;
    bsp2.owns.obj = malloc( sizeof( unsigned long int[4] ) );
    bsp2.owns.otc = malloc( sizeof( enum oType[4] ) );
    bsp2.owns.cnt = 4;

    bsp2.owns.obj[0] = (unsigned long int)&el0;
    bsp2.owns.obj[1] = (unsigned long int)&el1;
    bsp2.owns.obj[2] = (unsigned long int)&el2;
    bsp2.owns.obj[3] = (unsigned long int)&el3;//sp;
    bsp2.owns.otc[0] = ellipse;
    bsp2.owns.otc[1] = ellipse;
    bsp2.owns.otc[2] = ellipse;
    bsp2.owns.otc[3] = ellipse;//curve;

    el0.p0 = bsp2.pos; el0.size =300; el0.col = (Vec3d){148,176,0}; el0.tp = 1; el0.opq = 1.0; el0.ang = 0;
    el1.p0 = bsp2.pos; el1.size =420; el1.col = (Vec3d){165,99,150}; el1.tp = 2; el1.opq = 0.63; el1.ang = 0;
    el2.p0 = bsp2.pos; el2.size =600; el2.col = (Vec3d){95,56,160}; el2.tp = 3; el2.opq = 0.43; el2.ang = 0;
    el3.p0 = bsp2.pos; el3.size =690; el3.col = (Vec3d){129,255,129}; el3.tp = 4; el3.opq = 1.0; el3.ang = 0;

    sp.col = (Vec3d){256,0,0}; sp.opq = 1.0; sp.p0 = bsp2.pos;

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

    free( bsp0.owns.obj ); free( bsp0.owns.otc ); free ( p0 );
    free( bsp1.owns.obj ); free( bsp1.owns.otc ); free ( l0 );

    g_object_unref( app ); g_object_unref( pixels );

    printf("Exit\n"); fflush(stdout); return status;
}
