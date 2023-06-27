#include <math.h>
#include <stdio.h>
#include <gtk/gtk.h>

GdkPixbuf *pixels;
GtkWidget *sld_R, *sld_G, *sld_B, *sld_X, *sld_Y, *chk_X, *chk_Y, *chk_Z, *swch0, *swch1;
guchar *pix, oqins = 168;
double cx = 0.998629534755, sy = 0.052335956243, phi = 1.618033988749, dds = 0.0, df = 0.0,
       cx1 = 0.992546151641, sy1 = 0.121869343405, c90 = 0.0, s90 = 1.0, lln = 235.0,
       xshape = 1.0, scale = 900.00, scc = 0, scl = 300.0, rr = 1.0, gg = 0.1, bb = 0.5, rk, gk, bk,
       ca = 0, sa = 0, cb = 0, sb = 0, cq = 0, sq = 0, anX = 0, anY = 0, anZ = 0;
bool r_x = true, r_y = true, r_z = true, _x = true, _y = true, _z = true,
     _x_ = true, _y_ = true, _z_ = true, swst = true;

struct Vec2 { int mx,my,mz; bool is_visible; };
struct Vec3 { double mx,my,mz; unsigned char b,g,r; bool is_visible; };

struct Vec3 d3d[20], d3d3[20], d3dd[20], d3bd[20], d3cd[20], d3ed[20], fp3d[20];
struct Vec2 a2d[20], a2d3[20], a2dd[20], a2bd[20], a2cd[20], a2ed[20], fp2d[20];

struct Vec3 centerP( struct Vec3 a1, struct Vec3 a2, struct Vec3 b1, struct Vec3 b2, struct Vec3 c1)
{
    double p1 = a2.mx-a1.mx, m1 = a2.my-a1.my, l1 = a2.mz-a1.mz;

    double p2 = b2.mx-b1.mx, m2 = b2.my-b1.my, l2 = b2.mz-b1.mz;

    double a = ( ( p2 * ( b1.my - a1.my ) ) - ( m1 * ( b1.mx - a1.mx ) ) ) / ( ( p2 * m1 ) - ( p1 * m2 ) );

    if ( a != a ) a = ( ( m1 * ( b1.mz - a1.mz ) ) - ( l1 * ( b1.my - a1.my ) ) ) / ( ( m2 * l1 ) - ( m1 * l2 ) );

    double x1 = a1.mx + ( a * p1 ), y1 = a1.my + ( a * m1 ), z1 = a1.mz + ( a * l1 );

    double x2 = ( c1.mx + x1 ) / 2.0, y2 = ( c1.my + y1 ) / 2.0, z2 = ( c1.mz + z1 ) / 2.0;

    struct Vec3 res;

    res.mx = ( x2 + x1 ) / 2.0; res.my = ( y2 + y1 ) / 2.0; res.mz = ( z2 + z1 ) / 2.0;

    return res;
}

void jLine( struct Vec2 *a, struct Vec2 *b, guchar msk )
{
   int dx = abs( a->mx - b->mx ), dy = abs( a->my - b->my );

   int add = 0, slope = 9999, rmd = 6666, pos = 0, kx = 0, ky = 0;

   guchar cR = (double)(255) * ( 1 - rr ), cG = 255 - cR, cB = (double)(255) * ( 1 - bb );

   double sl = 0, ss = 0;

   if ( ( dx != 0 ) || ( dy != 0 ) ) {

   if ( b->mx >= a->mx ) kx = 1; else kx = -1; if ( b->my >= a->my ) ky = 1; else ky = -1;

   if ( dx >= dy ) { if ( dy != 0 ) { sl = (double)(dx) / (double)(dy); slope = dx / dy; ss = sl; rmd = slope; }

   for ( int i = 0; i <= dx; i++ ) { if ( i > slope )

   { add++; slope += rmd; sl += ss; if ( ( sl - slope ) > 1 ) slope++; }

   guchar *p = pix + ( ( a->my + ( ky * add ) ) * 700 + a->mx + ( kx * i ) ) * 4;

   if ( ( p[3] < 253 ) ) { p[0] = cR; p[1] = cG; p[2] = cB; p[3] = msk; } } }

   else { if ( dx != 0 ) { sl = (double)(dy) / (double)(dx); slope = dy / dx;  ss = sl; rmd = slope; }

   for ( int i = 0; i <= dy; i++ ) {

   if ( i > slope ) { add++; slope += rmd; sl += ss; if ( ( sl - slope ) > 1 ) slope++; }

   guchar *p = pix + ( ( a->my + ( ky * i ) ) * 700 + a->mx + ( kx * add ) ) * 4;

   if ( ( p[3] < 253 ) ) { p[0] = cR; p[1] = cG; p[2] = cB;  p[3] = msk; } } } }

}

void rotate3Dyt( struct Vec3 *p3d, double const cosa, double const sina, struct Vec2 *p2d )
{
    double ry = 0, rz = 0;
    ry = ( p3d->my * cosa ) - ( p3d->mz * sina ); rz = ( p3d->mz * cosa ) + ( p3d->my * sina );
    p3d->mz = rz; p3d->my = ry;
    p2d->mx = (int)( 350 + ( ( ( p3d->mx / ( -1.0 * ( rz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->my = (int)( 350 + ( ( ( p3d->my / ( -1.0 * ( rz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->mz = p3d->mz;
    if ( (p3d->mx/(-1.0*(rz-scale))*(scl-scc))<-350 || (p3d->mx/(-1.0*(rz-scale))*(scl-scc))>350 || p3d->is_visible == false ||
         (ry/(-1.0*(rz-scale))*(scl-scc))<-350 || (ry/(-1.0*(rz-scale))*(scl-scc))>350 ) p2d->is_visible=false;
    else p2d->is_visible=true;
}

void rotate3Dzt( struct Vec3 *p3d, double const cosa, double const sina, struct Vec2 *p2d )
{
    double rx = 0, ry = 0;
    rx = ( p3d->mx * cosa ) - ( p3d->my * sina ); ry = ( p3d->my * cosa ) + ( p3d->mx * sina );
    p3d->mx = rx; p3d->my = ry;
    p2d->mx = (int)( 350 + ( ( ( rx / ( -1.0 * ( p3d->mz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->my = (int)( 350 + ( ( ( ry / ( -1.0 * ( p3d->mz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->mz = p3d->mz;
    if ( (rx/(-1.0*(p3d->mz-scale))*(scl-scc))<-350 || (rx/(-1.0*(p3d->mz-scale))*(scl-scc))>350 || p3d->is_visible==false ||
         (ry/(-1.0*(p3d->mz-scale))*(scl-scc))<-350 || (ry/(-1.0*(p3d->mz-scale))*(scl-scc))>350 ) p2d->is_visible=false;
    else p2d->is_visible=true;
}

void rotate3Dxt( struct Vec3 *p3d, double const cosa, double const sina, struct Vec2 *p2d )
{
    double rx = 0, rz = 0;
    rx = ( p3d->mx * cosa ) + ( p3d->mz * sina ); rz = ( p3d->mz * cosa ) - ( p3d->mx * sina );
    p3d->mz = rz; p3d->mx = rx;
    p2d->mx = (int)( 350 + ( ( ( p3d->mx / ( -1.0 * ( rz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->my = (int)( 350 + ( ( ( p3d->my / ( -1.0 * ( rz - scale ) ) ) ) * ( scl - scc ) ) );
    p2d->mz = p3d->mz;
        if ( (rx/(-1.0*(rz-scale))*(scl-scc))<-350 || (rx/(-1.0*(rz-scale))*(scl-scc))>350 || p3d->is_visible==false ||
             (p3d->my/(-1.0*(rz-scale))*(scl-scc))<-350 || (p3d->my/(-1.0*(rz-scale))*(scl-scc))>350 ) p2d->is_visible=false;
        else p2d->is_visible=true;
}

static void chkVChng( GtkCheckButton* self )
{
    if ( (long int) self == (long int) chk_X ) {

        if ( _x_ ) { if ( _x ) { _x = !_x; } else {
            gtk_check_button_set_inconsistent( self, r_x ); r_x = !r_x;
            if ( r_x ) { _x = !_x; _x_ = !_x_; gtk_check_button_set_active( self, !_x ); }
        } } else { _x_ = !_x_; } return;

    }

    if ( (long int) self == (long int) chk_Y ) {

        if ( _y_ ) { if ( _y ) { _y = !_y; } else {
            gtk_check_button_set_inconsistent( self, r_y ); r_y = !r_y;
            if ( r_y ) { _y = !_y; _y_ = !_y_; gtk_check_button_set_active( self, !_y ); }
        } } else { _y_ = !_y_; } return;

    }

    if ( (long int) self == (long int) chk_Z )   {

        if ( _z_ ) { if ( _z ) { _z = !_z; } else {
            gtk_check_button_set_inconsistent( self, r_z ); r_z = !r_z;
            if ( r_z ) { _z = !_z; _z_ = !_z_; gtk_check_button_set_active( self, !_z ); }
        } } else { _z_ = !_z_; } return;

    }

}

void somelines( gdouble value )
{
    dds = 76 + ( value / 30 ); df = ( dds - ( 99 - dds ) * 3.3 ) / 100;

    d3bd[1].mx=((d3dd[1].mx-fp3d[0].mx)*df)+fp3d[0].mx; d3bd[1].my=((d3dd[1].my-fp3d[0].my)*df)+fp3d[0].my; d3bd[1].mz=((d3dd[1].mz-fp3d[0].mz)*df)+fp3d[0].mz;
    d3bd[0].mx=((d3dd[0].mx-fp3d[0].mx)*df)+fp3d[0].mx; d3bd[0].my=((d3dd[0].my-fp3d[0].my)*df)+fp3d[0].my; d3bd[0].mz=((d3dd[0].mz-fp3d[0].mz)*df)+fp3d[0].mz;
    d3bd[16].mx=((d3dd[16].mx-fp3d[0].mx)*df)+fp3d[0].mx; d3bd[16].my=((d3dd[16].my-fp3d[0].my)*df)+fp3d[0].my; d3bd[16].mz=((d3dd[16].mz-fp3d[0].mz)*df)+fp3d[0].mz;
    d3bd[12].mx=((d3dd[12].mx-fp3d[0].mx)*df)+fp3d[0].mx; d3bd[12].my=((d3dd[12].my-fp3d[0].my)*df)+fp3d[0].my; d3bd[12].mz=((d3dd[12].mz-fp3d[0].mz)*df)+fp3d[0].mz;
    d3bd[18].mx=((d3dd[18].mx-fp3d[0].mx)*df)+fp3d[0].mx; d3bd[18].my=((d3dd[18].my-fp3d[0].my)*df)+fp3d[0].my; d3bd[18].mz=((d3dd[18].mz-fp3d[0].mz)*df)+fp3d[0].mz;

    d3bd[4].mx=((d3dd[4].mx-fp3d[1].mx)*df)+fp3d[1].mx; d3bd[4].my=((d3dd[4].my-fp3d[1].my)*df)+fp3d[1].my; d3bd[4].mz=((d3dd[4].mz-fp3d[1].mz)*df)+fp3d[1].mz;
    d3bd[5].mx=((d3dd[5].mx-fp3d[1].mx)*df)+fp3d[1].mx; d3bd[5].my=((d3dd[5].my-fp3d[1].my)*df)+fp3d[1].my; d3bd[5].mz=((d3dd[5].mz-fp3d[1].mz)*df)+fp3d[1].mz;
    d3bd[14].mx=((d3dd[14].mx-fp3d[1].mx)*df)+fp3d[1].mx; d3bd[14].my=((d3dd[14].my-fp3d[1].my)*df)+fp3d[1].my; d3bd[14].mz=((d3dd[14].mz-fp3d[1].mz)*df)+fp3d[1].mz;
    d3bd[19].mx=((d3dd[19].mx-fp3d[1].mx)*df)+fp3d[1].mx; d3bd[19].my=((d3dd[19].my-fp3d[1].my)*df)+fp3d[1].my; d3bd[19].mz=((d3dd[19].mz-fp3d[1].mz)*df)+fp3d[1].mz;
    d3bd[17].mx=((d3dd[17].mx-fp3d[1].mx)*df)+fp3d[1].mx; d3bd[17].my=((d3dd[17].my-fp3d[1].my)*df)+fp3d[1].my; d3bd[17].mz=((d3dd[17].mz-fp3d[1].mz)*df)+fp3d[1].mz;

    d3bd[3].mx=((d3dd[3].mx-fp3d[2].mx)*df)+fp3d[2].mx; d3bd[3].my=((d3dd[3].my-fp3d[2].my)*df)+fp3d[2].my; d3bd[3].mz=((d3dd[3].mz-fp3d[2].mz)*df)+fp3d[2].mz;
    d3bd[2].mx=((d3dd[2].mx-fp3d[2].mx)*df)+fp3d[2].mx; d3bd[2].my=((d3dd[2].my-fp3d[2].my)*df)+fp3d[2].my; d3bd[2].mz=((d3dd[2].mz-fp3d[2].mz)*df)+fp3d[2].mz;
    d3bd[13].mx=((d3dd[13].mx-fp3d[2].mx)*df)+fp3d[2].mx; d3bd[13].my=((d3dd[13].my-fp3d[2].my)*df)+fp3d[2].my; d3bd[13].mz=((d3dd[13].mz-fp3d[2].mz)*df)+fp3d[2].mz;
    d3cd[16].mx=((d3dd[16].mx-fp3d[2].mx)*df)+fp3d[2].mx; d3cd[16].my=((d3dd[16].my-fp3d[2].my)*df)+fp3d[2].my; d3cd[16].mz=((d3dd[16].mz-fp3d[2].mz)*df)+fp3d[2].mz;
    d3cd[18].mx=((d3dd[18].mx-fp3d[2].mx)*df)+fp3d[2].mx; d3cd[18].my=((d3dd[18].my-fp3d[2].my)*df)+fp3d[2].my; d3cd[18].mz=((d3dd[18].mz-fp3d[2].mz)*df)+fp3d[2].mz;

    d3bd[7].mx=((d3dd[7].mx-fp3d[3].mx)*df)+fp3d[3].mx; d3bd[7].my=((d3dd[7].my-fp3d[3].my)*df)+fp3d[3].my; d3bd[7].mz=((d3dd[7].mz-fp3d[3].mz)*df)+fp3d[3].mz;
    d3bd[6].mx=((d3dd[6].mx-fp3d[3].mx)*df)+fp3d[3].mx; d3bd[6].my=((d3dd[6].my-fp3d[3].my)*df)+fp3d[3].my; d3bd[6].mz=((d3dd[6].mz-fp3d[3].mz)*df)+fp3d[3].mz;
    d3bd[15].mx=((d3dd[15].mx-fp3d[3].mx)*df)+fp3d[3].mx; d3bd[15].my=((d3dd[15].my-fp3d[3].my)*df)+fp3d[3].my; d3bd[15].mz=((d3dd[15].mz-fp3d[3].mz)*df)+fp3d[3].mz;
    d3cd[19].mx=((d3dd[19].mx-fp3d[3].mx)*df)+fp3d[3].mx; d3cd[19].my=((d3dd[19].my-fp3d[3].my)*df)+fp3d[3].my; d3cd[19].mz=((d3dd[19].mz-fp3d[3].mz)*df)+fp3d[3].mz;
    d3cd[17].mx=((d3dd[17].mx-fp3d[3].mx)*df)+fp3d[3].mx; d3cd[17].my=((d3dd[17].my-fp3d[3].my)*df)+fp3d[3].my; d3cd[17].mz=((d3dd[17].mz-fp3d[3].mz)*df)+fp3d[3].mz;


    d3cd[4].mx=((d3dd[4].mx-fp3d[4].mx)*df)+fp3d[4].mx; d3cd[4].my=((d3dd[4].my-fp3d[4].my)*df)+fp3d[4].my; d3cd[4].mz=((d3dd[4].mz-fp3d[4].mz)*df)+fp3d[4].mz;
    d3cd[0].mx=((d3dd[0].mx-fp3d[4].mx)*df)+fp3d[4].mx; d3cd[0].my=((d3dd[0].my-fp3d[4].my)*df)+fp3d[4].my; d3cd[0].mz=((d3dd[0].mz-fp3d[4].mz)*df)+fp3d[4].mz;
    d3bd[8].mx=((d3dd[8].mx-fp3d[4].mx)*df)+fp3d[4].mx; d3bd[8].my=((d3dd[8].my-fp3d[4].my)*df)+fp3d[4].my; d3bd[8].mz=((d3dd[8].mz-fp3d[4].mz)*df)+fp3d[4].mz;
    d3cd[12].mx=((d3dd[12].mx-fp3d[4].mx)*df)+fp3d[4].mx; d3cd[12].my=((d3dd[12].my-fp3d[4].my)*df)+fp3d[4].my; d3cd[12].mz=((d3dd[12].mz-fp3d[4].mz)*df)+fp3d[4].mz;
    d3cd[14].mx=((d3dd[14].mx-fp3d[4].mx)*df)+fp3d[4].mx; d3cd[14].my=((d3dd[14].my-fp3d[4].my)*df)+fp3d[4].my; d3cd[14].mz=((d3dd[14].mz-fp3d[4].mz)*df)+fp3d[4].mz;

    d3cd[6].mx=((d3dd[6].mx-fp3d[5].mx)*df)+fp3d[5].mx; d3cd[6].my=((d3dd[6].my-fp3d[5].my)*df)+fp3d[5].my; d3cd[6].mz=((d3dd[6].mz-fp3d[5].mz)*df)+fp3d[5].mz;
    d3cd[2].mx=((d3dd[2].mx-fp3d[5].mx)*df)+fp3d[5].mx; d3cd[2].my=((d3dd[2].my-fp3d[5].my)*df)+fp3d[5].my; d3cd[2].mz=((d3dd[2].mz-fp3d[5].mz)*df)+fp3d[5].mz;
    d3bd[10].mx=((d3dd[10].mx-fp3d[5].mx)*df)+fp3d[5].mx; d3bd[10].my=((d3dd[10].my-fp3d[5].my)*df)+fp3d[5].my; d3bd[10].mz=((d3dd[10].mz-fp3d[5].mz)*df)+fp3d[5].mz;
    d3cd[13].mx=((d3dd[13].mx-fp3d[5].mx)*df)+fp3d[5].mx; d3cd[13].my=((d3dd[13].my-fp3d[5].my)*df)+fp3d[5].my; d3cd[13].mz=((d3dd[13].mz-fp3d[5].mz)*df)+fp3d[5].mz;
    d3cd[15].mx=((d3dd[15].mx-fp3d[5].mx)*df)+fp3d[5].mx; d3cd[15].my=((d3dd[15].my-fp3d[5].my)*df)+fp3d[5].my; d3cd[15].mz=((d3dd[15].mz-fp3d[5].mz)*df)+fp3d[5].mz;

    d3cd[1].mx=((d3dd[1].mx-fp3d[6].mx)*df)+fp3d[6].mx; d3cd[1].my=((d3dd[1].my-fp3d[6].my)*df)+fp3d[6].my; d3cd[1].mz=((d3dd[1].mz-fp3d[6].mz)*df)+fp3d[6].mz;
    d3cd[5].mx=((d3dd[5].mx-fp3d[6].mx)*df)+fp3d[6].mx; d3cd[5].my=((d3dd[5].my-fp3d[6].my)*df)+fp3d[6].my; d3cd[5].mz=((d3dd[5].mz-fp3d[6].mz)*df)+fp3d[6].mz;
    d3bd[9].mx=((d3dd[9].mx-fp3d[6].mx)*df)+fp3d[6].mx; d3bd[9].my=((d3dd[9].my-fp3d[6].my)*df)+fp3d[6].my; d3bd[9].mz=((d3dd[9].mz-fp3d[6].mz)*df)+fp3d[6].mz;
    d3ed[12].mx=((d3dd[12].mx-fp3d[6].mx)*df)+fp3d[6].mx; d3ed[12].my=((d3dd[12].my-fp3d[6].my)*df)+fp3d[6].my; d3ed[12].mz=((d3dd[12].mz-fp3d[6].mz)*df)+fp3d[6].mz;
    d3ed[14].mx=((d3dd[14].mx-fp3d[6].mx)*df)+fp3d[6].mx; d3ed[14].my=((d3dd[14].my-fp3d[6].my)*df)+fp3d[6].my; d3ed[14].mz=((d3dd[14].mz-fp3d[6].mz)*df)+fp3d[6].mz;

    d3cd[7].mx=((d3dd[7].mx-fp3d[7].mx)*df)+fp3d[7].mx; d3cd[7].my=((d3dd[7].my-fp3d[7].my)*df)+fp3d[7].my; d3cd[7].mz=((d3dd[7].mz-fp3d[7].mz)*df)+fp3d[7].mz;
    d3cd[3].mx=((d3dd[3].mx-fp3d[7].mx)*df)+fp3d[7].mx; d3cd[3].my=((d3dd[3].my-fp3d[7].my)*df)+fp3d[7].my; d3cd[3].mz=((d3dd[3].mz-fp3d[7].mz)*df)+fp3d[7].mz;
    d3bd[11].mx=((d3dd[11].mx-fp3d[7].mx)*df)+fp3d[7].mx; d3bd[11].my=((d3dd[11].my-fp3d[7].my)*df)+fp3d[7].my; d3bd[11].mz=((d3dd[11].mz-fp3d[7].mz)*df)+fp3d[7].mz;
    d3ed[15].mx=((d3dd[15].mx-fp3d[7].mx)*df)+fp3d[7].mx; d3ed[15].my=((d3dd[15].my-fp3d[7].my)*df)+fp3d[7].my; d3ed[15].mz=((d3dd[15].mz-fp3d[7].mz)*df)+fp3d[7].mz;
    d3ed[13].mx=((d3dd[13].mx-fp3d[7].mx)*df)+fp3d[7].mx; d3ed[13].my=((d3dd[13].my-fp3d[7].my)*df)+fp3d[7].my; d3ed[13].mz=((d3dd[13].mz-fp3d[7].mz)*df)+fp3d[7].mz;

    d3cd[10].mx=((d3dd[10].mx-fp3d[8].mx)*df)+fp3d[8].mx; d3cd[10].my=((d3dd[10].my-fp3d[8].my)*df)+fp3d[8].my; d3cd[10].mz=((d3dd[10].mz-fp3d[8].mz)*df)+fp3d[8].mz;
    d3ed[16].mx=((d3dd[16].mx-fp3d[8].mx)*df)+fp3d[8].mx; d3ed[16].my=((d3dd[16].my-fp3d[8].my)*df)+fp3d[8].my; d3ed[16].mz=((d3dd[16].mz-fp3d[8].mz)*df)+fp3d[8].mz;
    d3ed[0].mx=((d3dd[0].mx-fp3d[8].mx)*df)+fp3d[8].mx; d3ed[0].my=((d3dd[0].my-fp3d[8].my)*df)+fp3d[8].my; d3ed[0].mz=((d3dd[0].mz-fp3d[8].mz)*df)+fp3d[8].mz;
    d3ed[2].mx=((d3dd[2].mx-fp3d[8].mx)*df)+fp3d[8].mx; d3ed[2].my=((d3dd[2].my-fp3d[8].my)*df)+fp3d[8].my; d3ed[2].mz=((d3dd[2].mz-fp3d[8].mz)*df)+fp3d[8].mz;
    d3cd[8].mx=((d3dd[8].mx-fp3d[8].mx)*df)+fp3d[8].mx; d3cd[8].my=((d3dd[8].my-fp3d[8].my)*df)+fp3d[8].my; d3cd[8].mz=((d3dd[8].mz-fp3d[8].mz)*df)+fp3d[8].mz;

    d3ed[4].mx=((d3dd[4].mx-fp3d[9].mx)*df)+fp3d[9].mx; d3ed[4].my=((d3dd[4].my-fp3d[9].my)*df)+fp3d[9].my; d3ed[4].mz=((d3dd[4].mz-fp3d[9].mz)*df)+fp3d[9].mz;
    d3ed[6].mx=((d3dd[6].mx-fp3d[9].mx)*df)+fp3d[9].mx; d3ed[6].my=((d3dd[6].my-fp3d[9].my)*df)+fp3d[9].my; d3ed[6].mz=((d3dd[6].mz-fp3d[9].mz)*df)+fp3d[9].mz;
    d3ed[8].mx=((d3dd[8].mx-fp3d[9].mx)*df)+fp3d[9].mx; d3ed[8].my=((d3dd[8].my-fp3d[9].my)*df)+fp3d[9].my; d3ed[8].mz=((d3dd[8].mz-fp3d[9].mz)*df)+fp3d[9].mz;
    d3ed[10].mx=((d3dd[10].mx-fp3d[9].mx)*df)+fp3d[9].mx; d3ed[10].my=((d3dd[10].my-fp3d[9].my)*df)+fp3d[9].my; d3ed[10].mz=((d3dd[10].mz-fp3d[9].mz)*df)+fp3d[9].mz;
    d3ed[17].mx=((d3dd[17].mx-fp3d[9].mx)*df)+fp3d[9].mx; d3ed[17].my=((d3dd[17].my-fp3d[9].my)*df)+fp3d[9].my; d3ed[17].mz=((d3dd[17].mz-fp3d[9].mz)*df)+fp3d[9].mz;

    d3ed[3].mx=((d3dd[3].mx-fp3d[10].mx)*df)+fp3d[10].mx; d3ed[3].my=((d3dd[3].my-fp3d[10].my)*df)+fp3d[10].my; d3ed[3].mz=((d3dd[3].mz-fp3d[10].mz)*df)+fp3d[10].mz;
    d3ed[1].mx=((d3dd[1].mx-fp3d[10].mx)*df)+fp3d[10].mx; d3ed[1].my=((d3dd[1].my-fp3d[10].my)*df)+fp3d[10].my; d3ed[1].mz=((d3dd[1].mz-fp3d[10].mz)*df)+fp3d[10].mz;
    d3cd[9].mx=((d3dd[9].mx-fp3d[10].mx)*df)+fp3d[10].mx; d3cd[9].my=((d3dd[9].my-fp3d[10].my)*df)+fp3d[10].my; d3cd[9].mz=((d3dd[9].mz-fp3d[10].mz)*df)+fp3d[10].mz;
    d3cd[11].mx=((d3dd[11].mx-fp3d[10].mx)*df)+fp3d[10].mx; d3cd[11].my=((d3dd[11].my-fp3d[10].my)*df)+fp3d[10].my; d3cd[11].mz=((d3dd[11].mz-fp3d[10].mz)*df)+fp3d[10].mz;
    d3ed[18].mx=((d3dd[18].mx-fp3d[10].mx)*df)+fp3d[10].mx; d3ed[18].my=((d3dd[18].my-fp3d[10].my)*df)+fp3d[10].my; d3ed[18].mz=((d3dd[18].mz-fp3d[10].mz)*df)+fp3d[10].mz;

    d3ed[5].mx=((d3dd[5].mx-fp3d[11].mx)*df)+fp3d[11].mx; d3ed[5].my=((d3dd[5].my-fp3d[11].my)*df)+fp3d[11].my; d3ed[5].mz=((d3dd[5].mz-fp3d[11].mz)*df)+fp3d[11].mz;
    d3ed[7].mx=((d3dd[7].mx-fp3d[11].mx)*df)+fp3d[11].mx; d3ed[7].my=((d3dd[7].my-fp3d[11].my)*df)+fp3d[11].my; d3ed[7].mz=((d3dd[7].mz-fp3d[11].mz)*df)+fp3d[11].mz;
    d3ed[9].mx=((d3dd[9].mx-fp3d[11].mx)*df)+fp3d[11].mx; d3ed[9].my=((d3dd[9].my-fp3d[11].my)*df)+fp3d[11].my; d3ed[9].mz=((d3dd[9].mz-fp3d[11].mz)*df)+fp3d[11].mz;
    d3ed[11].mx=((d3dd[11].mx-fp3d[11].mx)*df)+fp3d[11].mx; d3ed[11].my=((d3dd[11].my-fp3d[11].my)*df)+fp3d[11].my; d3ed[11].mz=((d3dd[11].mz-fp3d[11].mz)*df)+fp3d[11].mz;
    d3ed[19].mx=((d3dd[19].mx-fp3d[11].mx)*df)+fp3d[11].mx; d3ed[19].my=((d3dd[19].my-fp3d[11].my)*df)+fp3d[11].my; d3ed[19].mz=((d3dd[19].mz-fp3d[11].mz)*df)+fp3d[11].mz;
}

static gboolean sldVChng( GtkRange* self, gdouble value )
{

    if ( (long int) self == (long int) sld_R ) { rr = value / 100.0; return false; }
    if ( (long int) self == (long int) sld_G ) { gg = value / 100.0; return false; }
    if ( (long int) self == (long int) sld_B ) { bb = value / 100.0; return false; }
    if ( (long int) self == (long int) sld_Y ) { xshape = ( value ) / 100; return false; }
    if ( (long int) self == (long int) sld_X ) {

        if ( value == 690 ) gtk_widget_set_state_flags( GTK_WIDGET (sld_Y), GTK_STATE_FLAG_INSENSITIVE, 0);
        else gtk_widget_set_state_flags( GTK_WIDGET (sld_Y), GTK_STATE_FLAG_NORMAL, 1);

        somelines(value);

        return false; }

    return false;
}



static gboolean drawFrame( GtkWidget *widget, GdkFrameClock *fclock, gpointer udata )
{
    gboolean issw0 = gtk_switch_get_state ( GTK_SWITCH (swch0) ),
             issw1 = gtk_switch_get_state ( GTK_SWITCH (swch1) );

    if ( issw0 ) { if ( swst ) { if ( issw1 ) gtk_widget_set_state_flags( GTK_WIDGET (swch1), GTK_STATE_FLAG_CHECKED, 1);
                                 else gtk_widget_set_state_flags( GTK_WIDGET (swch1), GTK_STATE_FLAG_NORMAL, 1);
                                 swst = !swst; }
    if ( issw1 ) {
    for ( int x = 0; x < 700; x++ ) { for ( int y = 0; y < 700; y++ ) {
        guchar *p = pix + y * 2800 + x * 4;
        p[0] = (double)( p[0] + ( rand() & 0b00001111 ) ) * rr; //r0;
        p[1] = (double)( p[1] + ( rand() & 0b00001111 ) ) * gg; //g0;
        p[2] = (double)( p[2] + ( rand() & 0b00001111 ) ) * bb; //b0;
        p[3] = 252;
    } } } else { for ( int x = 0; x < 700; x++ ) { for ( int y = 0; y < 700; y++ ) {
                 guchar *p = pix + y * 2800 + x * 4;
                 p[0] = (double)( p[0] + ( rand() & 0b00001111 ) );
                 p[1] = (double)( p[1] + ( rand() & 0b00001111 ) );
                 p[2] = (double)( p[2] + ( rand() & 0b00001111 ) );
                 p[3] = 252;
             } } }
    }

    else {    if ( !swst ) { gtk_widget_set_state_flags( GTK_WIDGET (swch1), GTK_STATE_FLAG_INSENSITIVE, 0); swst = !swst; }

                for ( int x = 0; x < 700; x++ ) { for ( int y = 0; y < 700; y++ ) {

                guchar *p = pix + y * 2800 + x * 4;

                double xx = x - 350, yy = y - 350, dd = sqrt( (xx * xx) + (yy * yy) ), dk = tan(dd/350), dn;

                if ( dd > 256 ) dn = 0; else dn = tan ( dd / log(dk) ) * 0.0069;
                rk = ( dd < 101 ) ? 0 : ( dd < 303 ) ? ( dd - 101 ) / ( 303 - 101 ) : 1;
                gk = 1 - rk;
                bk = 1 - (rk + gk) * gk;

                double r0 = 169 * ( rr * dk) * rk; p[0] = r0;
                double g0 = 128 * ( gg * dn ) * gk; p[1] = g0;
                double b0 = 255 * bb * bk; p[2] = b0;
                p[3] = 252; } }

         }

    //printf("%i -\n", aa); fflush( stdout );

    ca = cos( anX * ( M_PI / 180 ) ); sa = sin( anX * ( M_PI / 180 ) );
    cb = cos( anY * ( M_PI / 180 ) ); sb = sin( anY * ( M_PI / 180 ) );
    cq = cos( anZ * ( M_PI / 180 ) ); sq = sin( anZ * ( M_PI / 180 ) );
    if ( r_x ) { if ( _x ) { anX += 1; if ( anX == 360 ) anX = 0; } else { anX -= 1; if ( anX == 0 ) anX = 360; } }
    if ( r_y ) { if ( _y ) { anY += 1; if ( anY == 360 ) anY = 0; } else { anY -= 1; if ( anY == 0 ) anY = 360; } }
    if ( r_z ) { if ( _z ) { anZ += 1; if ( anZ == 360 ) anZ = 0; } else { anZ -= 1; if ( anZ == 0 ) anZ = 360; } }

    for ( int c = 0; c <= 19; c++ ) {

        double x1 = d3d[c].mx,  y1 = d3d[c].my,  z1 = d3d[c].mz;
        double x2 = d3dd[c].mx, y2 = d3dd[c].my, z2 = d3dd[c].mz;
        double x3 = d3d3[c].mx, y3 = d3d3[c].my, z3 = d3d3[c].mz;
        double x4 = d3bd[c].mx, y4 = d3bd[c].my, z4 = d3bd[c].mz;
        double x5 = d3cd[c].mx, y5 = d3cd[c].my, z5 = d3cd[c].mz;
        double x6 = d3ed[c].mx, y6 = d3ed[c].my, z6 = d3ed[c].mz;
        double x7 = fp3d[c].mx, y7 = fp3d[c].my, z7 = fp3d[c].mz;

        d3d[c].mx  *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3d[c].my  *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3d[c].mz  *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3d3[c].mx *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3d3[c].my *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3d3[c].mz *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3dd[c].mx *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3dd[c].my *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;
        d3dd[c].mz *= ( dds - ( 99 - dds ) * 2 * xshape ) / 100.0;

          rotate3Dxt( &d3d[c], ca, sa, &a2d[c] );   rotate3Dxt( &d3dd[c], ca, sa, &a2dd[c] );
          rotate3Dxt( &d3d3[c], ca, sa, &a2d3[c] ); rotate3Dyt( &d3d[c], cb, sb, &a2d[c] );
          rotate3Dyt( &d3dd[c], cb, sb, &a2dd[c] ); rotate3Dyt( &d3d3[c], cb, sb, &a2d3[c] );
          rotate3Dzt( &d3d[c], cq, sq, &a2d[c] );   rotate3Dzt( &d3dd[c], cq, sq, &a2dd[c] );
          rotate3Dzt( &d3d3[c], cq, sq, &a2d3[c] ); rotate3Dxt( &d3bd[c], ca, sa, &a2bd[c] );
          rotate3Dxt( &d3cd[c], ca, sa, &a2cd[c] ); rotate3Dxt( &d3ed[c], ca, sa, &a2ed[c] );
          rotate3Dyt( &d3bd[c], cb, sb, &a2bd[c] ); rotate3Dyt( &d3cd[c], cb, sb, &a2cd[c] );
          rotate3Dyt( &d3ed[c], cb, sb, &a2ed[c] ); rotate3Dzt( &d3bd[c], cq, sq, &a2bd[c] );
          rotate3Dzt( &d3cd[c], cq, sq, &a2cd[c] ); rotate3Dzt( &d3ed[c], cq, sq, &a2ed[c] );
          rotate3Dxt( &fp3d[c], ca, sa, &fp2d[c] ); rotate3Dyt( &fp3d[c], cb, sb, &fp2d[c] );
          rotate3Dzt( &fp3d[c], cq, sq, &fp2d[c] );

        d3d[c].mx  = x1; d3d[c].my  = y1; d3d[c].mz  = z1;
        d3dd[c].mx = x2; d3dd[c].my = y2; d3dd[c].mz = z2;
        d3d3[c].mx = x3; d3d3[c].my = y3; d3d3[c].mz = z3;
        d3bd[c].mx = x4; d3bd[c].my = y4; d3bd[c].mz = z4;
        d3cd[c].mx = x5; d3cd[c].my = y5; d3cd[c].mz = z5;
        d3ed[c].mx = x6; d3ed[c].my = y6; d3ed[c].mz = z6;
        fp3d[c].mx = x7, fp3d[c].my = y7, fp3d[c].mz = z7;

    }


    jLine( &a2d[0], &a2d[8], 255 );  jLine( &a2d[0], &a2d[12], 255 ); jLine( &a2d[0], &a2d[16], 255 );
    jLine( &a2d[1], &a2d[9], 255 );  jLine( &a2d[1], &a2d[12], 255 ); jLine( &a2d[1], &a2d[18], 255 );
    jLine( &a2d[2], &a2d[10], 255);  jLine( &a2d[2], &a2d[13], 255 ); jLine( &a2d[2], &a2d[16], 255 );
    jLine( &a2d[3], &a2d[11], 255);  jLine( &a2d[3], &a2d[13], 255 ); jLine( &a2d[3], &a2d[18], 255 );
    jLine( &a2d[4], &a2d[14], 255 ); jLine( &a2d[4], &a2d[17], 255 ); jLine( &a2d[4], &a2d[8], 255  );
    jLine( &a2d[5], &a2d[19], 255 ); jLine( &a2d[5], &a2d[9], 255  ); jLine( &a2d[5], &a2d[14], 255 );
    jLine( &a2d[6], &a2d[15], 255 ); jLine( &a2d[6], &a2d[10], 255 ); jLine( &a2d[6], &a2d[17], 255 );
    jLine( &a2d[7], &a2d[11], 255 ); jLine( &a2d[7], &a2d[15], 255 ); jLine( &a2d[7], &a2d[19], 255 );
    jLine( &a2d[8], &a2d[10], 255 ); jLine( &a2d[9], &a2d[11], 255 ); jLine( &a2d[19], &a2d[17], 255 );
    jLine( &a2d[14],&a2d[12], 255 ); jLine( &a2d[15],&a2d[13], 255 ); jLine( &a2d[16], &a2d[18], 255 );


    jLine( &a2d3[0], &a2d3[8], 254 ); jLine( &a2d3[0], &a2d3[12], 254 ); jLine( &a2d3[0], &a2d3[16], 254 );
    jLine( &a2d3[1], &a2d3[9], 254 ); jLine( &a2d3[1], &a2d3[12], 254 ); jLine( &a2d3[1], &a2d3[18], 254 );
    jLine( &a2d3[2], &a2d3[10], 254); jLine( &a2d3[2], &a2d3[13], 254 ); jLine( &a2d3[2], &a2d3[16], 254 );
    jLine( &a2d3[3], &a2d3[11], 254); jLine( &a2d3[3], &a2d3[13], 254 ); jLine( &a2d3[3], &a2d3[18], 254 );
    jLine( &a2d3[4], &a2d3[14], 254); jLine( &a2d3[4], &a2d3[17], 254 ); jLine( &a2d3[4], &a2d3[8], 254 );
    jLine( &a2d3[5], &a2d3[19], 254); jLine( &a2d3[5], &a2d3[9], 254  ); jLine( &a2d3[5], &a2d3[14], 254);
    jLine( &a2d3[6], &a2d3[15], 254); jLine( &a2d3[6], &a2d3[10], 254 ); jLine( &a2d3[6], &a2d3[17], 254);
    jLine( &a2d3[7], &a2d3[11], 254); jLine( &a2d3[7], &a2d3[15], 254 ); jLine( &a2d3[7], &a2d3[19], 254);
    jLine( &a2d3[8], &a2d3[10], 254); jLine( &a2d3[9], &a2d3[11], 254 ); jLine( &a2d3[19], &a2d3[17], 254);
    jLine( &a2d3[14],&a2d3[12], 254); jLine( &a2d3[15],&a2d3[13], 254 ); jLine( &a2d3[16], &a2d3[18], 254);


    jLine( &a2dd[0], &a2dd[8], 253 ); jLine( &a2dd[0], &a2dd[12], 253 ); jLine( &a2dd[0], &a2dd[16], 253 );
    jLine( &a2dd[1], &a2dd[9], 253 ); jLine( &a2dd[1], &a2dd[12], 253 ); jLine( &a2dd[1], &a2dd[18], 253 );
    jLine( &a2dd[2], &a2dd[10], 253); jLine( &a2dd[2], &a2dd[13], 253 ); jLine( &a2dd[2], &a2dd[16], 253 );
    jLine( &a2dd[3], &a2dd[11], 253); jLine( &a2dd[3], &a2dd[13], 253 ); jLine( &a2dd[3], &a2dd[18], 253 );
    jLine( &a2dd[4], &a2dd[14], 253 ); jLine( &a2dd[4], &a2dd[17], 253 ); jLine( &a2dd[4], &a2dd[8], 253  );
    jLine( &a2dd[5], &a2dd[19], 253 ); jLine( &a2dd[5], &a2dd[9], 253  ); jLine( &a2dd[5], &a2dd[14], 253 );
    jLine( &a2dd[6], &a2dd[15], 253 ); jLine( &a2dd[6], &a2dd[10], 253 ); jLine( &a2dd[6], &a2dd[17], 253 );
    jLine( &a2dd[7], &a2dd[11], 253 ); jLine( &a2dd[7], &a2dd[15], 253 ); jLine( &a2dd[7], &a2dd[19], 253 );
    jLine( &a2dd[8], &a2dd[10], 253 ); jLine( &a2dd[9], &a2dd[11], 253 ); jLine( &a2dd[19], &a2dd[17], 253);
    jLine( &a2dd[14], &a2dd[12], 253); jLine( &a2dd[15],&a2dd[13], 253 ); jLine( &a2dd[16], &a2dd[18], 253);

    if ( dds != 99 ) {

        jLine( &a2dd[0], &a2bd[0], 252 );   jLine( &a2dd[12], &a2bd[12], 252 ); jLine( &a2dd[1], &a2bd[1], 252 );
        jLine( &a2dd[18], &a2bd[18], 252 ); jLine( &a2dd[16], &a2bd[16], 252 ); jLine( &a2dd[4], &a2bd[4], 252 );
        jLine( &a2dd[14], &a2bd[14], 252 ); jLine( &a2dd[5], &a2bd[5], 252 );   jLine( &a2dd[17], &a2bd[17], 252);
        jLine( &a2dd[19], &a2bd[19], 252 ); jLine( &a2dd[2], &a2bd[2], 252 );   jLine( &a2dd[3], &a2bd[3], 252 );
        jLine( &a2dd[13], &a2bd[13], 252 ); jLine( &a2dd[16], &a2cd[16], 252 ); jLine( &a2dd[18], &a2cd[18], 252 );
        jLine( &a2dd[6], &a2bd[6], 252 );   jLine( &a2dd[7], &a2bd[7], 252 );   jLine( &a2dd[15], &a2bd[15], 252 );
        jLine( &a2dd[17], &a2cd[17], 252 ); jLine( &a2dd[19], &a2cd[19], 252 ); jLine( &a2dd[8], &a2bd[8], 252 );
        jLine( &a2dd[4], &a2cd[4], 252 );   jLine( &a2dd[0], &a2cd[0], 252 );   jLine( &a2dd[12], &a2cd[12], 252 );
        jLine( &a2dd[14], &a2cd[14], 252 ); jLine( &a2dd[6], &a2cd[6], 252 );   jLine( &a2dd[10], &a2bd[10], 252 );
        jLine( &a2dd[2], &a2cd[2], 252 );   jLine( &a2dd[13], &a2cd[13], 252 ); jLine( &a2dd[15], &a2cd[15], 252 );
        jLine( &a2dd[5], &a2cd[5], 252 );  jLine( &a2dd[9], &a2bd[9], 252 );   jLine( &a2dd[1], &a2cd[1], 252 );
        jLine( &a2dd[12], &a2ed[12], 252 ); jLine( &a2dd[14], &a2ed[14], 252 ); jLine( &a2dd[7], &a2cd[7], 252 );
        jLine( &a2dd[11], &a2bd[11], 252 ); jLine( &a2dd[3], &a2cd[3], 252 );   jLine( &a2dd[13], &a2ed[13], 252 );
        jLine( &a2dd[15], &a2ed[15], 252 ); jLine( &a2dd[0], &a2ed[0], 252 );   jLine( &a2dd[16], &a2ed[16], 252 );
        jLine( &a2dd[2], &a2ed[2], 252 );   jLine( &a2dd[10], &a2cd[10], 252 ); jLine( &a2dd[8], &a2cd[8], 252 );
        jLine( &a2dd[4], &a2ed[4], 252 );   jLine( &a2dd[17], &a2ed[17], 252 ); jLine( &a2dd[6], &a2ed[6], 252 );
        jLine( &a2dd[10], &a2ed[10], 252 ); jLine( &a2dd[8], &a2ed[8], 252 );   jLine( &a2dd[1], &a2ed[1], 252 );
        jLine( &a2dd[3], &a2ed[3], 252 );   jLine( &a2dd[18], &a2ed[18], 252 ); jLine( &a2dd[11], &a2cd[11], 252 );
        jLine( &a2dd[9], &a2cd[9], 252 );   jLine( &a2dd[5], &a2ed[5], 252 );   jLine( &a2dd[7], &a2ed[7], 252 );
        jLine( &a2dd[9], &a2ed[9], 252 );   jLine( &a2dd[11], &a2ed[11], 252 ); jLine( &a2dd[19], &a2ed[19], 252 );
        jLine( &a2bd[0], &a2bd[12], 252 ); jLine( &a2bd[12], &a2bd[1], 252 ); jLine( &a2bd[1], &a2bd[18], 252 );
        jLine( &a2bd[18],&a2bd[16], 252 ); jLine( &a2bd[16], &a2bd[0], 252 ); jLine( &a2bd[4], &a2bd[14], 252 );
        jLine( &a2bd[14],&a2bd[5], 252 );  jLine( &a2bd[5], &a2bd[19], 252 ); jLine( &a2bd[19],&a2bd[17], 252 );
        jLine( &a2bd[17],&a2bd[4], 252 );  jLine( &a2bd[3], &a2bd[13], 252 ); jLine( &a2bd[13], &a2bd[2], 252 );
        jLine( &a2bd[2], &a2cd[16], 252 ); jLine( &a2cd[16], &a2cd[18], 252); jLine( &a2cd[18], &a2bd[3], 252 );
        jLine( &a2bd[7], &a2bd[15], 252 ); jLine( &a2bd[15], &a2bd[6], 252 ); jLine( &a2bd[6], &a2cd[17], 252 );
        jLine( &a2cd[17], &a2cd[19], 252); jLine( &a2cd[19], &a2bd[7], 252 ); jLine( &a2cd[4], &a2bd[8] , 252 );
        jLine( &a2bd[8],  &a2cd[0], 252 ); jLine( &a2cd[0], &a2cd[12], 252 ); jLine( &a2cd[12], &a2cd[14],252);
        jLine( &a2cd[14], &a2cd[4], 252 ); jLine( &a2cd[6], &a2bd[10], 252 ); jLine( &a2bd[10], &a2cd[2], 252 );
        jLine( &a2cd[2], &a2cd[13], 252 ); jLine( &a2cd[13],&a2cd[15], 252 ); jLine( &a2cd[15], &a2cd[6], 252 );
        jLine( &a2cd[5], &a2bd[9], 252  ); jLine( &a2bd[9],  &a2cd[1], 252 ); jLine( &a2cd[1], &a2ed[12], 252 );
        jLine( &a2ed[12], &a2ed[14], 252); jLine( &a2ed[14], &a2cd[5], 252 ); jLine( &a2cd[7], &a2bd[11], 252 );
        jLine( &a2bd[11], &a2cd[3], 252 ); jLine( &a2cd[3], &a2ed[13], 252 ); jLine( &a2ed[13], &a2ed[15],252);
        jLine( &a2ed[15], &a2cd[7], 252 ); jLine( &a2ed[0], &a2ed[16], 252 ); jLine( &a2ed[16], &a2ed[2], 252 );
        jLine( &a2ed[2], &a2cd[10], 252 ); jLine( &a2cd[10], &a2cd[8], 252 ); jLine( &a2cd[8], &a2ed[0] , 252 );
        jLine( &a2ed[4], &a2ed[17], 252 ); jLine( &a2ed[17], &a2ed[6], 252 ); jLine( &a2ed[6], &a2ed[10], 252 );
        jLine( &a2ed[10], &a2ed[8], 252); jLine( &a2ed[8],  &a2ed[4], 252 ); jLine( &a2ed[1], &a2ed[18], 252 );
        jLine( &a2ed[18], &a2ed[3], 252 ); jLine( &a2ed[3], &a2cd[11], 252 ); jLine( &a2cd[11], &a2cd[9], 252 );
        jLine( &a2cd[9], &a2ed[1], 252  ); jLine( &a2ed[5], &a2ed[19], 252); jLine( &a2ed[19], &a2ed[7], 252 );
        jLine( &a2ed[7], &a2ed[11], 252 ); jLine( &a2ed[11], &a2ed[9], 252 ); jLine( &a2ed[9], &a2ed[5], 252 );

    }


    gtk_picture_set_pixbuf( GTK_PICTURE (widget), pixels );

    gtk_widget_queue_draw( widget );

    //fillsub();

    return true;
}


static void activate( GtkApplication *app, gpointer udata )
{
    GtkWidget *window, *button, *fpos, *dwRect, *spr;

    window = gtk_application_window_new( app );
    gtk_window_set_title( GTK_WINDOW (window), "'" );
    gtk_window_set_default_size( GTK_WINDOW (window), 770, 835 );
    gtk_window_set_decorated( GTK_WINDOW (window), FALSE );
    gtk_window_set_resizable( GTK_WINDOW (window), FALSE );

    fpos = gtk_fixed_new(); gtk_window_set_child( GTK_WINDOW (window), fpos );

    button = gtk_button_new_with_label( "Q" );
    g_signal_connect_swapped( button, "clicked", G_CALLBACK (gtk_window_destroy), window );
    gtk_fixed_put( GTK_FIXED (fpos), button, 690, 30 );


    pixels = gdk_pixbuf_new( GDK_COLORSPACE_RGB, TRUE, 8, 700, 700 );
    pix = gdk_pixbuf_get_pixels( pixels );
    for ( int x = 0; x < 700; x++ ) { for ( int y = 0; y < 700; y++ ) {
    guchar *p = pix + y * 2800 + x * 4;
    double r0 = 169 * rr; double g0 = 128 * gg; double b0 = 255 * bb;
    p[0] = r0; p[1] = g0; p[2] = b0; p[3] = 252; } }


    dwRect = gtk_picture_new_for_pixbuf( pixels );
    gtk_picture_set_can_shrink( GTK_PICTURE (dwRect), FALSE );
    gtk_fixed_put( GTK_FIXED (fpos), dwRect, 35, 100 );
    gtk_widget_add_tick_callback( GTK_WIDGET (dwRect), drawFrame, NULL, NULL );

    spr = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    g_object_set( spr, "width-request", 770, NULL ); g_object_set( spr, "height-request", 6, NULL );
    gtk_fixed_put( GTK_FIXED (fpos), spr, 0, 0 );
    spr = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    g_object_set( spr, "width-request", 770, NULL ); g_object_set( spr, "height-request", 6, NULL );
    gtk_fixed_put( GTK_FIXED (fpos), spr, 0, 834 );
    spr = gtk_separator_new( GTK_ORIENTATION_VERTICAL );
    g_object_set( spr, "width-request", 6, NULL ); g_object_set( spr, "height-request", 834, NULL );
    gtk_fixed_put( GTK_FIXED (fpos), spr, 0, 0 );
    spr = gtk_separator_new( GTK_ORIENTATION_VERTICAL );
    g_object_set( spr, "width-request", 6, NULL ); g_object_set( spr, "height-request", 834, NULL );
    gtk_fixed_put( GTK_FIXED (fpos), spr, 770, 0 );

    spr = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    g_object_set( spr, "width-request", 706, NULL ); gtk_fixed_put( GTK_FIXED (fpos), spr, 32, 97 );
    spr = gtk_separator_new( GTK_ORIENTATION_HORIZONTAL );
    g_object_set( spr, "width-request", 706, NULL ); gtk_fixed_put( GTK_FIXED (fpos), spr, 32, 802 );
    spr = gtk_separator_new( GTK_ORIENTATION_VERTICAL );
    g_object_set( spr, "height-request", 706, NULL ); gtk_fixed_put( GTK_FIXED (fpos), spr, 32, 97 );
    spr = gtk_separator_new( GTK_ORIENTATION_VERTICAL );
    g_object_set( spr, "height-request", 706, NULL ); gtk_fixed_put( GTK_FIXED (fpos), spr, 737, 97 );

    sld_R = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 1, 100, 1 );
    g_object_set( sld_R, "width-request", 139, NULL );
    gtk_range_set_value ( GTK_RANGE (sld_R), 100 );
    gtk_fixed_put( GTK_FIXED (fpos), sld_R, 30, 10 );
    g_signal_connect( GTK_WIDGET (sld_R), "change-value", G_CALLBACK (sldVChng), NULL );

    sld_G = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 1, 42, 1 );
    g_object_set( sld_G, "width-request", 139, NULL );
    gtk_range_set_value ( GTK_RANGE (sld_G), 10 );
    gtk_fixed_put( GTK_FIXED (fpos), sld_G, 30, 35 );
    g_signal_connect( GTK_WIDGET (sld_G), "change-value", G_CALLBACK (sldVChng), NULL );

    sld_B = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 1, 100, 1 );
    g_object_set( sld_B, "width-request", 139, NULL );
    gtk_range_set_value ( GTK_RANGE (sld_B), 50 );
    gtk_fixed_put( GTK_FIXED (fpos), sld_B, 30, 60 );
    g_signal_connect( GTK_WIDGET (sld_B), "change-value", G_CALLBACK (sldVChng), NULL );

    sld_X = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 0, 690, 1 );
    g_object_set( sld_X, "width-request", 139, NULL );
    gtk_range_set_value ( GTK_RANGE (sld_X), 345 );
    somelines(345);
    gtk_fixed_put( GTK_FIXED (fpos), sld_X, 190, 22.5 );
    g_signal_connect( GTK_WIDGET (sld_X), "change-value", G_CALLBACK (sldVChng), NULL );

    sld_Y = gtk_scale_new_with_range( GTK_ORIENTATION_HORIZONTAL, 23, 100, 1 );
    g_object_set( sld_Y, "width-request", 139, NULL );
    gtk_range_set_value ( GTK_RANGE (sld_Y), 100 );
    //gtk_widget_set_state_flags( GTK_WIDGET (sld_Y), GTK_STATE_FLAG_INSENSITIVE, 0);
    gtk_fixed_put( GTK_FIXED (fpos), sld_Y, 190, 47.5 );
    g_signal_connect( GTK_WIDGET (sld_Y), "change-value", G_CALLBACK (sldVChng), NULL );

    chk_X = gtk_check_button_new_with_label( " x" );
    gtk_fixed_put( GTK_FIXED (fpos), chk_X, 365, 10 );
    g_signal_connect( GTK_WIDGET (chk_X), "toggled", G_CALLBACK (chkVChng), NULL );

    chk_Y = gtk_check_button_new_with_label( " y" );
    gtk_fixed_put( GTK_FIXED (fpos), chk_Y, 365, 35 );
    g_signal_connect( GTK_WIDGET (chk_Y), "toggled", G_CALLBACK (chkVChng), NULL );

    chk_Z = gtk_check_button_new_with_label( " z" );
    gtk_fixed_put( GTK_FIXED (fpos), chk_Z, 365, 60 );
    g_signal_connect( GTK_WIDGET (chk_Z), "toggled", G_CALLBACK (chkVChng), NULL );

    swch0 = gtk_switch_new();
    gtk_fixed_put( GTK_FIXED (fpos), swch0, 450, 20 );

    swch1 = gtk_switch_new();
    gtk_fixed_put( GTK_FIXED (fpos), swch1, 450, 55 );

    gtk_widget_set_visible( window, TRUE );

      d3d[0].mx = lln; d3d[0].my = lln; d3d[0].mz = lln; d3d[1].mx = lln; d3d[1].my =-lln; d3d[1].mz = lln;
      d3d[2].mx = lln; d3d[2].my = lln; d3d[2].mz =-lln; d3d[3].mx = lln; d3d[3].my =-lln; d3d[3].mz =-lln;
      d3d[4].mx =-lln; d3d[4].my = lln; d3d[4].mz = lln; d3d[5].mx =-lln; d3d[5].my =-lln; d3d[5].mz = lln;
      d3d[6].mx =-lln; d3d[6].my = lln; d3d[6].mz =-lln; d3d[7].mx =-lln; d3d[7].my =-lln; d3d[7].mz =-lln;

      d3d[8].mx = 0; d3d[8].my = phi * lln; d3d[8].mz = lln / phi; d3d[9].mx = 0; d3d[9].my = -phi * lln; d3d[9].mz = lln / phi;
      d3d[10].mx= 0; d3d[10].my= phi * lln; d3d[10].mz=-lln / phi; d3d[11].mx= 0; d3d[11].my= -phi * lln; d3d[11].mz=-lln / phi;

      d3d[12].mx = lln / phi; d3d[12].my = 0; d3d[12].mz = phi * lln; d3d[13].mx = lln / phi; d3d[13].my = 0; d3d[13].mz =-phi * lln;
      d3d[14].mx =-lln / phi; d3d[14].my = 0; d3d[14].mz = phi * lln; d3d[15].mx =-lln / phi; d3d[15].my = 0; d3d[15].mz =-phi * lln;

      d3d[16].mx = phi * lln; d3d[16].my = lln / phi; d3d[16].mz = 0;    d3d[17].mx =-phi * lln; d3d[17].my = lln / phi; d3d[17].mz = 0;
      d3d[18].mx = phi * lln; d3d[18].my =-lln / phi; d3d[18].mz = 0;    d3d[19].mx =-phi * lln; d3d[19].my =-lln / phi; d3d[19].mz = 0;

      for ( int i = 0; i < 20; i++ ) {

          d3d3[i].mx = d3d[i].mx * 1.3; d3d3[i].my = d3d[i].my * 1.3; d3d3[i].mz = d3d[i].mz * 1.3;
          d3dd[i].mx = d3d[i].mx * 1.6; d3dd[i].my = d3d[i].my * 1.6; d3dd[i].mz = d3d[i].mz * 1.6;

          d3bd[i].mx = d3dd[i].mx; d3bd[i].my = d3dd[i].my; d3bd[i].mz = d3dd[i].mz;
          d3cd[i].mx = d3dd[i].mx; d3cd[i].my = d3dd[i].my; d3cd[i].mz = d3dd[i].mz;
          d3ed[i].mx = d3dd[i].mx; d3ed[i].my = d3dd[i].my; d3ed[i].mz = d3dd[i].mz;

      }

      fp3d[0] = centerP( d3dd[1], d3dd[16], d3dd[0], d3dd[18], d3dd[12] );
      fp3d[1].mx =-fp3d[0].mx; fp3d[1].my = fp3d[0].my; fp3d[1].mz = fp3d[0].mz;
      fp3d[2].mx = fp3d[0].mx; fp3d[2].my = fp3d[0].my; fp3d[2].mz =-fp3d[0].mz;
      fp3d[3].mx =-fp3d[0].mx; fp3d[3].my = fp3d[0].my; fp3d[3].mz =-fp3d[0].mz;

      fp3d[4] = centerP( d3dd[4], d3dd[12], d3dd[0], d3dd[14], d3dd[8] );
      fp3d[5].mx = fp3d[4].mx; fp3d[5].my = fp3d[4].my; fp3d[5].mz =-fp3d[4].mz;
      fp3d[6].mx = fp3d[4].mx; fp3d[6].my =-fp3d[4].my; fp3d[6].mz = fp3d[4].mz;
      fp3d[7].mx = fp3d[4].mx; fp3d[7].my =-fp3d[4].my; fp3d[7].mz =-fp3d[4].mz;

      fp3d[8] = centerP( d3dd[0], d3dd[10], d3dd[2], d3dd[8], d3dd[16] );
      fp3d[9].mx =-fp3d[8].mx; fp3d[9].my = fp3d[8].my; fp3d[9].mz = fp3d[8].mz;
      fp3d[10].mx= fp3d[8].mx; fp3d[10].my=-fp3d[8].my; fp3d[10].mz= fp3d[8].mz;
      fp3d[11].mx=-fp3d[8].mx; fp3d[11].my=-fp3d[8].my; fp3d[11].mz= fp3d[8].mz;

}

int main()
{
    GtkApplication *app; int status;
    app = gtk_application_new( "org.gtk.example", G_APPLICATION_FLAGS_NONE );
    g_signal_connect( app, "activate", G_CALLBACK (activate), NULL );
    status = g_application_run( G_APPLICATION (app), 0, 0 );

    g_object_unref( app ); g_object_unref( pixels ); return status;
}
