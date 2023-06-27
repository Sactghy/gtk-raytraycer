#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double scwvec = 1.0, sccvec = 0.0;

struct { double x, y, z; } typedef Vec3d;
struct { double x, y; } typedef Vec2d;

enum oType { line, point, curve, ellipse };

struct { unsigned long int *obj; int cnt; enum oType *otc; } typedef Objects;

double norm( Vec3d *v )
{
    return v->x * v->x + v->y * v->y + v->z * v->z;
}

double length( Vec3d *v )
{
    return sqrt( v->x * v->x + v->y * v->y + v->z * v->z );
}

void normalize( Vec3d *v )
{
//  double len = norm( v ); if (len > 0) {
    double invLen = 1 / length( v );
    v->x *= invLen, v->y *= invLen, v->z *= invLen; //}
}

double dot( const Vec3d *v1, const Vec3d *v2 )
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

void cross( const Vec3d *v1, const Vec3d *v2, Vec3d *vres )
{
    vres->x = v1->y * v2->z - v1->z * v2->y;
    vres->y = v1->z * v2->x - v1->x * v2->z;
    vres->z = v1->x * v2->y - v1->y * v2->x;
}

void addvec( const Vec3d *v1, const Vec3d *v2, Vec3d *vres )
{
    vres->x = v1->x + v2->x;
    vres->y = v1->y + v2->y;
    vres->z = v1->z + v2->z;
}

void subvec( const Vec3d *v1, const Vec3d *v2, Vec3d *vres )
{
    vres->x = v1->x - v2->x;
    vres->y = v1->y - v2->y;
    vres->z = v1->z - v2->z;
}

void multscl( const Vec3d *v1, double k, Vec3d *vres )
{
    vres->x = v1->x * k;
    vres->y = v1->y * k;
    vres->z = v1->z * k;
}

void multvec( const Vec3d *v1, const Vec3d *v2, Vec3d *vres )
{
    vres->x = v1->x * v2->x;
    vres->y = v1->y * v2->y;
    vres->z = v1->z * v2->z;
}

struct { double m[4][4]; } typedef Matrix44;

Matrix44 rrx; double an0;

Matrix44 mxmult( Matrix44 *m, Matrix44 *rhs )
{
    Matrix44 mult;
    for ( int i = 0; i < 4; ++i ) { for ( int j = 0; j < 4; ++j ) {
    mult.m[i][j] = m->m[i][0]*rhs->m[0][j]+m->m[i][1]*rhs->m[1][j]+m->m[i][2]*rhs->m[2][j]+m->m[i][3]*rhs->m[3][j]; } }
    return mult;
}

void transpose( Matrix44 *m )
{
    *m = (Matrix44){ m->m[0][0], m->m[1][0], m->m[2][0], m->m[3][0],
                     m->m[0][1], m->m[1][1], m->m[2][1], m->m[3][1],
                     m->m[0][2], m->m[1][2], m->m[2][2], m->m[3][2],
                     m->m[0][3], m->m[1][3], m->m[2][3], m->m[3][3] };
}

void multVecMatrix( Matrix44 *m, Vec3d *src, Vec3d *dst )
{
        dst->x = src->x * m->m[0][0] + src->y * m->m[1][0] + src->z * m->m[2][0] + m->m[3][0];
        dst->y = src->x * m->m[0][1] + src->y * m->m[1][1] + src->z * m->m[2][1] + m->m[3][1];
        dst->z = src->x * m->m[0][2] + src->y * m->m[1][2] + src->z * m->m[2][2] + m->m[3][2];
        double w = src->x * m->m[0][3] + src->y * m->m[1][3] + src->z * m->m[2][3] + m->m[3][3];
        if ( w != 1 || w != 0 ) { dst->x /= w; dst->y /= w; dst->z /= w;  }
}

void multDirMatrix( Matrix44 *m, Vec3d *src, Vec3d *dst )
{
        double a, b, c;
        a = src->x * m->m[0][0] + src->y * m->m[1][0] + src->z * m->m[2][0];
        b = src->x * m->m[0][1] + src->y * m->m[1][1] + src->z * m->m[2][1];
        c = src->x * m->m[0][2] + src->y * m->m[1][2] + src->z * m->m[2][2];
        dst->x = a; dst->y = b; dst->z = c;
}

void inverse( Matrix44 *dst, Matrix44 *m )
{
    int i, j, k; Matrix44 s; Matrix44 t = *m; s = (Matrix44){0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    for ( i = 0; i < 3 ; i++ ) { int pivot = i; double pivotsize = t.m[i][i];
        if ( pivotsize < 0 ) pivotsize = -pivotsize;
        for ( j = i + 1; j < 4; j++ ) { double tmp = t.m[j][i];
            if ( tmp < 0 ) tmp = -tmp;
            if ( tmp > pivotsize ) { pivot = j; pivotsize = tmp; } }
        if ( pivotsize == 0 ) { printf("Can't compute the matrix!\n"); *dst = s; }
        if ( pivot != i ) { for ( j = 0; j < 4; j++ ) {
                double tmp = t.m[i][j];
                       t.m[i][j] = t.m[pivot][j];
                       t.m[pivot][j] = tmp;
                       tmp = s.m[i][j];
                       s.m[i][j] = s.m[pivot][j];
                       s.m[pivot][j] = tmp; } }
        for ( j = i + 1; j < 4; j++ ) { double f = t.m[j][i] / t.m[i][i];
            for ( k = 0; k < 4; k++ ) { t.m[j][k] -= f * t.m[i][k]; s.m[j][k] -= f * s.m[i][k]; } } }

    for ( i = 3; i >= 0; --i ) { double f;
        if ( ( f = t.m[i][i] ) == 0 ) { printf("Error ::: singular matrix!\n"); }
        for ( j = 0; j < 4; j++ ) { t.m[i][j] /= f; s.m[i][j] /= f; }
        for ( j = 0; j < i; j++ ) { f = t.m[j][i];
            for ( k = 0; k < 4; k++ ) { t.m[j][k] -= f * t.m[i][k]; s.m[j][k] -= f * s.m[i][k]; } }
    } *dst = s;
}

//const Matrix44<T>& invert() { *this = inverse(); return *this; }

struct { double opq, rad, rad2; int vis, isvis, type; Vec3d pos, col, o_pos; Objects owns; } typedef Sphere;

int intersectWSph( Sphere *wSph, Vec3d *orig, Vec3d *dir, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{
       Vec3d lpos;
       subvec( &wSph->pos, orig, &lpos );
       double tca = dot( &lpos, dir ), thc = sqrt( wSph->rad2 - ( dot( &lpos, &lpos ) - tca * tca )  );

       *dist1 = tca - thc; *dist2 = tca + thc;
       *dist1 = ( *dist1 < *dist2 ) ? *dist2 : *dist1;

       Vec3d Phit;
       Phit.x = ( orig->x + dir->x * *dist1 );// - wSph->pos.x;
       Phit.y = ( orig->y + dir->y * *dist1 );// - wSph->pos.x;
       Phit.z = ( orig->z + dir->z * *dist1 );// - wSph->pos.x;
       normalize( &Phit ); Phit.x *= scwvec; Phit.y *= scwvec; Phit.z *= scwvec;

       res_col->x = wSph->col.x + ( Phit.x / M_PI ) * 100;
       res_col->y = wSph->col.y - ( Phit.y / M_PI ) * 100;
       res_col->z = wSph->col.z + ( Phit.z / M_PI ) * 100;

       int rr, gg, bb;
       res_col->x = fabs(res_col->x); res_col->y = fabs(res_col->y); res_col->z = fabs(res_col->z);
       rr = res_col->x / 256; gg = res_col->y / 256; bb = res_col->z / 256;
       res_col->x = ( rr % 2 == 0 ) ? (double)((unsigned char)res_col->x) * wSph->opq : (255 - (double)((unsigned char)res_col->x)) * wSph->opq;
       res_col->y = ( gg % 2 == 0 ) ? (double)((unsigned char)res_col->y) * wSph->opq : (255 - (double)((unsigned char)res_col->y)) * wSph->opq;
       res_col->z = ( bb % 2 == 0 ) ? (double)((unsigned char)res_col->z) * wSph->opq : (255 - (double)((unsigned char)res_col->z)) * wSph->opq;

       *o2 = wSph->opq; return 1;
}

int intersectWSpc( Sphere *wSpc, Vec3d *orig, Vec3d *dir, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{
       Vec3d lpos;
       subvec( &wSpc->pos, orig, &lpos );
       double tca = dot( &lpos, dir ), thc = sqrt( wSpc->rad2 - ( dot( &lpos, &lpos ) - tca * tca )  );

       *dist1 = tca - thc; *dist2 = tca + thc;
       *dist1 = ( *dist1 < *dist2 ) ? *dist2 : *dist1;

       Vec3d Phit;
       Phit.x = ( orig->x + dir->x * *dist1 ) - wSpc->pos.x;
       Phit.y = ( orig->y + dir->y * *dist1 ) - wSpc->pos.x;
       Phit.z = ( orig->z + dir->z * *dist1 ) - wSpc->pos.x;
       normalize( &Phit ); Phit.x *= sccvec; Phit.y *= sccvec; Phit.z *= sccvec;

       unsigned int ais2 = 2;// (unsigned int) ( ( cos ( Phit.x ) + sin ( Phit.z * Phit.y ) ) ) ;


        double xx,// = atan ( sin(Phit.x * M_E) / ( cos ( Phit.y ) + 0.0001 ) ),
               yy = cos ( Phit.y * atan2 ( sin (Phit.z) , cos(Phit.x) ) ) * 35 ,
               zz = asin ( sin(Phit.z) / M_PI ) * 65;

        if ( ais2 % 2 == 0 ) { xx = yy; zz = xx; yy = zz; } else  { xx = zz ; zz = yy; yy = xx; }

       res_col->x = 0;//( wSpc->col.x + xx ) * 0.14;
       res_col->y = ( wSpc->col.y - yy ) * 0.69;
       res_col->z = ( wSpc->col.z + zz ) * 0.65;

       *o2 = wSpc->opq; return 1;
}

int intersectbSph( Sphere *bsp, Vec3d *o, Vec3d *d, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{
    Vec3d lpos; lpos.x = bsp->pos.x - o->x; lpos.y = bsp->pos.y - o->y; lpos.z = bsp->pos.z - o->z;
    double tca = dot( &lpos, d );
    if ( tca < 0 ) return 0;
    double d2 = dot( &lpos, &lpos ) - tca * tca;
    if ( d2 > bsp->rad2 ) return 0;
    double thc = sqrt( bsp->rad2 - d2 );

    *dist1 = tca - thc; *dist2 = tca + thc;
    *dist1 = ( *dist1 < *dist2 ) ? *dist2 : *dist1;

    *res_col = bsp->col; *o2 = bsp->opq; return 1;

 }

struct { double opq, rad; int vis, isvis, type; Vec3d pos, col, o_pos; enum oType otp; } typedef Point;

int intersectPnt( Point *pnt, Vec3d *o, Vec3d *d, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{
    Vec3d fvec, pvec, tvec, pps; pps = pnt->pos; double asc = 1;
    for ( int n = 0; n < 2; n++ ){
    subvec( &pps, o, &fvec );  pps.x *= 0.48;pps.y *= 0.48;pps.z *= 0.48;
    pvec = fvec; normalize( &pvec );
    cross( &pvec, d, &tvec );

    *dist1 = norm( &fvec );
    double pvar = norm( &tvec ) * *dist1 * asc; asc+=1.5;

    if ( pvar < pnt->rad ) { int xx = atan(pvar)*6;  if ( xx % 2 == 0 ) {

            *o2 = ( pnt->opq * ( pnt->rad - pvar ) * pvar);
            res_col->x = pnt->col.x * pvar / xx * 2;
            res_col->y = pnt->col.y * pvar / xx * 4;
            res_col->z = pnt->col.z * pvar / 2;

            return 1; }

    } } return 0;

}

struct { double opq; int vis, isvis, type; Vec3d n0, p0, p1, col; enum oType otp; } typedef Line;

int intersectLn( Line *ln, Vec3d *o, Vec3d *d, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{
        Vec3d u, r;  cross( &ln->n0, d, &u ); subvec( o, &ln->p1, &r );

        double is = dot ( &u, &r ), iss = is * is * 2;

        *dist1 = norm( &u ) * 1.36;

        if ( iss < *dist1 ) {

            Vec3d ppt0, ppt1, pp0, pp1; subvec( &ln->p0, o, &ppt0 );  cross( &ppt0, &u, &pp0 );
                                        subvec( &ln->p1, o, &ppt1 );  cross( &ppt1, &u, &pp1 );

//normalize( &pp0 ); normalize( &pp1 );

            double esc = 1.0, dd0 = dot( d, &pp0 ), dd1 = dot( d, &pp1 );

            if ( dd0 > 0.0 && dd1 < 0.0 && (int)(dd1+dd0)%2 == 0) {

            esc = ( dd0 > 1.0 ) ? ( dd1 < -1.0 ) ? 1.0 : -1.0 * dd1 : dd0;

            *o2 = ( *dist1 - fabs( is ) ) * esc * 2;
            if ( *o2 > 1 ) *o2 = 1; else if ( *o2 < 0 ) *o2 = 0;

            double ln1 = norm( &ppt0 ); normalize( &ppt0 );
            double dif1 = dot( &ppt0, d );
            double ln2 = norm( &ppt1 ); normalize( &ppt1 );
            double dif2 = dot( &ppt1, d );

            *dist1 = ( ( ln1 + ln2 ) * ( fabs( dif1 ) + fabs ( dif2 ) ) );

            res_col->y = ln->col.y + ( *o2 *69);
            res_col->x = ln->col.x * *o2 ;
            res_col->z = ln->col.z * ( 1 - *o2); *o2 *= 0.36;

            return 1; } else { Vec3d cln; cross( d, &ln->n0, &cln ); double scl = norm( &cln );

                if ( scl < 1 && scl > 0.01 ) {


                    Vec3d ppt0, ppt1, pp0, pp1; subvec( &ln->p0, o, &ppt0 );  cross( &ppt0, &u, &pp0 );
                                                subvec( &ln->p1, o, &ppt1 );  cross( &ppt1, &u, &pp1 );

 //normalize( &pp0 ); normalize( &pp1 );

                    double esc = 1.0, dd0 = dot( d, &pp0 ), dd1 = dot( d, &pp1 );

                    if ( dd0 > -0.5 && dd1 < 0.5 ) {
                    scl *= 1; *o2 = 1 - scl;

                    res_col->y = ln->col.y + ( *o2 *69);
                    res_col->x = ln->col.x * *o2 ;
                    res_col->z = ln->col.z * ( 1 - *o2); *o2 *= 0.36;
                    return 1; }

                }

                } return 0; }

    return 0;

}

int solveQ( double a, double b, double c, double *x0, double *x1 )
{
    double d = b * b - 4.0 * a * c;
    if ( d < 0 ) { return 0; }
    else if ( d == 0 ) { *x0 = *x1 = -0.5 * b / a; }
    else { double q = ( b > 0 ) ? -0.5 * ( b + sqrt( d ) ) : -0.5 * ( b - sqrt( d ) );
           *x0 = q / a ; *x1 = c / q ; }
    return 1;
}

struct { double opq, size, ang; int tp, vis, isvis; Vec3d p0, f0, col; enum oType otp; Matrix44 rmx; } typedef Ellipse;

int intersectEl( Ellipse *el, Vec3d *o, Vec3d *d, double *dist1, double *dist2, Vec3d *res_col1, Vec3d *res_col2, double *o1,  double *o2 )
{
    Vec3d f0; subvec( o, &el->p0, &f0 );


    Vec3d dv, fv;
    multDirMatrix( &el->rmx, d, &dv );
    multDirMatrix( &el->rmx, &f0, &fv );
    Vec3d dv0 = (Vec3d){ dv.x / 3, dv.y , dv.z / 2 };
    Vec3d fv0 = (Vec3d){ fv.x / 3, fv.y , fv.z / 2 };


    //multDirMatrix( &el->rmx, &fv0, &fv );
    //multDirMatrix( &el->rmx, &dv0, &dv );
     //printf("%f - %f\n", fv.x, fv0.x);

    double a, b, c;
           a = dot( &dv0, &dv0);
           b = 2.0 * dot( &dv0, &fv0 );
           c = dot( &fv0, &fv0 ) - el->size;

    if ( solveQ( a, b, c, dist1, dist2 ) ) { *o1 = *o2 = el->opq;

        if ( *dist1 <= 0 ) { *dist1 = INFINITY; } if ( *dist2 <= 0 ) { *dist2 = INFINITY; }

        if ( *dist1 > *dist2 ) { double tmp = *dist1; *dist1 = *dist2; *dist2 = tmp;
                                 if ( *dist1 == INFINITY ) return 0; }

        Vec3d Fhit, Phit, tt, tFhit, tPhit;
              multscl( d, *dist1, &tFhit ); addvec( o, &tFhit, &tt );
              subvec( &tt, &el->p0, &tFhit ); normalize( &tFhit );
              multscl( d, *dist2, &tPhit ); addvec( o, &tPhit, &tt );
              subvec( &tt, &el->p0, &tPhit ); normalize( &tPhit );
              //Fhit = ( ( o + d * dist1 ) - pos ) . normalize();
              //Phit = ( ( o + d * dist2 ) - pos ) . normalize();
              multDirMatrix( &el->rmx, &tFhit, &Fhit );
              multDirMatrix( &el->rmx, &tPhit, &Phit );
              multDirMatrix( &rrx, &Fhit, &tFhit ); Fhit = tFhit;
              multDirMatrix( &rrx, &Phit, &tPhit ); Phit = tPhit;

        double dotd1; dotd1 = ( dist1 < 0 ) ? dot( &Fhit, d ) : 1;
        double dotd2; dotd2 = ( dist2 < 0 ) ? dot( &Phit, d ) : 1;

        uint ais1, ais2;

        switch ( el->tp ) { case 0 : break;

        case 1 : ais1 = (uint) ( ( ( cos ( Fhit.z ) + atan ( cos ( Fhit.y ) - sin ( ( Fhit.x ) + 1 ) ) ) ) * 6 ); //scw );
                 ais2 = (uint) ( ( ( cos ( Phit.z ) + atan ( cos ( Phit.y ) - sin ( ( Phit.x ) + 1 ) ) ) ) * 6 ); // scw );
                 break;

        case 2 : ais1 = (uint) ( ( ( sin ( atan ( Fhit.z  ) * atan ( Fhit.y  ) ) + sin ( cos ( Fhit.y*1.619 ) - tan ( Fhit.x ) ) ) ) * 10 );
                 ais2 = (uint) ( ( ( sin ( atan ( Phit.z  ) * atan ( Phit.y  ) ) + sin ( cos ( Phit.y*1.619 ) - tan ( Phit.x ) ) ) ) * 10 );
                 break;

        case 3 : ais1 = (uint) ( sin( cos ( ( Fhit.x*1.619 - Fhit.z ) + tan ( -1*Fhit.z * Fhit.y ) ) )  * 6 );
                 ais2 = (uint) ( sin( cos ( ( Phit.x*1.619 - Fhit.z ) + tan ( -1*Phit.z * Phit.y ) ) )  * 6 );
                 break;

        }

        if ( ( el->tp != 4 && ais1 % 2 == 0 ) || dotd1 < 0 ) { *dist1 = INFINITY; } else {

        *res_col1 = (Vec3d){ el->col.x - ( Fhit.x / M_PI ) * 190.0 ,
                    el->col.y + ( Fhit.y / M_PI ) * 190.0 ,
                    el->col.z + ( Fhit.z / M_PI ) * 190.0 }; }

        if ( ( el->tp != 4 && ais2 % 2 == 0 ) || dotd2 < 0 ) { *dist2 = INFINITY; } else {

        *res_col2 = (Vec3d){ el->col.x - ( Phit.x / M_PI ) * 190.0 ,
                     el->col.y + ( Phit.y / M_PI ) * 190.0 ,
                     el->col.z + ( Phit.z / M_PI ) * 190.0 }; }

        if ( el->tp == 4 ) // <- != { *dist1 *= 0.0001; *dist2 *= 0.0002; }
        //else
        { double sca =fabs(*dist2-*dist1); if ( sca < 63 ) { *o1 *= sca/63; *o2 *= sca/63; *o1 *= *o1; *o2 *= *o2;} *dist1 *= 0.01; *dist2 *= 0.02;
            *dist1 *= 1000.0; *dist2 *= 1000.0;}
        *dist1 *= *dist1; *dist2 *= *dist2;

   return 1; } return 0;

}

struct { double opq; int vis, isvis; Vec3d p0, col; enum oType otp; } typedef Curve;

int intersectSp( Curve *sp, Vec3d *o, Vec3d *d, double *dist1, double *dist2, Vec3d *res_col, double *o2 )
{

    double a_=10,b_=2,c_=6,d_=14,e_=20,f_=4;
    double p = a_*(d->x*d->x) + 2*b_*d->x*d->y + c_*d->y;
    double q = (a_*o->x+b_*o->y+d_)*d->x + (b_*o->x+c_*o->y+e_)*d->y;
    double r = a_*(o->x*o->x) + 2*b_*o->x*o->y + c_*(o->y*o->y) + 2*d_*o->x + 2*e_*o->y + f_;

    if ( solveQ( p, q, r, dist1, dist2 ) ) {

        *o2 = sp->opq;
        *res_col = sp->col;

    return 1; }

    return 0;
}

void rotate( Vec3d *p0, Vec3d *ps, double cx, double sy, int tp )
{
    double tx, tz;

         if ( tp == 1 ) { p0->x = p0->x - ps->x; p0->y = p0->y - ps->y; p0->z = p0->z - ps->z;
                          tx=(p0->x*cx)-(p0->z*sy); tz=(p0->z*cx)+(p0->x*sy);
                          p0->z = tz; p0->x = tx; p0->x = p0->x + ps->x;
                          p0->y = p0->y + ps->y; p0->z = p0->z + ps->z; }

         if ( tp == 0 ) { p0->x = p0->x - ps->x; p0->y = p0->y - ps->y; p0->z = p0->z - ps->z;
                          tx=(p0->x*cx)-(p0->y*sy); tz=(p0->y*cx)+(p0->x*sy);
                          p0->y = tz; p0->x = tx; p0->x = p0->x + ps->x;
                          p0->y = p0->y + ps->y; p0->z = p0->z + ps->z; }

         if ( tp == 2 ) { p0->x = p0->x - ps->x; p0->y = p0->y - ps->y; p0->z = p0->z - ps->z;
                          tx=(p0->y*cx)-(p0->z*sy); tz=(p0->z*cx)+(p0->y*sy);
                          p0->z = tz; p0->y = tx; p0->x = p0->x + ps->x;
                          p0->y = p0->y + ps->y; p0->z = p0->z + ps->z; }
}

#endif // GEOMETRY_H
