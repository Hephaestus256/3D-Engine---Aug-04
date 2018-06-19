#ifndef INCLUDE_ZMAP
#define INCLUDE_ZMAP
#define ZMAP_WHOLE 65536 * 16384

#include "math3d.h"
#include "gfx2d.h"
#include "gfx3d.h"
#include "texclass.h"
//#include "/programs/include/fixed.h"

inline bool slice_x_axis
  (lin_relat* xz, point_3d p1, point_3d p2, point_3d p3);
inline bool slice_y_axis
  (lin_relat* yz, point_3d p1, point_3d p2, point_3d p3);
inline void zmap_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  double inv_z, double m_xy, double a);
inline void zmap_yz_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  double a, double b, int c, long diz);
inline void zmap_z_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  long iz);
inline void zmap_group
  (tri_group group, z_buffer_type* zmap, view_type view, int poly);
inline double calc_scan_m
  (view_type view, point_3d p1, point_3d p2, point_3d p3);
inline void zmap_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap,
  int poly, double d_inv_z, double m_xy, double a);
inline void zmap_yz_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap, int poly, double a,
  double b, int c, long diz);
inline void zmap_z_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap, int poly, long iz);
inline bool slice_xz (y_line_3d ref1, y_line_3d ref2, lin_relat* xz);
inline void clear_zmap (view_type view, z_buffer_type* zmap);


inline bool slice_x_axis
  (lin_relat* xz, point_3d p1, point_3d p2, point_3d p3)
{
//  if (approx_equal (p1.z, p2.z))
//    if (approx_equal (p2.z, p3.z))
//      return false;

  if (approx_equal (p1.y, p2.y))
    if (approx_equal (p2.y, p3.y)) // 1=2=3
      return false;
    else { // 1=2!3
      y_line_3d ref1 = calc_y_line_3d (p1, p3);
      y_line_3d ref2 = calc_y_line_3d (p2, p3);
      return slice_xz (ref1, ref2, xz);
    }
  else if (approx_equal (p2.y, p3.y)) { // 2=3!1
    y_line_3d ref1 = calc_y_line_3d (p3, p1);
    y_line_3d ref2 = calc_y_line_3d (p2, p1);
    return slice_xz (ref1, ref2, xz);
  }
  else if (approx_equal (p1.y, p3.y)) {
    y_line_3d ref1 = calc_y_line_3d (p1, p2);
    y_line_3d ref2 = calc_y_line_3d (p3, p2);
    return slice_xz (ref1, ref2, xz);
  }
  else {
    y_line_3d ref1 = calc_y_line_3d (p1, p2);
    y_line_3d ref2 = calc_y_line_3d (p2, p3);
    return slice_xz (ref1, ref2, xz);
  }
}


inline bool slice_xz (y_line_3d ref1, y_line_3d ref2, lin_relat* xz)
{
  point_3d pa, pb;
  double y;
  
  if (approx_equal (ref1.bz, ref2.bz))
    if (approx_equal (ref1.mz, ref2.mz))
      return false;
    else
      y = 1;
  else
    y = 0;
    
  pa.x = y * ref1.mx + ref1.bx;
  pa.z = y * ref1.mz + ref1.bz;
  pb.x = y * ref2.mx + ref2.bx;
  pb.z = y * ref2.mz + ref2.bz;
  *xz = calc_dxdz_line (pa, pb);

  return true;
}


inline double calc_scan_m
  (view_type view, point_3d p1, point_3d p2, point_3d p3)
{
  if (approx_equal (p1.z, p2.z)) {
    line_3d l3d = calc_line_3d (p1, p3);
    point_3d pa;
    pa.x = l3d.mx * p2.z + l3d.bx;
    pa.y = l3d.my * p2.z + l3d.by;
    return (view.zoom.x * (p2.x - pa.x)) / (view.zoom.y * (p2.y - pa.y));
  }
  else if (approx_equal (p1.z, p3.z)) { // 1=3!2
    line_3d l3d = calc_line_3d (p1, p2);
    point_3d pa;
    pa.x = l3d.mx * p3.z + l3d.bx;
    pa.y = l3d.my * p3.z + l3d.by;
    return (view.zoom.x * (p3.x - pa.x)) / (view.zoom.y * (p3.y - pa.y));
  }
  else if (approx_equal (p2.z, p3.z)) { // 2=3!1
    line_3d l3d = calc_line_3d (p1, p2);
    point_3d pa;
    pa.x = l3d.mx * p3.z + l3d.bx;
    pa.y = l3d.my * p3.z + l3d.by;
    return (view.zoom.x * (p3.x - pa.x)) / (view.zoom.y * (p3.y - pa.y));
  }
  else {
    line_3d l3d = calc_line_3d (p1, p3);
    point_3d pa;
    pa.x = l3d.mx * p2.z + l3d.bx;
    pa.y = l3d.my * p2.z + l3d.by;
    return (view.zoom.x * (p2.x - pa.x)) / (view.zoom.y * (p2.y - pa.y));
  }
}


inline bool slice_y_axis
  (lin_relat* yz, point_3d p1, point_3d p2, point_3d p3)
{
  if (approx_equal (p1.z, p2.z))
    if (approx_equal (p2.z, p3.z))
      return false;
    else
      *yz = calc_dydz_line (p2, p3);
  else
    *yz = calc_dydz_line (p1, p2);

  return true;
}


inline void zmap_group
  (tri_group group, z_buffer_type* zmap, view_type view, int poly)
{
  point_3d p1 = group.ref.orig->rel;
  point_3d p2 = group.ref.u->rel;
  point_3d p3 = group.ref.v->rel;
  unsigned long c = ZMAP_WHOLE;
  
  view_clip_class g;
  temp_tri result[1000];

  lin_relat xz, yz;
  if (slice_x_axis (&xz, p1, p2, p3)) {
    double m_xy = calc_scan_m (view, p1, p2, p3);
    double a = m_xy * view.center.y - view.center.x - view.zoom.x * xz.m;
    double d_inv_z = 1 / (view.zoom.x * xz.b);

    for (tri_type* t = group.first; t != NULL; t = t->next) {
      int ct = g.view_clip (t->t3d, &view, result);
      for (int i = 0; i < ct; i++)
        zmap_tri (result[i], view, zmap, poly, d_inv_z, m_xy, a);
    }
  }
  else if (slice_y_axis (&yz, p1, p2, p3)) {
    double a = -view.center.y - view.zoom.y * yz.m;
    double b = 1 / (view.zoom.y * yz.b);
    unsigned long diz = (unsigned long)(c * b);

    for (tri_type* t = group.first; t != NULL; t = t->next) {
      int ct = g.view_clip (t->t3d, &view, result);
      for (int i = 0; i < ct; i++)
        zmap_yz_tri (result[i], view, zmap, poly, a, b, c, diz);
    }
  }
  else {
    for (tri_type* t = group.first; t != NULL; t = t->next) {
      int ct = g.view_clip (t->t3d, &view, result);
      for (int i = 0; i < ct; i++)
        zmap_z_tri (result[i], view, zmap, poly,
        (unsigned long)(double(c) / p1.z));
    }
  }
}


inline void zmap_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap,
  int poly, double d_inv_z, double m_xy, double a)
{
  point_2d s1 = map_to_scrn (view, tri.p1);
  point_2d s2 = map_to_scrn (view, tri.p2);
  point_2d s3 = map_to_scrn (view, tri.p3);
  
  if (approx_lesser (s1.y, s2.y))
    if (approx_lesser (s3.y, s1.y)) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, d_inv_z, m_xy, a);
      zmap_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, d_inv_z, m_xy, a);
//        printf ("3, 1, 2");
    }
    else if (approx_greater (s3.y, s1.y))
      if (approx_greater (s3.y, s2.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, d_inv_z, m_xy, a);
        zmap_sect (s2.y, s3.y, line_31, line_23, zmap, poly, view, d_inv_z, m_xy, a);
//          printf ("1, 2, 3");
      }
      else if (approx_lesser (s3.y, s2.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_sect (s1.y, s3.y, line_31, line_12, zmap, poly, view, d_inv_z, m_xy, a);
        zmap_sect (s3.y, s2.y, line_23, line_12, zmap, poly, view, d_inv_z, m_xy, a);
//          printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, d_inv_z, m_xy, a);
//          printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      zmap_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, d_inv_z, m_xy, a);
//    printf ("1-3, 2");
    }
  else if (approx_greater (s1.y, s2.y))
    if (approx_lesser (s3.y, s2.y)) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_sect (s3.y, s2.y, line_23, line_31, zmap, poly, view, d_inv_z, m_xy, a);
      zmap_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, d_inv_z, m_xy, a);
//    printf ("3, 2, 1");
    }
    else if (approx_greater (s3.y, s2.y))
      if (approx_greater (s3.y, s1.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_sect (s2.y, s1.y, line_12, line_23, zmap, poly, view, d_inv_z, m_xy, a);
        zmap_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, d_inv_z, m_xy, a);
//      printf ("2, 1, 3");
      }
      else if (approx_lesser (s3.y, s1.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, d_inv_z, m_xy, a);
        zmap_sect (s3.y, s1.y, line_12, line_31, zmap, poly, view, d_inv_z, m_xy, a);
//      printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        zmap_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, d_inv_z, m_xy, a);
//      printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, d_inv_z, m_xy, a);
//    printf ("2-3, 1");
    }
  else
    if (approx_lesser (s3.y, s1.y)) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, d_inv_z, m_xy, a);
//    printf ("3, 1-2");
    }
    else if (approx_greater (s3.y, s1.y)) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, d_inv_z, m_xy, a);
//    printf ("1-2, 3");
    }
}


inline void zmap_yz_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap, int poly, double a,
  double b, int c, long diz)
{
  point_2d s1 = map_to_scrn (view, tri.p1);
  point_2d s2 = map_to_scrn (view, tri.p2);
  point_2d s3 = map_to_scrn (view, tri.p3);
  
//inline void zmap_yz_sect (double y1, double y2,
//  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
//  double a, double b, int c, long diz)

  if (approx_lesser (s1.y, s2.y))
    if (approx_lesser (s3.y, s1.y)) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_yz_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, a, b, c, diz);
      zmap_yz_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, a, b, c, diz);
//        printf ("3, 1, 2");
    }
    else if (approx_greater (s3.y, s1.y))
      if (approx_greater (s3.y, s2.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_yz_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, a, b, c, diz);
        zmap_yz_sect (s2.y, s3.y, line_31, line_23, zmap, poly, view, a, b, c, diz);
//          printf ("1, 2, 3");
      }
      else if (approx_lesser (s3.y, s2.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_yz_sect (s1.y, s3.y, line_31, line_12, zmap, poly, view, a, b, c, diz);
        zmap_yz_sect (s3.y, s2.y, line_23, line_12, zmap, poly, view, a, b, c, diz);
//          printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_yz_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, a, b, c, diz);
//          printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      zmap_yz_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, a, b, c, diz);
//    printf ("1-3, 2");
    }
  else if (approx_greater (s1.y, s2.y))
    if (approx_lesser (s3.y, s2.y)) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_yz_sect (s3.y, s2.y, line_23, line_31, zmap, poly, view, a, b, c, diz);
      zmap_yz_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, a, b, c, diz);
//    printf ("3, 2, 1");
    }
    else if (approx_greater (s3.y, s2.y))
      if (approx_greater (s3.y, s1.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_yz_sect (s2.y, s1.y, line_12, line_23, zmap, poly, view, a, b, c, diz);
        zmap_yz_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, a, b, c, diz);
//      printf ("2, 1, 3");
      }
      else if (approx_lesser (s3.y, s1.y)) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_yz_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, a, b, c, diz);
        zmap_yz_sect (s3.y, s1.y, line_12, line_31, zmap, poly, view, a, b, c, diz);
//      printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        zmap_yz_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, a, b, c, diz);
//      printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_yz_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, a, b, c, diz);
//    printf ("2-3, 1");
    }
  else
    if (approx_lesser (s3.y, s1.y)) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_yz_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, a, b, c, diz);
//    printf ("3, 1-2");
    }
    else if (approx_greater (s3.y, s1.y)) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_yz_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, a, b, c, diz);
//    printf ("1-2, 3");
    }
}


inline void zmap_z_tri
  (temp_tri tri, view_type view, z_buffer_type* zmap, int poly, long iz)
{
  point_2d s1 = map_to_scrn (view, tri.p1);
  point_2d s2 = map_to_scrn (view, tri.p2);
  point_2d s3 = map_to_scrn (view, tri.p3);
  
//inline void zmap_yz_sect (double y1, double y2,
//  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
//  double a, double b, int c, long diz)

  if (s1.y < s2.y)
    if (s3.y < s1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_z_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, iz);
      zmap_z_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, iz);
//        printf ("3, 1, 2");
    }
    else if (s3.y > s1.y)
      if (s3.y > s2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_z_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, iz);
        zmap_z_sect (s2.y, s3.y, line_31, line_23, zmap, poly, view, iz);
//          printf ("1, 2, 3");
      }
      else if (s3.y < s2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_z_sect (s1.y, s3.y, line_31, line_12, zmap, poly, view, iz);
        zmap_z_sect (s3.y, s2.y, line_23, line_12, zmap, poly, view, iz);
//          printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_z_sect (s1.y, s2.y, line_31, line_12, zmap, poly, view, iz);
//          printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      zmap_z_sect (s1.y, s2.y, line_23, line_12, zmap, poly, view, iz);
//    printf ("1-3, 2");
    }
  else if (s1.y > s2.y)
    if (s3.y < s2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_z_sect (s3.y, s2.y, line_23, line_31, zmap, poly, view, iz);
      zmap_z_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, iz);
//    printf ("3, 2, 1");
    }
    else if (s3.y > s2.y)
      if (s3.y > s1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_z_sect (s2.y, s1.y, line_12, line_23, zmap, poly, view, iz);
        zmap_z_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, iz);
//      printf ("2, 1, 3");
      }
      else if (s3.y < s1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
        zmap_z_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, iz);
        zmap_z_sect (s3.y, s1.y, line_12, line_31, zmap, poly, view, iz);
//      printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
        lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
        zmap_z_sect (s2.y, s3.y, line_12, line_23, zmap, poly, view, iz);
//      printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (s1, s2);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_z_sect (s2.y, s1.y, line_12, line_31, zmap, poly, view, iz);
//    printf ("2-3, 1");
    }
  else
    if (s3.y < s1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_z_sect (s3.y, s1.y, line_23, line_31, zmap, poly, view, iz);
//    printf ("3, 1-2");
    }
    else if (s3.y > s1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (s2, s3);
      lin_relat line_31 = calc_2d_dxdy_line (s3, s1);
      zmap_z_sect (s1.y, s3.y, line_31, line_23, zmap, poly, view, iz);
//    printf ("1-2, 3");
    }
}


inline void zmap_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  double d_inv_z, double m_xy, double a)
{
  static int c = ZMAP_WHOLE;
  double top = ifloor(y1 + .5);
  double bot = ifloor(y2 - .5);
  double l = lef.m * (top + .5) + lef.b;
  double r = rig.m * (top + .5) + rig.b;
  int y_offset = _screen.x_res * (int)top;
  unsigned long diz = (unsigned long)(c * d_inv_z);

  for (int y = int(top); y <= (int)bot; y++) {
    int x1 = (int)ifloor(l + .5);
    int x2 = (int)ifloor(r - .5);
    double inv_z = (a + (x1 + .5) - m_xy * y) * d_inv_z;
    unsigned long iz = (unsigned long)(c * inv_z);
    unsigned long start = (unsigned long)(&zmap[x1 + y_offset]);
    int len = x2 - x1;

    if (x1 <= x2)
    asm (
      "push %%esi \n"
      "movl %0, %%ebx \n"
      "movl %1, %%eax \n"
      "movl %2, %%edx \n"
      "movl %3, %%ecx \n"
      "movl %4, %%edi \n"
      "push %%ebp \n"
      "movl __pxl_ct, %%ebp \n"
      ".align 4 \n"
      "rep_zmap: \n"
      "  movl 0(%%ebx), %%esi \n"
      "  cmpl %%esi, %%eax \n"
      "  jle skip_zmap \n"
      "    movl %%eax, 0(%%ebx) \n"
      "    movl %%edi, 4(%%ebx) \n"
      "    incl %%ebp \n"
      "  skip_zmap: \n"
      "  addl %%edx, %%eax \n"
      "  addl $8, %%ebx \n"
      "decl %%ecx \n"
      "jge rep_zmap \n"
      "movl %%ebp, __pxl_ct \n"
      "pop %%ebp \n"
      "pop %%esi \n"
      :
      :"g" (start), "g" (iz), "g" (diz), "g" (len), "g" (poly)
      :"eax", "ebx", "ecx", "edx", "edi"
    );
//    for (int x = x1; x <= x2; x++) {
//      if (iz > zmap[x + y_offset].inv_z) {
//        zmap[x + y_offset].inv_z = iz;
//        zmap[x + y_offset].poly = poly;
//      }
//      iz += diz;
//    }
    l += lef.m;
    r += rig.m;
    y_offset += _screen.x_res;
  }
}


inline void zmap_yz_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  double a, double b, int c, long diz)
{
  double top = ifloor(y1 + .5);
  double bot = ifloor(y2 - .5);
  double l = lef.m * (top + .5) + lef.b;
  double r = rig.m * (top + .5) + rig.b;
  int y_offset = _screen.x_res * (int)top;
  long iz = (unsigned long)(c * (top + .5 + a) * b);
  
  for (int y = int(top); y <= (int)bot; y++) {
    int x1 = (int)ifloor(l + .5);
    int x2 = (int)ifloor(r - .5);
    unsigned long start = (unsigned long)(&zmap[x1 + y_offset]);
    int len = x2 - x1;

    if (x1 <= x2)
    asm (
      "push %%esi \n"
      "movl %0, %%ebx \n"
      "movl %1, %%eax \n"
      "movl %2, %%ecx \n"
      "movl %3, %%edi \n"
      "push %%ebp \n"
      "movl __pxl_ct, %%ebp \n"
      ".align 4 \n"
      "rep_zmap_yz: \n"
      "  movl 0(%%ebx), %%esi \n"
      "  cmpl %%esi, %%eax \n"
      "  jle skip_zmap_yz \n"
      "    movl %%eax, 0(%%ebx) \n"
      "    movl %%edi, 4(%%ebx) \n"
      "    incl %%ebp \n"      
      "  skip_zmap_yz: \n"
      "  addl $8, %%ebx \n"
      "decl %%ecx \n"
      "jge rep_zmap_yz \n"
      "movl %%ebp, __pxl_ct \n"
      "pop %%ebp \n"
      "pop %%esi \n"
      :
      :"g" (start), "g" (iz), "g" (len), "g" (poly)
      :"eax", "ebx", "ecx", "edx", "edi"
    );

    l += lef.m;
    r += rig.m;
    iz += diz;
    y_offset += _screen.x_res;
  }
}


inline void zmap_z_sect (double y1, double y2,
  lin_relat lef, lin_relat rig, z_buffer_type* zmap, int poly, view_type view,
  long iz)
{
  double top = ifloor(y1 + .5);
  double bot = ifloor(y2 - .5);
  double l = lef.m * (top + .5) + lef.b;
  double r = rig.m * (top + .5) + rig.b;
  int y_offset = _screen.x_res * (int)top;

  for (int y = int(top); y <= (int)bot; y++) {
    int x1 = (int)ifloor(l + .5);
    int x2 = (int)ifloor(r - .5);
    unsigned long start = (unsigned long)(&zmap[x1 + y_offset]);
    int len = x2 - x1;

    if (x1 <= x2)
    asm (
      "push %%esi \n"
      "movl %0, %%ebx \n"
      "movl %1, %%eax \n"
      "movl %2, %%ecx \n"
      "movl %3, %%edi \n"
      "push %%ebp \n"
      "movl __pxl_ct, %%ebp \n"
      ".align 4 \n"
      "rep_zmap_z: \n"
      "  movl 0(%%ebx), %%esi \n"
      "  cmpl %%esi, %%eax \n"
      "  jle skip_zmap_z \n"
      "    movl %%eax, 0(%%ebx) \n"
      "    movl %%edi, 4(%%ebx) \n"
      "    incl %%ebp \n"      
      "  skip_zmap_z: \n"
      "  addl $8, %%ebx \n"
      "decl %%ecx \n"
      "jge rep_zmap_z \n"
      "movl %%ebp, __pxl_ct \n"
      "pop %%ebp \n"
      "pop %%esi \n"      
      :
      :"g" (start), "g" (iz), "g" (len), "g" (poly)
      :"eax", "ebx", "ecx", "edx", "edi"
    );

    l += lef.m;
    r += rig.m;
    y_offset += _screen.x_res;
  }
}


inline void clear_zmap (view_type view, z_buffer_type* zmap)
{
  int len = _screen.x_res * _screen.y_res;
  
  asm (
    "movl %0, %%ecx \n"
    "movl %1, %%ebx \n"
    "xorl %%eax, %%eax \n"
    ".align 4 \n"
    "itscrap: \n"
    "  movl %%eax, 0(%%ebx, %%ecx, 8) \n"
    "decl %%ecx \n"
    "jnl itscrap \n"
    :
    :"m" (len), "m" (zmap)
    :"eax", "ebx", "ecx"
  );
}


#endif // !INCLUDE_ZMAP
