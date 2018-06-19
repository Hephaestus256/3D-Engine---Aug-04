#ifndef INCLUDE_MATH_3D
#define INCLUDE_MATH_3D

#include <math.h>
#include "math2d.h"

#ifndef CMP_PRECISION
#define CMP_PRECISION (1.0 / 1000000.0)
#endif

struct point_3d {
  double x, y, z;
};

struct line_3d {
  double mx, my, bx, by;
};

struct plane_type {
// z = m1 * x + b1
// b1 = m2 * y + b
  double m1, m2, b;
  double perp_m1, perp_m2;
  bool y_plane;
  bool m1_inf;
  double sec_m1;
  double csc_m1;
};

struct cylinder_3d {
  double radius;
  double height;
};

struct x_line_3d {
  double my, mz, by, bz;
};

struct y_line_3d {
  double mx, mz, bx, bz;
};


inline double dist_3d (point_3d p1, point_3d p2);
inline point_3d offset_3d (double dist, double angle_xz, double angle_yz);
inline point_3d p3d (double x, double y, double z);
inline line_3d calc_line_3d (point_3d p1, point_3d p2);
inline void calc_line_3d (line_3d* l3d, point_3d* p1, point_3d* p2);
inline point_3d intrapolate_3d (point_3d p1, point_3d p2, double z);
inline void calc_line_3d (line_3d* l3d, point_3d p1, point_3d p2);
inline point_3d intersect_line_plane_3d (plane_type plane, line_3d line);
inline x_line_3d calc_x_line_3d (point_3d p1, point_3d p2);
inline y_line_3d calc_y_line_3d (point_3d p1, point_3d p2);
inline void calc_y_line_3d (y_line_3d* line, point_3d p1, point_3d p2);
inline int point_rel_to_plane (plane_type plane, point_3d point, double rad);
inline int point_relative_to_plane (plane_type plane, point_3d point);
void calc_plane
  (plane_type* plane, point_3d point1, point_3d point2, point_3d point3);
inline bool is_within_xz_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3);
inline bool is_within_xy_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3);
inline bool is_within_zy_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3);
inline bool line_intersect_tri (line_3d line, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p);
inline bool line_intersect_tri_y (point_3d ref, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p);
inline bool line_intersect_tri_xy (lin_relat line, double z, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p);
inline y_line_3d parll_y_line_3d (y_line_3d line, point_3d point);
inline point_3d point_perp_line (y_line_3d line, point_3d point);
inline lin_relat calc_dz_dx_line (point_3d p1, point_3d p2);
inline double dz_dx (point_3d p1, point_3d p2);
inline point_3d offset_point_3d
  (point_3d point, double x, double y, double z);
inline line_3d calc_parll_line_3d (line_3d line, point_3d point);
inline line_equat calc_3d_dydx_line (point_3d p1, point_3d p2);
inline lin_relat calc_dydx_line (point_3d p1, point_3d p2);
inline lin_relat calc_dxdy_line (point_3d p1, point_3d p2);
inline lin_relat calc_dzdy_line (point_3d p1, point_3d p2);
inline lin_relat calc_parll_dzdy_line
  (point_3d p1, point_3d p2, double offset);
inline lin_relat calc_dydz_line (point_3d p1, point_3d p2);
inline lin_relat calc_z_parll_dydz_line
  (point_3d p1, point_3d p2, double offset);
inline double sqrd (double n);
inline bool closer (point_3d ref, point_3d p1, point_3d p2);
inline lin_relat calc_dxdz_line (point_3d p1, point_3d p2);
inline dbl_pair angles (point_3d p1, point_3d p2);
inline bool line_intersect_tri_x (point_3d ref, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p);
inline bool y_line_intersect_line_yz
  (point_3d p1a, point_3d p1b, point_3d p2a, point_3d p2b,
  lin_relat l, double z_offet, point_3d* result);
inline bool y_line_intersect_tri_yz
  (point_3d p1, point_3d p2, point_3d p3, point_3d pa, point_3d pb);
inline point_3d offset_3d
  (double r, point_3d p, double xz, double yz, double xy);
inline point_3d offset_3d_xz (double r, point_3d p, double a);
inline void get_x_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3);
inline void get_y_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3);
inline void get_z_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3);
inline int point_approx_to_plane (plane_type plane, point_3d point);


inline line_equat calc_3d_dydx_line (point_3d p1, point_3d p2)
{
  line_equat l;

  l.m = (p2.y - p1.y) / (p2.x - p1.x);
  l.b = p1.y - l.m * p1.x;
  
  return l;
}


inline line_3d calc_parll_line_3d (line_3d line, point_3d point)
{
  line_3d l3d;

  l3d.mx = line.mx;
  l3d.my = line.my;
  l3d.bx = point.x - point.z * line.mx;
  l3d.by = point.y - point.z * line.my;
  
  return l3d;
}


inline point_3d p3d (double x, double y, double z)
{
  point_3d temp;
  temp.x = x;
  temp.y = y;
  temp.z = z;
  return temp;
}


inline point_3d offset_point_3d
  (point_3d point, double x, double y, double z)
{
  point_3d temp;
  
  temp.x = point.x + x;
  temp.y = point.y + y;
  temp.z = point.z + z;

  return temp;
}


inline double dz_dx (point_3d p1, point_3d p2)
{
  return (p2.z - p1.z) / (p2.x - p1.x);
}


/*
inline lin_relat calc_dz_dx_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = dz_dx (p1, p2);
  line.b = p1.z - p1.x * line.m;
  
  return line;
}
*/

inline point_3d point_perp_line (y_line_3d line, point_3d point)
{
  point_3d temp;

  temp.y = (line.mx * point.x + line.mz * point.z + point.y -
           line.mx * line.bx - line.mz * line.bz) /
           (line.mx * line.mx + line.mz * line.mz + 1);
  temp.x = line.mx * temp.y + line.bx;
  temp.z = line.mz * temp.y + line.bz;

  return temp;
}


inline y_line_3d parll_y_line_3d (y_line_3d line, point_3d point)
{
  y_line_3d temp;

  temp.mx = line.mx;
  temp.mz = line.mz;
  temp.bx = point.x - point.y * temp.mx;
  temp.bz = point.z - point.y * temp.mz;

  return temp;
}


void calc_plane
  (plane_type* plane, point_3d point1, point_3d point2, point_3d point3)
{
// plane:
//     z = m1x + b1
//     b1 = m2y + b
//     z = m1x + m2y + b

  y_line_3d line;
  point_3d point_b;
  
  if (approx_equal (point1.y, point2.y))
    if (approx_equal (point1.y, point3.y)) { // 1 = 2 = 3
      plane->b = point1.y;
      plane->y_plane = true;
      plane->m1_inf = false;
    }
    else // (1 = 2) != 3
      if (!approx_equal (point1.x, point2.x)) {
        plane->m1 = (point2.z - point1.z) / (point2.x - point1.x);
        double b1 = point3.z - point3.x * plane->m1;
        double b2 = point1.z - point1.x * plane->m1;
        plane->m2 = (b2 - b1) / (point1.y - point3.y);
        plane->b = b1 - plane->m2 * point3.y;
        plane->y_plane = false;
        plane->m1_inf = false;
      }
      else {
        plane->m2 = (point3.x - point2.x) / (point3.y - point2.y);
        plane->perp_m2 = -(point3.y - point2.y) / (point3.x - point2.x);
        plane->b = point2.x - point2.y * plane->m2;
        plane->y_plane = false;
        plane->m1_inf = true;
      }
  else
    if (approx_equal (point1.y, point3.y)) { // (3 = 1) != 2
      calc_y_line_3d (&line, point2, point3);
      point_b.x = line.mx * point1.y + line.bx;
      point_b.z = line.mz * point1.y + line.bz;
      if (!approx_equal(point_b.x, point1.x)) {
        plane->m1 = (point1.z - point3.z) / (point1.x - point3.x);
        double b1 = point2.z - point2.x * plane->m1;
        double b2 = point3.z - point3.x * plane->m1;
        plane->m2 = (b2 - b1) / (point3.y - point2.y);
        plane->b = b1 - plane->m2 * point2.y;
        plane->y_plane = false;
        plane->m1_inf = false;
      }
      else { // possible bugs here
        line_3d z_line;
        if (!approx_equal (point1.z, point2.z)) {
          calc_line_3d (&z_line, point1, point2);
          double x = z_line.mx * point3.z + z_line.bx;
          double y = z_line.my * point3.z + z_line.by;
          plane->m2 = (point3.x - x) / (point3.y - y);
          plane->perp_m2 = -(point3.y - y) / (point3.x - x);
          plane->b = x - y * plane->m2;
        }
        else {
          calc_line_3d (&z_line, point2, point3);
          double x = z_line.mx * point1.z + z_line.bx;
          double y = z_line.my * point1.z + z_line.by;
          plane->m2 = (point1.x - x) / (point1.y - y);
          plane->perp_m2 = -(point1.y - y) / (point1.x - x);
          plane->b = x - y * plane->m2;
        }
        plane->m1_inf = true;
        plane->y_plane = false;
      }
    }
    else {
      calc_y_line_3d (&line, point1, point2);
      point_b.x = line.mx * point3.y + line.bx;
      point_b.z = line.mz * point3.y + line.bz;
      if (!approx_equal (point_b.x, point3.x)) {
        plane->m1 = (point_b.z - point3.z) / (point_b.x - point3.x);
        double b1 = point3.z - point3.x * plane->m1;
        double b2 = point1.z - point1.x * plane->m1;
        plane->m2 = (b2 - b1) / (point1.y - point3.y);
        plane->b = b1 - plane->m2 * point3.y;
        plane->y_plane = false;
        plane->m1_inf = false;
      }
      else {
        line_3d z_line;
        if (!approx_equal (point1.z, point2.z)) {
          calc_line_3d (&z_line, point1, point2);
          double x = z_line.mx * point3.z + z_line.bx;
          double y = z_line.my * point3.z + z_line.by;
          plane->m2 = (point3.x - x) / (point3.y - y);
          plane->b = x - y * plane->m2;
        }
        else {
          calc_line_3d (&z_line, point2, point3);
          double x = z_line.mx * point1.z + z_line.bx;
          double y = z_line.my * point1.z + z_line.by;
          plane->m2 = (point1.x - x) / (point1.y - y);
          plane->b = x - y * plane->m2;
        }
        plane->m1_inf = true;
        plane->y_plane = false;
      }
    }

  double am1 = atan(plane->m1);
  plane->sec_m1 = 1 / cos(am1);
  plane->csc_m1 = 1 / sin(am1);
  plane->perp_m1 = -1 / plane->m1;
  plane->perp_m2 = -1 / plane->m2;
}


inline int point_relative_to_plane (plane_type plane, point_3d point)
{
  if (plane.y_plane)
    if (point.y < plane.b)
      return -1;
    else if (point.y > plane.b)
      return 1;
    else
      return 0;
  else {
    if (plane.m1_inf) {
      double x = plane.m2 * point.y + plane.b;
      if (point.x < x)
        return -1;
      else if (point.x > x)
        return 1;
      else
        return 0;
    }
    else {
      double z = plane.m1 * point.x + plane.m2 * point.y + plane.b;
      if (point.z < z)
        return -1;
      else if (point.z > z)
        return 1;
      else
        return 0;
    }
  }
}


inline int point_approx_to_plane (plane_type plane, point_3d point)
{
  if (plane.y_plane) {
    if (approx_equal (point.y, plane.b))
      return 0;
    else if (point.y < plane.b)
      return -1;
    else
      return 1;
  }
  else {
    if (plane.m1_inf) {
       double x = plane.m2 * point.y + plane.b;
      if (approx_equal (point.x, x))
        return 0;
      else if (point.x < x)
        return -1;
      else
        return 1;
    }
    else {
      double z = plane.m1 * point.x + plane.m2 * point.y + plane.b;
      if (approx_equal (point.z, z))
        return 0;
      else if (point.z < z)
        return -1;
      else
        return 1;
    }
  }
}


inline int point_rel_to_plane (plane_type plane, point_3d point, double rad)
{
  if (plane.y_plane) {
    if (point.y < plane.b + rad) {
      return -1;
    }
    else if (point.y > plane.b + rad) {
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    if (plane.m1_inf) {
      double x = plane.m2 * point.y + plane.b;
      if (point.x < x)
        return -1;
      else if (point.x > x)
        return 1;
      else
        return 0;
    }
    else {
      double z = plane.m1 * point.x + plane.m2 * point.y + plane.b;
      if (point.z < z)
        return -1;
      else if (point.z > z)
        return 1;
      else
        return 0;
    }
  }
}



inline bool line_intersect_tri_x (point_3d ref, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p)
{
  if (plane.y_plane)
    return false;
  else if (plane.m1_inf) {
    p->x = plane.m2 * ref.y + plane.b;
    p->y = ref.y;
    p->z = ref.z;
    return is_within_zy_tri (*p, p1, p2, p3);
  }
  else
    if (is_approx_zero (plane.m1))
      return false;
    else {
      p->x = (ref.z - plane.m2 * ref.y - plane.b) / plane.m1;
      p->y = ref.y;
      p->z = ref.z;
      return is_within_zy_tri (*p, p1, p2, p3);
    }
}


inline bool line_intersect_tri_y (point_3d ref, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p)
{
  if (plane.y_plane) {
    p->x = ref.x;
    p->y = plane.b;
    p->z = ref.z;
    return is_within_xz_tri (*p, p1, p2, p3);
  }
  else if (plane.m1_inf)
    if (is_approx_zero(plane.m2))
      return false;
    else {
      p->x = ref.x;
      p->y = (ref.x - plane.b) / plane.m2;
      p->z = ref.z;
      return is_within_zy_tri (*p, p1, p2, p3);
    }
  else
    if (is_approx_zero(plane.m2))
      return false;
    else {
      p->x = ref.x;
      p->y = (ref.z - plane.b - plane.m1 * ref.x) / plane.m2;
      p->z = ref.z;
      return is_within_xy_tri (*p, p1, p2, p3);      
    }
}


inline bool line_intersect_tri_xy (line_3d line, double z, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p)
{
  if (plane.y_plane)
    if (is_approx_zero(line.mx))
      return false;
    else {
      p->x = (plane.b - line.bx) / line.mx;
      p->y = plane.b;
      p->z = z;
      return is_within_xz_tri (*p, p1, p2, p3);
    }
  else if (plane.m1_inf) {
    double denom = 1 - line.mx * plane.m2;
    if (is_approx_zero(denom))
      return false;
    else {
      p->y = (line.mx * plane.b + line.bx) / denom;
      p->x = plane.m2 * p->y + plane.b;
      p->z = z;
      return is_within_zy_tri (*p, p1, p2, p3);
    }
  }
  else
    if (is_approx_zero(plane.m2)) {
      p->x = plane.m1 * z + plane.b;
      p->y = line.mx * (plane.m1 * z + plane.b) + line.bx;
      p->z = z;
    }
    else {
      double denom = line.mx + plane.m1 / plane.m2;
      if (is_approx_zero(denom))
        return false;
      else {
        p->x = ((z - plane.b) / plane.m2 - line.bx) / denom;
        p->y = line.mx * p->x + line.bx;
        p->z = z;
        return is_within_xy_tri (*p, p1, p2, p3);
      }
    }
}


inline bool line_intersect_tri (line_3d line, plane_type plane,
  point_3d p1, point_3d p2, point_3d p3, point_3d* p)
{
  if (plane.y_plane)
    if (is_approx_zero(line.my))
      return false;
    else {
      if (is_approx_zero(line.mx)) {
        p->z = (plane.b - line.by) / line.my;
        p->y = plane.b;
        p->x = line.mx * p->z + line.bx;
      }
      else {
        p->x = line.mx * (plane.b - line.by) / line.my + line.bx;
        p->y = plane.b;
        p->z = (p->x - line.bx) / line.mx;
      }
      return is_within_xz_tri (*p, p1, p2, p3);
    }
  else if (plane.m1_inf) {
    double denom = plane.m2 * line.my - line.mx;
    if (is_approx_zero(denom))
      return false;
    else {
      p->z = (line.bx - plane.m2 * line.by - plane.b) / denom;
      p->x = line.mx * p->z + line.bx;
      p->y = line.my * p->z + line.by;
      return is_within_zy_tri (*p, p1, p2, p3);
    }
  }
  else {
    double denom = 1 - plane.m1 * line.mx - plane.m2 * line.my;
    if (is_approx_zero(denom))
      return false;
    else {
      p->z = (plane.m1 * line.bx + plane.m2 * line.by + plane.b) / denom;
      p->x = line.mx * p->z + line.bx;
      p->y = line.my * p->z + line.by;
      return is_within_xy_tri (*p, p1, p2, p3);
    }
  }
}


inline bool is_within_xz_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3)
{
  if (approx_equal(t1.z, t2.z))
    if (t3.z > t1.z) {
      lin_relat l = calc_dxdz_line (t2, t3);
      lin_relat r = calc_dxdz_line (t1, t3);
      if (p.z < t1.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l = calc_dxdz_line (t1, t3);
      lin_relat r = calc_dxdz_line (t2, t3);
      if (p.z > t1.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
  else if (t1.z < t2.z)
    if (approx_equal(t3.z, t1.z)) {
      lin_relat l = calc_dxdz_line (t2, t1);
      lin_relat r = calc_dxdz_line (t2, t3);
      if (p.z < t1.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.z, t2.z)) {
      lin_relat l = calc_dxdz_line (t1, t2);
      lin_relat r = calc_dxdz_line (t1, t3);
      if (p.z > t2.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
    else if (t3.z < t1.z) {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x < l12.m * p.z + l12.b)
        return false;
      else if (p.x > l23.m * p.z + l23.b)
        return false;
      else if (p.x < l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.z > t2.z) {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x < l12.m * p.z + l12.b)
        return false;
      else if (p.x < l23.m * p.z + l23.b)
        return false;
      else if (p.x > l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x < l12.m * p.z + l12.b)
        return false;
      else if (p.x > l23.m * p.z + l23.b)
        return false;
      else if (p.x > l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
  else // p1.z > p2.z
    if (approx_equal(t3.z, t1.z)) {
      lin_relat l = calc_dxdz_line (t2, t3);
      lin_relat r = calc_dxdz_line (t2, t1);
      if (p.z > t1.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.z, t2.z)) {
      lin_relat l = calc_dxdz_line (t1, t3);
      lin_relat r = calc_dxdz_line (t1, t2);
      if (p.z < t2.z)
        return false;
      else if (p.x < l.m * p.z + l.b)
        return false;
      else if (p.x > r.m * p.z + r.b)
        return false;
      else
        return true;
    }
    else if (t3.z < t2.z) {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x > l12.m * p.z + l12.b)
        return false;
      else if (p.x > l23.m * p.z + l23.b)
        return false;
      else if (p.x < l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.z > t1.z) {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x > l12.m * p.z + l12.b)
        return false;
      else if (p.x < l23.m * p.z + l23.b)
        return false;
      else if (p.x > l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dxdz_line (t1, t2);
      lin_relat l23 = calc_dxdz_line (t2, t3);
      lin_relat l31 = calc_dxdz_line (t3, t1);
      if (p.x > l12.m * p.z + l12.b)
        return false;
      else if (p.x < l23.m * p.z + l23.b)
        return false;
      else if (p.x < l31.m * p.z + l31.b)
        return false;
      else
        return true;
    }
}


inline bool is_within_xy_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3)
{
  if (approx_equal(t1.y, t2.y))
    if (t3.y > t1.y) {
      lin_relat l = calc_dxdy_line (t1, t3);
      lin_relat r = calc_dxdy_line (t2, t3);
      if (p.y < t1.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l = calc_dxdy_line (t2, t3);
      lin_relat r = calc_dxdy_line (t1, t3);
      if (p.y > t1.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
  else if (t1.y < t2.y)
    if (approx_equal(t3.y, t1.y)) {
      lin_relat l = calc_dxdy_line (t2, t3);
      lin_relat r = calc_dxdy_line (t2, t1);
      if (p.y < t1.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.y, t2.y)) {
      lin_relat l = calc_dxdy_line (t1, t3);
      lin_relat r = calc_dxdy_line (t1, t2);
      if (p.y > t2.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (t3.y < t1.y) {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x > l12.m * p.y + l12.b)
        return false;
      else if (p.x < l23.m * p.y + l23.b)
        return false;
      else if (p.x > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.y > t2.y) {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x > l12.m * p.y + l12.b)
        return false;
      else if (p.x > l23.m * p.y + l23.b)
        return false;
      else if (p.x < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x > l12.m * p.y + l12.b)
        return false;
      else if (p.x < l23.m * p.y + l23.b)
        return false;
      else if (p.x < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
  else // p1.y > p2.y
    if (approx_equal(t3.y, t1.y)) {
      lin_relat l = calc_dxdy_line (t2, t1);
      lin_relat r = calc_dxdy_line (t2, t3);
      if (p.y > t1.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.y, t2.y)) {
      lin_relat l = calc_dxdy_line (t1, t2);
      lin_relat r = calc_dxdy_line (t1, t3);
      if (p.y < t2.y)
        return false;
      else if (p.x < l.m * p.y + l.b)
        return false;
      else if (p.x > r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (t3.y < t2.y) {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x < l12.m * p.y + l12.b)
        return false;
      else if (p.x < l23.m * p.y + l23.b)
        return false;
      else if (p.x > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.y > t1.y) {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x < l12.m * p.y + l12.b)
        return false;
      else if (p.x > l23.m * p.y + l23.b)
        return false;
      else if (p.x < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dxdy_line (t1, t2);
      lin_relat l23 = calc_dxdy_line (t2, t3);
      lin_relat l31 = calc_dxdy_line (t3, t1);
      if (p.x < l12.m * p.y + l12.b)
        return false;
      else if (p.x > l23.m * p.y + l23.b)
        return false;
      else if (p.x > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
}


inline bool is_within_zy_tri
  (point_3d p, point_3d t1, point_3d t2, point_3d t3)
{
  if (approx_equal(t1.y, t2.y))
    if (t3.y > t1.y) {
      lin_relat l = calc_dzdy_line (t1, t3);
      lin_relat r = calc_dzdy_line (t2, t3);
      if (p.y < t1.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l = calc_dzdy_line (t2, t3);
      lin_relat r = calc_dzdy_line (t1, t3);
      if (p.y > t1.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
  else if (t1.y < t2.y)
    if (approx_equal(t3.y, t1.y)) {
      lin_relat l = calc_dzdy_line (t2, t3);
      lin_relat r = calc_dzdy_line (t2, t1);
      if (p.y < t1.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.y, t2.y)) {
      lin_relat l = calc_dzdy_line (t1, t3);
      lin_relat r = calc_dzdy_line (t1, t2);
      if (p.y > t2.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (t3.y < t1.y) {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z < l12.m * p.y + l12.b)
        return false;
      else if (p.z > l23.m * p.y + l23.b)
        return false;
      else if (p.z < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.y > t2.y) {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z < l12.m * p.y + l12.b)
        return false;
      else if (p.z < l23.m * p.y + l23.b)
        return false;
      else if (p.z > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z < l12.m * p.y + l12.b)
        return false;
      else if (p.z > l23.m * p.y + l23.b)
        return false;
      else if (p.z > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
  else // p1.y > p2.y
    if (approx_equal(t3.y, t1.y)) {
      lin_relat l = calc_dzdy_line (t2, t1);
      lin_relat r = calc_dzdy_line (t2, t3);
      if (p.y > t1.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (approx_equal(t3.y, t2.y)) {
      lin_relat l = calc_dzdy_line (t1, t2);
      lin_relat r = calc_dzdy_line (t1, t3);
      if (p.y < t2.y)
        return false;
      else if (p.z > l.m * p.y + l.b)
        return false;
      else if (p.z < r.m * p.y + r.b)
        return false;
      else
        return true;
    }
    else if (t3.y < t2.y) {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z > l12.m * p.y + l12.b)
        return false;
      else if (p.z > l23.m * p.y + l23.b)
        return false;
      else if (p.z < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else if (t3.y > t1.y) {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z > l12.m * p.y + l12.b)
        return false;
      else if (p.z < l23.m * p.y + l23.b)
        return false;
      else if (p.z > l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
    else {
      lin_relat l12 = calc_dzdy_line (t1, t2);
      lin_relat l23 = calc_dzdy_line (t2, t3);
      lin_relat l31 = calc_dzdy_line (t3, t1);
      if (p.z > l12.m * p.y + l12.b)
        return false;
      else if (p.z < l23.m * p.y + l23.b)
        return false;
      else if (p.z < l31.m * p.y + l31.b)
        return false;
      else
        return true;
    }
}


inline void calc_y_line_3d (y_line_3d* line, point_3d p1, point_3d p2)
{
//  x = m1y + b1
//  z = m2y + b2

  line->mx = (p2.x - p1.x) / (p2.y - p1.y);
  line->bx = p1.x - line->mx * p1.y;
  line->mz = (p2.z - p1.z) / (p2.y - p1.y);
  line->bz = p1.z - line->mz * p1.y;
}


inline y_line_3d calc_y_line_3d (point_3d p1, point_3d p2)
{
//  x = m1y + b1
//  z = m2y + b2
  y_line_3d line;
  
  line.mx = (p2.x - p1.x) / (p2.y - p1.y);
  line.bx = p1.x - line.mx * p1.y;
  line.mz = (p2.z - p1.z) / (p2.y - p1.y);
  line.bz = p1.z - line.mz * p1.y;

  return line;
}


inline x_line_3d calc_x_line_3d (point_3d p1, point_3d p2)
{
//  y = m1x + b1
//  z = m2x + b2
  x_line_3d line;
  
  line.my = (p2.y - p1.y) / (p2.x - p1.x);
  line.by = p1.y - line.my * p1.x;
  line.mz = (p2.z - p1.z) / (p2.x - p1.x);
  line.bz = p1.z - line.mz * p1.x;

  return line;
}


inline point_3d intersect_line_plane_3d (plane_type plane, line_3d line)
{
  point_3d p;

  p.z = (line.bx * plane.m1 + line.by * plane.m2 + plane.b) /
        (1 - plane.m1 * line.mx - plane.m2 * line.my);
  p.x = line.mx * p.z + line.bx;
  p.y = line.my * p.z + line.by;
  
  return p;
}


inline void calc_line_3d (line_3d* l3d, point_3d p1, point_3d p2)
{
  l3d->mx = (p2.x - p1.x) / (p2.z - p1.z);
  l3d->my = (p2.y - p1.y) / (p2.z - p1.z);
  l3d->bx = p1.x - l3d->mx * p1.z;
  l3d->by = p1.y - l3d->my * p1.z;
}


inline void calc_line_3d (line_3d* l3d, point_3d* p1, point_3d* p2)
{
  l3d->mx = (p2->x - p1->x) / (p2->z - p1->z);
  l3d->my = (p2->y - p1->y) / (p2->z - p1->z);
  l3d->bx = p1->x - l3d->mx * p1->z;
  l3d->by = p1->y - l3d->my * p1->z;
}


inline point_3d intrapolate_3d (point_3d p1, point_3d p2, double z)
{
  point_3d temp;
  double dx = (p2.x - p1.x) / (p2.z - p1.z);
  double dy = (p2.y - p1.y) / (p2.z - p1.z);

  temp.x = p1.x + dx * (z - p1.z);
  temp.y = p1.y + dy * (z - p1.z);
  temp.z = z;

  return temp;
}


inline line_3d calc_line_3d (point_3d p1, point_3d p2)
{
  line_3d l3d;
  
  l3d.mx = (p2.x - p1.x) / (p2.z - p1.z);
  l3d.my = (p2.y - p1.y) / (p2.z - p1.z);
  l3d.bx = p1.x - l3d.mx * p1.z;
  l3d.by = p1.y - l3d.my * p1.z;

  return l3d;
}


inline lin_relat calc_dxdz_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.x - p1.x) / (p2.z - p1.z);
  line.b = p1.x - p1.z * line.m;
  
  return line;
}


inline lin_relat calc_dzdx_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.z - p1.z) / (p2.x - p1.x);
  line.b = p1.z - p1.x * line.m;
  
  return line;
}


inline lin_relat calc_dzdy_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.z - p1.z) / (p2.y - p1.y);
  line.b = p1.z - p1.y * line.m;
  
  return line;
}


inline lin_relat calc_parll_dzdy_line
  (point_3d p1, point_3d p2, double offset)
{
  lin_relat line;

  line.m = (p2.z - p1.z) / (p2.y - p1.y);
  line.b = p1.z - p1.y * line.m + offset;
  
  return line;
}


inline lin_relat calc_dydz_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.y - p1.y) / (p2.z - p1.z);
  line.b = p1.y - p1.z * line.m;
  
  return line;
}


inline lin_relat calc_z_parll_dydz_line
  (point_3d p1, point_3d p2, double offset)
{
  lin_relat line;

  line.m = (p2.y - p1.y) / (p2.z - p1.z);
  line.b = p1.y - (p1.z + offset) * line.m;
  
  return line;
}


inline lin_relat calc_dydx_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.y - p1.y) / (p2.x - p1.x);
  line.b = p1.y - p1.x * line.m;
  
  return line;
}


inline lin_relat calc_dxdy_line (point_3d p1, point_3d p2)
{
  lin_relat line;

  line.m = (p2.x - p1.x) / (p2.y - p1.y);
  line.b = p1.x - p1.y * line.m;
  
  return line;
}


inline point_3d offset_3d (double dist, double angle_xz, double angle_yz)
{
  point_3d p3d;
  
  p3d.x = sin(angle_xz) * cos(angle_yz) * dist;
  p3d.y = -sin(angle_yz) * dist;
  p3d.z = cos(angle_xz) * cos(angle_yz) * dist;

  return p3d;
}


inline double sqrd (double n)
{
  return n * n;
}


inline bool closer (point_3d ref, point_3d p1, point_3d p2)
{
  double d1 =
    sqrt(sqrd(ref.x - p1.x) + sqrd(ref.y - p1.y) + sqrd(ref.z - p1.z));
  double d2 =
    sqrt(sqrd(ref.x - p2.x) + sqrd(ref.y - p2.y) + sqrd(ref.z - p2.z));
    
  if (d1 < d2)
    return true;
  else
    return false;
}


inline dbl_pair angles (point_3d p1, point_3d p2)
{
  dbl_pair a;

  a.x = atan2 (p2.z - p1.z, p2.x - p1.x);
  a.y = atan2 (p2.y - p1.y, p2.z - p1.z);

  return a;
}


inline double dist_3d (point_3d p1, point_3d p2)
{
  return sqrt(sqrd(p2.x - p1.x) + sqrd(p2.y - p1.y) + sqrd(p2.z - p1.z));
}


inline bool y_line_intersect_line_yz
  (point_3d p1a, point_3d p1b, point_3d p2a, point_3d p2b,
  lin_relat l, double z_offset, point_3d* result)
{
  // 1: y-line
  // 2: y = mz + b    m = dy/dz
  result->x = p1a.x;
  result->y = l.m * p1a.z + l.b;
  result->z = p1a.z;
  
  if (approx_between (result->y, p1a.y, p1b.y))
    if (approx_between (p1a.z, p2a.z + z_offset, p2b.z + z_offset))
      return true;

  return false;
}


inline point_3d offset_3d_xz (double r, point_3d p, double a)
{
  double dx = r * cos(a);
  double dz = r * sin(a);
  point_3d result;
  
  result.x = p.x + dx;
  result.y = p.y;
  result.z = p.z + dz;

  return result;
}


inline void get_x_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.x <= p2.x)
    if (p3.x <= p1.x) { // 3,1,2
      *min = p3.x - radius;
      *max = p2.x + radius;
    }
    else if (p3.x >= p2.x) { // 1,2,3
      *min = p1.x - radius;
      *max = p3.x + radius;
    }
    else { // 1,3,2
      *min = p1.x - radius;
      *max = p2.x + radius;
    }
  else
    if (p3.x <= p2.x) { // 3,2,1
      *min = p3.x - radius;
      *max = p1.x + radius;
    }
    else if (p3.x >= p1.x) { // 2,1,3
      *min = p2.x - radius;
      *max = p3.x + radius;
    }
    else { // 2,3,1
      *min = p2.x - radius;
      *max = p1.x + radius;
    }
}


inline void get_y_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.y <= p2.y)
    if (p3.y <= p1.y) { // 3,1,2
      *min = p3.y - radius;
      *max = p2.y + radius;
    }
    else if (p3.y >= p2.y) { // 1,2,3
      *min = p1.y - radius;
      *max = p3.y + radius;
    }
    else { // 1,3,2
      *min = p1.y - radius;
      *max = p2.y + radius;
    }
  else
    if (p3.y <= p2.y) { // 3,2,1
      *min = p3.y - radius;
      *max = p1.y + radius;
    }
    else if (p3.y >= p1.y) { // 2,1,3
      *min = p2.y - radius;
      *max = p3.y + radius;
    }
    else { // 2,3,1
      *min = p2.y - radius;
      *max = p1.y + radius;
    }
}


inline void get_z_range_3d
  (double* min, double* max, double radius,
  point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.z <= p2.z)
    if (p3.z <= p1.z) { // 3,1,2
      *min = p3.z - radius;
      *max = p2.z + radius;
    }
    else if (p3.z >= p2.z) { // 1,2,3
      *min = p1.z - radius;
      *max = p3.z + radius;
    }
    else { // 1,3,2       `
      *min = p1.z - radius;
      *max = p2.z + radius;
    }
  else
    if (p3.z <= p2.z) { // 3,2,1
      *min = p3.z - radius;
      *max = p1.z + radius;
    }
    else if (p3.z >= p1.z) { // 2,1,3
      *min = p2.z - radius;
      *max = p3.z + radius;
    }
    else { // 2,3,1
      *min = p2.z - radius;
      *max = p1.z + radius;
    }
}


#endif // !INCLUDE_MATH_3D
