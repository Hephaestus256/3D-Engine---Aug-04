#include <stdio.h>
#include <math.h>
#include "math2d.h"
#include "math3d.h"
#include "texclass.h"

#define INFRONT 1
#define BEHIND 0
#define INTERSECT -1
#define CLEAR -2
#define PARALELL -3
#define DONTKNOW -4

bool and_angle
  (double t1, double b1, double t2, double b2, double* t, double* b);
inline void get_view_range_xz (tri_group g, double* b, double* t);
inline void get_view_range_yz (tri_group g, double* b, double* t);
inline void get_view_range_xy (tri_group g, double* b, double* t);
inline bool get_common_range_xz
  (tri_group g1, tri_group g2, double* b, double* t);
inline bool get_common_range_yz
  (tri_group g1, tri_group g2, double* b, double* t);
inline int get_point_case_xz
  (tri_group g, double a, double b, double inv, point_type* p, double angle);
inline int get_point_case_yz
  (tri_group g, double a, double b, double inv, point_type* p, double angle);
inline int get_point_case_yz_pn
  (tri_group gp, tri_group gn, point_type* p, double angle);
inline int get_point_case_yz_in
  (tri_group gi, tri_group gn, point_type* p, double angle);
inline int get_point_case_xy_pn
  (tri_group gp, tri_group gn, point_type* p, double angle);
inline int get_point_case_xy_in
  (tri_group gi, tri_group gn, point_type* p, double angle);
inline int get_point_case_xy_nn
  (tri_group g1, tri_group g2, point_type* p, double angle);
inline int get_tri_case_xz
  (tri_group g, tri_group ref, tri_type tri, double angle);
inline int get_tri_case_xz_in
  (tri_group g, tri_group ref, tri_type tri, double angle);
inline int get_tri_case_yz
  (tri_group g, tri_group ref, tri_type tri, double angle);
inline int get_tri_case_yz_pn
  (tri_group gp, tri_group gn, tri_type tri, double angle);
inline int get_tri_case_yz_in
  (tri_group gp, tri_group gn, tri_type tri, double angle);
inline int get_tri_case_xy_pn
  (tri_group gp, tri_group gn, tri_type tri, double angle);
inline int get_tri_case_xy_ip
  (tri_group gi, tri_group gp, tri_type tri, double angle);
inline int get_tri_case_xy_ii
  (tri_group g1, tri_group g2, tri_type tri, double angle);
inline int get_tri_case_xy_in
  (tri_group gp, tri_group gn, tri_type tri, double angle);
inline int get_tri_case_xy_nn
  (tri_group g1, tri_group g2, tri_type tri, double angle);
inline int tri_rel_tri_xz
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back);
inline int tri_rel_tri_yz
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back);
inline int tri_rel_tri_yz_pn
  (tri_group gp, tri_group gn, tri_type front, tri_type back);
inline int tri_rel_tri_yz_np
  (tri_group gn, tri_group gp, tri_type front, tri_type back);
inline int tri_rel_tri_xy_ip
  (tri_group front_i, tri_group back_p, tri_type front, tri_type back);
inline int tri_rel_tri_xy_pi
  (tri_group front_p, tri_group back_i, tri_type front, tri_type back);
inline int tri_rel_tri_xy_ii
  (tri_group g1, tri_group g2, tri_type front, tri_type back);
inline int tri_rel_tri
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back);
inline int parse_point_cases (int c1, int c2, int c3);
void show_case (int c, int l);


inline bool get_common_range_xz
  (tri_group g1, tri_group g2, double* a1, double* a2)
{
  double b1, t1;
  double b2, t2;
  
  get_view_range_xz (g1, &b1, &t1);
  get_view_range_xz (g2, &b2, &t2);

  bool boo1 = and_angle (b1, t1, b2 - 2 * PI, t2 - 2 * PI, a1, a2);
  bool boo2 = and_angle (b1 - 2 * PI, t1 - 2 * PI, b2, t2, a1, a2);
  bool boo3 = and_angle (b1, t1, b2, t2, a1, a2);

  return boo1 || boo2 || boo3;
}


inline bool get_common_range_yz
  (tri_group g1, tri_group g2, double* a1, double* a2)
{
  double b1, t1;
  double b2, t2;
  
  get_view_range_yz (g1, &b1, &t1);
  get_view_range_yz (g2, &b2, &t2);

  bool boo1 = and_angle (b1, t1, b2 - 2 * PI, t2 - 2 * PI, a1, a2);
  bool boo2 = and_angle (b1 - 2 * PI, t1 - 2 * PI, b2, t2, a1, a2);
  bool boo3 = and_angle (b1, t1, b2, t2, a1, a2);

  return boo1 || boo2 || boo3;
}


inline bool get_common_range_xy
  (tri_group g1, tri_group g2, double* a1, double* a2)
{
  double b1, t1;
  double b2, t2;
  
  get_view_range_xy (g1, &b1, &t1);
  get_view_range_xy (g2, &b2, &t2);

  bool boo1 = and_angle (b1, t1, b2 - 2 * PI, t2 - 2 * PI, a1, a2);
  bool boo2 = and_angle (b1 - 2 * PI, t1 - 2 * PI, b2, t2, a1, a2);
  bool boo3 = and_angle (b1, t1, b2, t2, a1, a2);

  return boo1 || boo2 || boo3;
}


inline void get_view_range_xz (tri_group g, double* b, double* t)
{
  if (g.plane.m1_inf)
    if (g.vis_side == 1) {
      *b = PI;
      *t = 2 * PI;
    }
    else {
      *b = 0;
      *t = PI;
    }
  else
    if (g.vis_side == 1) {
      double a = atan(g.plane.m1);
      if (a < 0)
        *b = a + 2 * PI;
      else
        *b = a;
      *t = *b + PI;
    }
    else {
      *b = atan(g.plane.m1) + PI;
      *t = *b + PI;
    }
}


inline void get_view_range_yz (tri_group g, double* b, double* t)
{
  if (g.plane.y_plane)
    if (g.vis_side == -1) {
      *b = 0.5 * PI;
      *t = 1.5 * PI;
    }
    else {
      *b = 1.5 * PI;
      *t = 2.5 * PI;
    }
  else if (g.plane.m1_inf)
    if (g.vis_side == -1) {
      *b = 0.5 * PI;
      *t = 1.5 * PI;
    }
    else {
      *b = 1.5 * PI;
      *t = 2.5 * PI;
    }
  else
    if (g.vis_side == 1) {
      double a = atan(g.plane.m2);
      if (a < 0)
        *b = a + 2 * PI;
      else
        *b = a;
      *t = *b + PI;
    }
    else {
      *b = atan(g.plane.m2) + PI;
      *t = *b + PI;
    }
}


inline void get_view_range_xy (tri_group g, double* b, double* t)
{
  if (g.plane.y_plane)
    if (g.vis_side == 1) {
      *b = PI;
      *t = 2 * PI;
    }
    else {
      *b = 0;
      *t = PI;
    }
  else if (g.plane.m1_inf)
    if (is_approx_zero (g.plane.m2))
      if (g.vis_side == -1) {
        *b = 0.5 * PI;
        *t = 1.5 * PI;
      }
      else {
        *b = 1.5 * PI;
        *t = 2.5 * PI;
      }
    else if (g.plane.m2 < 0)
      if (g.vis_side == -1) {
        double a = atan(-1 / g.plane.m2);
        if (a < 0)
          *b = a + 2 * PI;
        else
          *b = a;
        *t = *b + PI;
      }
      else {
        *b = atan(-1 / g.plane.m2) + PI;
        *t = *b + PI;
      }
    else
      if (g.vis_side == -1) {
        *b = atan(-1 / g.plane.m2) + PI;
        *t = *b + PI;
      }
      else {
        double a = atan(-1 / g.plane.m2);
        if (a < 0)
          *b = a + 2 * PI;
        else
          *b = a;
        *t = *b + PI;
      }
  else
    if (is_approx_zero (g.plane.m2))
      if (g.plane.m1 > 0)
        if (g.vis_side == 1) {
          *b = 0.5 * PI;
          *t = 1.5 * PI;
        }
        else {
          *b = 1.5 * PI;
          *t = 2.5 * PI;
        }
      else
        if (g.vis_side == -1) {
          *b = 0.5 * PI;
          *t = 1.5 * PI;
        }
        else {
          *b = 1.5 * PI;
          *t = 2.5 * PI;
        }
    else if (g.plane.m2 < 0)
      if (g.vis_side == -1) {
        double a = atan(g.plane.m1 / g.plane.m2);
        if (a < 0)
          *b = a + 2 * PI;
        else
          *b = a;
        *t = *b + PI;
      }
      else {
        *b = atan(g.plane.m1 / g.plane.m2) + PI;
        *t = *b + PI;
      }
    else
      if (g.vis_side == -1) {
        *b = atan(g.plane.m1 / g.plane.m2) + PI;
        *t = *b + PI;
      }
      else {
        double a = atan(g.plane.m1 / g.plane.m2);
        if (a < 0)
          *b = a + 2 * PI;
        else
          *b = a;
        *t = *b + PI;
      }
}


inline bool and_angle
  (double b1, double t1, double b2, double t2, double* a1, double* a2)
{
  if (b1 <= b2)
    if (t1 <= b2)
      return false;
    else {
      *a1 = t1;
      *a2 = b2;
      return true;
    }
  else
    if (b1 < t2) {
      *a1 = b1;
      *a2 = t2;
      return true;
    }
    else
      return false;
}


inline int get_point_case_xz
  (tri_group g, double a, double b, double inv, point_type* p, double angle)
{
  double cx = (p->abs.y * a + b) * inv;
  double cz = g.plane.m1 * cx + g.plane.m2 * p->abs.y + g.plane.b;

  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.z, cz))
    return 0;
  else if (equivalent_angles (atan2 (p->abs.z - cz, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xz_in
  (tri_group gi, tri_group gn, point_type* p, double angle)
{
  double cx = gi.plane.m2 * p->abs.y + gi.plane.b;
  double cz = gn.plane.m1 * cx + gn.plane.m2 * p->abs.y + gn.plane.b;

  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.z, cz))
    return 0;
  else if (equivalent_angles (atan2 (p->abs.z - cz, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_yz
  (tri_group g, double a, double b, double inv, point_type* p, double angle)
{
  double cy = (p->abs.x * a + b) * inv;
  double cz = g.plane.m1 * p->abs.x + g.plane.m2 * cy + g.plane.b;

  if (approx_equal (p->abs.y, cy) && approx_equal (p->abs.z, cz))
    return 0;
  else if (equivalent_angles (atan2 (p->abs.z - cz, p->abs.y - cy), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_yz_pn
  (tri_group gp, tri_group gn, point_type* p, double angle)
{
  double cy = gp.plane.b;
  double cz = gn.plane.m1 * p->abs.x + gn.plane.m2 * cy + gn.plane.b;

  if (approx_equal (p->abs.y, cy) && approx_equal (p->abs.z, cz))
    return 0;
  else if (equivalent_angles (atan2 (p->abs.z - cz, p->abs.y - cy), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_yz_in
  (tri_group gi, tri_group gn, point_type* p, double angle)
{
  double cy = (p->abs.x - gi.plane.b) / gi.plane.m2;
  double cz = gn.plane.m1 * p->abs.x + gn.plane.m2 * cy + gn.plane.b;

  if (approx_equal (p->abs.y, cy) && approx_equal (p->abs.z, cz))
    return 0;
  else if (equivalent_angles (atan2 (p->abs.z - cz, p->abs.y - cy), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xy_ip
  (tri_group gi, tri_group gp, point_type* p, double angle)
{
  double cy = gp.plane.b;
  double cx = gi.plane.m2 * cy + gi.plane.b;
  
  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.y, cy))
    return 0;
  else if (equivalent_angles (atan2 (cy - p->abs.y, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xy_pn
  (tri_group gp, tri_group gn, point_type* p, double angle)
{
  double cy = gp.plane.b;
  double cx = (gn.plane.m2 * cy + gn.plane.b - p->abs.z) /
    gn.plane.m1;

  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.y, cy))
    return 0;
  else if (equivalent_angles (atan2 (cy - p->abs.y, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xy_ii
  (tri_group g1, tri_group g2, point_type* p, double angle)
{
  double cy = (g2.plane.b - g1.plane.b) / (g1.plane.m2 - g2.plane.m2);
  double cx = g1.plane.m2 * cy + g1.plane.b;
  
  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.y, cy))
    return 0;
  else if (equivalent_angles (atan2 (cy - p->abs.y, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xy_in
  (tri_group gi, tri_group gn, point_type* p, double angle)
{
  double cy = (p->abs.z - gn.plane.b + gi.plane.b * gn.plane.m1) /
    (gi.plane.m2 * gn.plane.m1 + gn.plane.m2);
  double cx = gi.plane.m2 * cy + gi.plane.b;
  
  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.y, cy))
    return 0;
  else if (equivalent_angles (atan2 (cy - p->abs.y, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_point_case_xy_nn
  (tri_group g1, tri_group g2, point_type* p, double angle)
{
  double cx, cy = (g2.plane.m1 * p->abs.z - g2.plane.m1 * g1.plane.b -
    g1.plane.m1 * p->abs.z + g1.plane.m1 * g2.plane.b) /
    (g1.plane.m2 * g2.plane.m1 - g1.plane.m1 * g2.plane.m2);

  if (is_approx_zero (g1.plane.m1))
    cx = (p->abs.z - g2.plane.m2 * cy - g2.plane.b) / g2.plane.m1;
  else
    cx = (p->abs.z - g1.plane.m2 * cy - g1.plane.b) / g1.plane.m1;

  if (approx_equal (p->abs.x, cx) && approx_equal (p->abs.y, cy))
    return 0;
  else if (equivalent_angles (atan2 (cy - p->abs.y, p->abs.x - cx), angle))
    return 1;
  else
    return -1;
}


inline int get_tri_case_xz
  (tri_group g, tri_group ref, tri_type tri, double angle)
{
  double inv = 1 / (g.plane.m1 - ref.plane.m1);
  double a = ref.plane.m2 - g.plane.m2;
  double b = ref.plane.b - g.plane.b;

  int pc1 = get_point_case_xz (g, a, b, inv, tri.t3d.p1, angle);
  int pc2 = get_point_case_xz (g, a, b, inv, tri.t3d.p2, angle);
  int pc3 = get_point_case_xz (g, a, b, inv, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xz_in
  (tri_group g, tri_group ref, tri_type tri, double angle)
{
  int pc1 = get_point_case_xz_in (g, ref, tri.t3d.p1, angle);
  int pc2 = get_point_case_xz_in (g, ref, tri.t3d.p2, angle);
  int pc3 = get_point_case_xz_in (g, ref, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_yz
  (tri_group g, tri_group ref, tri_type tri, double angle)
{
  double inv = 1 / (g.plane.m2 - ref.plane.m2);
  double a = ref.plane.m1 - g.plane.m1;
  double b = ref.plane.b - g.plane.b;

  int pc1 = get_point_case_yz (g, a, b, inv, tri.t3d.p1, angle);
  int pc2 = get_point_case_yz (g, a, b, inv, tri.t3d.p2, angle);
  int pc3 = get_point_case_yz (g, a, b, inv, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_yz_pn
  (tri_group gp, tri_group gn, tri_type tri, double angle)
{
  int pc1 = get_point_case_yz_pn (gp, gn, tri.t3d.p1, angle);
  int pc2 = get_point_case_yz_pn (gp, gn, tri.t3d.p2, angle);
  int pc3 = get_point_case_yz_pn (gp, gn, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_yz_in
  (tri_group gi, tri_group gn, tri_type tri, double angle)
{
  int pc1 = get_point_case_yz_in (gi, gn, tri.t3d.p1, angle);
  int pc2 = get_point_case_yz_in (gi, gn, tri.t3d.p2, angle);
  int pc3 = get_point_case_yz_in (gi, gn, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xy_pn
  (tri_group gp, tri_group gn, tri_type tri, double angle)
{
  int pc1 = get_point_case_xy_pn (gp, gn, tri.t3d.p1, angle);
  int pc2 = get_point_case_xy_pn (gp, gn, tri.t3d.p2, angle);
  int pc3 = get_point_case_xy_pn (gp, gn, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xy_nn
  (tri_group g1, tri_group g2, tri_type tri, double angle)
{
  int pc1 = get_point_case_xy_nn (g1, g2, tri.t3d.p1, angle);
  int pc2 = get_point_case_xy_nn (g1, g2, tri.t3d.p2, angle);
  int pc3 = get_point_case_xy_nn (g1, g2, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xy_ip
  (tri_group gi, tri_group gp, tri_type tri, double angle)
{
  int pc1 = get_point_case_xy_ip (gi, gp, tri.t3d.p1, angle);
  int pc2 = get_point_case_xy_ip (gi, gp, tri.t3d.p2, angle);
  int pc3 = get_point_case_xy_ip (gi, gp, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xy_ii
  (tri_group g1, tri_group g2, tri_type tri, double angle)
{
  int pc1 = get_point_case_xy_ii (g1, g2, tri.t3d.p1, angle);
  int pc2 = get_point_case_xy_ii (g1, g2, tri.t3d.p2, angle);
  int pc3 = get_point_case_xy_ii (g1, g2, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int get_tri_case_xy_in
  (tri_group gi, tri_group gn, tri_type tri, double angle)
{
  int pc1 = get_point_case_xy_in (gi, gn, tri.t3d.p1, angle);
  int pc2 = get_point_case_xy_in (gi, gn, tri.t3d.p2, angle);
  int pc3 = get_point_case_xy_in (gi, gn, tri.t3d.p3, angle);

  return parse_point_cases (pc1, pc2, pc3);
}


inline int tri_rel_tri_xz
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back)
{
  int front_c, back_c;
  
  if (front_g.plane.y_plane) // p_
    if (back_g.plane.y_plane) // pp
      if (front_g.vis_side != back_g.vis_side)
        return CLEAR;
      else if (front_g.vis_side == -1)
        if (approx_equal (front_g.plane.b, back_g.plane.b))
          return CLEAR;
        else if (front_g.plane.b < back_g.plane.b)
          return INFRONT;
        else
          return BEHIND;
      else
        if (approx_equal (front_g.plane.b, back_g.plane.b))
          return CLEAR;
        else if (front_g.plane.b > back_g.plane.b)
          return INFRONT;
        else
          return BEHIND;
    else
      return DONTKNOW;
  else if (front_g.plane.m1_inf) // i_
    if (back_g.plane.y_plane) // ip
      return DONTKNOW;
    else if (back_g.plane.m1_inf) // ii
      return DONTKNOW;
    else { // in
      double front_a, back_a;
//      get_common_range_xz (front_g, back_g, &front_a, &back_a);
      if (!get_common_range_xz (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_xz_in (front_g, back_g, front, front_a);
      back_c = get_tri_case_xz_in (front_g, back_g, back, back_a);
    }
  else // n_
    if (back_g.plane.y_plane) // np
      return DONTKNOW;
    else if (back_g.plane.m1_inf) { // ni
      double front_a, back_a;
//      get_common_range_xz (front_g, back_g, &front_a, &back_a);
      if (!get_common_range_xz (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_xz_in (back_g, front_g, front, front_a);
      back_c = get_tri_case_xz_in (back_g, front_g, back, back_a);
    }
    else // nn
      if (approx_equal (front_g.plane.m1, back_g.plane.m1))
        return PARALELL;
      else {
        double front_a, back_a;
//        get_common_range_xz (front_g, back_g, &front_a, &back_a);
        if (!get_common_range_xz (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xz (front_g, back_g, front, front_a);
        back_c = get_tri_case_xz (back_g, front_g, back, back_a);
      }

  if ((front_c == 0) && (back_c == 0))
    return INTERSECT;
  else if (front_c < back_c)
    return BEHIND;
  else if (front_c == back_c)
    return CLEAR;
  else
    return INFRONT;
}


inline int tri_rel_tri_yz
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back)
{
  int front_c, back_c;
  
  if (front_g.plane.y_plane) // p_
    if (back_g.plane.y_plane) // pp
      return DONTKNOW;
    else if (back_g.plane.m1_inf) // pi
      return DONTKNOW;    
    else { // pn
      double front_a, back_a;
//      get_common_range_yz (front_g, back_g, &front_a, &back_a);
      if (!get_common_range_yz (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_yz_pn (front_g, back_g, front, front_a);
      back_c = get_tri_case_yz_pn (front_g, back_g, back, back_a);
    }
  else if (front_g.plane.m1_inf)
    if (back_g.plane.y_plane) // ip
      return DONTKNOW;        
    else if (back_g.plane.m1_inf) // ii
      return DONTKNOW;        
    else
      if (is_approx_zero (front_g.plane.m2))
        return DONTKNOW;
      else {
        double front_a, back_a;
//        get_common_range_yz (front_g, back_g, &front_a, &back_a);
        if (!get_common_range_yz (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_yz_in (front_g, back_g, front, front_a);
        back_c = get_tri_case_yz_in (front_g, back_g, back, back_a);
      }
  else
    if (back_g.plane.y_plane) { // np
      double front_a, back_a;
//      get_common_range_yz (front_g, back_g, &front_a, &back_a);
      if (!get_common_range_yz (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_yz_pn (back_g, front_g, front, front_a);
      back_c = get_tri_case_yz_pn (back_g, front_g, back, back_a);
    }
    else if (back_g.plane.m1_inf) // ni
      if (is_approx_zero (front_g.plane.m2))
        return DONTKNOW;
      else {
        double front_a, back_a;
//        get_common_range_yz (front_g, back_g, &front_a, &back_a);
        if (!get_common_range_yz (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_yz_in (back_g, front_g, front, front_a);
        back_c = get_tri_case_yz_in (back_g, front_g, back, back_a);
      }
    else // nn
      if (approx_equal (front_g.plane.m2, back_g.plane.m2))
        return PARALELL;
      else {
        double front_a, back_a;
//        get_common_range_yz (front_g, back_g, &front_a, &back_a);
        if (!get_common_range_yz (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_yz (front_g, back_g, front, front_a);
        back_c = get_tri_case_yz (back_g, front_g, back, back_a);
      }

  if ((front_c == 0) && (back_c == 0))
    return INTERSECT;
  else if (front_c < back_c)
    return BEHIND;
  else if (front_c == back_c)
    return CLEAR;
  else
    return INFRONT;
}


inline int tri_rel_tri_xy
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back)
{
  int front_c, back_c;
  
  if (front_g.plane.y_plane) // p_
    if (back_g.plane.y_plane) // pp
      return DONTKNOW;
    else if (back_g.plane.m1_inf) { // pi
      double front_a, back_a;
      if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_xy_ip (back_g, front_g, front, front_a);
      back_c = get_tri_case_xy_ip (back_g, front_g, back, back_a);
    }
    else // pn
      if (is_approx_zero (back_g.plane.m1))
        return PARALELL;
      else {
        double front_a, back_a;
        if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xy_pn (front_g, back_g, front, front_a);
        back_c = get_tri_case_xy_pn (front_g, back_g, back, back_a);
      }
  else if (front_g.plane.m1_inf)
    if (back_g.plane.y_plane) { // ip
      double front_a, back_a;
      if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      front_c = get_tri_case_xy_ip (front_g, back_g, front, front_a);
      back_c = get_tri_case_xy_ip (front_g, back_g, back, back_a);
    }
    else if (back_g.plane.m1_inf) { // ii
      double front_a, back_a;
      if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
        return CLEAR;
      if (approx_equal (front_g.plane.m2, back_g.plane.m2))
        return PARALELL;
      else {
        front_c = get_tri_case_xy_ii (front_g, back_g, front, front_a);
        back_c = get_tri_case_xy_ii (front_g, back_g, back, back_a);
      }
    }
    else // in
      if (is_approx_zero (front_g.plane.m1 * back_g.plane.m1 + back_g.plane.m2))
        return PARALELL;
      else {
        double front_a = 0, back_a = 0;
        if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xy_in (front_g, back_g, front, front_a);
        back_c = get_tri_case_xy_in (front_g, back_g, back, back_a);
      }
  else
    if (back_g.plane.y_plane) // np
      if (is_approx_zero (back_g.plane.m1))
        return PARALELL;
      else {
        double front_a, back_a;
        if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xy_pn (back_g, front_g, front, front_a);
        back_c = get_tri_case_xy_pn (back_g, front_g, back, back_a);
      }
    else if (back_g.plane.m1_inf) // ni
      if (is_approx_zero (front_g.plane.m1 * back_g.plane.m1 + back_g.plane.m2))
        return PARALELL;
      else {
        double front_a, back_a;
        if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xy_in (back_g, front_g, front, front_a);
        back_c = get_tri_case_xy_in (back_g, front_g, back, back_a);
      }
    else // nn
      if (approx_equal (front_g.plane.m2 * back_g.plane.m1, front_g.plane.m1 * back_g.plane.m2))
        return PARALELL;
      else {
        double front_a, back_a;
        if (!get_common_range_xy (front_g, back_g, &front_a, &back_a))
          return CLEAR;
        front_c = get_tri_case_xy_nn (front_g, back_g, front, front_a);
        back_c = get_tri_case_xy_nn (front_g, back_g, back, back_a);
      }

  if ((front_c == 0) && (back_c == 0))
    return INTERSECT;
  else if (front_c < back_c)
    return BEHIND;
  else if (front_c == back_c)
    return CLEAR;
  else
    return INFRONT;
}


inline int parse_point_cases (int pc1, int pc2, int pc3)
{
  if (pc1 == -1)
    if (pc2 == -1)
      if (pc3 == -1)
        return -1; // -1, -1, -1
      else if (pc3 == 0)
        return -1; // -1, -1, 0
      else
        return 0; // -1, -1, 1
    else if (pc2 == 0)
      if (pc3 == -1)
        return -1; // -1, 0, -1
      else if (pc3 == 0)
        return -1; // -1, 0, 0
      else
        return -1; // -1, 0, 1
    else
      if (pc3 == -1)
        return 0; // -1, 1, -1
      else if (pc3 == 0)
        return 0; // -1, 1, 0
      else
        return 0; // -1, 1, 1
  else if (pc1 == 0)
    if (pc2 == -1)
      if (pc3 == -1)
        return -1; // 0, -1, -1
      else if (pc3 == 0)
        return -1; // 0, -1, 0
      else
        return 0; // 0, -1, 1
    else if (pc2 == 0)
      if (pc3 == -1)
        return -1; // 0, 0, -1
      else if (pc3 == 0)
        return 0; // 0, 0, 0
      else
        return 1; // 0, 0, 1
    else
      if (pc3 == -1)
        return 0; // 0, 1, -1
      else if (pc3 == 0)
        return 1; // 0, 1, 0
      else
        return 1; // 0, 1, 1
  else
    if (pc2 == -1)
      if (pc3 == -1)
        return 0; // 1, -1, -1
      else if (pc3 == 0)
        return 0; // 1, -1, 0
      else
        return 0; // 1, -1, 1
    else if (pc2 == 0)
      if (pc3 == -1)
        return -1; // 1, 0, -1
      else if (pc3 == 0)
        return 1; // 1, 0, 0
      else
        return 1; // 1, 0, 1
    else
      if (pc3 == -1)
        return 0; // 1, 1, -1
      else if (pc3 == 0)
        return 1; // 1, 1, 0
      else
        return 1; // 1, 1, 1
}


inline int tri_rel_tri
  (tri_group front_g, tri_group back_g, tri_type front, tri_type back)
{
  int c;

  c = tri_rel_tri_xz (front_g, back_g, front, back);
  if (c >= -1)
    return c;

  c = tri_rel_tri_yz (front_g, back_g, front, back);
  if (c >= -1)
    return c;

  c = tri_rel_tri_xy (front_g, back_g, front, back);
  if (c == PARALELL)
    if (front_g.vis_side == back_g.vis_side)
      if (approx_equal (front_g.plane.b, back_g.plane.b))
        return CLEAR;
      else if (front_g.plane.b > back_g.plane.b)
        if (front_g.vis_side == -1)
          return BEHIND;
        else
          return INFRONT;
      else
        if (front_g.vis_side == -1)
          return INFRONT;        
        else
          return BEHIND;        
    else
      return CLEAR;
  else
    return c;
}


void show_case (int c, int l)
{
  if (c == INFRONT)
    bprint ("In front", l);
  else if (c == BEHIND)
    bprint ("Behind", l);
  else if (c == INTERSECT)
    bprint ("Intersect", l);
  else if (c == CLEAR)
    bprint ("Clear", l);
  else if (c == PARALELL)
    bprint ("Paralell", l);
  else if (c == DONTKNOW)
    bprint ("Dont know", l);
  else
    bprint ("Error", l);
}


