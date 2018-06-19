#ifndef INCLUDE_COLIS
#define INCLUDE_COLOS

#include "gfx3d.h"
#include "math3d.h"
#include "texclass.h"

point_type _test_p[2];

struct path_type {
  point_3d strt, stp;
};

inline bool colis_path_gen (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool check_colis_group (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* closest,
  point_3d curr_top, point_3d tent_top);
inline bool exclude_colis_cube (point_3d min, point_3d max,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_group_range_flat
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_group_range
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_path_flat (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_path_inf (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_tri_range
  (tri_type* group, double clear,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_tri_range_flat
  (tri_type* group, double clear,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool exclude_colis_perim_x
  (tri_group group, cylinder_3d cyl, point_3d curr_top, point_3d tent_top);
inline bool exclude_colis_perim_xy
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d tent_top, point_3d curr_bot, point_3d tent_bot);
inline bool path_int_plane_x (tri_group group, cylinder_3d cyl,
  point_3d* result,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool path_int_plane_xy (tri_group group, cylinder_3d cyl,
  point_3d* result_a, point_3d* result_b,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool path_int_plane_xyz (tri_group group, cylinder_3d cyl,
  point_3d* int_a, point_3d* int_b,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline point_3d path_pen_xy
  (plane_type plane, point_3d pa, point_3d pb, double offset);
inline point_3d line_pen_xy
  (plane_type plane, point_3d p1, double offset);
inline int sort_intersections_xy_neg
  (plane_type plane, point_3d t, point_3d b, point_3d* first, point_3d* last);
inline int sort_intersections_xy_pos
  (plane_type plane, point_3d t, point_3d b, point_3d* first, point_3d* last);
inline point_3d calc_new_loc
  (cylinder_3d cyl, point_3d pa, point_3d pb, point_3d ntersect);
inline void sort_intersections_xyz
  (plane_type plane, point_3d curr_t, point_3d curr_b, point_3d t, point_3d b,
  point_3d* first, point_3d* last);
inline point_3d line_pen_xyz
  (plane_type plane, point_3d p1, double offset);
inline point_3d path_pen_xyz
  (plane_type plane, point_3d pa, point_3d pb, double offset);
inline bool colis_xz_element_up (point_3d* result, point_3d* closest,
  point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r);
inline bool colis_xz_element_dn (point_3d* result, point_3d* closest,
  point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r);
inline bool colis_xz_tri_up (point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r);
inline bool colis_xz_tri_dn (point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r);
inline point_3d int_edge_and_z (point_3d mid, point_3d p1, point_3d p2);
inline bool colis_xz_plane
  (point_3d* result, point_3d* closest, double offset,
  point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3);
inline bool colis_xz_tris
  (tri_group group, point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb);
inline bool colis_xz_tris_rev
  (tri_group group, point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb);
inline point_3d int_yz_edge_and_y (point_3d p1, point_3d mid, point_3d p2);
inline bool colis_yz_tri_up (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
  point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_tri_dn (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
  point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_tri
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_group
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_group_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_element_up (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r);
inline bool colis_yz_element_dn (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r);
inline bool path_int_xz_horz
  (point_3d pa, point_3d pb,
  point_3d l, point_3d r, point_3d* result);
inline bool path_int_yz_horz
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d l, point_3d r, bool* slide, plane_type* res_plane, point_3d* result);
inline bool path_int_xz_edge
  (point_3d pa, point_3d pb,
   point_3d p1, point_3d p2, point_3d* result);
inline bool path_int_yz_edge
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d p1, point_3d p2, bool* slide, plane_type* res_plane, point_3d* result);
inline void store_closer (point_3d* closest, point_3d result, point_3d orig,
  plane_type plane, bool* slide, plane_type* res_plane);
inline point_3d path_pen_xz (double y, point_3d pa, point_3d pb);
inline point_3d line_pen_xz (double y, point_3d p);
inline bool path_int_plane_xz_pos (double y, bool* slide, point_3d* result,
  point_3d curr_top, point_3d tent_top);
inline bool path_int_plane_xz_neg (double y, bool* slide, point_3d* result,
  point_3d curr_bot, point_3d tent_bot);
inline bool colis_xz_quad_neg
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d cc, point_3d cw, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_xz_quad_neg_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d cc, point_3d cw, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_xz_quad_pos
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_xz_quad_pos_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool xz_quad_up_top
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d t, point_3d l, point_3d r, double offset,
  tri_group group, point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool xz_quad_up_bot
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d t, point_3d l, point_3d r, double offset,
  tri_group group, point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool xz_quad_dn_top
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d b, point_3d l, point_3d r, double offset,
  tri_group group, point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool xz_quad_dn_bot
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d b, point_3d l, point_3d r, double offset,
  tri_group group, point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_tri_tb
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   tri_group group, double offset,
   point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_group_x
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_yz_group_x_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_element_up (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r);
inline bool colis_xy_element_dn (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r);
inline bool path_int_xy_horz
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d l, point_3d r, bool* slide, plane_type* res_plane, point_3d* result);
inline bool path_int_xy_edge
  (point_3d pa, point_3d pb, plane_type plane, double offset,
   point_3d p1, point_3d p2, bool* slide, plane_type* res_plane, point_3d* result);
inline bool colis_xy_group
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_group_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline point_3d int_xy_edge_and_y (point_3d p1, point_3d mid, point_3d p2);
inline bool colis_xy_tri_up
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_tri_dn
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_tri
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_tri_tb
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   tri_group group, double offset,
   point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_group_flat
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline bool colis_xy_group_flat_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot);
inline plane_type edge_plane (point_3d p1, point_3d p2);
point_3d deflect_path
  (plane_type plane, cylinder_3d cyl, point_3d ntersect, point_3d curr, point_3d tent);
inline point_3d deflect_xyz
  (plane_type plane, point_3d ref, point_3d p);
inline point_3d deflect_xy
  (plane_type plane, point_3d ref, point_3d p);
void move_player (view_type view, camera_type* cursor, cylinder_3d cyl,
  point_type** first_point, tri_group** first_group);
bool check_colis_cyl (tri_group* first_group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* closest, point_3d curr, point_3d tent);
int coplanar (plane_type plane, point_3d p, point_3d ref);
void disp_plane (plane_type plane);


bool check_colis_cyl (tri_group* first_group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* closest, point_3d curr, point_3d tent)
{
  if (same_point (tent, curr))
    return false;
  else {
    bool hit = false;
    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (check_colis_group (*g, cyl, slide, res_plane, closest, curr, tent))
        hit = true;
    return hit;
  }
}


inline bool check_colis_group (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* closest,
  point_3d curr_top, point_3d tent_top)
{
  point_3d result;
  point_3d curr_bot = offset_point_3d (curr_top, 0, cyl.height, 0);
  point_3d tent_bot = offset_point_3d (tent_top, 0, cyl.height, 0);
  
//  _spec[0].abs = curr_top;
//  _spec[1].abs = curr_bot;
//  _spec[0].select = true;

  if (group.plane.y_plane)
    return colis_path_flat (group, cyl, slide, res_plane, &result, closest, curr_top, curr_bot, tent_top, tent_bot);
  else if (group.plane.m1_inf)
    return colis_path_inf (group, cyl, slide, res_plane, &result, closest, curr_top, curr_bot, tent_top, tent_bot);
  else
    return colis_path_gen (group, cyl, slide, res_plane, &result, closest, curr_top, curr_bot, tent_top, tent_bot);
}


inline bool colis_path_flat (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d res, pa, pb;
  
  if (group.vis_side == -1) {
    if (path_int_plane_xz_neg (group.plane.b, slide, &pa, curr_bot, tent_bot))
      if (colis_group_range_flat (group, cyl, curr_top, curr_bot, tent_top, tent_bot))
        if (colis_xz_tris (group, result, &res, cyl.radius, pa, pb)) {
          res = calc_new_loc (cyl, curr_top, tent_top, res);
          if (dist_3d (res, curr_top) < dist_3d (*closest, curr_top)) {
            *result = *closest = res;
            *res_plane = group.plane;
            return true;
          }
          else
            return false;
        }
  }
  else {
    if (path_int_plane_xz_pos (group.plane.b, slide, &pa, curr_top, tent_top))
      if (colis_group_range_flat (group, cyl, curr_top, curr_bot, tent_top, tent_bot))
        if (colis_xz_tris_rev (group, result, &res, cyl.radius, pa, pb)) {
          res = calc_new_loc (cyl, curr_top, tent_top, res);
          if (dist_3d (res, curr_top) < dist_3d (*closest, curr_top)) {
            *result = *closest = res;
            *res_plane = group.plane;            
            return true;
          }
          else
            return false;
        }
  }

  return false;
/*
  else {
    if (approx_greater_inc (curr_top.y, group.plane.b))
      if (approx_lesser_inc (tent_top.y, group.plane.b))
        if (approx_equal (curr_top.y, tent_top.y))
          return false;
        else {
          if (path_int_plane_xz_neg (group.plane.b, &pa, &pb, curr_top, curr_bot, tent_top, tent_bot))
            if (colis_group_range_flat (group, cyl, curr_top, curr_bot, tent_top, tent_bot)) {
              y_line_3d path = calc_y_line_3d (curr_top, tent_top);
              result->x = path.mx * group.plane.b + path.bx;
              result->y = group.plane.b;
              result->z = path.mz * group.plane.b + path.bz;
              return colis_xz_tris (group, result, cyl.radius, pa, pb);
            }
        }
  }
*/
}


inline bool colis_path_inf (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  if (approx_zero (group.plane.m2)) {
    if (exclude_colis_perim_x (group, cyl, curr_top, tent_top))
      if (colis_group_range (group, cyl, curr_top, curr_bot, tent_top, tent_bot)) {
//        if (path_int_plane_x (group, cyl, result, curr_top, curr_bot, tent_top, tent_bot)) {
          point_3d result_bot = offset_point_3d (*result, 0, cyl.height - 1, 0);
          if (group.vis_side == -1) {
            if (colis_yz_group_x_rev (slide, res_plane, result, closest, group, cyl, *result, result_bot, curr_top, curr_bot, tent_top, tent_bot)) {
              *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
              return true;
            }
          }
          else {
            if (colis_yz_group_x (slide, res_plane, result, closest, group, cyl, *result, result_bot, curr_top, curr_bot, tent_top, tent_bot)) {
              *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
              return true;
            }
          }
      }
//        if (path_int_plane_x (group, cyl, result, curr_top, curr_bot, tent_top, tent_bot)) {
//          *closest = calc_new_loc (cyl, curr_top, tent_top, *result);
//          return true;
//        }
  }
  else {
    point_3d result_a, result_b;
//    if (exclude_colis_perim_xy (group, cyl, curr_top, curr_bot, tent_top, tent_bot)) {
//      if (colis_group_range (group, cyl, curr_top, curr_bot, tent_top, tent_bot))
          if (group.vis_side == -1) {
            if (colis_yz_group_rev (slide, res_plane, result, closest, group, cyl, result_a, result_b, curr_top, curr_bot, tent_top, tent_bot)) {
              *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
              return true;
            }
          }
          else {
            if (colis_yz_group (slide, res_plane, result, closest, group, cyl, result_a, result_b, curr_top, curr_bot, tent_top, tent_bot)) {
              *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
              return true;
            }
          }
//    }
  }

  return false;
}


inline bool path_int_plane_x (tri_group group, cylinder_3d cyl,
  point_3d* result,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  if (approx_equal (curr_top.x, tent_top.x))
    return false;

  if (group.vis_side == -1) {
    if (approx_lesser_inc (curr_top.x, group.plane.b - cyl.radius))
      if (approx_greater_inc (tent_top.x, group.plane.b - cyl.radius)) {
        x_line_3d path = calc_x_line_3d (curr_top, tent_top);
        result->x = group.plane.b - cyl.radius;
        result->y = path.my * result->x + path.by;
        result->z = path.mz * result->x + path.bz;
        return true;
      }
  }
  else {
    if (approx_greater_inc (curr_top.x, group.plane.b + cyl.radius))
      if (approx_lesser_inc (tent_top.x, group.plane.b + cyl.radius)) {
        x_line_3d path = calc_x_line_3d (curr_top, tent_top);
        result->x = group.plane.b + cyl.radius;
        result->y = path.my * result->x + path.by;
        result->z = path.mz * result->x + path.bz;
        return true;
      }
  }

  return false;
}


inline bool path_int_plane_xz_neg (double y, bool* slide, point_3d* result,
  point_3d curr_bot, point_3d tent_bot)
{
  *slide = false;
  
  if (approx_greater (curr_bot.y, y))
    return false;
  else if (approx_lesser (curr_bot.y, y))
    if (approx_greater_inc (tent_bot.y, y)) {
      *result = path_pen_xz (y, curr_bot, tent_bot);
      return true;
    }
    else
      return false;
  else
    if (approx_greater (tent_bot.y, y)) {
      *result = path_pen_xz (y, curr_bot, tent_bot);
      return true;
    }
    else if (approx_lesser (tent_bot.y, y))
      return false;
    else {
      *slide = true;
      return false;
    }

/*
  if (approx_lesser_inc (curr_bot.y, y))
    if (approx_lesser_inc (tent_bot.y, y))
      return false;
    else if (approx_lesser (tent_top.y, y)) {
      *pa = path_pen_xz (y, curr_bot, tent_bot);
      *pb = line_pen_xz (y, tent_top);
    }
    else {
      *pa = path_pen_xz (y, curr_bot, tent_bot);
      *pb = path_pen_xz (y, curr_top, tent_top);
    }
  else if (approx_lesser (curr_top.y, y))
    if (approx_lesser_inc (tent_bot.y, y)) {
      *pa = line_pen_xz (y, curr_top);
      *pb = path_pen_xz (y, curr_bot, tent_bot);
    }
    else if (approx_lesser (tent_top.y, y)) {
      *pa = line_pen_xz (y, curr_top);
      *pb = line_pen_xz (y, tent_top);
    }
    else {
      *pa = line_pen_xz (y, curr_top);
      *pb = path_pen_xz (y, curr_top, tent_top);
    }
  else
    return false;
  */
}


inline bool path_int_plane_xz_pos (double y, bool* slide, point_3d* result,
  point_3d curr_top, point_3d tent_top)
{
  *slide = false;
  
  if (approx_lesser (curr_top.y, y))
    return false;
  else if (approx_greater (curr_top.y, y))
    if (approx_lesser_inc (tent_top.y, y)) {
      *result = path_pen_xz (y, curr_top, tent_top);
      return true;
    }
    else
      return false;
  else
    if (approx_lesser (tent_top.y, y)) {
      *result = path_pen_xz (y, curr_top, tent_top);
      return true;
    }
    else if (approx_greater (tent_top.y, y))
      return false;
    else {
      *slide = true;
      return false;
    }
}


inline bool path_int_plane_xy (tri_group group, cylinder_3d cyl,
  point_3d* int_a, point_3d* int_b,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  if (group.vis_side == -1) {
    if (approx_lesser_inc (curr_top.x, group.plane.m2 * curr_top.y + group.plane.b - cyl.radius))
      if (approx_lesser_inc (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b - cyl.radius))
        if (approx_lesser_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b - cyl.radius))
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) // 1
            return false;
          else { // 2
            *int_a = path_pen_xy (group.plane, curr_bot, tent_bot, -cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, -cyl.radius);
            bprint ("2", 7);
            return true;
          }
        else
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 3
            *int_a = path_pen_xy (group.plane, curr_top, tent_top, -cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, -cyl.radius);
            bprint ("3", 7);
            return true;
          }
          else { // 4
            bprint ("4", 7);
            sort_intersections_xy_neg (group.plane,
              path_pen_xy (group.plane, curr_top, tent_top, -cyl.radius),
              path_pen_xy (group.plane, curr_bot, tent_bot, -cyl.radius), int_a, int_b);
            return true;
          }
      else
        if (approx_lesser_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b - cyl.radius))
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 5
            return false;
//            *int_a = line_pen_xy (group.plane, curr_top);
//            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot);
          }
          else { // 6
            *int_a = line_pen_xy (group.plane, curr_top, -cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, -cyl.radius);
            bprint ("6a", 7);
            return true;
          }
        else
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) // 7
            return false; // N/A
          else { // 8
            bprint ("8a", 7);
            *int_a = line_pen_xy (group.plane, curr_top, -cyl.radius);
            *int_b = path_pen_xy (group.plane, curr_top, tent_top, -cyl.radius);
            return true;
          }
    else
      if (approx_lesser_inc (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b - cyl.radius))
        if (approx_lesser_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b - cyl.radius))
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 9
            return false;
//            *int_a = line_pen_xy (group.plane, curr_top);
//            *int_b = path_pen_xy (group.plane, curr_top, tent_top);
          }
          else // 10
            false; // N/A
        else
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 11
            bprint ("11", 7);
            *int_a = line_pen_xy (group.plane, curr_top, -cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, -cyl.radius);
            return true;
          }
          else { // 12
            bprint ("12", 7);
            *int_a = line_pen_xy (group.plane, curr_top, -cyl.radius);
            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot, -cyl.radius);
            return true;
          }
      else
        if (approx_lesser_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b - cyl.radius))
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 13
            return false;
//            *int_a = path_pen_xy (group.plane, curr_top, tent_top);
//            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot);
          }
          else { // 14
            return false;
//            *int_a = path_pen_xy (group.plane, curr_top, tent_top);
//            *int_b = line_pen_xy (group.plane, tent_top);
          }
        else
          if (approx_lesser_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius)) { // 15
            return false;
//            *int_a = path_pen_xy (group.plane, curr_bot, tent_bot);
//            *int_b = line_pen_xy (group.plane, tent_top);
          }
          else // 16
            return false;
  }
  else {
    if (approx_greater_inc (curr_top.x, group.plane.m2 * curr_top.y + group.plane.b + cyl.radius))
      if (approx_greater_inc (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b + cyl.radius))
        if (approx_greater_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b + cyl.radius))
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) // 1
            return false;
          else { // 2
            *int_a = path_pen_xy (group.plane, curr_bot, tent_bot, cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, cyl.radius);
            bprint ("2", 7);
            return true;
          }
        else
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 3
            *int_a = path_pen_xy (group.plane, curr_top, tent_top, cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, cyl.radius);
            bprint ("3", 7);
            return true;
          }
          else { // 4
            bprint ("4", 7);
            sort_intersections_xy_pos (group.plane,
              path_pen_xy (group.plane, curr_top, tent_top, cyl.radius),
              path_pen_xy (group.plane, curr_bot, tent_bot, cyl.radius), int_a, int_b);
            return true;
          }
      else
        if (approx_greater_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b + cyl.radius))
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 5
            return false;
//            *int_a = line_pen_xy (group.plane, curr_top);
//            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot);
          }
          else { // 6
            *int_a = line_pen_xy (group.plane, curr_top, cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, cyl.radius);
            bprint ("6", 7);
            return true;
          }
        else
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) // 7
            return false; // N/A
          else { // 8
            bprint ("8", 7);
            *int_a = line_pen_xy (group.plane, curr_top, cyl.radius);
            *int_b = path_pen_xy (group.plane, curr_top, tent_top, cyl.radius);
            return true;
          }
    else
      if (approx_greater_inc (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b + cyl.radius))
        if (approx_greater_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b + cyl.radius))
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 9
            return false;
//            *int_a = line_pen_xy (group.plane, curr_top);
//            *int_b = path_pen_xy (group.plane, curr_top, tent_top);
          }
          else // 10
            false; // N/A
        else
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 11
            bprint ("11", 7);
            *int_a = line_pen_xy (group.plane, curr_top, cyl.radius);
            *int_b = line_pen_xy (group.plane, tent_top, cyl.radius);
            return true;
          }
          else { // 12
            bprint ("12", 7);
            *int_a = line_pen_xy (group.plane, curr_top, cyl.radius);
            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot, cyl.radius);
            return true;
          }
      else
        if (approx_greater_inc (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b + cyl.radius))
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 13
            return false;
//            *int_a = path_pen_xy (group.plane, curr_top, tent_top);
//            *int_b = path_pen_xy (group.plane, curr_bot, tent_bot);
          }
          else { // 14
            return false;
//            *int_a = path_pen_xy (group.plane, curr_top, tent_top);
//            *int_b = line_pen_xy (group.plane, tent_top);
          }
        else
          if (approx_greater_inc (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius)) { // 15
            return false;
//            *int_a = path_pen_xy (group.plane, curr_bot, tent_bot);
//            *int_b = line_pen_xy (group.plane, tent_top);
          }
          else // 16
            return false;
  }
}


inline point_3d path_pen_xz (double y, point_3d pa, point_3d pb)
{
  point_3d r;
  
  y_line_3d path = calc_y_line_3d (pa, pb);
  r.x = path.mx * y + path.bx;
  r.y = y;
  r.z = path.mz * y + path.bz;
  
  return r;
}


inline point_3d line_pen_xz (double y, point_3d p)
{
  point_3d r;

  r.x = p.x;
  r.y = y;
  r.z = p.z;
  
  return r;
}


inline point_3d path_pen_xy
  (plane_type plane, point_3d pa, point_3d pb, double offset)
{
  point_3d p;

  if (approx_equal (pa.y, pb.y)) {
    lin_relat l = calc_dzdx_line (pa, pb);
    p.x = plane.m2 * pa.y + plane.b + offset;
    p.y = pa.y;
    p.z = l.m * p.x + l.b;
  }
  else {
    y_line_3d l = calc_y_line_3d (pa, pb);
    p.y = ((plane.b + offset) - l.bx) / (l.mx - plane.m2);
    p.x = plane.m2 * p.y + plane.b + offset;
    p.z = l.mz * p.y + l.bz;
  }

  return p;
}


inline point_3d line_pen_xy
  (plane_type plane, point_3d p1, double offset)
{
  point_3d p;

  p.x = p1.x;
  p.y = (p1.x - (plane.b + offset)) / plane.m2;
  p.z = p1.z;

  return p;
}


inline int sort_intersections_xy_neg
  (plane_type plane, point_3d t, point_3d b, point_3d* first, point_3d* last)
{
  if (plane.m2 > 0) {
    *first = t;
    *last = b;
    return -1;
  }
  else {
    *first = b;
    *last = t;
    return 1;
  }
}


inline int sort_intersections_xy_pos
  (plane_type plane, point_3d t, point_3d b, point_3d* first, point_3d* last)
{
  if (plane.m2 < 0) {
    *first = t;
    *last = b;
    return -1;
  }
  else {
    *first = b;
    *last = t;
    return 1;
  }
}


inline point_3d path_pen_xyz
  (plane_type plane, point_3d pa, point_3d pb, double offset)
{
  point_3d p;

  if (approx_equal (pa.z, pb.z)) {
    lin_relat l = calc_dydx_line (pa, pb);
    p.x = (pa.z - (plane.b + offset) - plane.m2 * l.b) / (plane.m1 + l.m * plane.m2);
    p.y = l.m * p.x + l.b;
    p.z = pa.z;
  }
  else {
    line_3d l = calc_line_3d (pa, pb);
    p.z = (plane.m1 * l.bx + plane.m2 * l.by + (plane.b + offset)) /
      (1 - plane.m1 * l.mx - plane.m2 * l.my);
    p.x = l.mx * p.z + l.bx;
    p.y = l.my * p.z + l.by;
  }

  return p;
}


inline point_3d line_pen_xyz
  (plane_type plane, point_3d p1, double offset)
{
  point_3d p;

  p.x = p1.x;
  p.y = (p1.z - plane.m1 * p1.x - (plane.b + offset)) / plane.m2;
  p.z = p1.z;

  return p;
}


inline void sort_intersections_xyz
  (plane_type plane, point_3d curr_t, point_3d curr_b, point_3d t, point_3d b,
  point_3d* first, point_3d* last)
{
  double d1 = dist_3d (curr_t, t);
  double d2 = dist_3d (curr_b, b);

  if (d1 < d2) {
    *first = t;
    *last = b;
  }
  else {
    *first = b;
    *last = t;
  }
}


inline bool path_int_plane_xz (tri_group group, cylinder_3d cyl,
  point_3d* result,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  if (approx_equal (curr_top.z, group.plane.m1 * curr_top.x + group.plane.b - group.plane.sec_m1 * cyl.radius))
    if (approx_equal (tent_top.z, group.plane.m1 * tent_top.x + group.plane.b - group.plane.sec_m1 * cyl.radius))
      return false;
      
  if (group.vis_side == -1) {
    if (approx_lesser_inc (curr_top.z, group.plane.m1 * curr_top.x + group.plane.b - group.plane.sec_m1 * cyl.radius))
      if (approx_greater_inc (tent_top.z, group.plane.m1 * tent_top.x + group.plane.b - group.plane.sec_m1 * cyl.radius))
        if (approx_equal (curr_top.x, tent_top.x)) {
          if (approx_equal (curr_top.z, tent_top.z))
            return false;
          else {
            lin_relat path = calc_dydz_line (curr_top, tent_top);
            result->z = group.plane.m1 * curr_top.x + group.plane.b - group.plane.sec_m1 * cyl.radius;
            result->y = path.m * result->z + path.b;
            result->x = curr_top.x;
            return true;
          }
        }
        else {
          x_line_3d path = calc_x_line_3d (curr_top, tent_top);
          result->x = (group.plane.b - group.plane.sec_m1 * cyl.radius - path.bz) / (path.mz - group.plane.m1);
          result->y = path.my * result->x + path.by;
          result->z = path.mz * result->x + path.bz;
          return true;
        }
  }
  else {
    if (approx_greater_inc (curr_top.z, group.plane.m1 * curr_top.x + group.plane.b + group.plane.sec_m1 * cyl.radius))
      if (approx_lesser_inc (tent_top.z, group.plane.m1 * tent_top.x + group.plane.b + group.plane.sec_m1 * cyl.radius)) {
        if (approx_equal (curr_top.x, tent_top.x)) {
          if (approx_equal (curr_top.z, tent_top.z))
            return false;
          else {
            lin_relat path = calc_dydz_line (curr_top, tent_top);
            result->z = group.plane.m1 * curr_top.x + group.plane.b + group.plane.sec_m1 * cyl.radius ;
            result->y = path.m * result->z + path.b;
            result->x = curr_top.x;
            return true;
          }
        }
        else {
          x_line_3d path = calc_x_line_3d (curr_top, tent_top);
          result->x = (group.plane.b + group.plane.sec_m1 * cyl.radius - path.bz) / (path.mz - group.plane.m1);
          result->y = path.my * result->x + path.by;
          result->z = path.mz * result->x + path.bz;
          return true;
        }
      }
  }

  return false;
}


inline bool exclude_colis_perim_x
  (tri_group group, cylinder_3d cyl, point_3d curr_top, point_3d tent_top)
{
  if (approx_lesser (curr_top.x, group.plane.b - cyl.radius))
    if (approx_lesser (tent_top.x, group.plane.b - cyl.radius))
      return false;
  if (approx_greater (curr_top.x, group.plane.b + cyl.radius))
    if (approx_greater (tent_top.x, group.plane.b + cyl.radius))
      return false;

  return true;
}


inline bool exclude_colis_perim_xy
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d curr_bot, point_3d tent_top, point_3d tent_bot)
{
  if (group.plane.m2 > 0) {
    if (approx_lesser (curr_top.x, group.plane.m2 * curr_top.y + group.plane.b - cyl.radius))
      if (approx_lesser (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b - cyl.radius))
        return false;
    if (approx_greater (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b + cyl.radius))
      if (approx_greater (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b + cyl.radius))
        return false;
  }
  else {
    if (approx_lesser (curr_bot.x, group.plane.m2 * curr_bot.y + group.plane.b - cyl.radius))
      if (approx_lesser (tent_bot.x, group.plane.m2 * tent_bot.y + group.plane.b - cyl.radius))
        return false;
    if (approx_greater (curr_top.x, group.plane.m2 * curr_top.y + group.plane.b + cyl.radius))
      if (approx_greater (tent_top.x, group.plane.m2 * tent_top.y + group.plane.b + cyl.radius))
        return false;
  }
  
  return true;
}


inline bool colis_group_range_flat
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d min, max;
  double clear = cyl.radius + .001;
  
  for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
    if (colis_tri_range_flat (tri, clear, curr_top, curr_bot, tent_top, tent_bot))
      return true;

  return false;
}


inline bool colis_tri_range_flat
  (tri_type* tri, double clear,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d min, max;
  
  get_z_range_3d (&min.z, &max.z, clear,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (curr_top.z < min.z)
    if (tent_top.z < min.z)
      return false;
  if (curr_top.z > max.z)
    if (tent_top.z > max.z)
      return false;

  get_x_range_3d (&min.x, &max.x, clear,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (curr_top.x < min.x)
    if (tent_top.x < min.x)
      return false;
  if (curr_top.x > max.x)
    if (tent_top.x > max.x)
      return false;

  return true;
}


inline bool colis_group_range
  (tri_group group, cylinder_3d cyl,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d min, max;
  double clear = cyl.radius + .001;
  
  for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
    if (colis_tri_range (tri, clear, curr_top, curr_bot, tent_top, tent_bot))
      return true;

  return false;
}


inline bool colis_tri_range
  (tri_type* tri, double clear,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d min, max;
  
  get_z_range_3d (&min.z, &max.z, clear,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (curr_top.z < min.z)
    if (tent_top.z < min.z)
      return false;
  if (curr_top.z > max.z)
    if (tent_top.z > max.z)
      return false;

  get_x_range_3d (&min.x, &max.x, clear,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (curr_top.x < min.x)
    if (tent_top.x < min.x)
      return false;
  if (curr_top.x > max.x)
    if (tent_top.x > max.x)
      return false;

  get_y_range_3d (&min.y, &max.y, .001,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (curr_bot.y < min.y)
    if (tent_bot.y < min.y)
      return false;
  if (curr_top.y > max.y)
    if (tent_top.y > max.y)
      return false;

  return true;
}


inline bool colis_path_gen (tri_group group, cylinder_3d cyl,
  bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d int_a, int_b;
  
  if (approx_zero (group.plane.m2)) {
//    if (colis_group_range (group, cyl, curr_top, curr_bot, tent_top, tent_bot))
    if (group.vis_side == 1) {
      if (colis_xy_group_flat_rev (slide, res_plane, result, closest, group, cyl,
        int_a, int_b, curr_top, curr_bot, tent_top, tent_bot)) {
//        *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
        return true;
      }
    }
    else
      if (colis_xy_group_flat (slide, res_plane, result, closest, group, cyl,
        int_a, int_b, curr_top, curr_bot, tent_top, tent_bot)) {
//        *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
        return true;
      }
  }
  else {
//    if (path_int_plane_xyz (group, cyl, &int_a, &int_b, curr_top, curr_bot, tent_top, tent_bot))
//      if (colis_group_range (group, cyl, curr_top, curr_bot, tent_top, tent_bot))
    if (group.vis_side == 1) {
      if (colis_xy_group_rev (slide, res_plane, result, closest, group, cyl,
        int_a, int_b, curr_top, curr_bot, tent_top, tent_bot)) {
//        *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
        return true;
      }
    }
    else {
      if (colis_xy_group (slide, res_plane, result, closest, group, cyl,
        int_a, int_b, curr_top, curr_bot, tent_top, tent_bot)) {
//        *closest = calc_new_loc (cyl, curr_top, tent_top, *closest);
        return true;
      }
    }
  }
  
  return false;
}


inline bool path_int_plane_xyz (tri_group group, cylinder_3d cyl,
  point_3d* int_a, point_3d* int_b,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  double offset = group.plane.sec_m1 * cyl.radius;

  if (approx_equal (curr_top.z, group.plane.m1 * curr_top.x + group.plane.m2 * curr_top.y + group.plane.b - offset))
    if (approx_equal (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b - offset))
      return false;
  if (approx_equal (curr_bot.z, group.plane.m1 * curr_bot.x + group.plane.m2 * curr_bot.y + group.plane.b - offset))
    if (approx_equal (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset))
      return false;
      
  if (group.vis_side == -1) {
    if (approx_lesser (curr_top.z, group.plane.m1 * curr_top.x + group.plane.m2 * curr_top.y + group.plane.b - offset))
      if (approx_lesser (curr_bot.z, group.plane.m1 * curr_bot.x + group.plane.m2 * curr_bot.y + group.plane.b - offset))
        if (approx_lesser (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b - offset))
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) // 1
            return false;
          else { // 2
            *int_a = path_pen_xyz (group.plane, curr_bot, tent_bot, -offset);
            *int_b = line_pen_xyz (group.plane, tent_top, -offset);
            bprint ("2", 7);
            return true;
          }
        else
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 3
            *int_a = path_pen_xyz (group.plane, curr_top, tent_top, -offset);
            *int_b = line_pen_xyz (group.plane, tent_top, -offset);
            bprint ("3", 7);
            return true;
          }
          else { // 4
            bprint ("4", 7);
            sort_intersections_xyz (group.plane, curr_top, curr_bot,
              path_pen_xyz (group.plane, curr_top, tent_top, -offset),
              path_pen_xyz (group.plane, curr_bot, tent_bot, -offset), int_a, int_b);
            return true;
          }
      else
        if (approx_lesser (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b - offset))
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 5
            return false;
//            *int_a = line_pen_xyz (group.plane, curr_top);
//            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot);
          }
          else { // 6
            *int_a = line_pen_xyz (group.plane, curr_top, -offset);
            *int_b = line_pen_xyz (group.plane, tent_top, -offset);
            bprint ("6", 7);
            return true;
          }
        else
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) // 7
            return false; // N/A
          else { // 8
            bprint ("8", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, -offset);
            *int_b = path_pen_xyz (group.plane, curr_top, tent_top, -offset);
            return true;
          }
    else
      if (approx_lesser (curr_bot.z, group.plane.m1 * curr_bot.x + group.plane.m2 * curr_bot.y + group.plane.b - offset))
        if (approx_lesser (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b - offset))
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 9
            return false;
//            *int_a = line_pen_xyz (group.plane, curr_top);
//            *int_b = path_pen_xyz (group.plane, curr_top, tent_top);
          }
          else // 10
            return false; // N/A
        else
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 11
            bprint ("11", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, -offset);
            *int_b = line_pen_xyz (group.plane, tent_top, -offset);
            return true;
          }
          else { // 12
            bprint ("12", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, -offset);
            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot, -offset);
            return true;
          }
      else
        if (approx_lesser (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b - offset))
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 13
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_top, tent_top);
//            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot);
          }
          else { // 14
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_top, tent_top);
//            *int_b = line_pen_xyz (group.plane, tent_top);
          }
        else
          if (approx_lesser (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b - offset)) { // 15
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_bot, tent_bot);
//            *int_b = line_pen_xyz (group.plane, tent_top);
          }
          else // 16
            return false;
  }
  else {
    if (approx_greater (curr_top.z, group.plane.m1 * curr_top.x + group.plane.m2 * curr_top.y + group.plane.b + offset))
      if (approx_greater (curr_bot.z, group.plane.m1 * curr_bot.x + group.plane.m2 * curr_bot.y + group.plane.b + offset))
        if (approx_greater (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b + offset))
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) // 1
            return false;
          else { // 2
            *int_a = path_pen_xyz (group.plane, curr_bot, tent_bot, offset);
            *int_b = line_pen_xyz (group.plane, tent_top, offset);
            bprint ("2", 7);
            return true;
          }
        else
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 3
            *int_a = path_pen_xyz (group.plane, curr_top, tent_top, offset);
            *int_b = line_pen_xyz (group.plane, tent_top, offset);
            bprint ("3", 7);
            return true;
          }
          else { // 4
            bprint ("4", 7);
            sort_intersections_xyz (group.plane, curr_top, curr_bot,
              path_pen_xyz (group.plane, curr_top, tent_top, offset),
              path_pen_xyz (group.plane, curr_bot, tent_bot, offset), int_a, int_b);
            return true;
          }
      else
        if (approx_greater (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b + offset))
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 5
            return false;
//            *int_a = line_pen_xyz (group.plane, curr_top);
//            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot);
          }
          else { // 6
            *int_a = line_pen_xyz (group.plane, curr_top, offset);
            *int_b = line_pen_xyz (group.plane, tent_top, offset);
            bprint ("6", 7);
            return true;
          }
        else
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) // 7
            return false; // N/A
          else { // 8
            bprint ("8", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, offset);
            *int_b = path_pen_xyz (group.plane, curr_top, tent_top, offset);
            return true;
          }
    else
      if (approx_greater (curr_bot.z, group.plane.m1 * curr_bot.x + group.plane.m2 * curr_bot.y + group.plane.b + offset))
        if (approx_greater (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b + offset))
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 9
            return false;
//            *int_a = line_pen_xyz (group.plane, curr_top);
//            *int_b = path_pen_xyz (group.plane, curr_top, tent_top);
          }
          else // 10
            return false; // N/A
        else
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 11
            bprint ("11", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, offset);
            *int_b = line_pen_xyz (group.plane, tent_top, offset);
            return true;
          }
          else { // 12
            bprint ("12", 7);
            *int_a = line_pen_xyz (group.plane, curr_top, offset);
            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot, offset);
            return true;
          }
      else
        if (approx_greater (tent_top.z, group.plane.m1 * tent_top.x + group.plane.m2 * tent_top.y + group.plane.b + offset))
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 13
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_top, tent_top);
//            *int_b = path_pen_xyz (group.plane, curr_bot, tent_bot);
          }
          else { // 14
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_top, tent_top);
//            *int_b = line_pen_xyz (group.plane, tent_top);
          }
        else
          if (approx_greater (tent_bot.z, group.plane.m1 * tent_bot.x + group.plane.m2 * tent_bot.y + group.plane.b + offset)) { // 15
            return false;
//            *int_a = path_pen_xyz (group.plane, curr_bot, tent_bot);
//            *int_b = line_pen_xyz (group.plane, tent_top);
          }
          else // 16
            return false;
  }
}


inline point_3d calc_new_loc
  (cylinder_3d cyl, point_3d pa, point_3d pb, point_3d ntersect)
{
  point_3d p;
  
  if (approx_equal (pa.z, pb.z))
    if (approx_equal (pa.x, pb.x)) {
      point_3d p_bot = offset_point_3d (pa, 0, cyl.height, 0);
      if (approx_lesser_inc (ntersect.y, pa.y))
        p = ntersect;
      else if (approx_greater_inc (ntersect.y, p_bot.y))
        p = offset_point_3d (ntersect, 0, -cyl.height, 0);
      else
        p = pb;
    }
    else {
      x_line_3d path = calc_x_line_3d (pa, pb);
      p.x = ntersect.x;
      p.y = path.my * ntersect.x + path.by;
      p.z = ntersect.z;//path.mz * ntersect.x + path.bz;
    }
  else {
    line_3d path = calc_line_3d (pa, pb);
    p.x = ntersect.x;//path.mx * ntersect.z + path.bx;
    p.y = path.my * ntersect.z + path.by;
    p.z = ntersect.z;
  }
  
  return p;
}


inline bool colis_yz_group
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane = *res_plane;
  
  if (path_int_plane_xy (group, cyl, &result_a, &result_b, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_b;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri (slide, &plane, result, &closer, group, cyl.radius, result_a, result_b, tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = plane;
    }
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_yz_group_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane = *res_plane;
  
  if (path_int_plane_xy (group, cyl, &result_a, &result_b, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_b;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri (slide, &plane, result, &closer, group, -cyl.radius, result_a, result_b, tri->t3d.p3->abs, tri->t3d.p2->abs, tri->t3d.p1->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = plane;
    }      
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_yz_group_x
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane;
  
  if (path_int_plane_x (group, cyl, &result_a, curr_top, curr_bot, tent_top, tent_bot)) {
    result_b = offset_point_3d (result_a, 0, cyl.height - 1, 0);
    closer = result_b;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri (slide, &plane, result, &closer, group, cyl.radius, result_a, result_b, tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = group.plane;
    }      
    return true;
  }
  else {
    return false;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_yz_group_x_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane;
  
  if (path_int_plane_x (group, cyl, &result_a, curr_top, curr_bot, tent_top, tent_bot)) {
    result_b = offset_point_3d (result_a, 0, cyl.height - 1, 0);
    closer = result_b;    
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri (slide, &plane, result, &closer, group, -cyl.radius, result_a, result_b, tri->t3d.p3->abs, tri->t3d.p2->abs, tri->t3d.p1->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = group.plane;
    }
    return true;
  }
  else {
    return false;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_xz_tris
  (tri_group group, point_3d* result, point_3d* closest,
   double offset, point_3d pa, point_3d pb)
{
  for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
    if (colis_xz_plane (result, closest, offset, pa, pb,
    tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs))
      return true;

  return false;
}


inline bool colis_xz_tris_rev
  (tri_group group, point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb)
{
  for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
    if (colis_xz_plane (result, closest, offset, pa, pb,
    tri->t3d.p3->abs, tri->t3d.p2->abs, tri->t3d.p1->abs))
      return true;

  return false;
}


inline bool colis_yz_tri
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  
  if (approx_lesser (p1.y, p2.y))
    if (approx_lesser (p3.y, p1.y)) { // 3,1,2
      point_3d mid = int_yz_edge_and_y (p3, p1, p2);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, mid, p1, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, mid, p1, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_greater (p3.y, p2.y)) { // 1,2,3
      point_3d mid = int_yz_edge_and_y (p1, p2, p3);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, mid, p2, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, mid, p2, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_equal (p3.y, p1.y)) // 1-3,2
      return colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, p3, p1, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_equal (p3.y, p2.y)) // 1,2-3
      return colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, p3, p2, curr_top, curr_bot, tent_top, tent_bot);
    else { // 1,3,2
      point_3d mid = int_yz_edge_and_y (p1, p3, p2);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, p3, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, p3, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
  else if (approx_greater (p1.y, p2.y))
    if (approx_lesser (p3.y, p2.y)) { // 3,2,1
      point_3d mid = int_yz_edge_and_y (p3, p2, p1);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, p2, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, p2, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 2,1,3
      point_3d mid = int_yz_edge_and_y (p2, p1, p3);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, p1, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, p1, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_equal (p3.y, p2.y)) // 2-3,1
      return colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, p2, p3, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_equal (p3.y, p1.y)) // 2,1-3
      return colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, p1, p3, curr_top, curr_bot, tent_top, tent_bot);
    else { // 2,3,1
      point_3d mid = int_yz_edge_and_y (p2, p3, p1);
      if (colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, mid, p3, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, mid, p3, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
  else
    if (approx_lesser (p3.y, p1.y)) // 3,1-2
      return colis_yz_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, p2, p1, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_greater (p3.y, p1.y)) // 1-2,3
      return colis_yz_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, p1, p2, curr_top, curr_bot, tent_top, tent_bot);
    else // 1-2-3
      return false;

  return colis;
}


inline bool colis_yz_tri_tb
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   tri_group group, double offset,
   point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  
  if (approx_lesser (p1.y, p2.y))
    if (approx_lesser (p3.y, p1.y)) { // 3,1,2
      point_3d mid = int_yz_edge_and_y (p3, p1, p2);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, mid, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p2, mid, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p2.y)) { // 1,2,3
      point_3d mid = int_yz_edge_and_y (p1, p2, p3);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, mid, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, mid, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p1.y)) { // 1-3,2
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p2, p3, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p2, p3, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p2.y)) { // 1,2-3
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, p3, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p1, p3, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else { // 1,3,2
      point_3d mid = int_yz_edge_and_y (p1, p3, p2);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, p3, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p3, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
  else if (approx_greater (p1.y, p2.y))
    if (approx_lesser (p3.y, p2.y)) { // 3,2,1
      point_3d mid = int_yz_edge_and_y (p3, p2, p1);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, p2, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p2, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 2,1,3
      point_3d mid = int_yz_edge_and_y (p2, p1, p3);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, p1, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, p1, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p2.y)) { // 2-3,1
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p1, p2, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p2, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p1.y)) { // 2,1-3
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, p1, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p2, p1, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else { // 2,3,1
      point_3d mid = int_yz_edge_and_y (p2, p3, p1);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, mid, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, mid, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
  else
    if (approx_lesser (p3.y, p1.y)) { // 3,1-2
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, p2, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p3, p2, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 1-2,3
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p3, p1, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, p1, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else // 1-2-3
      return false;

  return colis;
}


inline point_3d int_edge_and_z (point_3d p1, point_3d mid, point_3d p2)
{
  point_3d result;
  lin_relat edge = calc_dxdz_line (p1, p2);
  
  result.z = mid.z;
  result.x = edge.m * mid.z + edge.b;
  
  return result;
}


inline point_3d int_yz_edge_and_y (point_3d p1, point_3d mid, point_3d p2)
{
  point_3d result;
  lin_relat edge = calc_dzdy_line (p1, p2);
  
  result.y = mid.y;
  result.z = edge.m * mid.y + edge.b;
  
  return result;
}


inline bool colis_yz_tri_up
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  point_3d tl, tr, bl, br;
  bool colis = false;
  plane_type plane = *res_plane;
  
  if (group.vis_side == 1) {
    tl = offset_point_3d (t, offset, 0, -offset);
    tr = offset_point_3d (t, offset, 0, offset);
    bl = offset_point_3d (l, offset, 0, -offset);
    br = offset_point_3d (r, offset, 0, offset);
  }
  else {
    tl = offset_point_3d (t, offset, 0, offset);
    tr = offset_point_3d (t, offset, 0, -offset);
    bl = offset_point_3d (l, offset, 0, offset);
    br = offset_point_3d (r, offset, 0, -offset);
  }
//  _spec[0].abs = tl;
//  _spec[1].abs = tr;
//  _spec[1].select = true;
//  _spec[2].abs = bl;
//  _spec[3].abs = br;

  if (colis_yz_element_up (slide, &plane, result, group.plane, offset, pa, pb, tl, bl, br)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;    
  }
  if (colis_yz_element_dn (slide, &plane, result, group.plane, offset, pa, pb, br, tl, tr)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;
  }

  return colis;
}


inline bool colis_yz_tri_dn
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  point_3d tl, tr, bl, br;
  bool colis = false;
  plane_type plane = *res_plane;
  
  if (group.vis_side == 1) {
    tl = offset_point_3d (l, offset, 0, -offset);
    tr = offset_point_3d (r, offset, 0, offset);
    bl = offset_point_3d (b, offset, 0, -offset);
    br = offset_point_3d (b, offset, 0, offset);
//    return xz_quad_dn_top
//     (closest, b, l, r, offset, group, curr_top, curr_bot, tent_top, tent_bot);
  }
  else {
    tl = offset_point_3d (l, offset, 0, offset);
    tr = offset_point_3d (r, offset, 0, -offset);
    bl = offset_point_3d (b, offset, 0, offset);
    br = offset_point_3d (b, offset, 0, -offset);
//    return xz_quad_dn_bot
//     (closest, b, l, r, offset, group, curr_top, curr_bot, tent_top, tent_bot);
  }
  
  if (colis_yz_element_up (slide, &plane, result, group.plane, offset, pa, pb, tl, bl, br)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;
  }
  if (colis_yz_element_dn (slide, &plane, result, group.plane, offset, pa, pb, br, tl, tr)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;    
  }

  return colis;
}


inline bool colis_yz_element_up (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r)
{
  bool a_in = true;
  
  if (approx_greater_inc (pa.y, l.y)) {
    a_in = false;
    if (approx_greater (pb.y, l.y))
      return false;
    else
      if (path_int_yz_horz (pa, pb, plane, offset, l, r, slide, res_plane, result))
        return true;
  }

  lin_relat le = calc_dzdy_line (t, l);
  if (approx_lesser_inc (pa.z, le.m * pa.y + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.z, le.m * pb.y + le.b))
      return false;
    else
      if (path_int_yz_edge (pa, pb, plane, offset, t, l, slide, res_plane, result))
        return true;
  }

  lin_relat re = calc_dzdy_line (t, r);
  if (approx_greater_inc (pa.z, re.m * pa.y + re.b)) {
    a_in = false;  
    if (approx_greater (pb.z, re.m * pb.y + re.b))
      return false;
    else
      if (path_int_yz_edge (pa, pb, plane, offset, t, r, slide, res_plane, result))
        return true;
  }

  if (a_in) {
    *result = pa;
//    _spec[41].abs = *result;
    return true;
  }
  else
    return false;
}


inline bool colis_yz_element_dn (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r)
{
  bool a_in = true;
  
  if (approx_lesser_inc (pa.y, l.y)) {
    a_in = false;
    if (approx_lesser (pb.y, l.y))
      return false;
    else
      if (path_int_yz_horz (pa, pb, plane, offset, l, r, slide, res_plane, result))
        return true;
  }

  lin_relat le = calc_dzdy_line (b, l);
  if (approx_lesser_inc (pa.z, le.m * pa.y + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.z, le.m * pb.y + le.b))
      return false;
    else
      if (path_int_yz_edge (pa, pb, plane, offset, l, b, slide, res_plane, result))
        return true;
  }

  lin_relat re = calc_dzdy_line (b, r);
  if (approx_greater_inc (pa.z, re.m * pa.y + re.b)) {
    a_in = false;  
    if (approx_greater (pb.z, re.m * pb.y + re.b))
      return false;
    else
      if (path_int_yz_edge (pa, pb, plane, offset, r, b, slide, res_plane, result))
        return true;
  }

  if (a_in) {
    *result = pa;
    *res_plane = plane;
//    _spec[41].abs = *result;
    return true;
  }
  else
    return false;
}


inline bool path_int_xz_horz
  (point_3d pa, point_3d pb,
  point_3d l, point_3d r, bool* slide, plane_type* res_plane, point_3d* result)
{
  if (approx_equal (pa.z, pb.z)) {
    *result = pa;
    return true;
  }
  
  lin_relat path = calc_dxdz_line (pa, pb);
  double int_x = path.m * l.z + path.b;
  
  if (approx_greater_inc (int_x, l.x))
    if (approx_lesser_inc (int_x, r.x)) {
      result->y = l.y;
      result->z = l.z;
      result->x = int_x;
//      _spec[42].abs = *result;
      return true;
    }

  return false;
}


inline bool path_int_yz_horz
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d l, point_3d r, bool* slide, plane_type* res_plane, point_3d* result)
{
  if (approx_equal (pa.y, pb.y)) {
    *result = pa;
    return true;
  }
  
  lin_relat path = calc_dzdy_line (pa, pb);
  double int_z = path.m * l.y + path.b;
  
  if (approx_greater_inc (int_z, l.z))
    if (approx_lesser_inc (int_z, r.z)) {
      result->y = l.y;
      result->z = int_z;
      result->x = plane.m2 * l.y + plane.b + offset;
      *res_plane = edge_plane (l, r);
//      _spec[42].abs = *result;
      return true;
    }

  return false;
}


inline bool path_int_yz_edge
  (point_3d pa, point_3d pb, plane_type plane, double offset,
   point_3d p1, point_3d p2, bool* slide, plane_type* res_plane, point_3d* result)
{
  lin_relat path;
  lin_relat edge = calc_dzdy_line (p1, p2);

  if (approx_equal (pa.y, pb.y))
    result->y = pa.y;
  else {
    path = calc_dzdy_line (pa, pb);
    result->y = (path.b - edge.b) / (edge.m - path.m);
  }
    
  if (approx_greater_inc (result->y, p1.y))
    if (approx_lesser_inc (result->y, p2.y)) {
      result->z = edge.m * result->y + edge.b;
      result->x = plane.m2 * result->y + plane.b + offset;
      *res_plane = edge_plane (p1, p2);
      return true;
    }

  return false;
}


inline bool path_int_xz_edge
  (point_3d pa, point_3d pb,
   point_3d p1, point_3d p2, bool* slide, plane_type* res_plane, point_3d* result)
{
  lin_relat edge = calc_dxdz_line (p1, p2);
  
  if (approx_equal (pa.z, pb.z))
    result->z = pa.z;
  else {
    lin_relat path = calc_dxdz_line (pa, pb);
    result->z = (path.b - edge.b) / (edge.m - path.m);
  }
    
  if (approx_lesser_inc (result->z, p1.z))
    if (approx_greater_inc (result->z, p2.z)) {
      result->y = p1.y;
      result->x = edge.m * result->z + edge.b;
      return true;
    }

  return false;
}


/*
inline bool path_int_yz_r_edge
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d p1, point_3d p2, point_3d* result)
{
  lin_relat path = calc_dzdy_line (pa, pb);
  lin_relat edge = calc_dzdy_line (p1, p2);

  if (approx_equal (pa.y, pb.y))
    result->y = pa.y;
  else
    result->y = (path.b - edge.b) / (edge.m - path.m);
    
  if (approx_greater_inc (result->y, p1.y))
    if (approx_lesser_inc (result->y, p2.y)) {
      result->z = edge.m * result->y + edge.b;
      result->x = plane.m2 * result->y + plane.b + offset;
      return true;
    }

  return false;
}
*/


inline bool colis_xz_plane
  (point_3d* result, point_3d* closest, double offset,
  point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3)
{
  if (approx_greater (p1.z, p2.z))
    if (approx_greater (p3.z, p1.z)) { // 3,1,2
      point_3d mid = int_edge_and_z (p3, p1, p2);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p3, mid, p1))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p2, mid, p1))
        return true;
    }
    else if (approx_lesser (p3.z, p2.z)) { // 1,2,3
      point_3d mid = int_edge_and_z (p1, p2, p3);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p1, mid, p2))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p3, mid, p2))
        return true;
    }
    else if (approx_equal (p3.z, p1.z)) // 1-3,2
      return colis_xz_tri_dn (result, closest, offset, pa, pb, p2, p3, p1);
    else if (approx_equal (p3.z, p2.z)) // 1,2-3
      return colis_xz_tri_up (result, closest, offset, pa, pb, p1, p3, p2);
    else { // 1,3,2
      point_3d mid = int_edge_and_z (p1, p3, p2);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p1, p3, mid))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p2, p3, mid))
        return true;
    }
  else if (approx_lesser (p1.z, p2.z))
    if (approx_greater (p3.z, p2.z)) { // 3,2,1
      point_3d mid = int_edge_and_z (p3, p2, p1);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p3, p2, mid))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p1, p2, mid))
        return true;
    }
    else if (approx_lesser (p3.z, p1.z)) { // 2,1,3
      point_3d mid = int_edge_and_z (p2, p1, p3);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p2, p1, mid))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p3, p1, mid))
        return true;
    }
    else if (approx_equal (p3.z, p2.z)) // 2-3,1
      return colis_xz_tri_dn (result, closest, offset, pa, pb, p1, p2, p3);
    else if (approx_equal (p3.z, p1.z)) // 2,1-3
      return colis_xz_tri_up (result, closest, offset, pa, pb, p2, p1, p3);
    else { // 2,3,1
      point_3d mid = int_edge_and_z (p2, p3, p1);
      if (colis_xz_tri_up (result, closest, offset, pa, pb, p2, mid, p3))
        return true;
      if (colis_xz_tri_dn (result, closest, offset, pa, pb, p1, mid, p3))
        return true;
    }
  else
    if (approx_greater (p3.z, p1.z)) // 3,1-2
      return colis_xz_tri_up (result, closest, offset, pa, pb, p3, p2, p1);
    else if (approx_lesser (p3.z, p1.z)) // 1-2,3
      return colis_xz_tri_dn (result, closest, offset, pa, pb, p3, p1, p2);
    else // 1-2-3
      return false;

  return false;
}


inline bool colis_xz_tri_up (point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r)
{
  return colis_xz_element_up (result, closest, pa, pb, t, l, r);
}


inline bool colis_xz_tri_dn (point_3d* result, point_3d* closest,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r)
{
  return colis_xz_element_dn (result, closest, pa, pb, b, l, r);
}


inline bool colis_xz_element_up (point_3d* result, point_3d* closest,
  point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r)
{
  if (approx_lesser (pa.z, l.z))
    return false;

  lin_relat le = calc_dxdz_line (t, l);
  if (approx_lesser (pa.x, le.m * pa.z + le.b))
    return false;

  lin_relat re = calc_dxdz_line (t, r);
  if (approx_greater (pa.x, re.m * pa.z + re.b))
    return false;

  *closest = *result = pa;
  return true;
/*
  bool a_in = true;

  if (approx_lesser_inc (pa.z, l.z)) {
    a_in = false;
    if (approx_lesser (pb.z, l.z))
      return false;
    else
      if (path_int_xz_horz (pa, pb, l, r, result))
        return true;
  }

  lin_relat le = calc_dxdz_line (t, l);
  if (approx_lesser_inc (pa.x, le.m * pa.z + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.x, le.m * pb.z + le.b))
      return false;
    else
      if (path_int_xz_edge (pa, pb, t, l, result))
        return true;
  }

  lin_relat re = calc_dxdz_line (t, r);
  if (approx_greater_inc (pa.x, re.m * pa.z + re.b)) {
    a_in = false;  
    if (approx_greater (pb.x, re.m * pb.z + re.b))
      return false;
    else
      if (path_int_xz_edge (pa, pb, t, r, result))
        return true;
  }

  if (a_in) {
    *result = pa;
    return true;
  }
  else
    return false;
  */
}


inline bool colis_xz_element_dn (point_3d* result, point_3d* closest,
  point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r)
{
  if (approx_greater (pa.z, l.z))
    return false;

  lin_relat le = calc_dxdz_line (l, b);
  if (approx_lesser (pa.x, le.m * pa.z + le.b))
    return false;

  lin_relat re = calc_dxdz_line (r, b);
  if (approx_greater (pa.x, re.m * pa.z + re.b))
    return false;

  *closest = *result = pa;
  return true;

  /*
  bool a_in = true;

  if (approx_greater_inc (pa.z, l.z)) {
    a_in = false;
    if (approx_greater (pb.z, l.z))
      return false;
    else
      if (path_int_xz_horz (pa, pb, l, r, result))
        return true;
  }

  lin_relat le = calc_dxdz_line (b, l);
  if (approx_lesser_inc (pa.x, le.m * pa.z + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.x, le.m * pb.z + le.b))
      return false;
    else
      if (path_int_xz_edge (pa, pb, l, b, result))
        return true;
  }

  lin_relat re = calc_dxdz_line (b, r);
  if (approx_greater_inc (pa.x, re.m * pa.z + re.b)) {
    a_in = false;  
    if (approx_greater (pb.x, re.m * pb.z + re.b))
      return false;
    else
      if (path_int_xz_edge (pa, pb, r, b, result))
        return true;
  }

  if (a_in) {
    *result = pa;
    return true;
  }
  else
    return false;
  */
}


inline void store_closer
  (point_3d* closest, point_3d result, point_3d orig,
   plane_type plane, bool* slide, plane_type* res_plane)
{
  if (closer (orig, result, *closest)) {
    *closest = result;
    *res_plane = plane;
  }
}


inline bool colis_xz_quad_neg
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d closer, pa, pb;

  if (path_int_plane_xz_neg (p1.y, slide, &pa, curr_bot, tent_bot)) {
    double angle = atan2 (ref2.z - ref1.z, ref2.x - ref1.x);
    bool colis = false;

    double mul = sqrt(2);
    point_3d p1a = offset_3d_xz (radius, p1, angle - PI);
    point_3d p1b = offset_3d_xz (radius * mul, p1, angle - 3 * PI / 4);
    point_3d p2a = offset_3d_xz (radius, p2, angle);
    point_3d p2b = offset_3d_xz (radius * mul, p2, angle - PI / 4);

    if (colis_xz_plane (result, &closer, radius, pa, pb, p2b, p2a, p1b))
      colis = true;
    if (colis_xz_plane (result, &closer, radius, pa, pb, p2a, p1a, p1b))
      colis = true;

    if (colis) {
      closer = calc_new_loc (cyl, curr_top, tent_top, closer);
      if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top))
        *closest = closer;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}


inline bool colis_xz_quad_neg_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d closer, pa, pb;
  
  if (path_int_plane_xz_neg (p1.y, slide, &pa, curr_bot, tent_bot)) {
    double angle = atan2 (ref2.z - ref1.z, ref2.x - ref1.x);
    bool colis = false;

    double mul = sqrt(2);
    point_3d p1a = offset_3d_xz (radius, p1, angle - PI);
    point_3d p1b = offset_3d_xz (radius * mul, p1, angle - 3 * PI / 4);
    point_3d p2a = offset_3d_xz (radius, p2, angle);
    point_3d p2b = offset_3d_xz (radius * mul, p2, angle - PI / 4);

//    _spec[20].abs = p1a;
//    _spec[21].abs = p1b;
//    _spec[21].select = true;
//    _spec[22].abs = p2a;
//    _spec[23].abs = p2b;
    
    if (colis_xz_plane (result, &closer, radius, pa, pb, p2a, p2b, p1b))
      colis = true;
    if (colis_xz_plane (result, &closer, radius, pa, pb, p2a, p1b, p1a))
      colis = true;

    if (colis) {
      closer = calc_new_loc (cyl, curr_top, tent_top, closer);
      if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top))
        *closest = closer;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}


inline bool colis_xz_quad_pos
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d closer, pa, pb;
  
  if (path_int_plane_xz_pos (p1.y, slide, &pa, curr_top, tent_top)) {
    double angle = atan2 (ref2.z - ref1.z, ref2.x - ref1.x);
    bool colis = false;

    double mul = sqrt(2);
    point_3d p1a = offset_3d_xz (radius, p1, angle - PI);
    point_3d p1b = offset_3d_xz (radius * mul, p1, angle - 3 * PI / 4);
    point_3d p2a = offset_3d_xz (radius, p2, angle);
    point_3d p2b = offset_3d_xz (radius * mul, p2, angle - PI / 4);

//    _spec[20].abs = p1a;
//    _spec[21].abs = p1b;
//    _spec[21].select = true;
//    _spec[22].abs = p2a;
//    _spec[23].abs = p2b;
    
//    if (colis_xz_plane (result, radius, pa, pb, p2b, p2a, p1b))
//      colis = true;
//    if (colis_xz_plane (result, radius, pa, pb, p2a, p1a, p1b))
//      colis = true;
    if (colis_xz_plane (result, &closer, radius, pa, pb, p1b, p2a, p2b))
      colis = true;
    if (colis_xz_plane (result, &closer, radius, pa, pb, p1b, p1a, p2a))
      colis = true;

    if (colis) {
      closer = calc_new_loc (cyl, curr_top, tent_top, closer);
      if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top))
        *closest = closer;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}


inline bool colis_xz_quad_pos_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
  point_3d p1, point_3d p2, point_3d ref1, point_3d ref2, double radius,
  point_3d curr_top, point_3d curr_bot,
  point_3d tent_top, point_3d tent_bot)
{
  point_3d closer, ntersect, pb;
  
  if (path_int_plane_xz_pos (p1.y, slide, &ntersect, curr_top, tent_top)) {
    double angle = atan2 (ref2.z - ref1.z, ref2.x - ref1.x);
    bool colis = false;

    double mul = sqrt(2);
    point_3d p1a = offset_3d_xz (radius, p1, angle - PI);
    point_3d p1b = offset_3d_xz (radius * mul, p1, angle - 3 * PI / 4);
    point_3d p2a = offset_3d_xz (radius, p2, angle);
    point_3d p2b = offset_3d_xz (radius * mul, p2, angle - PI / 4);

//    _spec[20].abs = p1a;
//    _spec[21].abs = p1b;
//    _spec[21].select = true;
//    _spec[22].abs = p2a;
//    _spec[23].abs = p2b;
    
    if (colis_xz_plane (result, &closer, radius, ntersect, pb, p1b, p1a, p2b))
      colis = true;
    if (colis_xz_plane (result, &closer, radius, ntersect, pb, p1a, p2a, p2b))
      colis = true;

//    if (colis)
//      closer = calc_new_loc (cyl, curr_top, tent_top, closer);
      
//    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top))
    if (colis) {
      closer = calc_new_loc (cyl, curr_top, tent_top, closer);
      if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top))
        *closest = closer;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}


inline bool xz_quad_up_top
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   point_3d t, point_3d l, point_3d r, double offset,
   tri_group group, point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  if (group.vis_side == 1) 
    return colis_xz_quad_neg_rev (slide, res_plane, result, closest, cyl, t, t, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
  else
    return colis_xz_quad_neg_rev (slide, res_plane, result, closest, cyl, t, t, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
}


inline bool xz_quad_up_bot
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   point_3d t, point_3d l, point_3d r, double offset,
   tri_group group, point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  if (group.vis_side == 1)
    return colis_xz_quad_pos_rev (slide, res_plane, result, closest, cyl, l, r, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
  else
    return colis_xz_quad_pos_rev (slide, res_plane, result, closest, cyl, l, r, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
}


inline bool xz_quad_dn_top
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   point_3d b, point_3d l, point_3d r, double offset,
   tri_group group, point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  if (group.vis_side == 1) 
    return colis_xz_quad_neg_rev (slide, res_plane, result, closest, cyl, l, r, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
  else
    return colis_xz_quad_neg_rev (slide, res_plane, result, closest, cyl, l, r, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
}


inline bool xz_quad_dn_bot
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   point_3d b, point_3d l, point_3d r, double offset,
   tri_group group, point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  if (group.vis_side == 1) 
    return colis_xz_quad_pos_rev (slide, res_plane, result, closest, cyl, b, b, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
  else
    return colis_xz_quad_pos_rev (slide, res_plane, result, closest, cyl, b, b, l, r, offset,
    curr_top, curr_bot, tent_top, tent_bot);
}


inline bool colis_xy_element_up (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r)
{
  bool a_in = true;
  *res_plane = plane;
//  exit(123);
  
  if (approx_greater_inc (pa.y, l.y)) {
    a_in = false;
    if (approx_greater (pb.y, l.y))
      return false;
    else
      if (path_int_xy_horz (pa, pb, plane, offset, l, r, slide, res_plane, result))
        return true;
  }

  lin_relat le = calc_dxdy_line (t, l);
  if (approx_lesser_inc (pa.x, le.m * pa.y + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.x, le.m * pb.y + le.b))
      return false;
    else
      if (path_int_xy_edge (pa, pb, plane, offset, t, l, slide, res_plane, result))
        return true;
  }

  lin_relat re = calc_dxdy_line (t, r);
  if (approx_greater_inc (pa.x, re.m * pa.y + re.b)) {
    a_in = false;  
    if (approx_greater (pb.x, re.m * pb.y + re.b))
      return false;
    else
      if (path_int_xy_edge (pa, pb, plane, offset, t, r, slide, res_plane, result))
        return true;
  }

  if (a_in) {
    *result = pa;
    *res_plane = plane;
    return true;
  }
  else
    return false;
}


inline bool colis_xy_element_dn (bool* slide, plane_type* res_plane, point_3d* result, plane_type plane,
  double offset, point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r)
{
  bool a_in = true;
  *res_plane = plane;
//  exit(123);
  
  if (approx_lesser_inc (pa.y, l.y)) {
    a_in = false;
    if (approx_lesser (pb.y, l.y))
      return false;
    else
      if (path_int_xy_horz (pa, pb, plane, offset, l, r, slide, res_plane, result))
        return true;
  }

  lin_relat le = calc_dxdy_line (b, l);
  if (approx_lesser_inc (pa.x, le.m * pa.y + le.b)) {
    a_in = false;  
    if (approx_lesser (pb.x, le.m * pb.y + le.b))
      return false;
    else
      if (path_int_xy_edge (pa, pb, plane, offset, l, b, slide, res_plane, result))
        return true;
  }

  lin_relat re = calc_dxdy_line (b, r);
  if (approx_greater_inc (pa.x, re.m * pa.y + re.b)) {
    a_in = false;  
    if (approx_greater (pb.x, re.m * pb.y + re.b))
      return false;
    else
      if (path_int_xy_edge (pa, pb, plane, offset, r, b, slide, res_plane, result))
        return true;
  }

  if (a_in) {
    *result = pa;
    *res_plane = plane;    
    return true;
  }
  else
    return false;
}


inline bool path_int_xy_horz
  (point_3d pa, point_3d pb, plane_type plane, double offset,
  point_3d l, point_3d r, bool* slide, plane_type* res_plane, point_3d* result)
{
  if (approx_equal (pa.y, pb.y)) {
    *result = pa;
    *res_plane = edge_plane (pa, pb);
    return true;
  }
  
  lin_relat path = calc_dxdy_line (pa, pb);
  double int_x = path.m * l.y + path.b;
  
  if (approx_greater_inc (int_x, l.x))
    if (approx_lesser_inc (int_x, r.x)) {
      result->y = l.y;
      result->z = plane.m1 * int_x + plane.m2 * l.y + plane.b + offset;
      result->x = int_x;
      *res_plane = edge_plane (pa, pb);
      return true;
    }

  return false;
}


inline bool path_int_xy_edge
  (point_3d pa, point_3d pb, plane_type plane, double offset,
   point_3d p1, point_3d p2, bool* slide, plane_type* res_plane, point_3d* result)
{
  lin_relat path;
  lin_relat edge = calc_dxdy_line (p1, p2);

  if (approx_equal (pa.y, pb.y))
    result->y = pa.y;
  else {
    path = calc_dxdy_line (pa, pb);
    result->y = (path.b - edge.b) / (edge.m - path.m);
  }
    
  if (approx_greater_inc (result->y, p1.y))
    if (approx_lesser_inc (result->y, p2.y)) {
      result->x = edge.m * result->y + edge.b;
      result->z = plane.m1 * result->x + plane.m2 * result->y + plane.b + offset;
      *res_plane = edge_plane (p1, p2);
      return true;
    }

  return false;
}


inline bool colis_xy_group
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane = *res_plane;
  
  if (path_int_plane_xyz (group, cyl, &result_a, &result_b, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_b;
//    _spec[0].abs = result_a;
//    _spec[1].abs = result_b;
//    _spec[0].select = true;
    double offset = group.plane.sec_m1 * -cyl.radius;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_xy_tri (slide, &plane, result, &closer, group, offset, result_a, result_b, tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs, curr_top, curr_bot, tent_top, tent_bot)) {
        colis = true;
      }
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = plane;
    }
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_xy_group_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane = *res_plane;
  
  if (path_int_plane_xyz (group, cyl, &result_a, &result_b, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_b;
    double offset = group.plane.sec_m1 * cyl.radius;
//    _spec[0].abs = result_a;
//    _spec[1].abs = result_b;
//    _spec[0].select = true;
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_xy_tri (slide, &plane, result, &closer, group, offset, result_a, result_b, tri->t3d.p3->abs, tri->t3d.p2->abs, tri->t3d.p1->abs, curr_top, curr_bot, tent_top, tent_bot)) {
        colis = true;
      }
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = plane;
    }
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_xy_group_flat
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane;
  
  if (path_int_plane_xz (group, cyl, &result_a, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_a;
    result_b = offset_point_3d (result_a, 0, cyl.height, 0);
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_xy_tri (slide, &plane, result, &closer, group, -cyl.radius, result_a, result_b, tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = group.plane;
    }
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_xy_group_flat_rev
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, cylinder_3d cyl,
   point_3d pa, point_3d pb,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  point_3d closer, result_a, result_b;
  plane_type plane;
  
  if (path_int_plane_xz (group, cyl, &result_a, curr_top, curr_bot, tent_top, tent_bot)) {
    closer = result_a;
    result_b = offset_point_3d (result_a, 0, cyl.height, 0);
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_xy_tri (slide, &plane, result, &closer, group, cyl.radius, result_a, result_b, tri->t3d.p3->abs, tri->t3d.p2->abs, tri->t3d.p1->abs, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
  }

  if (colis) {
    closer = calc_new_loc (cyl, curr_top, tent_top, closer);
    if (dist_3d (closer, curr_top) < dist_3d (*closest, curr_top)) {
      *closest = closer;
      *res_plane = group.plane;
    }
    return true;
  }
  else {
    for (tri_type* tri = group.first; tri != NULL; tri = tri->next)
      if (colis_yz_tri_tb (slide, res_plane, result, closest, cyl, group, cyl.radius,
      tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs,
      curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    return colis;
  }
}


inline bool colis_xy_tri
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  
  if (approx_lesser (p1.y, p2.y))
    if (approx_lesser (p3.y, p1.y)) { // 3,1,2
      point_3d mid = int_xy_edge_and_y (p3, p1, p2);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, mid, p1, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, mid, p1, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_greater (p3.y, p2.y)) { // 1,2,3
      point_3d mid = int_xy_edge_and_y (p1, p2, p3);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, mid, p2, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, mid, p2, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_equal (p3.y, p1.y)) // 1-3,2
      return colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, p3, p1, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_equal (p3.y, p2.y)) // 1,2-3
      return colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, p3, p2, curr_top, curr_bot, tent_top, tent_bot);
    else { // 1,3,2
      point_3d mid = int_xy_edge_and_y (p1, p3, p2);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p1, p3, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p2, p3, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
  else if (approx_greater (p1.y, p2.y))
    if (approx_lesser (p3.y, p2.y)) { // 3,2,1
      point_3d mid = int_xy_edge_and_y (p3, p2, p1);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, p2, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, p2, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 2,1,3
      point_3d mid = int_xy_edge_and_y (p2, p1, p3);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, p1, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, p1, mid, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
    else if (approx_equal (p3.y, p2.y)) // 2-3,1
      return colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, p2, p3, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_equal (p3.y, p1.y)) // 2,1-3
      return colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, p1, p3, curr_top, curr_bot, tent_top, tent_bot);
    else { // 2,3,1
      point_3d mid = int_xy_edge_and_y (p2, p3, p1);
      if (colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p2, mid, p3, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
      if (colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p1, mid, p3, curr_top, curr_bot, tent_top, tent_bot))
        colis = true;
    }
  else
    if (approx_lesser (p3.y, p1.y)) // 3,1-2
      return colis_xy_tri_up (slide, res_plane, result, closest, group, offset, pa, pb, p3, p2, p1, curr_top, curr_bot, tent_top, tent_bot);
    else if (approx_greater (p3.y, p1.y)) // 1-2,3
      return colis_xy_tri_dn (slide, res_plane, result, closest, group, offset, pa, pb, p3, p1, p2, curr_top, curr_bot, tent_top, tent_bot);
    else // 1-2-3
      return false;

  return colis;
}


inline bool colis_xy_tri_up
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d t, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  point_3d tl, tr, bl, br;
  bool colis = false;
  double angle = atan2(r.z - l.z, r.x - l.x);
  double dist = fabs(offset) * sqrt(2);
  plane_type plane = *res_plane;
  
  if (group.vis_side == 1) {
    tl = offset_3d_xz (dist, t, angle + 3 * PI / 4);
    tr = offset_3d_xz (dist, t, angle + PI / 4);
    bl = offset_3d_xz (dist, l, angle + 3 * PI / 4);
    br = offset_3d_xz (dist, r, angle + PI / 4);
  }
  else {
    tl = offset_3d_xz (dist, t, angle + 5 * PI / 4);
    tr = offset_3d_xz (dist, t, angle - PI / 4);
    bl = offset_3d_xz (dist, l, angle + 5 * PI / 4);
    br = offset_3d_xz (dist, r, angle - PI / 4);
  }

  if (colis_xy_element_up (slide, &plane, result, group.plane, offset, pa, pb, tl, bl, br)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;    
  }
  if (colis_xy_element_dn (slide, &plane, result, group.plane, offset, pa, pb, br, tl, tr)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;
  }

  return colis;
}


inline bool colis_xy_tri_dn
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, tri_group group, double offset,
   point_3d pa, point_3d pb, point_3d b, point_3d l, point_3d r,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  point_3d tl, tr, bl, br;
  bool colis = false;
  double angle = atan2(r.z - l.z, r.x - l.x);
  double dist = fabs(offset) * sqrt(2);
  plane_type plane = *res_plane;
  
  if (group.vis_side == 1) {
    tl = offset_3d_xz (dist, l, angle + 3 * PI / 4);
    tr = offset_3d_xz (dist, r, angle + PI / 4);
    bl = offset_3d_xz (dist, b, angle + 3 * PI / 4);
    br = offset_3d_xz (dist, b, angle + PI / 4);
  }
  else {
    tl = offset_3d_xz (dist, l, angle + 5 * PI / 4);
    tr = offset_3d_xz (dist, r, angle - PI / 4);
    bl = offset_3d_xz (dist, b, angle + 5 * PI / 4);
    br = offset_3d_xz (dist, b, angle - PI / 4);
  }

//  _spec[4].abs = tl;
//  _spec[5].abs = tr;
//  _spec[5].select = true;
//  _spec[6].abs = bl;
//  _spec[7].abs = br;

  if (colis_xy_element_up (slide, &plane, result, group.plane, offset, pa, pb, tl, bl, br)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;
  }
  if (colis_xy_element_dn (slide, &plane, result, group.plane, offset, pa, pb, br, tl, tr)) {
    store_closer (closest, *result, pa, plane, slide, res_plane);
    colis = true;    
  }

//  *res_plane = group.plane;
  return colis;
}


inline point_3d int_xy_edge_and_y (point_3d p1, point_3d mid, point_3d p2)
{
  point_3d result;
  lin_relat edge = calc_dxdy_line (p1, p2);
  
  result.y = mid.y;
  result.x = edge.m * mid.y + edge.b;
  
  return result;
}


inline bool colis_xy_tri_tb
  (bool* slide, plane_type* res_plane, point_3d* result, point_3d* closest, cylinder_3d cyl,
   tri_group group, double offset,
   point_3d p1, point_3d p2, point_3d p3,
   point_3d curr_top, point_3d curr_bot,
   point_3d tent_top, point_3d tent_bot)
{
  bool colis = false;
  
  if (approx_lesser (p1.y, p2.y))
    if (approx_lesser (p3.y, p1.y)) { // 3,1,2
      point_3d mid = int_xy_edge_and_y (p3, p1, p2);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, mid, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p2, mid, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p2.y)) { // 1,2,3
      point_3d mid = int_xy_edge_and_y (p1, p2, p3);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, mid, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, mid, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p1.y)) { // 1-3,2
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p2, p3, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p2, p3, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p2.y)) { // 1,2-3
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, p3, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p1, p3, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else { // 1,3,2
      point_3d mid = int_xy_edge_and_y (p1, p3, p2);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p1, p3, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p3, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
  else if (approx_greater (p1.y, p2.y))
    if (approx_lesser (p3.y, p2.y)) { // 3,2,1
      point_3d mid = int_xy_edge_and_y (p3, p2, p1);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, p2, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p2, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 2,1,3
      point_3d mid = int_xy_edge_and_y (p2, p1, p3);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, p1, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, p1, mid, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p2.y)) { // 2-3,1
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p1, p2, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, p2, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_equal (p3.y, p1.y)) { // 2,1-3
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, p1, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p2, p1, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else { // 2,3,1
      point_3d mid = int_xy_edge_and_y (p2, p3, p1);
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p2, mid, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p1, mid, p3, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
  else
    if (approx_lesser (p3.y, p1.y)) { // 3,1-2
      if (xz_quad_up_top (slide, res_plane, result, closest, cyl, p3, p2, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_up_bot (slide, res_plane, result, closest, cyl, p3, p2, p1, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else if (approx_greater (p3.y, p1.y)) { // 1-2,3
      if (xz_quad_dn_top (slide, res_plane, result, closest, cyl, p3, p1, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
      if (xz_quad_dn_bot (slide, res_plane, result, closest, cyl, p3, p1, p2, offset,
        group, curr_top, curr_bot, tent_top, tent_bot))
        return true;
    }
    else // 1-2-3
      return false;

  return colis;
}


inline plane_type edge_plane (point_3d p1, point_3d p2)
{
  plane_type plane;
  
  if (approx_equal (p1.x, p2.x)) {
    plane.y_plane = false;
    plane.m1_inf = true;
    plane.b = p1.x;
    plane.m2 = 0;
  }
  else {
    lin_relat biv = calc_dzdx_line (p1, p2);
    plane.y_plane = false;
    plane.m1_inf = false;
    plane.b = biv.b;
    plane.m1 = biv.m;
    plane.m2 = 0;
  }

  return plane;
}


point_3d deflect_path
  (plane_type plane, cylinder_3d cyl,
  point_3d ntersect, point_3d curr, point_3d tent)
{
  point_3d new_loc;
  
  if (plane.y_plane) {
    new_loc.x = tent.x;
    new_loc.y = ntersect.y;
    new_loc.z = tent.z;
  }
  else if (plane.m1_inf) {
    point_3d imag = tent;
    new_loc = deflect_xy (plane, ntersect, imag);
  }
  else {
    point_3d imag = tent;
    new_loc = deflect_xyz (plane, ntersect, imag);
  }

  return new_loc;
}


inline point_3d deflect_xy
  (plane_type plane, point_3d ref, point_3d p)
{
  point_3d result;

  if (approx_zero (plane.m2)) {
    result.x = ref.x;
    result.y = p.y;
    result.z = p.z;
//    _spec[10].abs = result;
//    _spec[10].select = true;
  }
  else {
    double b1 = ref.x - plane.m2 * ref.y;
    double b2 = p.x + 1 / plane.m2 * p.y;
    result.y = (plane.m2 * (b2 - b1)) / (plane.m2 * plane.m2 + 1);
    result.x = plane.m2 * result.y + b1;
    result.z = p.z;
  }

  return result;
}


inline point_3d deflect_xyz
  (plane_type plane, point_3d ref, point_3d p)
{
  point_3d result;
  
  double imag_b = ref.z - plane.m1 * ref.x - plane.m2 * ref.y;
  double b1 = p.x + p.z * plane.m1;
  double b2 = p.y + p.z * plane.m2;
  result.z = (plane.m1 * b1 + plane.m2 * b2 + imag_b) /
    (1 + plane.m1 * plane.m1 + plane.m2 * plane.m2);
  result.x = b1 - result.z * plane.m1;
  result.y = b2 - result.z * plane.m2;

  return result;
}


void move_player (view_type view, camera_type* cursor, cylinder_3d cyl,
  point_type** first_point, tri_group** first_group)
{
  point_3d result;
  point_3d curr = offset_point_3d (cursor->pos, 0, -1, 0);
  check_controls__make__
    (cursor, view, first_point, first_group);
  point_3d tent = offset_point_3d (cursor->pos, 0, -1, 0);
  point_3d closest = tent;
  plane_type res_plane;
  bool hit;
  bool slide = false;
  res_plane.y_plane = false;
  res_plane.m1_inf = true;
//    do {
  hit = check_colis_cyl (*first_group, cyl, &slide, &res_plane, &closest, curr, tent);
//  if (hit)
//    tent = deflect_path (res_plane, cyl, closest, curr, tent);
//        curr = closest;
//      }
//    } while (!same_point (tent, curr) && hit);
  bprint ("hit: ", double(hit), 15);
  disp_plane (res_plane);
//  bprint ("diff: ", curr.z - closest.z, 16);

/*
  if (hit)
    if (coplanar (res_plane, curr, closest) && coplanar (res_plane, tent, closest))
      bprint ("slide", 16);
    else
      bprint ("no slide", 16);
  else
    bprint ("no hit", 16);
*/
//  bprint ("slide: ", double(slide), 16);
  
//  hit = check_colis_cyl (*first_group, cyl, &res_plane, &closest, curr, tent);
/*
  if (_collis) {
    hit = check_colis_cyl (*first_group, cyl, &slide, &res_plane, &closest, curr, tent);
    while (hit) {
      tent = deflect_path (res_plane, cyl, closest, curr, tent);
      curr = closest;
      hit = check_colis_cyl (*first_group, cyl, &slide, &res_plane, &closest, curr, tent);
    }
  }
*/
///*
  if (_collis) {
    if (hit)
      tent = closest;//deflect_path (res_plane, cyl, closest, curr, tent);
    cursor->pos = offset_point_3d (tent, 0, 1, 0);
  }
//*/
  cursor->pos = offset_point_3d (tent, 0, 1, 0);

  _spec[0].abs = closest;
  _spec[1].abs = tent;
  _spec[0].select = true;
}


void disp_plane (plane_type plane)
{
  if (plane.y_plane)
    bprint ("y_plane", 17);
  else if (plane.m1_inf)
    bprint ("m1_inf", 17);
  else
    bprint ("general", 17);
}


#endif
