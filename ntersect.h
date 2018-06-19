#include "math3d.h"
#include "gfx3d.h"

inline void intersect_result (tri_group g, ray_type ray, point_3d ref,
  point_3d inter, tri_type* t, point_3d* p3d_result, tri_type** tri_result);
tri_type* ray_int_tri (ray_type ray, tri_group* first_group,
  point_3d* p3d_result);
  

inline void intersect_result (tri_group g, ray_type ray, point_3d ref,
  point_3d inter, tri_type* t, point_3d* p3d_result, tri_type** tri_result)
{
  if (same_direction (ray, ref, inter))
    if (point_relative_to_plane (g.plane, ref) == g.vis_side) {
      if (*tri_result == NULL) {
        *p3d_result = inter;
        *tri_result = t;
      }
      else
        if (closer(ray.pos, inter, *p3d_result)) {
          *p3d_result = inter;
          *tri_result = t;
        }
    }
}


tri_type* ray_int_tri (ray_type ray, tri_group* first_group,
  point_3d* p3d_result)
{
  line_3d line;
  point_3d inter, ref;
  tri_type* tri_result = NULL;

  int line_case = calc_aim_line (ray, &line, &ref);
  
  if (line_case == -1) {
    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (g->vis_side == 1) {
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_x (ray.pos, g->plane,
          t->t3d.p2->abs, t->t3d.p1->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      }
      else {
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_x (ray.pos, g->plane,
          t->t3d.p1->abs, t->t3d.p2->abs, t->t3d.p3->abs, &inter))
          intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      }
  }
  else if (line_case == 0) {
    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (g->vis_side == 1) {
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_y (ray.pos, g->plane,
          t->t3d.p2->abs, t->t3d.p1->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      }
      else {
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_y (ray.pos, g->plane,
          t->t3d.p1->abs, t->t3d.p2->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      }
  }
  else if (line_case == 1) {
    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (g->vis_side == 1)
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_xy (line, ray.pos.z, g->plane,
          t->t3d.p2->abs, t->t3d.p1->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      else
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri_xy (line, ray.pos.z, g->plane,
          t->t3d.p1->abs, t->t3d.p2->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
  }
  else
    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (g->vis_side == 1)
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri (line, g->plane,
          t->t3d.p2->abs, t->t3d.p1->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;
      else
        for (tri_type* t = g->first; t != NULL; t = t->next)
          if (line_intersect_tri (line, g->plane,
          t->t3d.p1->abs, t->t3d.p2->abs, t->t3d.p3->abs, &inter))
            intersect_result (*g, ray, ref, inter, t, p3d_result, &tri_result);
          else
            ;

  return tri_result;
}


