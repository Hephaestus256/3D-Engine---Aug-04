#ifndef INCLUDE_TEXTRI
#define INCLUDE_TEXTRI

#include "math2d.h"
#include "math3d.h"
#include "gfx2d.h"
#include "gfx3d.h"
#include "/programs/include/fixed.h"
#include "texclass.h"
#include "texline.h"
#include "triclip.h"
#include <time.h>

inline long long curr_time();
inline bool is_front_side_of_tri (tri_type* tri);
inline bool point_left_of_edge (point_3d pn, point_3d p1, point_3d p2);
inline bool point_right_of_edge (point_3d pn, point_3d p1, point_3d p2);
uv_acc_type* create_uv_acc();
inline void draw_tri_group
  (tri_group g, view_type* v, uv_acc_type* acc, z_buffer_type* zmap, int num);
inline bool calc_horz_bmp_edge
  (line_3d l3d, line_2d* ref_2d, double scan_m, view_type view);
inline bool calc_vert_bmp_edge
  (line_3d l3d, line_2d* ref_2d, double scan_m, view_type view);

inline void draw_tri_group
  (tri_group g, view_type* v, uv_acc_type* acc, z_buffer_type* zmap, int num)
{
  if (point_relative_to_plane (g.plane, v->camera.pos) == g.vis_side) {
    tri_group_class group;
    _texdat = (int)acc;
//    bmp_type* b =
    _bmp = g.ref.bmp->dat;
    group.init_tex_plane (g, v, zmap, acc);
    group.group_num = num;
#ifdef MODE_MAKE
    for (group.curr = g.first; group.curr != NULL; group.curr = group.curr->next)
      group.draw_curr_tri();
#else
    for (group.curr = g.first; group.curr <= g.last; group.curr++)
      group.draw_curr_tri();
#endif
  }
}


void tri_group_class::draw_curr_tri()
{
  tri_type tri = *curr;
  #ifdef MODE_MAKE
  if (tri.select) {
    long m = long(CLOCK_GRAN * .5);
    long t = curr_time() % m;
    double rat;
    if (t < m / 2) {
      rat = t / double(m / 2);
      _bk_light[2] = int(256 * rat) << 8;
      _bk_light[1] = int(256 * rat) << 8;
      _bk_light[0] = int(256 * rat) << 8;
    }
    else {
      _bk_light[2] = 256 << 8;
      _bk_light[1] = 256 << 8;
      _bk_light[0] = 256 << 8;
    }
  }
  else {
    _bk_light[0] = tri.bk_light[0];
    _bk_light[1] = tri.bk_light[1];
    _bk_light[2] = tri.bk_light[2];
  }
  #else
  _bk_light[0] = tri.bk_light[0];
  _bk_light[1] = tri.bk_light[1];
  _bk_light[2] = tri.bk_light[2];
  #endif

  _exp_shad =
    mk_rgb_light (_bk_light[2] >> 8, _bk_light[1] >> 8, _bk_light[0] >> 8);

  view_clip_class t3d;
  temp_tri result[64];
  int ct = t3d.view_clip (curr->t3d, view, result);
  
  for (int i = 0; i < ct; i++) {
    temp_tri t_in[2], t_out[2];
    int in_ct, out_ct;
    z_cut (result[i], view->dither_dist, t_in, t_out, &in_ct, &out_ct);

    for (int j = 0; j < out_ct; j++)
      draw_tri_no_dither (t_out[j].p1, t_out[j].p2, t_out[j].p3);
    for (int j = 0; j < in_ct; j++)
      draw_tri_dither (t_in[j].p1, t_in[j].p2, t_in[j].p3);
  }
}


inline void tri_group_class::draw_tri_dither
  (point_3d p1, point_3d p2, point_3d p3)
{
//  int rot;

//  init_tex_plane (,
//    &dat, &rot, &dat.bmp, *view);
//  plane_type plane;

 // plane_type pl = _pl;
 // tri_cut tc = _tc;
  
//  tri_2d ins[100], outs[100];

//  int sm = 0;
//  int c = 255;
  
//      tc.disect_tri (map_to_scrn (*view, p1),
//                     map_to_scrn (*view, p2),
//                     map_to_scrn (*view, p3), ins, outs);

//  if (0 && (pl.b > 0) && (fabs(pl.m1) < 1000000) && (is_real_num(pl.m1)) && !pl.y_plane)

  if (tex_init.rot) {
    init_vert_tri (map_to_scrn(*view, p1),
                   map_to_scrn(*view, p2),
                   map_to_scrn(*view, p3));
    if (tex_init.scan_m > 0)
      draw_vert_pos_tex_tri_dither (map_to_scrn(*view, p1),
                             map_to_scrn(*view, p2),
                             map_to_scrn(*view, p3));
    else
      draw_vert_neg_tex_tri_dither (map_to_scrn(*view, p1),
                             map_to_scrn(*view, p2),
                             map_to_scrn(*view, p3));
  }
  else {
    init_horz_tri (map_to_scrn(*view, p1),
                   map_to_scrn(*view, p2),
                   map_to_scrn(*view, p3));
    draw_horz_tex_tri_dither (map_to_scrn(*view, p1),
                       map_to_scrn(*view, p2),
                       map_to_scrn(*view, p3));
  }
/*
    if (dat.rot) {
      _seg_mode = 0;
      init_horz_tri (map_to_scrn (*view, p1),
                     map_to_scrn (*view, p2),
                     map_to_scrn (*view, p3), &_uv_init, dat.scan_m);
      draw_horz_tex_tri (map_to_scrn (*view, p1),
                         map_to_scrn (*view, p2),
                         map_to_scrn (*view, p3), dat, c, *view);
    }
    else {
      _seg_mode = 0;
      init_vert_tri (map_to_scrn (*view, p1),
                     map_to_scrn (*view, p2),
                     map_to_scrn (*view, p3), &_uv_init, dat.scan_m);
      draw_vert_tex_tri (map_to_scrn (*view, p1),
                         map_to_scrn (*view, p2),
                         map_to_scrn (*view, p3), dat, c, *view);
    }
*/
}


inline void tri_group_class::draw_tri_no_dither
  (point_3d p1, point_3d p2, point_3d p3)
{
//  int rot;

//  init_tex_plane (,
//    &dat, &rot, &dat.bmp, *view);
//  plane_type plane;

 // plane_type pl = _pl;
 // tri_cut tc = _tc;
  
//  tri_2d ins[100], outs[100];

//  int sm = 0;
//  int c = 255;
  
//      tc.disect_tri (map_to_scrn (*view, p1),
//                     map_to_scrn (*view, p2),
//                     map_to_scrn (*view, p3), ins, outs);

//  if (0 && (pl.b > 0) && (fabs(pl.m1) < 1000000) && (is_real_num(pl.m1)) && !pl.y_plane)

  if (tex_init.rot) {
    init_vert_tri (map_to_scrn(*view, p1),
                   map_to_scrn(*view, p2),
                   map_to_scrn(*view, p3));
    if (tex_init.scan_m > 0) {
      draw_vert_pos_tex_tri (map_to_scrn(*view, p1),
                             map_to_scrn(*view, p2),
                             map_to_scrn(*view, p3));
    }
    else {
      draw_vert_neg_tex_tri (map_to_scrn(*view, p1),
                             map_to_scrn(*view, p2),
                             map_to_scrn(*view, p3));
    }
  }
  else {
    init_horz_tri (map_to_scrn(*view, p1),
                   map_to_scrn(*view, p2),
                   map_to_scrn(*view, p3));
    draw_horz_tex_tri (map_to_scrn(*view, p1),
                       map_to_scrn(*view, p2),
                       map_to_scrn(*view, p3));
  }
/*
    if (dat.rot) {
      _seg_mode = 0;
      init_horz_tri (map_to_scrn (*view, p1),
                     map_to_scrn (*view, p2),
                     map_to_scrn (*view, p3), &_uv_init, dat.scan_m);
      draw_horz_tex_tri (map_to_scrn (*view, p1),
                         map_to_scrn (*view, p2),
                         map_to_scrn (*view, p3), dat, c, *view);
    }
    else {
      _seg_mode = 0;
      init_vert_tri (map_to_scrn (*view, p1),
                     map_to_scrn (*view, p2),
                     map_to_scrn (*view, p3), &_uv_init, dat.scan_m);
      draw_vert_tex_tri (map_to_scrn (*view, p1),
                         map_to_scrn (*view, p2),
                         map_to_scrn (*view, p3), dat, c, *view);
    }
*/
}


void tri_group_class::draw_solid_horz_line (int x1, int x2, int y, int c)
{
  for (int x = x1; x <= x2; x++)
    pxl (x, y, 255, c, c, view->dest);
}


void tri_group_class::draw_horz_tex_tri
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect (p3.y, p1.y, line_23, line_31);
      draw_horz_sect (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect (p1.y, p2.y, line_31, line_12);
        draw_horz_sect (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect (p1.y, p3.y, line_31, line_12);
        draw_horz_sect (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_horz_sect (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect (p3.y, p2.y, line_23, line_31);
      draw_horz_sect (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect (p2.y, p1.y, line_12, line_23);
        draw_horz_sect (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect (p2.y, p3.y, line_12, line_23);
        draw_horz_sect (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_horz_sect (p2.y, p1.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


void tri_group_class::draw_vert_pos_tex_tri
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect (p3.y, p1.y, line_23, line_31);
      draw_vert_pos_sect (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect (p1.y, p2.y, line_31, line_12);
        draw_vert_pos_sect (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect (p1.y, p3.y, line_31, line_12);
        draw_vert_pos_sect (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_vert_pos_sect (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect (p3.y, p2.y, line_23, line_31);
      draw_vert_pos_sect (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect (p2.y, p1.y, line_12, line_23);
        draw_vert_pos_sect (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect (p2.y, p3.y, line_12, line_23);
        draw_vert_pos_sect (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_vert_pos_sect (p2.y, p3.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


void tri_group_class::draw_vert_neg_tex_tri
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect (p3.y, p1.y, line_23, line_31);
      draw_vert_neg_sect (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect (p1.y, p2.y, line_31, line_12);
        draw_vert_neg_sect (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect (p1.y, p3.y, line_31, line_12);
        draw_vert_neg_sect (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_vert_neg_sect (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect (p3.y, p2.y, line_23, line_31);
      draw_vert_neg_sect (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect (p2.y, p1.y, line_12, line_23);
        draw_vert_neg_sect (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect (p2.y, p3.y, line_12, line_23);
        draw_vert_neg_sect (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_vert_neg_sect (p2.y, p3.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


void tri_group_class::draw_horz_tex_tri_dither
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect_dither (p3.y, p1.y, line_23, line_31);
      draw_horz_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect_dither (p1.y, p2.y, line_31, line_12);
        draw_horz_sect_dither (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect_dither (p1.y, p3.y, line_31, line_12);
        draw_horz_sect_dither (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect_dither (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_horz_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect_dither (p3.y, p2.y, line_23, line_31);
      draw_horz_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect_dither (p2.y, p1.y, line_12, line_23);
        draw_horz_sect_dither (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_horz_sect_dither (p2.y, p3.y, line_12, line_23);
        draw_horz_sect_dither (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_horz_sect_dither (p2.y, p1.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect_dither (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_horz_sect_dither (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


void tri_group_class::draw_vert_pos_tex_tri_dither
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect_dither (p3.y, p1.y, line_23, line_31);
      draw_vert_pos_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect_dither (p1.y, p2.y, line_31, line_12);
        draw_vert_pos_sect_dither (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect_dither (p1.y, p3.y, line_31, line_12);
        draw_vert_pos_sect_dither (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect_dither (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_vert_pos_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect_dither (p3.y, p2.y, line_23, line_31);
      draw_vert_pos_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect_dither (p2.y, p1.y, line_12, line_23);
        draw_vert_pos_sect_dither (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_pos_sect_dither (p2.y, p3.y, line_12, line_23);
        draw_vert_pos_sect_dither (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_vert_pos_sect_dither (p2.y, p3.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect_dither (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_pos_sect_dither (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


void tri_group_class::draw_vert_neg_tex_tri_dither
  (point_2d p1, point_2d p2, point_2d p3)
{
  if (p1.y < p2.y)
    if (p3.y < p1.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect_dither (p3.y, p1.y, line_23, line_31);
      draw_vert_neg_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("3, 1, 2");
    }
    else if (p3.y > p1.y)
      if (p3.y > p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect_dither (p1.y, p2.y, line_31, line_12);
        draw_vert_neg_sect_dither (p2.y, p3.y, line_31, line_23);
//        printf ("1, 2, 3");
      }
      else if (p3.y < p2.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect_dither (p1.y, p3.y, line_31, line_12);
        draw_vert_neg_sect_dither (p3.y, p2.y, line_23, line_12);
//        printf ("1, 3, 2");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect_dither (p1.y, p2.y, line_31, line_12);
//        printf ("1, 2-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      draw_vert_neg_sect_dither (p1.y, p2.y, line_23, line_12);
//      printf ("1-3, 2");
    }
  else if (p1.y > p2.y)
    if (p3.y < p2.y) {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect_dither (p3.y, p2.y, line_23, line_31);
      draw_vert_neg_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("3, 2, 1");
    }
    else if (p3.y > p2.y)
      if (p3.y > p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect_dither (p2.y, p1.y, line_12, line_23);
        draw_vert_neg_sect_dither (p1.y, p3.y, line_31, line_23);
//        printf ("2, 1, 3");
      }
      else if (p3.y < p1.y) {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
        draw_vert_neg_sect_dither (p2.y, p3.y, line_12, line_23);
        draw_vert_neg_sect_dither (p3.y, p1.y, line_12, line_31);
//        printf ("2, 3, 1");
      }
      else {
        lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
        lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
        draw_vert_neg_sect_dither (p2.y, p3.y, line_12, line_23);
//        printf ("2, 1-3");
      }
    else {
      lin_relat line_12 = calc_2d_dxdy_line (p1, p2);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect_dither (p2.y, p1.y, line_12, line_31);
//      printf ("2-3, 1");
    }
  else
    if (p3.y < p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect_dither (p3.y, p1.y, line_23, line_31);
//      printf ("3, 1-2");
    }
    else if (p3.y > p1.y) {
      lin_relat line_23 = calc_2d_dxdy_line (p2, p3);
      lin_relat line_31 = calc_2d_dxdy_line (p3, p1);
      draw_vert_neg_sect_dither (p1.y, p3.y, line_31, line_23);
//      printf ("1-2, 3");
    }
}


/*
inline void tri_group_class::draw_left_sect_pos
  (double top, double bot, lin_relat left_edge, lin_relat right_edge)
{
  *l = init_edge (left_edge, *y);

  for (; *y <= sect_bot (bot); (*y)++) {
//    draw_solid_horz_line (line_start(*l), line_end(*r), (int)*y, _color);
    gfx_h_line ((int)*y, line_start(*l), line_end(*r), tex_init.scan_m, 0, *view);
    *l += left_edge.m;
    *r += right_edge.m;
  }
}


inline void tri_group_class::draw_right_sect_pos
  (double bot, int* y,
   line_equat left_edge, line_equat right_edge, dbl* l, dbl* r)
{
  *r = init_edge (right_edge, *y);

  for (; *y <= sect_bot (bot); (*y)++) {
//    draw_solid_horz_line (line_start(*l), line_end(*r), (int)*y, _color);
    gfx_h_line ((int)*y, line_start(*l), line_end(*r), tex_init.scan_m, 0, *view);
    *l += left_edge.m;
    *r += right_edge.m;
  }
}


inline void tri_group_class::draw_first_sect_neg
  (double bot, double top, int* y,
   line_equat left_edge, line_equat right_edge, dbl* l, dbl* r)
{
  *y = sect_bot (bot);
  *l = init_edge (left_edge, *y);
  *r = init_edge (right_edge, *y);

  for (; *y >= sect_top (top); (*y)--) {
//    draw_solid_horz_line (line_start(*l), line_end(*r), (int)*y, _color);
//    gfx_v_line_neg ((int)*y, line_start(*l), line_end(*r), scan_m, 0, *view);
    *l -= left_edge.m;
    *r -= right_edge.m;
  }
}


inline void tri_group_class::draw_left_sect_neg
  (double top, int* y,
   line_equat left_edge, line_equat right_edge, dbl* l, dbl* r)
{
  *l = init_edge (left_edge, *y);

  for (; *y >= sect_top (top); (*y)--) {
//    draw_solid_horz_line (line_start(*l), line_end(*r), (int)*y, _color);
//    gfx_v_line_neg ((int)*y, line_start(*l), line_end(*r), scan_m, 0, *view);
    *l -= left_edge.m;
    *r -= right_edge.m;
  }
}


inline void tri_group_class::draw_right_sect_neg
  (double top, int* y,
   line_equat left_edge, line_equat right_edge, dbl* l, dbl* r)
{
  *r = init_edge (right_edge, *y);

  for (; *y >= sect_top (top); (*y)--) {
//    draw_solid_horz_line (line_start(*l), line_end(*r), (int)*y, _color);
//    gfx_h_line ((int)*y, line_start(*l), line_end(*r), scan_m, 0, *view);
    *l -= left_edge.m;
    *r -= right_edge.m;
  }
}
*/

void tri_group_class::init_tex_plane
  (tri_group group, view_type* v, z_buffer_type* zm, uv_acc_type* acc)
{
  int dydx = 1;
  int dxdy = 0;
  view = v;
  zmap = zm;
  tex_init.dest = (int)v->dest;
  tex_init.acc = acc;
  tex_init.bmp = group.ref.bmp->dat;
  tex_init.orig = group.ref.orig->rel.z;
  tex_init.bk_light[0] = 128 << 8;
  tex_init.bk_light[1] = 128 << 8;
  tex_init.bk_light[2] = 128 << 8;
  
  if (approx_equal (group.ref.orig->rel.z, group.ref.u->rel.z))
    if (approx_equal (group.ref.orig->rel.z, group.ref.v->rel.z))
      if (fabs(group.ref.orig->rel.x - group.ref.u->rel.x) < fabs(group.ref.orig->rel.y - group.ref.u->rel.y)) { // flat
        tex_init.rot = dxdy;
        point_2d p0 = map_to_scrn (*view, group.ref.orig);
        point_2d pu = map_to_scrn (*view, group.ref.u);
        point_2d pv = map_to_scrn (*view, group.ref.v);
        tex_init.scan_m = (pu.x - p0.x) / (pu.y - p0.y);
        tex_init.s1 = p0.x - p0.y * tex_init.scan_m;
        tex_init.s2 = pv.x - pv.y * tex_init.scan_m;
        line_equat top;
        calc_line (&top, p0, pv);
        double bot_b = pu.y - top.m * pu.x;
        tex_init.edge1.m = tex_init.edge2.m = top.m / (1 - tex_init.scan_m * top.m);
        tex_init.edge1.b = top.b / (1 - tex_init.scan_m * top.m);
        tex_init.edge2.b = bot_b / (1 - tex_init.scan_m * top.m);
        tex_init.jump = 1;
        tex_init.inv_diff_1 = 1 / (tex_init.s2 - tex_init.s1);
      }
      else { // flat
        tex_init.rot = dxdy;
        point_2d p0 = map_to_scrn (*view, group.ref.orig);
        point_2d pu = map_to_scrn (*view, group.ref.u);
        point_2d pv = map_to_scrn (*view, group.ref.v);
        tex_init.scan_m = (pv.x - p0.x) / (pv.y - p0.y);
        tex_init.s1 = p0.x - p0.y * tex_init.scan_m;
        tex_init.s2 = pu.x - pu.y * tex_init.scan_m;
        line_equat top;
        calc_line (&top, p0, pu);
        double bot_b = pv.y - top.m * pv.x;
        tex_init.edge1.m =
        tex_init.edge2.m = top.m / (1 - tex_init.scan_m * top.m);
        tex_init.edge1.b = top.b / (1 - tex_init.scan_m * top.m);
        tex_init.edge2.b = bot_b / (1 - tex_init.scan_m * top.m);
        tex_init.jump = 2;
        tex_init.inv_diff_1 = 1 / (tex_init.s2 - tex_init.s1);
      }
    else { // u-0, v
      double dx = (group.ref.u->rel.x - group.ref.orig->rel.x) * view->zoom.x;
      double dy = (group.ref.u->rel.y - group.ref.orig->rel.y) * view->zoom.y;
      tex_init.inv_diff_1 = 1 / (group.ref.v->rel.z - group.ref.orig->rel.z);
      tex_init.jump = 3;
      if (fabs(dx) > fabs(dy)) {
        tex_init.scan_m = dy / dx;
        tex_init.rot = dydx;
        line_3d l3d_a = calc_line_3d (group.ref.orig, group.ref.v);
        line_3d l3d_b = calc_parll_line_3d (l3d_a, group.ref.u->rel);
        calc_horz_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_horz_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
        else {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
      }
      else {
        tex_init.scan_m = dx / dy;
        tex_init.rot = dxdy;
        line_3d l3d_a = calc_line_3d (group.ref.orig, group.ref.v);
        line_3d l3d_b = calc_parll_line_3d (l3d_a, group.ref.u->rel);
        calc_vert_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_vert_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
        else {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
      }
    }
  else
    if (approx_equal (group.ref.orig->rel.z, group.ref.v->rel.z)) {
      double dx = (group.ref.v->rel.x - group.ref.orig->rel.x) * view->zoom.x;
      double dy = (group.ref.v->rel.y - group.ref.orig->rel.y) * view->zoom.y;
      tex_init.jump = 4;
      tex_init.inv_diff_1 = 1 / (group.ref.u->rel.z - group.ref.orig->rel.z);
      if (fabs(dx) > fabs(dy)) {
        tex_init.scan_m = dy / dx;
        tex_init.rot = dydx;
        line_3d l3d_a = calc_line_3d (group.ref.orig, group.ref.u);
        line_3d l3d_b = calc_parll_line_3d (l3d_a, group.ref.v->rel);
        calc_horz_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_horz_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
        else {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
      }
      else {
        tex_init.scan_m = dx / dy;
        tex_init.rot = dxdy;
        line_3d l3d_a = calc_line_3d (group.ref.orig, group.ref.u);
        line_3d l3d_b = calc_parll_line_3d (l3d_a, group.ref.v->rel);
        calc_vert_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_vert_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
        else {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
      }
    }
    else {
      line_3d l3d_a = calc_line_3d (group.ref.orig, group.ref.u);
      line_3d l3d_b = calc_line_3d (group.ref.orig, group.ref.v);
      double x2 = l3d_a.mx * group.ref.v->rel.z + l3d_a.bx;
      double y2 = l3d_a.my * group.ref.v->rel.z + l3d_a.by;
      double dx = view->zoom.x * (x2 - group.ref.v->rel.x);
      double dy = view->zoom.y * (y2 - group.ref.v->rel.y);
      tex_init.inv_diff_1 = 1 / (group.ref.u->rel.z - group.ref.orig->rel.z);
      tex_init.inv_diff_2 = 1 / (group.ref.v->rel.z - group.ref.orig->rel.z);
      tex_init.jump = 5;
      
      if (fabs(dx) > fabs(dy)) {
        tex_init.scan_m = dy / dx;
        tex_init.rot = dydx;
        calc_horz_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_horz_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
        else {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
      }
      else {
        tex_init.scan_m = dx / dy;
        tex_init.rot = dxdy;
        calc_vert_bmp_edge (l3d_b, &tex_init.edge2, tex_init.scan_m, *view);
        if (calc_vert_bmp_edge (l3d_a, &tex_init.edge1, tex_init.scan_m, *view)) {
          tex_init.zoom = view->zoom.x;
          tex_init.center = view->center.x;
          tex_init.ref_2d.m = 1;
          tex_init.ref_2d.b = tex_init.scan_m * tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.mx;
          tex_init.ref_3d.b = l3d_a.bx;
        }
        else {
          tex_init.zoom = view->zoom.y;
          tex_init.center = view->center.y;
          tex_init.ref_2d.m = tex_init.edge1.m;
          tex_init.ref_2d.b = tex_init.edge1.b;
          tex_init.ref_3d.m = l3d_a.my;
          tex_init.ref_3d.b = l3d_a.by;
        }
      }
    }
}

/*
void update_uv ()
  //eax = uv_calc_group.ref* init, ecx = uv_acc_type* dat, edx = int s
{
  uv_calc_ref* init;
  uv_acc_type* dat;
  int s;
  
  asm (
    "movl %%eax, %0 \n"
    "movl %%ebx, %1 \n"
    "movl %%ecx, %2 \n"
  :"=g" (init), "=g" (dat), "=g" (s)
  :
  :"memory"
  );
  
  init = &_uv_init;
  dat = &_texdat[s + DAT_MID];
  
  if (init->jump == 1) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    dat->u1 = 0;
    dat->u2 = init->bmp->width;
    dat->v1 =
    dat->v2 = init->bmp->height * (s - init->s1) * init->inv_diff_1;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 2) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    dat->u1 =
    dat->u2 = init->bmp->width * (s - init->s1) * init->inv_diff_1;
    dat->v1 = 0;
    dat->v2 = init->bmp->height;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 3) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 = 0;
    dat->u2 = init->bmp->width;
    dat->v1 =
    dat->v2 = (z - init->orig) * init->inv_diff_1 * init->bmp->height;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 4) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 =
    dat->u2 = (z - init->orig) * init->inv_diff_1 * init->bmp->width;
    dat->v1 = 0;
    dat->v2 = init->bmp->height;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 5) {
    dat->ind1 = init->edge2.m * s + init->edge2.b;
    dat->ind2 = init->edge1.m * s + init->edge1.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 = 0;
    dat->u2 = (z - init->orig) * init->inv_diff_1 * init->bmp->width;
    dat->v1 = (z - init->orig) * init->inv_diff_2 * init->bmp->height;
    dat->v2 = 0;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }

  if (s < _dtemp[20])
    _dtemp[20] = s;
  if (s > _dtemp[21])
    _dtemp[21] = s;

  _dtemp[22]++;
}
*/

inline void update_uv
  (uv_calc_ref* init, uv_acc_type* dat, int s, int num)
{
  dat->init = num;
/*
  dat->last_ind = -99;
  dat->u1 = 0;
  dat->u2 = 0;
  dat->v1 = 0;
  dat->v2 = 0;
  dat->du.valu = 0;
  dat->dv.valu = 0;

  return;
*/
  if (init->jump == 1) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    dat->u1 = 0;
    dat->u2 = init->bmp.side;
    dat->v1 =
    dat->v2 = init->bmp.side * (s - init->s1) * init->inv_diff_1;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 2) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    dat->u1 =
    dat->u2 = init->bmp.side * (s - init->s1) * init->inv_diff_1;
    dat->v1 = 0;
    dat->v2 = init->bmp.side;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 3) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 = 0;
    dat->u2 = init->bmp.side;
    dat->v1 =
    dat->v2 = (z - init->orig) * init->inv_diff_1 * init->bmp.side;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 4) {
    dat->ind1 = init->edge1.m * s + init->edge1.b;
    dat->ind2 = init->edge2.m * s + init->edge2.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 =
    dat->u2 = (z - init->orig) * init->inv_diff_1 * init->bmp.side;
    dat->v1 = 0;
    dat->v2 = init->bmp.side;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
  else if (init->jump == 5) {
    dat->ind1 = init->edge2.m * s + init->edge2.b;
    dat->ind2 = init->edge1.m * s + init->edge1.b;
    dat->inv_len = 1 / (dat->ind2 - dat->ind1);
    dat->last_ind = -99;
    double p = init->ref_2d.m * s + init->ref_2d.b;
    double z = (init->zoom * init->ref_3d.b) /
      (p - init->center - init->zoom * init->ref_3d.m);
    dat->u1 = 0;
    dat->u2 = (z - init->orig) * init->inv_diff_1 * init->bmp.side;
    dat->v1 = (z - init->orig) * init->inv_diff_2 * init->bmp.side;
    dat->v2 = 0;
    dat->du = cv_sng((dat->u2 - dat->u1) * dat->inv_len);
    dat->dv = cv_sng((dat->v2 - dat->v1) * dat->inv_len);
  }
//  dat->init = 1;
}


inline void tri_group_class::init_horz_tri
  (point_2d p1, point_2d p2, point_2d p3)
{
  double s1 = calc_dxdy_b (p1, tex_init.scan_m);
  double s2 = calc_dxdy_b (p2, tex_init.scan_m);
  double s3 = calc_dxdy_b (p3, tex_init.scan_m);
  
  double min = s1;
  if (s2 < min)
    min = s2;
  if (s3 < min)
    min = s3;
    
  double max = s1;
  if (s2 > max)
    max = s2;
  if (s3 > max)
    max = s3;

  for (int s = (int)floor(min) - 1; s <= (int)ceil(max) + 1; s++)
    if (tex_init.acc[s].init != group_num)
      update_uv (&tex_init, &tex_init.acc[s], s, group_num);
}


inline void tri_group_class::init_vert_tri
  (point_2d p1, point_2d p2, point_2d p3)
{
  double s1 = calc_dydx_b (p1, tex_init.scan_m);
  double s2 = calc_dydx_b (p2, tex_init.scan_m);
  double s3 = calc_dydx_b (p3, tex_init.scan_m);
  
  double min = s1;
  if (s2 < min)
    min = s2;
  if (s3 < min)
    min = s3;
    
  double max = s1;
  if (s2 > max)
    max = s2;
  if (s3 > max)
    max = s3;

  for (int s = (int)floor(min) - 1; s <= (int)ceil(max) + 1; s++)
    if (tex_init.acc[s].init != group_num)
      update_uv (&tex_init, &tex_init.acc[s], s, group_num);
}


/*
void gfx_h_seg2 (int x, int y1, int y2, double scan_m, int color,
  view_type view)
{
  asm volatile (
    "push %%edi \n"
    "movl __x2_buff, %%edi \n"
    "addl __scrn, %%edi \n"
    "movl __f, %%esi \n"
    "push %%ebp \n"
    "movl %%esp, __temp \n"
    "leal __texdat, %%ebp \n"
    "addl __temp4, %%ebp \n"
//    "movl __poly, %%edx \n"
    "movl __temp8, %%ebx \n"

  "l_gfx_h_seg21: \n"

//    "movl (%%ebx), %%eax \n"
//    "cmpl %%eax, __poly \n"
//    "jne s_gfx_h_seg21 \n"

    "movl 0 + 16(%%ebp), %%eax \n"
    "movl __temp11, %%ecx \n"
    "cmpl %%ecx, %%eax \n"
    "je s_gfx_hseg22 \n"
      "fldl __dtemp + 8 * 1 \n" // load y
      "fldl 0 + 52(%%ebp) \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fldl 0 + 68(%%ebp) \n"
      "fmulp %%st(0), %%st(1) \n"
      "fstpl __dtemp + 8 * 0 \n"

      "fldl 0 + 20(%%ebp) \n"
      "fldl 0 + 36(%%ebp) \n"
      "fsubp %%st(0), %%st(1) \n"
      "fldl __dtemp + 0 * 8 \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl 0 + 20(%%ebp) \n"
      "faddp %%st(0), %%st(1) \n"
      "fldl __f2sng \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl __half \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fistpl 0 + 0(%%ebp) \n"

      "fldl 0 + 28(%%ebp) \n"
      "fldl 0 + 44(%%ebp) \n"
      "fsubp %%st(0), %%st(1) \n"
      "fldl __dtemp + 0 * 8 \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl 0 + 28(%%ebp) \n"
      "faddp %%st(0), %%st(1) \n"
      "fldl __f2sng \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl __half \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fistpl 0 + 4(%%ebp) \n"
      "movl %%ecx, 0 + 16(%%ebp) \n"
    "s_gfx_hseg22: \n"

    "movl 88 + 16(%%ebp), %%eax \n"
    "cmpl %%ecx, %%eax \n"
    "je s_gfx_hseg23 \n"

      "fldl __dtemp + 8 * 1 \n" // load y
      "fldl 88 + 52(%%ebp) \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fldl 88 + 68(%%ebp) \n"
      "fmulp %%st(0), %%st(1) \n"
      "fstpl __dtemp + 8 * 0 \n"

      "fldl 88 + 20(%%ebp) \n"
      "fldl 88 + 36(%%ebp) \n"
      "fsubp %%st(0), %%st(1) \n"
      "fldl __dtemp + 0 * 8 \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl 88 + 20(%%ebp) \n"
      "faddp %%st(0), %%st(1) \n"
      "fldl __f2sng \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl __half \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fistpl 88 + 0(%%ebp) \n"

      "fldl 88 + 28(%%ebp) \n"
      "fldl 88 + 44(%%ebp) \n"
      "fsubp %%st(0), %%st(1) \n"
      "fldl __dtemp + 0 * 8 \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl 88 + 28(%%ebp) \n"
      "faddp %%st(0), %%st(1) \n"
      "fldl __f2sng \n"
      "fmulp %%st(0), %%st(1) \n"
      "fldl __half \n"
      "fsubrp %%st(0), %%st(1) \n"
      "fistpl 88 + 4(%%ebp) \n"
      "movl %%ecx, 88 + 16(%%ebp) \n"
    "s_gfx_hseg23: \n"

    "movl %%ebx, __temp7 \n"
    "movl 0(%%ebp), %%ecx \n"
    "movl 4(%%ebp), %%esp \n"
    "movl 8(%%ebp), %%eax \n"
    "movl 12(%%ebp), %%ebx \n"
    "addl %%ecx, %%eax \n"
    "addl %%esp, %%ebx \n"
    "movl %%eax, 0(%%ebp) \n"
    "movl %%ebx, 4(%%ebp) \n"
    "movl 88 + 0(%%ebp), %%eax \n"
    "movl 88 + 4(%%ebp), %%ebx \n"
    "subl %%ecx, %%eax \n"
    "subl %%esp, %%ebx \n"

    "imull %%esi \n"
    "shll $16, %%edx \n"
    "shrl $16, %%eax \n"
    "addl %%eax, %%edx \n"
    "addl %%edx, %%ecx \n" // ecx = u (16.16)

    "movl %%ebx, %%eax \n"
    "imull %%esi \n"
    "shll $16, %%edx \n"
    "shrl $16, %%eax \n"
    "addl %%eax, %%edx \n"
    "addl %%esp, %%edx \n" // edx = v (16.16)

    "movl %%ecx, %%esp \n"
    "movl %%edx, %%eax \n"
    "addl __temp21, %%esp \n"
    "addl __temp22, %%eax \n"
    "andl $127 << 16, %%esp \n"
    "andl $127 << 16, %%eax \n"
    "shrl $16, %%esp \n"
    "shrl $16 - 7, %%eax \n"
    "addl %%esp, %%eax \n"
    "leal (%%eax, %%eax, 2), %%eax \n"
    "addl __temp19, %%eax \n"

    "shrl $16, %%ecx \n"
//    "movl %%ecx, %%esp \n"
//    "andl $127, %%esp \n"
//    "movl %%edx, %%eax \n"
//    "andl $127 << 16, %%eax \n"
//    "shrl $16 - 7, %%eax \n"
//    "addl %%esp, %%eax \n"
//    "leal (%%eax, %%eax, 2), %%eax \n"
//    "addl __temp19, %%eax \n"

    "andl $63 << 16, %%edx \n"
    "andl $63 << 0, %%ecx \n"
    "shrl $16 - 6, %%edx \n"

    "addl %%ecx, %%edx \n"
    "movl __temp5, %%ebx \n"
    "leal (%%edx, %%edx, 2), %%edx \n"
    "movl __temp3, %%esp \n"

//    "xorl %%ecx, %%ecx \n"

//    "movl (%%eax), %%eax \n"
//    "addl __temp19, %%eax \n"
    
    "movl __temp16, %%ecx \n"
    "movl 0(%%eax), %%esp \n"
    "andl $255, %%esp \n"
    "movb 0(%%ebx, %%edx), %%cl \n"
    "shll $8, %%esp \n"
    "addl %%esp, %%ecx \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 0(%%edi) \n"
    
    "movl __temp17, %%ecx \n"
    "movl 1(%%eax), %%esp \n"
    "andl $255, %%esp \n"
    "movb 1(%%ebx, %%edx), %%cl \n"
    "shll $8, %%esp \n"
    "addl %%esp, %%ecx \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 1(%%edi) \n"
    
    "movl __temp18, %%ecx \n"
    "movl 2(%%eax), %%esp \n"
    "andl $255, %%esp \n"    
    "movb 2(%%ebx, %%edx), %%cl \n"
    "shll $8, %%esp \n"
    "addl %%esp, %%ecx \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 2(%%edi) \n"

    "movl __temp17, %%esp \n"
    "movl %%esp, 16(%%ebp) \n"
    "movl __temp7, %%ebx \n"

    "s_gfx_h_seg21: \n"

    "addl $8, %%ebx \n"
    "addl $3, %%edi \n"
    "addl $88, %%ebp \n"
  "decl __temp2 \n"
  "jge l_gfx_h_seg21 \n"

    "movl __temp, %%esp \n"
    "pop %%ebp \n"
    "pop %%edi \n"
  :
  :"eax" (_temp4), "g" (_bmp.texel)
  :"memory", "ebx", "ecx", "edx", "esi"
  );
}


void draw_gfx_v_seg (long start, long array)
{
  asm volatile (
    "draw_gfx_v_seg_start_1: \n"
    "push %%ebp \n"
    "movl %0, %%edi \n"
    "movl %1, %%ebp \n"
    "movl %%esp, __esp \n"
    "movl $0, %%esi \n"

    "draw_gfx_v_seg_rep_1: \n"
    "  movl 0(%%ebp), %%ecx \n"
    "  movl 4(%%ebp), %%esp \n"
    "  movl 8(%%ebp), %%eax \n"
    "  movl 12(%%ebp), %%ebx \n"
    "  addl %%ecx, %%eax \n"
    "  addl %%esp, %%ebx \n"
    "  movl %%eax, 0(%%ebp) \n"
    "  movl %%ebx, 4(%%ebp) \n"
    "  movl 88 + 0(%%ebp), %%eax \n"
    "  movl 88 + 4(%%ebp), %%ebx \n"
    "  subl %%ecx, %%eax \n"
    "  subl %%esp, %%ebx \n"

    "  imull %%esi \n"
    "  shll $16, %%edx \n"
    "  shrl $16, %%eax \n"
    "  addl %%eax, %%edx \n"
    "  addl %%edx, %%ecx \n"

    "  movl %%ebx, %%eax \n"
    "  imull %%esi \n"
    "  shll $16, %%edx \n"
    "  shrl $16, %%eax \n"
    "  addl %%eax, %%edx \n"
    "  addl %%esp, %%edx \n"

    "  andl $63 << 16, %%ecx \n"
    "  andl $63 << 16, %%edx \n"
    "  shrl $16, %%ecx \n"
    "  shrl $16 - 6, %%edx \n"
    "  addl %%ecx, %%edx \n"
    "  movl __tri_data + 28, %%ebx \n"
    "  leal (%%edx, %%edx, 2), %%edx \n"

    "  movl __tri_data + 16, %%ecx \n"
    "  movb 0(%%ebx, %%edx), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 0(%%edi) \n"
    
    "  movl __tri_data + 20, %%ecx \n"
    "  movb 1(%%ebx, %%edx), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 1(%%edi) \n"
    
    "  movl __tri_data + 24, %%ecx \n"
    "  movb 2(%%ebx, %%edx), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 2(%%edi) \n"
    
    "  addl $3, %%edi \n"
//    "  addl $88, %%ebp \n"
    "decl __tri_data + 0 \n"
    "jge draw_gfx_v_seg_rep_1 \n"

    "movl __esp, %%esp \n"
    "pop %%ebp \n"
    :
    :"g" (start), "g" (array)
    :"eax", "ebx", "ecx", "edx", "esi", "edi"
  );
}


void gfx_v_line_pos (int y, int x1, int x2, double scan_m, int color,
  view_type view)
{
  if (x1 > x2)
    return;

  _tri_data.count = x2 - x1;
//  _tri_data.back_r = _temp15;
//  _tri_data.back_g = _temp16;
//  _tri_data.back_b = _temp17;
  _tri_data.bmp.texel = _bmp.texel;

  sng m = cv_sng(scan_m);
  long s = int(65536 * (double(y + .5) - double(x1 + .5) * scan_m));
  int scrn = (unsigned int)_x2_buff + (unsigned int)3 * ((unsigned int)x1 + (unsigned int)y * (unsigned int)320);
  _temp3 = x1;
  _temp4 = x2 + 1;//int(&_texdat[DAT_MID + s_scan.whole]);
  _temp5 = m.valu;
  _temp6 = s & 65535;
  _temp7 = int(&_texdat[DAT_MID + (s >> 16)]);
  
//  while (_temp3 <= x2) {
//    int scan = s_scan.whole;//(int)floor(f_scan);

    asm (
      "push %%ebp \n"
      "push %%esi \n"
      "movl %%esp, __esp \n"
      "movl %1, %%esi \n"
      "movl %3, %%edi \n"
      "movl %0, %%ebp \n"
      "gfx_v_line_pos_rep1: \n"
      
      // start update u & v
      "movl 16(%%ebp), %%eax \n"
      "movl __temp3, %%edx \n"
      "cmp %%eax, %%edx \n"
      "je gfx_v_line_pos_skip1 \n"
      "  fild __temp3 \n" // load x
      "  fldl __half \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl 0 + 52(%%ebp) \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fldl 0 + 68(%%ebp) \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fstpl __dtemp + 8 * 0 \n"

      "  fldl 0 + 20(%%ebp) \n"
      "  fldl 0 + 36(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 0 + 20(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 0 + 0(%%ebp) \n"

      "  fldl 0 + 28(%%ebp) \n"
      "  fldl 0 + 44(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 0 + 28(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 0 + 4(%%ebp) \n"

      "  fild __temp3 \n" // load x
      "  fldl __half \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl 88 + 52(%%ebp) \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fldl 88 + 68(%%ebp) \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fstpl __dtemp + 8 * 0 \n"

      "  fldl 88 + 20(%%ebp) \n"
      "  fldl 88 + 36(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 88 + 20(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 76(%%ebp) \n"

      "  fldl 88 + 28(%%ebp) \n"
      "  fldl 88 + 44(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 88 + 28(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 80(%%ebp) \n"
      "gfx_v_line_pos_skip1: \n"
      // end update u & v
      
      "movl 0(%%ebp), %%ebx \n"
      "movl 8(%%ebp), %%ecx \n"
      "movl 76(%%ebp), %%eax \n"
      "movl 8+88(%%ebp), %%esp \n"
      "addl %%ebx, %%ecx \n"
      "addl %%eax, %%esp \n"
      "movl %%ecx, 0(%%ebp) \n"
      "movl %%esp, 76(%%ebp) \n"

      "subl %%ebx, %%eax \n"
      "imull %%esi \n"
      "shll $16, %%edx \n"
      "shrl $16, %%eax \n"
      "addl %%eax, %%edx \n"
      "addl %%edx, %%ebx \n"
      "shrl $16, %%ebx \n"
      "andl $63, %%ebx \n"
      "movl %%ebx, __temp1 \n"

      "movl 4(%%ebp), %%ebx \n"
      "movl 12(%%ebp), %%ecx \n"
      "movl 80(%%ebp), %%eax \n"
      "movl 12+88(%%ebp), %%esp \n"
      "addl %%ebx, %%ecx \n"
      "addl %%eax, %%esp \n"
      "movl %%ecx, 4(%%ebp) \n"
      "movl %%esp, 80(%%ebp) \n"

      "subl %%ebx, %%eax \n"
      "imull %%esi \n"
      "shll $16, %%edx \n"
      "shrl $16, %%eax \n"
      "addl %%eax, %%edx \n"
      "addl %%edx, %%ebx \n"
      "andl $63 << 16, %%ebx \n"
      "shrl $16 - 6, %%ebx \n"
      "movl __temp1, %%eax \n"
      "addl %%ebx, %%eax \n"
      "leal (%%eax, %%eax, 2), %%eax \n"
      "movl __tri_data + 28, %%ebx \n"

//      "movb 0(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 0(%%edi) \n"
//      "movb 1(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 1(%%edi) \n"
//      "movb 2(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 2(%%edi) \n"

      "addl $3, %%edi \n"

      "movl __temp3, %%eax \n"
      "inc %%eax \n"
      "movl %%eax, 16(%%ebp) \n"
      "movl %%eax, __temp3 \n"
      
      "subw __temp5, %%si \n"
      "jnc gfx_v_line_pos_skip2 \n"
      "  subl $88, %%ebp \n"
      "gfx_v_line_pos_skip2: \n"

      "cmp __temp4, %%eax \n"
      "jl gfx_v_line_pos_rep1 \n"
      
      "movl __esp, %%esp \n"
      "pop %%esi \n"
      "pop %%ebp \n"
    :
    :"g" (_temp7), "g" (_temp6), "g" (_bmp.texel), "g" (scrn)
    :"eax", "ebx", "ecx", "edx", "edi", "esi", "ebp", "memory"
    );

//    _texdat[DAT_MID + scan].last_ind = x + 1;
//    s -= m.valu;
//    scrn += 3;
//  }

}


void gfx_v_line_neg (int y, int x1, int x2, double scan_m, int color,
  view_type view)
{
  if (x1 > x2)
    return;

  if (_flash) {
    for (int x = x1; x <= x2; x++)
      pxl (x, y, 0, 128, 0);
    return;
  }

  _tri_data.count = x2 - x1;
//  _tri_data.back_r = _temp16;//128 << 8;
//  _tri_data.back_g = _temp17;//128 << 8;
//  _tri_data.back_b = _temp18;//128 << 8;
  _tri_data.bmp.texel = _bmp.texel;

  sng m = cv_sng(-scan_m);
  long s = int(65536 * (double(y + .5) - double(x1 + .5) * scan_m));
  int scrn = (unsigned int)_x2_buff + (unsigned int)3 * ((unsigned int)x1 + (unsigned int)y * (unsigned int)320);
  _temp3 = x1;
  _temp4 = x2 + 1;//int(&_texdat[DAT_MID + s_scan.whole]);
  _temp5 = m.valu;
  _temp6 = s & 65535;
  _temp7 = int(&_texdat[DAT_MID + (s >> 16)]);
  
    asm (
      "push %%ebp \n"
      "push %%esi \n"
      "movl %%esp, __esp \n"
      "movl %1, %%esi \n"
      "movl %3, %%edi \n"
      "movl %0, %%ebp \n"

      // start update u & v
      "gfx_v_line_neg_rep1: \n"
      "movl 16(%%ebp), %%eax \n"
      "movl __temp3, %%edx \n"
      "cmp %%eax, %%edx \n"
      "je gfx_v_line_neg_skip1 \n"
      "  fild __temp3 \n" // load x
      "  fldl __half \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl 0 + 52(%%ebp) \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fldl 0 + 68(%%ebp) \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fstpl __dtemp + 8 * 0 \n"

      "  fldl 0 + 20(%%ebp) \n"
      "  fldl 0 + 36(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 0 + 20(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 0 + 0(%%ebp) \n"

      "  fldl 0 + 28(%%ebp) \n"
      "  fldl 0 + 44(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 0 + 28(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 0 + 4(%%ebp) \n"

      "  fild __temp3 \n" // load x
      "  fldl __half \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl 88 + 52(%%ebp) \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fldl 88 + 68(%%ebp) \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fstpl __dtemp + 8 * 0 \n"

      "  fldl 88 + 20(%%ebp) \n"
      "  fldl 88 + 36(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 88 + 20(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 76(%%ebp) \n"

      "  fldl 88 + 28(%%ebp) \n"
      "  fldl 88 + 44(%%ebp) \n"
      "  fsubp %%st(0), %%st(1) \n"
      "  fldl __dtemp + 0 * 8 \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl 88 + 28(%%ebp) \n"
      "  faddp %%st(0), %%st(1) \n"
      "  fldl __f2sng \n"
      "  fmulp %%st(0), %%st(1) \n"
      "  fldl __half \n"
      "  fsubrp %%st(0), %%st(1) \n"
      "  fistpl 80(%%ebp) \n"
      "gfx_v_line_neg_skip1: \n"
      // end update u & v
      
      "movl 0(%%ebp), %%ebx \n"
      "movl 8(%%ebp), %%ecx \n"
      "movl 76(%%ebp), %%eax \n"
      "movl 8+88(%%ebp), %%esp \n"
      "addl %%ebx, %%ecx \n"
      "addl %%eax, %%esp \n"
      "movl %%ecx, 0(%%ebp) \n"
      "movl %%esp, 76(%%ebp) \n"

      "subl %%ebx, %%eax \n"
      "imull %%esi \n"
      "shll $16, %%edx \n"
      "shrl $16, %%eax \n"
      "addl %%eax, %%edx \n"
      "addl %%edx, %%ebx \n"
      "shrl $16, %%ebx \n"
      "andl $63, %%ebx \n"
      "movl %%ebx, __temp1 \n"

      "movl 4(%%ebp), %%ebx \n"
      "movl 12(%%ebp), %%ecx \n"
      "movl 80(%%ebp), %%eax \n"
      "movl 12+88(%%ebp), %%esp \n"
      "addl %%ebx, %%ecx \n"
      "addl %%eax, %%esp \n"
      "movl %%ecx, 4(%%ebp) \n"
      "movl %%esp, 80(%%ebp) \n"

      "subl %%ebx, %%eax \n"
      "imull %%esi \n"
      "shll $16, %%edx \n"
      "shrl $16, %%eax \n"
      "addl %%eax, %%edx \n"
      "addl %%edx, %%ebx \n"
      "andl $63 << 16, %%ebx \n"
      "shrl $16 - 6, %%ebx \n"
      "movl __temp1, %%eax \n"
      "addl %%ebx, %%eax \n"
      "leal (%%eax, %%eax, 2), %%eax \n"
      "movl __tri_data + 28, %%ebx \n"

    "  movl __tri_data + 16, %%ecx \n"
    "  movb 0(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 0(%%edi) \n"
    
    "  movl __tri_data + 20, %%ecx \n"
    "  movb 1(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 1(%%edi) \n"
    
    "  movl __tri_data + 24, %%ecx \n"
    "  movb 2(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 2(%%edi) \n"

//      "movb 0(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 0(%%edi) \n"
//      "movb 1(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 1(%%edi) \n"
//      "movb 2(%%eax, %%ebx), %%dl \n"
//      "movb %%dl, 2(%%edi) \n"

      "addl $3, %%edi \n"

      "movl __temp3, %%eax \n"
      "inc %%eax \n"
      "movl %%eax, 16(%%ebp) \n"
      "movl %%eax, __temp3 \n"
      
      "addw __temp5, %%si \n"
      "jnc gfx_v_line_neg_skip2 \n"
      "  addl $88, %%ebp \n"
      "gfx_v_line_neg_skip2: \n"

      "cmp __temp4, %%eax \n"
      "jl gfx_v_line_neg_rep1 \n"
      
      "movl __esp, %%esp \n"
      "pop %%esi \n"
      "pop %%ebp \n"
    :
    :"g" (_temp7), "g" (_temp6), "g" (_bmp.texel), "g" (scrn)
    :"eax", "ebx", "ecx", "edx", "edi", "esi", "ebp", "memory"
    );
}


inline void find_n_and_dn
  (double* n, double* dn, double m1, double b1, double m2, double b2)
{
  *dn = 1 / (m2 - m1);
  *n = (b1 - b2) * *dn;
}


inline void find_n_and_dn2
  (double* n, double* dn, double m1, double b1, double m2, double b2)
{
  *n = (b2 * m1 + b1) / (1 - m1 * m2);
  *dn = 1 / (1 - m1 * m2);
}


void calc_scan_m
  (double* dx, double* dy, point_type ref1, point_type p1, point_type p2,
   view_type view)
{
  line_3d l3d;

  calc_line_3d (&l3d, p1.rel, p2.rel);
  double x2 = l3d.mx * ref1.rel.z + l3d.bx;
  double y2 = l3d.my * ref1.rel.z + l3d.by;
  *dx = view.zoom.x * (x2 - ref1.rel.x);
  *dy = view.zoom.y * (y2 - ref1.rel.y);
}


void calc_bmp_edge
  (line_equat* edge, point_type p_near, point_type p_far, view_type view)
{
  if (p_near.rel.z >= view.z_cutoff)
    calc_line (edge, p_near.scr, p_far.scr);
  else {
    point_type p3d;
    p3d.rel = intrapolate_3d (p_near, p_far, view.z_cutoff);
    rel_to_scr (p3d.rel, &p3d.scr, view);
    calc_line (edge, p3d.scr, p_far.scr);
  }
}


inline bool calc_scrn_edge
  (line_equat* line, point_3d p1, point_3d p2, point_2d* pa, point_2d* pb,
   view_type view)
{
  if (p1.z > p2.z)
    if (p2.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, intrapolate_3d (p1, p2, view.z_cutoff)));
    else
      return false;
  else if (p1.z < p2.z)
    if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else if (p2.z > view.z_cutoff)
      calc_2d_dydx_line (line, *pa = map_to_scrn (view, intrapolate_3d (p1, p2, view.z_cutoff)), *pb = map_to_scrn (view, p2));
    else
      return false;
  else
    if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else
      return false;

  return true;
}


inline bool calc_scrn_edge2
  (line_equat* line, point_3d p1, point_3d p2, point_2d* pa, point_2d* pb,
   view_type view)
{
//  calc_2d_dydx_line
//    (line, map_to_scrn (view, p1), map_to_scrn (view, p2));

//  return true;
  
  if (p1.z > p2.z)
    if (p2.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, intrapolate_3d (p1, p2, view.z_cutoff)));
  else if (p1.z < p2.z)
    if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else if (p2.z > view.z_cutoff)
      calc_2d_dydx_line (line, *pa = map_to_scrn (view, intrapolate_3d (p1, p2, view.z_cutoff)), *pb = map_to_scrn (view, p2));
    else
      return false;
  else
    if (p1.z > view.z_cutoff)
      calc_2d_dydx_line
        (line, *pa = map_to_scrn (view, p1), *pb = map_to_scrn (view, p2));
    else
      return false;

//  return true;
}
*/

inline bool calc_horz_bmp_edge
  (line_3d l3d, line_2d* ref_2d, double scan_m, view_type view)
{
  point_2d p2d_1, p2d_2;
  point_3d p3d_1, p3d_2;

  p3d_1.x = l3d.mx * 1 + l3d.bx;
  p3d_1.y = l3d.my * 1 + l3d.by;
  p3d_1.z = 1;
  p3d_2.x = l3d.mx * 2 + l3d.bx;
  p3d_2.y = l3d.my * 2 + l3d.by;
  p3d_2.z = 2;

  p2d_1 = map_to_scrn (view, p3d_1);
  p2d_2 = map_to_scrn (view, p3d_2);

  double dx = p2d_2.x - p2d_1.x;
  double dy = p2d_2.y - p2d_1.y;

  if (is_approx_zero (dx)) {
    ref_2d->m = 0;
    ref_2d->b = p2d_1.x;
    return true;
  }
  else {
    line_equat l2d;
    l2d.m = dy / dx;
    l2d.b = p2d_1.y - p2d_1.x * l2d.m;
    ref_2d->m = 1 / (l2d.m - scan_m);
    ref_2d->b = -l2d.b * ref_2d->m;
    return false;
  }
}


inline bool calc_vert_bmp_edge
  (line_3d l3d, line_2d* ref_2d, double scan_m, view_type view)
{
  point_2d p2d_1, p2d_2;
  point_3d p3d_1, p3d_2;

  p3d_1.x = l3d.mx * 1 + l3d.bx;
  p3d_1.y = l3d.my * 1 + l3d.by;
  p3d_1.z = 1;
  p3d_2.x = l3d.mx * 2 + l3d.bx;
  p3d_2.y = l3d.my * 2 + l3d.by;
  p3d_2.z = 2;

  p2d_1 = map_to_scrn (view, p3d_1);
  p2d_2 = map_to_scrn (view, p3d_2);

  double dx = p2d_2.x - p2d_1.x;
  double dy = p2d_2.y - p2d_1.y;

  if (is_approx_zero (dy)) {
    ref_2d->m = 0;
    ref_2d->b = p2d_1.y;
    return true;
  }
  else {
    line_equat l2d;
    l2d.m = dx / dy;
    l2d.b = p2d_1.x - p2d_1.y * l2d.m;
    ref_2d->m = 1 / (l2d.m - scan_m);
    ref_2d->b = -l2d.b * ref_2d->m;
    return false;
  }
}


uv_acc_type* create_uv_acc()
{
  uv_acc_type* acc_start = (uv_acc_type*) malloc (
    (_screen.y_res + _screen.x_res * 2 + 40) * sizeof(uv_acc_type));
  return acc_start + _screen.x_res + 20;
}


inline bool is_front_side_of_tri (tri_type* tri)
{
  point_2d p1, p2, p3;
  p1.x = tri->t3d.p1->rel.x;// / tri->t3d.p1->rel.z;
  p1.y = tri->t3d.p1->rel.y;// / tri->t3d.p1->rel.z;
  p2.x = tri->t3d.p2->rel.x;// / tri->t3d.p2->rel.z;
  p2.y = tri->t3d.p2->rel.y;// / tri->t3d.p2->rel.z;
  p3.x = tri->t3d.p3->rel.x;// / tri->t3d.p3->rel.z;
  p3.y = tri->t3d.p3->rel.y;// / tri->t3d.p3->rel.z;

  if (approx_equal (p1.y, p2.y))
    if (approx_equal (p1.y, p3.y)) { // 1-2-3
      bprint ("1-2-3", 8);
      return false;
    }
    else if (p3.y > p1.y) { // 1-2,3
      bprint ("1-2,3", 8);
      if (p1.x < p2.x)
        return true;
      else
        return false;
    }
    else { // 3,1-2
      bprint ("3,1-2", 8);
      if (p1.x > p2.x)
        return true;
      else
        return false;
    }
  else if (p1.y > p2.y)
    if (approx_equal (p3.y, p1.y)) { // 2,1-3
      bprint ("2,1-3", 8);
      if (p1.x < p3.x)
        return true;
      else
        return false;
    }
    else if (approx_equal (p3.y, p2.y)) { // 2-3,1
      bprint ("2-3,1", 8);
      if (p2.x < p3.x)
        return true;
      else
        return false;
    }
    else if (p3.y > p1.y) { // 2,1,3
      bprint ("2,1,3", 8);
      return point_left_of_edge
        (tri->t3d.p1->rel, tri->t3d.p2->rel, tri->t3d.p3->rel);
    }
    else if (p3.y < p2.y) { // 3,2,1
      bprint ("3,2,1", 8);
      return point_left_of_edge
        (tri->t3d.p2->rel, tri->t3d.p3->rel, tri->t3d.p1->rel);
    }
    else { // 2,3,1
      bprint ("2,3,1", 8);
      return point_right_of_edge
        (tri->t3d.p3->rel, tri->t3d.p2->rel, tri->t3d.p1->rel);
    }
  else
    if (approx_equal (p3.y, p1.y)) { // 1-3,2
      bprint ("1-3,2", 8);
      if (p3.x < p1.x)
        return true;
      else
        return false;
    }
    else if (approx_equal (p3.y, p2.y)) { // 1,2-3
      bprint ("1,2-3", 8);
      if (p3.x < p2.x)
        return true;
      else
        return false;
    }
    else if (p3.y > p2.y) { // 1,2,3
      bprint ("1,2,3", 8);
      return point_right_of_edge
        (tri->t3d.p2->rel, tri->t3d.p1->rel, tri->t3d.p3->rel);
    }
    else if (p3.y < p1.y) { // 3,1,2
      bprint ("3,1,2", 8);
      return point_right_of_edge
        (tri->t3d.p1->rel, tri->t3d.p3->rel, tri->t3d.p2->rel);
    }
    else { // 1,3,2
      bprint ("1,3,2", 8);
      return point_left_of_edge
        (tri->t3d.p3->rel, tri->t3d.p1->rel, tri->t3d.p2->rel);
    }
}


inline bool point_left_of_edge (point_3d pn, point_3d p1, point_3d p2)
{
  lin_relat l = calc_dxdy_line (p1, p2);
  if (pn.x < l.m * pn.y + l.b)
    return true;
  else
    return false;
}


inline bool point_right_of_edge (point_3d pn, point_3d p1, point_3d p2)
{
  lin_relat l = calc_dxdy_line (p1, p2);
  if (pn.x > l.m * pn.y + l.b)
    return true;
  else
    return false;
}


inline long long curr_time()
{
  return (long long)uclock() / ((long long)UCLOCKS_PER_SEC / CLOCK_GRAN);
}


/*
inline void find_x_and_dx
  (double* n, double* dn, line_equat edge, double scan_m, double scan_b)
{
  if (edge.inv_m != 0) {
    *dn = 1 / (edge.m - scan_m);
    *n = (scan_b - edge.b) * *dn;
  }
  else {
    *dn = 0;
    *n = edge.inv_b;
  }
}


inline void find_y_and_dy
  (double* n, double* dn, line_equat edge, double scan_m, double scan_b)
{
  if (edge.m != 0) {
    *dn = 1 / (edge.inv_m - scan_m);
    *n = (scan_b - edge.inv_b) * *dn;
  }
  else {
    *dn = 0;
    *n = edge.b;
  }
}
*/


#endif //!INCLUDE_TEXTRI
