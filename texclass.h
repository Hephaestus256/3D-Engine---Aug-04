#ifndef INCLUDE_TEXCLASS
#define INCLUDE_TEXCLASS

#include "gfx2d.h"
#include "gfx3d.h"
#include "\programs\include\fixed.h"

struct z_buffer_type {
  unsigned long inv_z;
  int poly;
};

struct tri_3d {
  point_type* p1;
  point_type* p2;
  point_type* p3;
};

struct uv_acc_type {
  sng u, v; // 0, 4
  sng du, dv; // 8, 12
  int last_ind; // 16
  double u1, v1; // 20, 28
  double u2, v2; // 36, 44
  double ind1, ind2; // 52, 60
  double inv_len; // 68
  sng next_u, next_v; // 76, 80
  int init; // 84
}; // 88 bytes long

struct uv_calc_ref {
  uv_acc_type* acc;
  long dest;
  int bk_light[3];
  bmp_type bmp;
  double inv_diff_1;
  double inv_diff_2;
  double scan_m;
  int rot;
  int jump;
  double s1, s2;
  line_2d edge1, edge2;
  lin_relat ref_2d;
  lin_relat ref_3d;
  double orig;
  double center;
  double zoom;
};

struct tex_ref {
  point_type* orig;
  point_type* u;
  point_type* v;
  #ifdef MODE_MAKE
  mk_bmp* bmp;
  #else
  bmp_type* bmp;
  #endif
};

#ifdef MODE_MAKE

struct tri_type {
  tri_3d t3d;
  unsigned long bk_light[3];
  bool select;
  tri_type* next;
};

struct tri_group {
  tri_type* first;
  tex_ref ref;
  plane_type plane;
  int vis_side;
  bool select;  
  tri_group* next;
};

#else

struct tri_type {
  tri_3d t3d;
  unsigned long bk_light[3];
};

struct tri_group {
  tri_type* first;
  tri_type* last;
  tex_ref ref;
  plane_type plane;
  int vis_side;
};

#endif

class tri_group_class {
public:
  z_buffer_type* zmap;
  int group_num;
  tri_type* curr;
  view_type* view;
  void draw_curr_tri();
  void init_tex_plane
    (tri_group group, view_type* v, z_buffer_type* zm, uv_acc_type* acc);
private:
  uv_calc_ref tex_init;
  void tri_group_class::gfx_v_line_pos (int y, int x1, int x2, double scan_m);
  void tri_group_class::gfx_v_line_neg (int y, int x1, int x2, double scan_m);
  void tri_group_class::draw_horz_tex_tri
    (point_2d p1, point_2d p2, point_2d p3);
  void tri_group_class::draw_vert_pos_tex_tri
    (point_2d p1, point_2d p2, point_2d p3);
  void tri_group_class::draw_vert_neg_tex_tri
    (point_2d p1, point_2d p2, point_2d p3);
  void tri_group_class::draw_horz_tex_tri_dither
    (point_2d p1, point_2d p2, point_2d p3);
  void tri_group_class::draw_vert_pos_tex_tri_dither
    (point_2d p1, point_2d p2, point_2d p3);
  void tri_group_class::draw_vert_neg_tex_tri_dither
    (point_2d p1, point_2d p2, point_2d p3);
  inline void tri_group_class::draw_horz_sect
    (double t, double b, lin_relat l, lin_relat r);
  inline void tri_group_class::draw_vert_pos_sect
    (double t, double b, lin_relat l, lin_relat r);
  inline void tri_group_class::draw_vert_neg_sect
    (double t, double b, lin_relat l, lin_relat r);
  inline void tri_group_class::draw_horz_sect_dither
    (double t, double b, lin_relat l, lin_relat r);
  inline void tri_group_class::draw_vert_pos_sect_dither
    (double t, double b, lin_relat l, lin_relat r);
  inline void tri_group_class::draw_vert_neg_sect_dither
    (double t, double b, lin_relat l, lin_relat r);
  void gfx_h_line (int y, int x1, int x2, double scan_m, int color,
    view_type view);
  void gfx_h_seg1 ();
  void gfx_v_line_pos (int y, int x1, int x2, double scan_m, int color,
    view_type view);
  void gfx_v_line_neg (int y, int x1, int x2, double scan_m, int color,
    view_type view);
  inline void init_horz_tri (point_2d p1, point_2d p2, point_2d p3);
  inline void init_vert_tri (point_2d p1, point_2d p2, point_2d p3);
  void draw_tri_dither (point_3d p1, point_3d p2, point_3d p3);
  void draw_tri_no_dither (point_3d p1, point_3d p2, point_3d p3);
  void dither_cut (point_3d p1, point_3d p2, point_3d p3);
  inline point_3d tri_group_class::clip_edge_lr
    (point_3d p1, point_3d p2, double view_m);
  inline point_3d tri_group_class::clip_edge_tb
    (point_3d p1, point_3d p2, double view_m);
  inline void edge_clip_left (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_right (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_top (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_bot (point_3d p1, point_3d p2, point_3d p3);
  void draw_solid_horz_line (int x1, int x2, int y, int c);
  void draw_tex_tri_pos (point_2d p1, point_2d p2, point_2d p3);
  void draw_tex_tri_neg (point_2d p1, point_2d p2, point_2d p3);
  inline void draw_first_sect_pos
    (double top, double bot, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
  inline void draw_left_sect_pos
    (double bot, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
  inline void draw_right_sect_pos
    (double bot, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
  inline void draw_first_sect_neg
    (double bot, double top, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
  inline void draw_left_sect_neg
    (double top, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
  inline void draw_right_sect_neg
    (double top, double* y,
     line_equat left_edge, line_equat right_edge, double* l, double* r);
};

#endif //!TEXCLASS
