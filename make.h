#ifndef INCLUDE_MAKE
#define INCLUDE_MAKE
#ifndef MODE_MAKE
#define MODE_MAKE
#endif
#define CLOCK_GRAN (long long)10000
#include "\programs\zmap\gfx2d.h"

struct mk_bmp {
  char filename[13];
  bmp_type dat;
  mk_bmp* next;
  int n;
};
mk_bmp* _first_mkb;
mk_bmp* _bmp_pal[10];
mk_bmp* curr_bmp;

#include "\programs\zmap\control.h"
#include "\programs\zmap\math2d.h"
#include "\programs\zmap\math3d.h"
#include "\programs\zmap\gfx3d.h"
#include "\programs\zmap\texclass.h"
#include "\programs\zmap\textri.h"
//#include "trisort.h"
#include "ntersect.h"
#include <time.h>
#include <string.h>

struct mk_dat_type {
  int cursor_mode;
  int plane_mode;
  double mouse_rot;
  double view_rad;
  bool angle_snap;
  double angle_snap_grad;
  int new_group;
  tri_group* working_group;
} _mk_dat;

point_type _spec[50];

void show_plane_points (point_type* first, tri_group* g, view_type view);
void delete_point (point_type** first, point_type* del);
void delete_group (tri_group** first, tri_group* del);
void resize_refs (tri_group* first_group);
inline point_type* get_point_in_scope
  (point_type* first_point, view_type view);
void parse_right_tri
  (point_type* p1, point_type* p2, point_type* p3,
   point_type** orig, point_type** pa, point_type** pb);
point_3d extrapolate_quad
  (point_type* p1, point_type* p2, point_type* p3,
   point_type** orig, point_type** pa, point_type** pb);
point_3d extrapolate_quad (point_type* p1, point_type* p2, point_type* p3);
mk_bmp* load_mk_bmp (const char* filename);
double get_valu_buff (char* buff);
void ginput (int* inp_pos, char* keyb_inp);
point_type* find_last_point (point_type* p);
tri_group* find_last_group (tri_group* g);
point_type* read_point (int handle);
tri_type* read_tri (int handle, point_type** point_dat);
tri_group* read_group (int handle, point_type** point_dat, mk_bmp** bmp_dat);
inline void select_point (view_type view, point_type* first_point);
void free_all (point_type** first_point, tri_group** first_group);
tri_type* create_tri (tri_group* g,
  point_type* p1, point_type* p2, point_type* p3, view_type view);
inline void sort_clockwise
  (point_type** p1, point_type** p2, point_type** p3, view_type view);
inline void swap_p3d (point_type** p1, point_type** p2);
void check_controls__make__
  (camera_type* camera, view_type view,
   point_type** first_point, tri_group* first_group);
inline void show_cursor (view_type view);
inline void calc_view (camera_type cursor, camera_type* view);
tri_group* create_tex_ref (tri_group first_group, point_type** first_point,
  bmp_type* curr_bmp, tri_type* tri, view_type view);
dbl_pair perp_to_plane (plane_type plane, tri_type t, point_3d ref);
inline tri_group* find_tri_group
  (tri_group* first_group, tri_type* tri);
inline tri_type* find_last_tri (tri_group g);
void save_scene (const char* filename, ray_type* cursor,
  int mouse_x, int mouse_y, point_type* first_point, tri_group* first_group);
void load_scene (const char* filename, ray_type* cursor,
   int* mouse_x, int* mouse_y,
   point_type** first_point, tri_group** first_group);
int get_point_selection
  (point_type* first_point, int ct, point_type** sel);
point_type* create_point (point_3d abs, point_type** first_point);
point_type* create_corner_point (point_3d abs, point_type** first_point);
point_type* create_ref_point (point_3d abs, point_type** first_point);
void clear_point_selection (point_type* first_point);
void free_all (point_type** first_point, tri_group** first_group);
tri_group* create_group
  (tri_group** first_group, point_type** first_point, tri_type* tri,
  mk_bmp* bmp, point_3d orig, point_3d u, point_3d v, point_3d cam);
mk_bmp* find_last_bmp (mk_bmp* first);
int count_tris (tri_group* g);
inline bool same_point (point_3d p1, point_3d p2);
bool point_is_used (tri_group* first_g, point_type* p);
inline void show_corner_special (point_type corner, view_type view);
point_3d point_follow_cursor
  (camera_type camera, camera_type prev_cam, point_3d p);
void point_follow_cursor
  (camera_type curs, camera_type prev_curs, bool alt, point_type* p);
int all_p_sel_group (tri_group* g);
void get_vis_angle (tri_group group, double* ax, double* ay, double* axy);
dbl_pair perp_to_plane (tri_group group);
int get_vis_side (tri_group prev_group, double dax, double day);
int get_vis_side (tri_group group);
inline int point_grt_line
  (double x0, double y0, double x1, double y1, double x2, double y2);
inline int point_les_line
  (double x0, double y0, double x1, double y1, double x2, double y2);
void update_refs
  (tri_group prev_group, tri_group* curr_group,
   tri_type curr_tri, point_type fixed);
point_type* get_unsel_point (tri_group& group);
void get_rel_angles (plane_type plane, double* xz, double* yz);
void get_rel_angles
  (tri_group& g, point_3d p1, point_3d p2, point_3d p3,
   double* ax, double* ay, double* axy);
bool coplanar (tri_group g);
point_3d* copy_sel_points (point_type* first_point);
bool is_crap (tri_group* g);
void kill_crap (tri_group** first);
void correct_refs (tri_group* g, point_type fixed);
double alt_fall_a (plane_type plane);
void get_trans_angles
  (tri_group& ga, point_3d p1a, point_3d p2a, point_3d p3a,
   tri_group& gb, point_3d p1b, point_3d p2b, point_3d p3b,
   double* axa, double* aya, double* axya,
   double* axb, double* ayb, double* axyb);
void snap_angle (double* a, double grad);


inline void show_cursor (view_type view)
{
  point_2d p;
  p.x = _screen.center.x;
  p.y = _screen.center.y;
 
  if (_mk_dat.cursor_mode) {
    int r = int(7 * (_screen.x_res / 320.0));
    for (int x = -r; x <= r; x++) {
      pxl (p.x + x, p.y - r, 255, 0, 0, view.dest);
      pxl (p.x + x, p.y + r, 255, 0, 0, view.dest);
    }
    for (int y = -r; y <= r; y++) {
      pxl (p.x - r, p.y + y, 255, 0, 0, view.dest);
      pxl (p.x + r, p.y + y, 255, 0, 0, view.dest);
    }
  }
  else {
    int s = int(5 * (_screen.x_res / 320.0));
    for (int y = -s; y < s; y++)
      pxl (p.x, p.y + y, 255, 0, 0, view.dest);
    for (int x = -s; x < s; x++)
      pxl (p.x + x, p.y, 255, 0, 0, view.dest);
  }
}


inline void show_corner (point_type corner, view_type view)
{
  int r1 = int(3 * (_screen.x_res / 320.0));
  int r2 = int(9 * (_screen.x_res / 320.0));
  
  if (corner.visible)
//  if (!corner.ref)
  if (corner.rel.z >= .00001) {
    point_2d p = map_to_scrn (view, corner.rel);
    int r = int(r1 * (view.zoom.x / corner.rel.z));
    int x, y;
    
    if (r < r1)
      r = r1;
    else if (r > r2)
      r = r2;
      
    if (corner.select) {
      for (y = -r, x = -r; y < r; y++, x++)
        pxl (p.x + x, p.y + y, 0, 0, 255, view.dest);
      for (y = -r, x = r; y < r; y++, x--)
        pxl (p.x + x, p.y + y, 0, 0, 255, view.dest);
    }
    else {
      if (!corner.ref) {
        for (y = -r, x = -r; y < r; y++, x++)
          pxl (p.x + x, p.y + y, 255, 255, 255, view.dest);
        for (y = -r, x = r; y < r; y++, x--)
          pxl (p.x + x, p.y + y, 255, 255, 255, view.dest);
      }
      else {
        for (y = -r, x = -r; y < r; y++, x++)
          pxl (p.x + x, p.y + y, 0, 255, 0, view.dest);
        for (y = -r, x = r; y < r; y++, x--)
          pxl (p.x + x, p.y + y, 0, 255, 0, view.dest);
      }
    }

    pxl (p.x + 0, p.y + 0, 0, 255, 255, view.dest);
    pxl (p.x + 1, p.y + 0, 0, 255, 255, view.dest);
    pxl (p.x + 0, p.y + 1, 0, 255, 255, view.dest);
    pxl (p.x + 1, p.y + 1, 0, 255, 255, view.dest);
  }
}


inline void show_corner_special (point_type corner, view_type view)
{
  int r1 = int(3 * (_screen.x_res / 320.0));
  int r2 = int(9 * (_screen.x_res / 320.0));
  
  if (corner.visible)
  if (!corner.ref)
  if (corner.rel.z >= .00001) {
    point_2d p = map_to_scrn (view, corner.rel);
    int r = int(r1 * (view.zoom.x / corner.rel.z));
    int x, y;
    
    if (r < r1)
      r = r1;
    else if (r > r2)
      r = r2;
      
    if (corner.select) {
      for (y = -r, x = -r; y < r; y++, x++)
        pxl (p.x + x, p.y + y, 0, 0, 255, view.dest);
      for (y = -r, x = r; y < r; y++, x--)
        pxl (p.x + x, p.y + y, 0, 0, 255, view.dest);
    }
    else {
      if (!corner.ref) {
        for (y = -r, x = -r; y < r; y++, x++)
          pxl (p.x + x, p.y + y, 255, 0, 0, view.dest);
        for (y = -r, x = r; y < r; y++, x--)
          pxl (p.x + x, p.y + y, 255, 0, 0, view.dest);
      }
      else {
        for (y = -r, x = -r; y < r; y++, x++)
          pxl (p.x + x, p.y + y, 0, 255, 0, view.dest);
        for (y = -r, x = r; y < r; y++, x--)
          pxl (p.x + x, p.y + y, 0, 255, 0, view.dest);
      }
    }

    pxl (p.x + 0, p.y + 0, 0, 255, 255, view.dest);
    pxl (p.x + 1, p.y + 0, 0, 255, 255, view.dest);
    pxl (p.x + 0, p.y + 1, 0, 255, 255, view.dest);
    pxl (p.x + 1, p.y + 1, 0, 255, 255, view.dest);
  }
}


inline point_2d map_abs_point (point_3d p, view_type view)
{
  rotate (&p.x, &p.z, view.camera.pos.x, view.camera.pos.z,
    view.camera.angle_xz);
  rotate (&p.z, &p.y, view.camera.pos.z, view.camera.pos.y,
    view.camera.angle_yz);
  rotate (&p.x, &p.y, view.camera.pos.x, view.camera.pos.y,
    view.camera.angle_xy);
  
  p.x = p.x - view.camera.pos.x;
  p.y = p.y - view.camera.pos.y;
  p.z = p.z - view.camera.pos.z;

  return map_to_scrn (view, p);
}


inline void calc_view (camera_type cursor, camera_type* camera)
{
  if (_mk_dat.cursor_mode & 1) {
    camera->angle_yz = cursor.angle_yz;
    camera->angle_xz = cursor.angle_xz;
    camera->angle_xy = cursor.angle_xy;
    camera->pos = cursor.pos;
  }
  else {
    camera->angle_yz = cursor.angle_yz;
    camera->angle_xz = cursor.angle_xz;
    camera->angle_xy = cursor.angle_xy;
    camera->pos = offset_point_3d (cursor.pos, 0, 0, -_mk_dat.view_rad);
    rotate (&camera->pos.z, &camera->pos.y,
            cursor.pos.z, cursor.pos.y, -cursor.angle_yz);
    rotate (&camera->pos.x, &camera->pos.z,
            cursor.pos.x, cursor.pos.z, -cursor.angle_xz);
  }
}


void cleanup_and_exit()
{
  key_delete();
  set_vesa_mode(0x3);
  exit(0);
}


inline void select_point (view_type view, point_type* first_point)
{
  double r = 7 * (_screen.x_res / 320.0);
  
  for (point_type* p = first_point; p != NULL; p = p->next) {
    point_2d s = map_to_scrn (view, *p);
    if (s.x > _screen.center.x - r)
      if (s.x < _screen.center.x + r)
        if (s.y > _screen.center.y - r)
          if (s.y < _screen.center.y + r)
            if (p->visible)
              if (!p->ref) {
                p->select = !p->select;
                return;
              }
  }
}


inline point_type* get_point_in_scope
  (point_type* first_point, view_type view)
{
  double r = 7 * (_screen.x_res / 320.0);
  
  for (point_type* p = first_point; p != NULL; p = p->next) {
    point_2d s = map_to_scrn (view, *p);
    if (p->visible)
    if (s.x > _screen.center.x - r)
      if (s.x < _screen.center.x + r)
        if (s.y > _screen.center.y - r)
          if (s.y < _screen.center.y + r)
            if (!p->ref)
              return p;
  }

  return NULL;
}


tri_group* create_tex_ref (
  tri_group** first_group, point_type** first_point, mk_bmp* curr_bmp,
  tri_type* tri, view_type view)
{
//#define TEXELS_PER_FOOT 16
//#define FEET_PER_UNIT 1
//#define TEXELS_PER_UNIT TEXELS_PER_FOOT * FEET_PER_UNIT

  plane_type plane;
  double rotat;
  point_3d temp[3];
  double w = curr_bmp->dat.side / TEXELS_PER_UNIT;
  double h = w;
  rotat = -PI / 6;//-5 * (PI / 180);

  point_2d center;
  center.x = (tri->t3d.p1->abs.x);// + tri->t3d.p2->abs.x + tri->t3d.p3->abs.x) / 3;
  center.y = (tri->t3d.p1->abs.y);// + tri->t3d.p2->abs.y + tri->t3d.p3->abs.y) / 3;
  point_2d offset;  
  offset.x = 0;//center.x - tri->t3d.p1->abs.x;
  offset.y = 0;//center.y - tri->t3d.p1->abs.y;

  double ax, ay;
  calc_plane (&plane, tri->t3d.p1->abs, tri->t3d.p2->abs, tri->t3d.p3->abs);
  if (plane.y_plane) {
    ax = 0;
    ay = -PI / 2;
    if (point_relative_to_plane (plane, view.camera.pos) == -1)
      ay += PI;
  }
  else if (plane.m1_inf) {
    ax = -PI / 2;
    ay = -atan(plane.m2);
    if (point_relative_to_plane (plane, view.camera.pos) == 1) {
      ax += PI;
      ay = -ay;
    }
  }
  else if (approx_zero(plane.m1)) {
    ax = 0;
    if (approx_equal(tri->t3d.p1->abs.z, tri->t3d.p2->abs.z))
      if (approx_equal(tri->t3d.p1->abs.y, tri->t3d.p2->abs.y))
        ay = -atan2(tri->t3d.p2->abs.z - tri->t3d.p3->abs.z, tri->t3d.p2->abs.y - tri->t3d.p3->abs.y);
      else
        ay = -atan2(tri->t3d.p2->abs.z - tri->t3d.p1->abs.z, tri->t3d.p2->abs.y - tri->t3d.p1->abs.y);
    else
      ay = -atan2(tri->t3d.p2->abs.z - tri->t3d.p1->abs.z, tri->t3d.p2->abs.y - tri->t3d.p1->abs.y);
    if (point_relative_to_plane (plane, view.camera.pos) == 1) {
      ax += PI;
      ay = -ay;
    }
  }
  else
    if (approx_equal(tri->t3d.p2->abs.x, tri->t3d.p3->abs.x)) {
      ax = atan(plane.m1);
      double m2 = -1 / plane.m1;
      double b2 = tri->t3d.p3->abs.z - m2 * tri->t3d.p3->abs.x;
      x_line_3d ref = calc_x_line_3d (tri->t3d.p1->abs, tri->t3d.p2->abs);
      point_3d pb;
      pb.x = (b2 - ref.bz) / (ref.mz - m2);
      pb.y = ref.my * pb.x + ref.by;
      pb.z = ref.mz * pb.x + ref.bz;

      if (pb.z > tri->t3d.p3->abs.z)
        ay = atan(edge_fall_m (pb, tri->t3d.p3->abs));
      else
        ay = atan(-edge_fall_m (pb, tri->t3d.p3->abs));

      if (point_relative_to_plane (plane, view.camera.pos) == 1) {
        ax += PI;
        ay = -ay;
      }
    }
    else {
      ax = atan(plane.m1);
      double m2 = -1 / plane.m1;
      double b2 = tri->t3d.p1->abs.z - m2 * tri->t3d.p1->abs.x;
      x_line_3d ref = calc_x_line_3d (tri->t3d.p2->abs, tri->t3d.p3->abs);
      point_3d pb;
      pb.x = (b2 - ref.bz) / (ref.mz - m2);
      pb.y = ref.my * pb.x + ref.by;
      pb.z = ref.mz * pb.x + ref.bz;

      if (pb.z > tri->t3d.p1->abs.z)
        ay = atan(edge_fall_m (pb, tri->t3d.p1->abs));
      else
        ay = atan(-edge_fall_m (pb, tri->t3d.p1->abs));

      if (point_relative_to_plane (plane, view.camera.pos) == 1) {
        ax += PI;
        ay = -ay;
      }
    }

  temp[0].x = tri->t3d.p1->abs.x - w / 2 + offset.x;
  temp[0].y = tri->t3d.p1->abs.y - h / 2 + offset.y;
  temp[0].z = tri->t3d.p1->abs.z;
  temp[1].x = tri->t3d.p1->abs.x + w / 2 + offset.x;
  temp[1].y = tri->t3d.p1->abs.y - h / 2 + offset.y;
  temp[1].z = tri->t3d.p1->abs.z;
  temp[2].x = tri->t3d.p1->abs.x - w / 2 + offset.x;
  temp[2].y = tri->t3d.p1->abs.y + h / 2 + offset.y;
  temp[2].z = tri->t3d.p1->abs.z;

  rotate (&temp[0].x, &temp[0].y, center.x, center.y, rotat);
  rotate (&temp[1].x, &temp[1].y, center.x, center.y, rotat);
  rotate (&temp[2].x, &temp[2].y, center.x, center.y, rotat);
  
  rotate (&temp[0].z, &temp[0].y, tri->t3d.p1->abs.z, tri->t3d.p1->abs.y, ay);
  rotate (&temp[0].x, &temp[0].z, tri->t3d.p1->abs.x, tri->t3d.p1->abs.z, ax);
  rotate (&temp[1].z, &temp[1].y, tri->t3d.p1->abs.z, tri->t3d.p1->abs.y, ay);
  rotate (&temp[1].x, &temp[1].z, tri->t3d.p1->abs.x, tri->t3d.p1->abs.z, ax);
  rotate (&temp[2].z, &temp[2].y, tri->t3d.p1->abs.z, tri->t3d.p1->abs.y, ay);
  rotate (&temp[2].x, &temp[2].z, tri->t3d.p1->abs.x, tri->t3d.p1->abs.z, ax);

  return create_group (first_group, first_point, tri, curr_bmp,
    temp[0], temp[1], temp[2], view.camera.pos);
}


tri_type* create_tri (tri_group* g,
  point_type* p1, point_type* p2, point_type* p3, view_type view)
{
  tri_type* t = (tri_type*)malloc(sizeof(tri_type));

  if (g != NULL) {
    tri_type* prev_last = find_last_tri (*g);
    prev_last->next = t;
  }

  t->bk_light[0] = 256 << 8;
  t->bk_light[1] = 256 << 8;
  t->bk_light[2] = 256 << 8;
  t->select = false;
  t->next = NULL;

  point_2d center, s1, s2, s3;

  rel_to_scr(p1->rel, &s1, view);
  rel_to_scr(p2->rel, &s2, view);
  rel_to_scr(p3->rel, &s3, view);
  
  center.x = (s1.x + s2.x + s3.x) / 3;
  center.y = (s1.y + s2.y + s3.y) / 3;

  double a1 = atan2(center.y - s1.y, s1.x - center.x);
  double a2 = atan2(center.y - s2.y, s2.x - center.x);
  double a3 = atan2(center.y - s3.y, s3.x - center.x);

  t->t3d.p1 = p1;
  t->t3d.p2 = p2;
  t->t3d.p3 = p3;
  
  if (!is_front_side_of_tri (t)) {
    t->t3d.p1 = p3;
    t->t3d.p2 = p2;
    t->t3d.p3 = p1;
  }
  
/*
  if (a1 < a2)
    if (a3 > a2) {
//      printf ("3, 2, 1"); // reverse
      t->t3d.p1 = p3;
      t->t3d.p2 = p2;
      t->t3d.p3 = p1;
    }
    else if (a3 < a1) {
//      printf ("2, 1, 3"); // reverse
      t->t3d.p1 = p3;
      t->t3d.p2 = p2;
      t->t3d.p3 = p1;
    }
    else {
//      printf ("2, 3, 1"); // forward
      t->t3d.p1 = p1;
      t->t3d.p2 = p2;
      t->t3d.p3 = p3;
    }
  else
    if (a3 > a1) {
//      printf ("3, 1, 2"); // forward
      t->t3d.p1 = p1;
      t->t3d.p2 = p2;
      t->t3d.p3 = p3;
    }
    else if (a3 < a2) {
//      printf ("1, 2, 3"); // forward
      t->t3d.p1 = p1;
      t->t3d.p2 = p2;
      t->t3d.p3 = p3;
    }
    else {
//      printf ("1, 3, 2"); // reverse
      t->t3d.p1 = p3;
      t->t3d.p2 = p2;
      t->t3d.p3 = p1;
    }
*/
  return t;
}


inline void swap_p3d (point_type** p1, point_type** p2)
{
  point_type* temp = *p1;
  **p1 = **p2;
  **p2 = *temp;
}


void check_controls__make__
  (camera_type* camera, view_type view,
  point_type** first_point, tri_group** first_group)
{
  double rec_inc = 2;
  double snap_grad = 15;
  static bool type_point = false;
  static char keyb_inp[20] = "";
  static int inp_pos = 0;
  static int inp_n = 0;
  static point_3d new_pt;
  static double offset_xz = 0;
  static double offset_yz = 0;
  static int prev_mx = 512, prev_my = 512;
  static camera_type prev_curs;
  static bool move_mode = false;
  static int move_ct;
  static point_3d* move_p3d;
  static bool move_alt;
  prev_curs = *camera;
  
  get_mouse_stat_game();
  
  if (!type_point) {
    if (_keyb_stat[0x1D].press) {
      tri_type* t = ray_int_tri (view.camera, *first_group, NULL);
      if (t != NULL) {
        tri_group* g = find_tri_group (*first_group, t);
        for (int k = 0x2; k <= 0xA; k++)
          if (_keyb_stat[k].change) {
            g->ref.bmp = _bmp_pal[k - 1];
            resize_refs (g);
          }
        if (_keyb_stat[0xB].change) {
          g->ref.bmp = _bmp_pal[0];
          resize_refs (g);
        }
      }
    }
    else {
      for (int k = 0x2; k <= 0xA; k++)
        if (_keyb_stat[k].change)
          curr_bmp = _bmp_pal[k - 1];
      if (_keyb_stat[0xB].change)
        curr_bmp = _bmp_pal[0];
    }

  if (_mk_dat.plane_mode) {
//    camera->angle_xy = _mouse_status.x;
    if (_mk_dat.angle_snap)
      camera->angle_xy = deg_to_rad (snap_grad * int(rad_to_deg (
        camera->angle_xy + deg_to_rad(snap_grad / 2)) / snap_grad));
    if (_keyb_stat[0x16].press)
      if (_mk_dat.new_group)
        if (_keyb_stat[0x2A].press) {
          point_type* sel[10];
          int sel_ct = get_point_selection (*first_point, 3, sel);
          if (sel_ct == 3)
            create_corner_point
              (extrapolate_quad (sel[0], sel[1], sel[2]), first_point);
        }
        else {
          point_type* sel[10];
          int sel_ct = get_point_selection (*first_point, 3, sel);
          point_type* orig;
          point_type* pa;
          point_type* pb;
          if (sel_ct == 3) {
            _mk_dat.new_group = false;
            calc_view (*camera, &view.camera);
            point_type* pn = create_corner_point
              (extrapolate_quad (sel[0], sel[1], sel[2], &orig, &pa, &pb),
              first_point);
            calc_rel_point (view.camera.pos, view.camera.angle_xz,
              view.camera.angle_yz, view.camera.angle_xy, pn);

            tri_type* t1 = create_tri (NULL, orig, pa, pb, view);
            tri_group* g =
              create_tex_ref (first_group, first_point, curr_bmp, t1, view);
            create_tri (g, pn, pa, pb, view);
            clear_point_selection (*first_point);
          }
        }
      else
        if (_keyb_stat[0x2A].press) {
          point_type* sel[10];
          int sel_ct = get_point_selection (*first_point, 3, sel);
          if (sel_ct == 3)
            create_corner_point
              (extrapolate_quad (sel[0], sel[1], sel[2]), first_point);
        }
        else {
          point_type* sel[10];
          int sel_ct = get_point_selection (*first_point, 3, sel);
          point_type* orig;
          point_type* pa;
          point_type* pb;
          if (sel_ct == 3) {
            calc_view (*camera, &view.camera);
            point_type* pn = create_corner_point
              (extrapolate_quad (sel[0], sel[1], sel[2], &orig, &pa, &pb),
              first_point);
            calc_rel_point (view.camera.pos, view.camera.angle_xz,
              view.camera.angle_yz, view.camera.angle_xy, pn);
            create_tri (_mk_dat.working_group, orig, pa, pb, view);
            create_tri (_mk_dat.working_group, pn, pa, pb, view);
            clear_point_selection (*first_point);
          }
        }
    if (_keyb_stat[0x1D].press) {
      if (_keyb_stat['K'].press) {
        camera->pos.x -= sin(camera->angle_xz + PI / 2) * rec_inc / 12.0;
        camera->pos.z -= cos(camera->angle_xz + PI / 2) * rec_inc / 12.0;
      }
      if (_keyb_stat['M'].press) {
        camera->pos.x += sin(camera->angle_xz + PI / 2) * rec_inc / 12.0;
        camera->pos.z += cos(camera->angle_xz + PI / 2) * rec_inc / 12.0;
      }
      if (_keyb_stat[72].press) {
        camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        camera->pos.y -= sin(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
      }
      if (_keyb_stat[80].press) {
        camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        camera->pos.y += sin(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
      }
    }
    else {
      if (_keyb_stat['K'].press) {
        camera->pos.x -= sin(camera->angle_xz + PI / 2) * rec_inc;
        camera->pos.z -= cos(camera->angle_xz + PI / 2) * rec_inc;
      }
      if (_keyb_stat['M'].press) {
        camera->pos.x += sin(camera->angle_xz + PI / 2) * rec_inc;
        camera->pos.z += cos(camera->angle_xz + PI / 2) * rec_inc;
      }
      if (_keyb_stat[72].press) {
        camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
        camera->pos.y -= sin(camera->angle_yz + PI / 2) * rec_inc;
        camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
      }
      if (_keyb_stat[80].press) {
        camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
        camera->pos.y += sin(camera->angle_yz + PI / 2) * rec_inc;
        camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
      }
    }
    if (_keyb_stat[0x14].change)
      if (_mk_dat.new_group) {
        point_type* sel[10];
        int pct = get_point_selection (*first_point, 3, sel);
        if (pct == 3) {
          tri_type* t = create_tri (NULL, sel[0], sel[1], sel[2], view);
          _mk_dat.new_group = false;
          _mk_dat.working_group =
            create_tex_ref (first_group, first_point, curr_bmp, t, view);
          clear_point_selection (*first_point);
        }
      }
      else {
        point_type* sel[10];
        int pct = get_point_selection (*first_point, 3, sel);
        if (pct == 3) {
          create_tri (_mk_dat.working_group, sel[0], sel[1], sel[2], view);
          clear_point_selection (*first_point);
        }
      }
  }
  else {
    if (_keyb_stat[0x18].change) {
      move_mode = !move_mode;
      if (_keyb_stat[0x36].press)
        move_alt = true;
      else
        move_alt = false;
    }

    if (_keyb_stat[0x1D].press) {
      camera->angle_xz = _mouse_status.x;
      camera->angle_yz = -_mouse_status.y;
    }
    else {
      camera->angle_xz = _mouse_status.x;
      camera->angle_yz = -_mouse_status.y;
    }

    if (_mk_dat.angle_snap) {
      snap_angle (&camera->angle_xz, snap_grad);
      snap_angle (&camera->angle_yz, snap_grad);
//      camera->angle_xz = deg_to_rad (snap_grad * int(rad_to_deg (
//        camera->angle_xz + deg_to_rad(snap_grad / 2)) / snap_grad));
//      camera->angle_yz = deg_to_rad (snap_grad * int(rad_to_deg (
//        camera->angle_yz + deg_to_rad(snap_grad / 2)) / snap_grad));
    }

    if (_keyb_stat[0x38].press) {
      if (_keyb_stat['K'].change)
        offset_xz -= deg_to_rad(15);
      if (_keyb_stat['M'].change)
        offset_xz += deg_to_rad(15);
      if (_keyb_stat[72].change)
        offset_yz += deg_to_rad(15);
      if (_keyb_stat[80].change)
        offset_yz -= deg_to_rad(15);
    }
    else {
      if (_keyb_stat[0x1D].press) {
        if (_keyb_stat['K'].press) {
          camera->pos.x -= sin(camera->angle_xz + PI / 2) * rec_inc / 12.0;
          camera->pos.z -= cos(camera->angle_xz + PI / 2) * rec_inc / 12.0;
        }
        if (_keyb_stat['M'].press) {
          camera->pos.x += sin(camera->angle_xz + PI / 2) * rec_inc / 12.0;
          camera->pos.z += cos(camera->angle_xz + PI / 2) * rec_inc / 12.0;
        }
        if (_keyb_stat[72].press) {
          camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
          camera->pos.y -= sin(camera->angle_yz + PI / 2) * rec_inc / 12.0;
          camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        }
        if (_keyb_stat[80].press) {
          camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
          camera->pos.y += sin(camera->angle_yz + PI / 2) * rec_inc / 12.0;
          camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc / 12.0;
        }
      }
      else {
        if (_keyb_stat['K'].press) {
          camera->pos.x -= sin(camera->angle_xz + PI / 2) * rec_inc;
          camera->pos.z -= cos(camera->angle_xz + PI / 2) * rec_inc;
        }
        if (_keyb_stat['M'].press) {
          camera->pos.x += sin(camera->angle_xz + PI / 2) * rec_inc;
          camera->pos.z += cos(camera->angle_xz + PI / 2) * rec_inc;
        }
        if (_keyb_stat[72].press) {
          camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
          camera->pos.y -= sin(camera->angle_yz + PI / 2) * rec_inc;
          camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
        }
        if (_keyb_stat[80].press) {
          camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
          camera->pos.y += sin(camera->angle_yz + PI / 2) * rec_inc;
          camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz + PI / 2) * rec_inc;
        }
      }
    }
    if (_keyb_stat[16].press) {
      camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz) * rec_inc;
      camera->pos.y -= sin(camera->angle_yz) * rec_inc;
      camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz) * rec_inc;
    }
    if (_keyb_stat[30].press) {
      camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz) * rec_inc;
      camera->pos.y += sin(camera->angle_yz) * rec_inc;
      camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz) * rec_inc;
    }
    if (_keyb_stat[0x14].change) {
      point_type* sel[10];
      get_point_selection (*first_point, 3, sel);
      tri_type* t = create_tri (NULL, sel[0], sel[1], sel[2], view);
      create_tex_ref (first_group, first_point, curr_bmp, t, view);
      clear_point_selection (*first_point);
    }
    if (_keyb_stat[0x16].change)
      if (_keyb_stat[0x2A].press) {
        point_type* sel[10];
        int sel_ct = get_point_selection (*first_point, 3, sel);
        if (sel_ct == 3)
          create_corner_point
            (extrapolate_quad (sel[0], sel[1], sel[2]), first_point);
      }
      else {
        point_type* sel[10];
        int sel_ct = get_point_selection (*first_point, 3, sel);
        point_type* orig;
        point_type* pa;
        point_type* pb;
        if (sel_ct == 3) {
          calc_view (*camera, &view.camera);
          point_type* pn = create_corner_point
            (extrapolate_quad (sel[0], sel[1], sel[2], &orig, &pa, &pb),
            first_point);
          calc_rel_point (view.camera.pos, view.camera.angle_xz,
            view.camera.angle_yz, view.camera.angle_xy, pn);

          tri_type* t1 = create_tri (NULL, orig, pa, pb, view);
          tri_group* g =
            create_tex_ref (first_group, first_point, curr_bmp, t1, view);
          create_tri (g, pn, pa, pb, view);
          clear_point_selection (*first_point);
        }
        else
          bprint ("ct: ", sel_ct, 9);
      }
  }
//***************************************************************************
  if (_keyb_stat[0x3F].change)
    _collis = !_collis;
  if (_keyb_stat[0x13].change) {
    camera->pos.x = double(int(camera->pos.x));
    camera->pos.y = double(int(camera->pos.y));
    camera->pos.z = double(int(camera->pos.z));
  }
  if (_keyb_stat[0x2E].change)
    for (point_type* p = *first_point; p != NULL; p = p->next)
      if (!point_is_used (*first_group, p))
        delete_point (first_point, p);
  if (_keyb_stat[0x53].change) {
    for (point_type* p = *first_point; p != NULL; p = p->next)
      if (p->select)
        if (!point_is_used (*first_group, p))
          delete_point (first_point, p);
    for (tri_group* g = *first_group; g != NULL; g = g->next)
      if (g->select)
        delete_group (first_group, g);
  }

  if (_keyb_stat[0x32].change) {
    //point_type* sel[10];
    //get_point_selection (*first_point, 1, sel);
    point_type* p = get_point_in_scope (*first_point, view);
    if (p != NULL) {
      camera->pos = p->abs;
      _mk_dat.cursor_mode = false;
    }
  }

  if (_keyb_stat[0x26].press) {
    tri_type* i = ray_int_tri (view.camera, *first_group, NULL);
    if (i != NULL) {
      tri_group* g = find_tri_group (*first_group, i);
      if (_keyb_stat[0x0C].press)
        for (tri_type* t = g->first; t != NULL; t = t->next) {
          t->bk_light[0] -= 1 << 8;
          t->bk_light[1] -= 1 << 8;
          t->bk_light[2] -= 1 << 8;
        }
      if (_keyb_stat[0x0D].press)
        for (tri_type* t = g->first; t != NULL; t = t->next) {
          t->bk_light[0] += 1 << 8;
          t->bk_light[1] += 1 << 8;
          t->bk_light[2] += 1 << 8;
        }
    }
  }
  else {
    if (_keyb_stat[0x0C].press)
    _mk_dat.view_rad -= 2.5;
    if (_keyb_stat[0x0D].press)
      _mk_dat.view_rad += 2.5;
    if (_keyb_stat[0x0E].press)
      _mk_dat.view_rad = view.dither_dist / DITHER_ZOOM;
  }
  if (_keyb_stat[0x2C].change) {
    offset_xz = offset_yz = 0;
    set_mouse_coord (MOUSE_MAX / 2, MOUSE_MAX / 2);
  }
  
//  if (_keyb_stat[0x29].press)
//    type_point = true;
  if (_keyb_stat[0x31].press) {
    set_mouse_coord (MOUSE_MAX / 2, MOUSE_MAX / 2);
    camera->pos.x = 0;
    camera->pos.y = 0;
    camera->pos.z = 0;
    free_all(first_point, first_group);
    _bmp_pal[0] = load_mk_bmp ("brick1.gfx");
    _bmp_pal[1] = load_mk_bmp ("tex2.gfx");
    _bmp_pal[2] = load_mk_bmp ("plank1.gfx");
    _bmp_pal[3] = load_mk_bmp ("wood1.gfx");
    _bmp_pal[4] = load_mk_bmp ("wood2.gfx");
    _bmp_pal[5] = load_mk_bmp ("concret1.gfx");
    _bmp_pal[6] = load_mk_bmp ("concret2.gfx");
    _bmp_pal[7] = load_mk_bmp ("concret3.gfx");
    _bmp_pal[8] = load_mk_bmp ("bbrick.gfx");
    _bmp_pal[9] = load_mk_bmp ("xzplane.gfx");
    curr_bmp = _first_mkb;
  }
  if (_keyb_stat[0x3C].change) {
    int mx, my;
    get_mouse_stat_raw (&mx, &my, NULL);
    save_scene("3dfx.dat", camera, mx, my, *first_point, *first_group);
  }
  if (_keyb_stat[0x3D].change) {
    int mx, my;  
    load_scene("3dfx.dat", camera, &mx, &my, first_point, first_group);
    set_mouse_coord (mx, my);
  }
  if (_keyb_stat[0x12].change) {
    point_3d pt;
    set_mouse_coord (MOUSE_MAX / 2, MOUSE_MAX / 2);
    tri_type* t = ray_int_tri (view.camera, *first_group, &pt);
    if (t != NULL) {
      clear_point_selection (*first_point);
      tri_group* g = find_tri_group (*first_group, t);
      plane_type plane;
      calc_plane (&plane,
        (t->t3d.p1->abs), (t->t3d.p2->abs), (t->t3d.p3->abs));
      dbl_pair a = perp_to_plane (*g);//plane, *t, view.camera.pos);
      camera->pos = pt;
      camera->angle_xz = a.x;
      camera->angle_yz = a.y;
      _mk_dat.plane_mode = true;
      _mk_dat.new_group = false;
      _mk_dat.working_group = g;
      show_plane_points (*first_point, g, view);
    }
  }
  if (_keyb_stat[0x21].change) {
    point_3d pt;
    tri_type* t = ray_int_tri (view.camera, *first_group, &pt);
    if (t != NULL) {
      plane_type plane;
      calc_plane (&plane,
        (t->t3d.p1->abs), (t->t3d.p2->abs), (t->t3d.p3->abs));
      dbl_pair a = perp_to_plane (plane, *t, view.camera.pos);
      camera->pos = pt;
      camera->angle_xz = a.x;
      camera->angle_yz = a.y;
//      _mk_dat.plane_mode = true;
    }
  }
  if (_keyb_stat[0x17].change) {
    point_3d pt;
    tri_type* t = ray_int_tri (view.camera, *first_group, &pt);
    if (t != NULL)
      camera->pos = pt;
  }
  if (_keyb_stat[0x1F].change) {
    _mk_dat.angle_snap = !_mk_dat.angle_snap;
  }
  if (_keyb_stat[0x0F].change)
    if (_mk_dat.cursor_mode) {
      _mk_dat.cursor_mode = false;
      camera->pos.x += sin(camera->angle_xz) * cos(camera->angle_yz) *
        _mk_dat.view_rad;
      camera->pos.y -= sin(camera->angle_yz) * _mk_dat.view_rad;
      camera->pos.z += cos(camera->angle_xz) * cos(camera->angle_yz) *
        _mk_dat.view_rad;
    }
    else {
      _mk_dat.cursor_mode = true;
      camera->pos.x -= sin(camera->angle_xz) * cos(camera->angle_yz) *
        _mk_dat.view_rad;
      camera->pos.y += sin(camera->angle_yz) * _mk_dat.view_rad;
      camera->pos.z -= cos(camera->angle_xz) * cos(camera->angle_yz) *
        _mk_dat.view_rad;
    }
  if (_keyb_stat[0x01].change) {
    clear_point_selection (*first_point);
    for (tri_group* g = *first_group; g != NULL; g = g->next)
      for (tri_type* t = g->first; t != NULL; t = t->next)
        t->select = false;
  }
  else if (_keyb_stat[0x19].change) { // P
    for (point_type* p = *first_point; p != NULL; p = p->next)
      p->visible = true;
    _mk_dat.plane_mode = !_mk_dat.plane_mode;
    _mk_dat.new_group = true;
    camera->angle_xy = 0;
  }
  if (_keyb_stat[0x2D].change) {
    free_all (first_point, first_group);
    cleanup_and_exit();
  }
  if (_keyb_stat[0x1C].change)
    if (_mk_dat.cursor_mode)
      select_point (view, *first_point);
    else
      create_corner_point (camera->pos, first_point);
  if (_mouse_status.button[0].change)
    if (_mk_dat.cursor_mode)
      select_point (view, *first_point);
    else
      create_corner_point (camera->pos, first_point);

  if (_mouse_status.button[1].change) {
    tri_type* t = ray_int_tri (view.camera, *first_group, NULL);
    tri_group* g = find_tri_group (*first_group, t);
    if (g != NULL) {
      g->select = !g->select;
      for (tri_type* tn = g->first; tn != NULL; tn = tn->next)
        tn->select = !tn->select;
    }
  }
  }
  
//    grab:
  if (move_mode) {
    point_3d* p_backup = copy_sel_points (*first_point);
    bool all_coplanar = true;
    
    for (point_type* p = *first_point; p != NULL; p = p->next)
      if (p->select)
        point_follow_cursor (*camera, prev_curs, move_alt, p);

    for (tri_group* g = *first_group; g != NULL; g = g->next)
      if (!coplanar(*g))
        all_coplanar = false;
        
    if (all_coplanar)
      for (tri_group* g = *first_group; g != NULL; g = g->next) {
        int sct = all_p_sel_group (g);
        if (sct == 2) {
          point_follow_cursor (*camera, prev_curs, move_alt, g->ref.u);
          point_follow_cursor (*camera, prev_curs, move_alt, g->ref.v);
          point_follow_cursor (*camera, prev_curs, move_alt, g->ref.orig);
          calc_plane (&g->plane, g->ref.u->abs, g->ref.v->abs, g->ref.orig->abs);
          g->vis_side = get_vis_side (g->first->t3d.p1->abs, g->first->t3d.p2->abs, g->first->t3d.p3->abs);
          //get_vis_side (*g);
        }
        else if (sct == 1) {
          point_type* fixed = get_unsel_point(*g);
          tri_group pg = *g;
          calc_plane (&g->plane, g->first->t3d.p1->abs, g->first->t3d.p2->abs, g->first->t3d.p3->abs);
          update_refs (pg, g, *g->first, *fixed);
          calc_plane (&g->plane, g->ref.u->abs, g->ref.v->abs, g->ref.orig->abs);
          g->vis_side = get_vis_side (g->first->t3d.p1->abs, g->first->t3d.p2->abs, g->first->t3d.p3->abs);
          //get_vis_side (*g);
        }
      }
    else {
      int pct = 0;
      for (point_type* p = *first_point; p != NULL; p = p->next)
        if (p->select) {
          p->abs = p_backup[pct];
          pct++;
        }
    }
    free(p_backup);
  }

  for (int i = 0; i < 128; i++)
    _keyb_stat[i].change = false;
}


dbl_pair perp_to_plane (tri_group group)
{
  dbl_pair a;
  double ax, ay, axy;
  
  get_vis_angle (group, &ax, &ay, &axy);
  if (group.plane.y_plane)
    a.y = simplify_angle (ay + PI);
  else if (group.plane.m1_inf) {
    a.x = simplify_angle (PI / 2 - ax + PI);
    a.y = simplify_angle (axy);
  }
  else {
    a.x = simplify_angle (PI / 2 - ax + PI);
    a.y = atan(cos(ax) * tan(ay));
  }
  
  return a;
}


dbl_pair perp_to_plane (plane_type plane, tri_type t, point_3d ref)
{
  dbl_pair a;
  point_3d p1 = t.t3d.p1->abs;
  point_3d p2 = t.t3d.p2->abs;
  point_3d p3 = t.t3d.p3->abs;
  
  if (plane.y_plane) {
    a.x = 0;
    if (ref.y > plane.b)
      a.y = PI / 2;
    else
      a.y = -PI / 2;
  }
  else if (plane.m1_inf)
    if (ref.x > plane.m2 * ref.y + plane.b) {
      a.x = -PI / 2;
      a.y = -atan(plane.m2);
    }
    else {
      a.x = PI / 2;
      a.y = atan(plane.m2);      
    }
  else {
    if (ref.z < plane.b)
      a.x = -atan(plane.m1);
    else
      a.x = PI - atan(plane.m1);
    if (is_approx_zero (plane.m1))
      a.y = atan(plane.m2);
    else {
      double m2 = -1 / plane.m1;
      double b2 = p1.z - m2 * p1.x;
      x_line_3d ref = calc_x_line_3d (p2, p3);
      point_3d pb;
      pb.x = (b2 - ref.bz) / (ref.mz - m2);
      pb.y = ref.my * pb.x + ref.by;
      pb.z = ref.mz * pb.x + ref.bz;

      if (pb.z > p1.z)
        a.y = -PI / 2 + atan(1 / edge_fall_m (pb, p1));
      else
        a.y = -PI / 2 + atan(-1 / edge_fall_m (pb, p1));
    }
  }

  return a;
}


inline tri_group* find_tri_group (tri_group* first_group, tri_type* tri)
{
  for (tri_group* g = first_group; g != NULL; g = g->next)
    for (tri_type* t = g->first; t != NULL; t = t->next)
      if (t == tri)
        return g;

  return NULL;
}


inline tri_type* find_last_tri (tri_group g)
{
  for (tri_type* t = g.first; t != NULL; t = t->next)
    if (t->next == NULL)
      return t;
}


int get_point_selection
  (point_type* first_point, int ct, point_type** sel)
{
  int n = 0;
  
  for (point_type* i = first_point; i != NULL; i = i->next)
    if (i->visible)
      if (i->select) {
        if (n < ct)
          sel[n] = i;
        n++;
      }

  return n;
}


point_type* create_point (point_3d abs, point_type** first_point)
{
  point_type* temp = (point_type*)malloc(sizeof(point_type));
    
  temp->abs = abs;
  temp->next = NULL;
  temp->visible = true;
  
  if (*first_point == NULL)
    *first_point = temp;
  else {
    point_type* last = find_last_point (*first_point);
    last->next = temp;
  }

  return temp;
}


inline point_type* create_corner_point
  (point_3d abs, point_type** first_point)
{
  for (point_type* p = *first_point; p != NULL; p = p->next)
    if (!p->ref)
      if (same_point (p->abs, abs)) {
        p->select = true;
        p->visible = true;
        return p;
      }

  point_type* temp = create_point (abs, first_point);
  temp->select = true;
  temp->ref = false;

  return temp;
}


point_type* create_ref_point (point_3d abs, point_type** first_point)
{
  point_type* temp = create_point (abs, first_point);
  temp->select = false;
  temp->ref = true;

  return temp;
}


void clear_point_selection (point_type* first_point)
{
  for (point_type* p = first_point; p != NULL; p = p->next)
    p->select = false;
}


void save_scene (const char* filename, ray_type* cursor,
   int mouse_x, int mouse_y, point_type* first_point, tri_group* first_group)
{
  int handle;
  int point_ct = 0;
  int group_ct = 0;
  int bmp_ct = 0;
  
  for (point_type* p = first_point; p != NULL; p = p->next) {
    p->n = point_ct;
    point_ct++;    
  }

  for (tri_group* g = first_group; g != NULL; g = g->next)
    group_ct++;

  for (mk_bmp* b = _first_mkb; b != NULL; b = b->next, bmp_ct++)
    b->n = bmp_ct;

  if (_dos_open(filename, 0x02, &handle))
    exit(5);

  int version = 0;
  _dos_write(handle, &version, sizeof(int), NULL);
  _dos_write(handle, &_mk_dat, sizeof(mk_dat_type), NULL);
  _dos_write(handle, cursor, sizeof(ray_type), NULL);
  _dos_write(handle, &mouse_x, sizeof(int), NULL);
  _dos_write(handle, &mouse_y, sizeof(int), NULL);
  _dos_write(handle, &point_ct, sizeof(int), NULL);
  _dos_write(handle, &group_ct, sizeof(int), NULL);
  _dos_write(handle, &bmp_ct, sizeof(int), NULL);

  // Write bmps
  for (mk_bmp* b = _first_mkb; b != NULL; b = b->next)
    _dos_write(handle, b->filename, 13, NULL);
  // Write working texture
  _dos_write(handle, &curr_bmp->n, sizeof(int), NULL);
  // Write bmp pallette
  for (int i = 0; i < 10; i++)
    _dos_write(handle, &_bmp_pal[i]->n, sizeof(int), NULL);
  // Write points
  for (point_type* p = first_point; p != NULL; p = p->next)
    _dos_write(handle, p, sizeof(point_type), NULL);

  for (tri_group* g = first_group; g != NULL; g = g->next) {
    int tri_ct = 0;
    for (tri_type* t = g->first; t != NULL; t = t->next)
      tri_ct++;
    _dos_write(handle, &tri_ct, sizeof(int), NULL);
    _dos_write(handle, &g->ref.bmp->n, sizeof(int), NULL);
    _dos_write(handle, &g->ref.orig->n, sizeof(int), NULL);
    _dos_write(handle, &g->ref.u->n, sizeof(int), NULL);
    _dos_write(handle, &g->ref.v->n, sizeof(int), NULL);
    _dos_write(handle, &g->vis_side, sizeof(int), NULL);
    for (tri_type* t = g->first; t != NULL; t = t->next) {
      _dos_write(handle, &t->t3d.p1->n, sizeof(int), NULL);
      _dos_write(handle, &t->t3d.p2->n, sizeof(int), NULL);
      _dos_write(handle, &t->t3d.p3->n, sizeof(int), NULL);
      _dos_write(handle, &t->select, sizeof(bool), NULL);
      _dos_write(handle, t->bk_light, 3 * sizeof(int), NULL);
    }
  }

  _dos_close(handle);
}


void load_scene (const char* filename, ray_type* cursor,
   int* mouse_x, int* mouse_y,
   point_type** first_point, tri_group** first_group)
{
  int handle, point_ct, group_ct, bmp_ct, n;
 
  free_all (first_point, first_group);
  
  if (_dos_open(filename, 0x02, &handle))
    exit(5);
  
  int version;
  _dos_read (handle, &version, sizeof(int), NULL);
  _dos_read (handle, &_mk_dat, sizeof(mk_dat_type), NULL);
  _dos_read (handle, cursor, sizeof(ray_type), NULL);
  _dos_read (handle, mouse_x, sizeof(int), NULL);
  _dos_read (handle, mouse_y, sizeof(int), NULL);
  _dos_read (handle, &point_ct, sizeof(int), NULL);
  _dos_read (handle, &group_ct, sizeof(int), NULL);
  _dos_read (handle, &bmp_ct, sizeof(int), NULL);

  mk_bmp** bmp_dat = (mk_bmp**)malloc(sizeof(mk_bmp*) * bmp_ct);

  if (bmp_ct)
    for (int i = 0; i < bmp_ct; i++) {
      char filename[13];
      _dos_read(handle, filename, 13, NULL);
      load_mk_bmp(filename);
    }
  else
    _first_mkb = NULL;
    
  n = 0;
  for (mk_bmp* b = _first_mkb; b != NULL; b = b->next, n++)
    bmp_dat[n] = b;
  
  _dos_read (handle, &n, sizeof(int), NULL);
  curr_bmp = bmp_dat[n];
  
  for (int i = 0; i < 10; i++) {
    _dos_read(handle, &n, sizeof(int), NULL);
    _bmp_pal[i] = bmp_dat[n];
  }

  if (point_ct) {
    point_type* last = *first_point = read_point (handle);
    for (int i = 2; i <= point_ct; i++)
      last = last->next = read_point (handle);
    last->next = NULL;
  }
  else
    *first_point = NULL;
    
  point_type** point_dat =
    (point_type**)malloc(sizeof(point_type*) * point_ct);
  n = 0;
  for (point_type* p = *first_point; p != NULL; p = p->next, n++)
    point_dat[n] = p;

  if (group_ct) {
    tri_group* last = *first_group = read_group (handle, point_dat, bmp_dat);
    for (int i = 2; i <= group_ct; i++)
      last = last->next = read_group (handle, point_dat, bmp_dat);
    last->next = NULL;
  }
  else
    *first_group = NULL;

  free(point_dat);

  _dos_close(handle);

  for (tri_group* g = *first_group; g != NULL; g = g->next)
    g->select = false;

  for (point_type* p = *first_point; p != NULL; p = p->next)
    p->visible = true;
//    for (tri_type* t = g->first; t != NULL; t = t->next) {
//      t->bk_light[0] =
//    }
//  resize_refs (*first_group);
}


point_type* read_point (int handle)
{
  point_type* p = (point_type*)malloc (sizeof(point_type));
  _dos_read (handle, p, sizeof(point_type), NULL);
  return p;
}


tri_type* read_tri (int handle, point_type** point_dat)
{
  int pt;
  tri_type* t = (tri_type*)malloc (sizeof(tri_type));
  
  _dos_read(handle, &pt, sizeof(int), NULL);
  t->t3d.p1 = point_dat[pt];
  _dos_read(handle, &pt, sizeof(int), NULL);
  t->t3d.p2 = point_dat[pt];
  _dos_read(handle, &pt, sizeof(int), NULL);
  t->t3d.p3 = point_dat[pt];
  _dos_read(handle, &t->select, sizeof(bool), NULL);
  _dos_read(handle, t->bk_light, 3 * sizeof(int), NULL);

  return t;
}


tri_group* read_group (int handle, point_type** point_dat, mk_bmp** bmp_dat)
{
  tri_group* g = (tri_group*)malloc (sizeof(tri_group));
  int tri_ct, bmp, pt;
  
  _dos_read(handle, &tri_ct, sizeof(int), NULL);
  _dos_read(handle, &bmp, sizeof(int), NULL);
  g->ref.bmp = bmp_dat[bmp];
  _dos_read(handle, &pt, sizeof(int), NULL);
  g->ref.orig = point_dat[pt];
  _dos_read(handle, &pt, sizeof(int), NULL);
  g->ref.u = point_dat[pt];
  _dos_read(handle, &pt, sizeof(int), NULL);
  g->ref.v = point_dat[pt];
  _dos_read(handle, &g->vis_side, sizeof(int), NULL);
  calc_plane (&g->plane, g->ref.u->abs, g->ref.v->abs, g->ref.orig->abs);
  
  tri_type* last = g->first = read_tri (handle, point_dat);
  for (int i = 2; i <= tri_ct; i++) {
    tri_type* t = read_tri (handle, point_dat);
    last->next = t;
    last = t;
  }
  last->next = NULL;

  return g;
}


void free_all (point_type** first_point, tri_group** first_group)
{
  for (mk_bmp* b = _first_mkb; b != NULL; ) {
    mk_bmp* fr = b;
    b = b->next;
    free(fr);
  }
  
  for (point_type* p = *first_point; p != NULL; ) {
    point_type* fr = p;
    p = p->next;
    free(fr);
  }

  for (tri_group* g = *first_group; g != NULL; g = g->next)
    for (tri_type* t = g->first; t != NULL; ) {
      tri_type* fr = t;
      t = t->next;
      free(fr);
    }

  for (tri_group* g = *first_group; g != NULL; ) {
    tri_group* fr = g;
    g = g->next;
    free(fr);
  }

  point_type* pn = NULL;
  
  *first_point = NULL;
  *first_group = NULL;
  _first_mkb = NULL;
}


point_type* find_last_point (point_type* p)
{
  while (p->next != NULL)
    p = p->next;
  return p;
}


tri_group* find_last_group (tri_group* g)
{
  while (g->next != NULL)
    g = g->next;
  return g;
}


tri_group* create_group
  (tri_group** first_group, point_type** first_point, tri_type* tri,
   mk_bmp* bmp, point_3d orig, point_3d u, point_3d v, point_3d cam)
{
  tri_group* g = (tri_group*)malloc (sizeof(tri_group));

  g->first = tri;
  g->ref.bmp = bmp;
  g->ref.orig = create_ref_point (orig, first_point);
  g->ref.u = create_ref_point (u, first_point);
  g->ref.v = create_ref_point (v, first_point);
  g->next = NULL;
  g->select = false;
  
  calc_plane (&g->plane, orig, u, v);
  g->vis_side = point_relative_to_plane (g->plane, cam);
  
  if (*first_group == NULL)
    *first_group = g;
  else {
    tri_group* last = find_last_group (*first_group);
    last->next = g;
  }
  
  return g;
}


void ginput (int* inp_pos, char* keyb_inp)
{
  for (int k = 0x2; k <= 0xA; k++)
    if (_keyb_stat[k].change) {
      keyb_inp[*inp_pos] = '1' + k - 0x2;
      keyb_inp[*inp_pos + 1] = '\0';
      (*inp_pos)++;
    }

  if (_keyb_stat[0xB].change) {
    keyb_inp[*inp_pos] = '0';
    keyb_inp[*inp_pos + 1] = '\0';
    (*inp_pos)++;
  }
  else if (_keyb_stat[0xC].change) {
    keyb_inp[*inp_pos] = '-';
    keyb_inp[*inp_pos + 1] = '\0';
    (*inp_pos)++;
  }
  else if (_keyb_stat[0x34].change) {
    keyb_inp[*inp_pos] = '.';
    keyb_inp[*inp_pos + 1] = '\0';
    (*inp_pos)++;    
  }
}


double get_valu_buff (char* buff)
{
  int dig = 0;
  int dec = -1;
  int start;
  
  while (buff[dig] != '\0') {
    if (buff[dig] == '.')
      dec = dig;
    dig++;
  }

  if (buff[0] == '-')
    start = 1;
  else
    start = 0;

  double result = 0;
  
  if (dec < 0) {
    double mul = 1;
    for (int i = dig - 1; i >= start; i--) {
      result += double(buff[i] - '0') * mul;
      mul *= 10;
    }
  }
  else {
    double mul = .1;
    double f = 0;
    for (int i = dec + 1; i < dig; i++) {
      f += (double)(buff[i] - '0') * mul;
      mul *= .1;
    }
    mul = 1;
    for (int i = dec - 1; i >= start; i--) {
      result += (double)(buff[i] - '0') * mul;
      mul *= 10;
    }
    result += f;
  }

  if (buff[0] == '-')
    return -result;
  else
    return result;
}


mk_bmp* load_mk_bmp (const char* filename)
{
  mk_bmp* bmp = (mk_bmp*)malloc (sizeof(mk_bmp));
  load_bmp (filename, &bmp->dat);
  strcpy (bmp->filename, filename);
  mk_bmp* p_last = find_last_bmp (_first_mkb);
  if (p_last == NULL)
    _first_mkb = bmp;
  else
    p_last->next = bmp;
  bmp->next = NULL;

  return bmp;
}


mk_bmp* find_last_bmp (mk_bmp* first)
{
  mk_bmp* b;
  
  if (first == NULL)
    return NULL;
  else {
    for (b = first; b->next != NULL; b = b->next)
      ;
    return b;
  }
}


point_3d extrapolate_quad (point_type* p1, point_type* p2, point_type* p3)
{
  point_type* orig;
  point_type* pa;
  point_type* pb;
  point_3d pnew;
  parse_right_tri (p1, p2, p3, &orig, &pa, &pb);

  double dx1 = pa->abs.x - orig->abs.x;
  double dy1 = pa->abs.y - orig->abs.y;
  double dz1 = pa->abs.z - orig->abs.z;
  double dx2 = pb->abs.x - orig->abs.x;
  double dy2 = pb->abs.y - orig->abs.y;
  double dz2 = pb->abs.z - orig->abs.z;

  pnew.x = pa->abs.x + dx2;
  pnew.y = pa->abs.y + dy2;
  pnew.z = pa->abs.z + dz2;

  return pnew;
}


point_3d extrapolate_quad
  (point_type* p1, point_type* p2, point_type* p3,
   point_type** orig, point_type** pa, point_type** pb)
{
  point_3d pnew;

  parse_right_tri (p1, p2, p3, orig, pa, pb);
  double dx1 = (*pa)->abs.x - (*orig)->abs.x;
  double dy1 = (*pa)->abs.y - (*orig)->abs.y;
  double dz1 = (*pa)->abs.z - (*orig)->abs.z;
  double dx2 = (*pb)->abs.x - (*orig)->abs.x;
  double dy2 = (*pb)->abs.y - (*orig)->abs.y;
  double dz2 = (*pb)->abs.z - (*orig)->abs.z;
  pnew.x = (*pa)->abs.x + dx2;
  pnew.y = (*pa)->abs.y + dy2;
  pnew.z = (*pa)->abs.z + dz2;

  return pnew;
}


void parse_right_tri
  (point_type* p1, point_type* p2, point_type* p3,
   point_type** orig, point_type** pa, point_type** pb)
{
  double dist_01 = dist_3d (p1->abs, p2->abs);
  double dist_12 = dist_3d (p2->abs, p3->abs);
  double dist_20 = dist_3d (p3->abs, p1->abs);

  if (dist_01 > dist_12)
    if (dist_01 > dist_20) { // 01
      *orig = p3;
      *pa = p1;
      *pb = p2;
    }
    else { // 20
      *orig = p2;
      *pa = p1;
      *pb = p3;
    }
  else
    if (dist_12 > dist_20) { // 12
      *orig = p1;
      *pa = p2;
      *pb = p3;
    }
    else { // 20
      *orig = p2;
      *pa = p1;
      *pb = p3;
    }
}


int count_corners (point_type* c)
{
  int ct = 0;
  
  for (; c != NULL; c = c->next)
    if (!c->ref)
      ct++;

  return ct;
}


int count_ref_pts (point_type* first)
{
  int ct = 0;
  
  for (; first != NULL; first = first->next)
    if (first->ref)
      ct++;

  return ct;
}


int count_tris (tri_group* g)
{
  int ct = 0;
  
  for (; g != NULL; g = g->next)
    for (tri_type* t = g->first; t != NULL; t = t->next)
      ct++;

  return ct;
}


void resize_refs (tri_group* g)
{
  double side = g->ref.bmp->dat.side;
  double ux = g->ref.u->abs.x - g->ref.orig->abs.x;
  double uy = g->ref.u->abs.y - g->ref.orig->abs.y;
  double uz = g->ref.u->abs.z - g->ref.orig->abs.z;
  double len = sqrt(ux * ux + uy * uy + uz * uz);
  
  double nux = (ux / len) / TEXELS_PER_UNIT * side;
  double nuy = (uy / len) / TEXELS_PER_UNIT * side;
  double nuz = (uz / len) / TEXELS_PER_UNIT * side;

  double vx = g->ref.v->abs.x - g->ref.orig->abs.x;
  double vy = g->ref.v->abs.y - g->ref.orig->abs.y;
  double vz = g->ref.v->abs.z - g->ref.orig->abs.z;
  double nvx = (vx / len) / TEXELS_PER_UNIT * side;
  double nvy = (vy / len) / TEXELS_PER_UNIT * side;
  double nvz = (vz / len) / TEXELS_PER_UNIT * side;

  g->ref.u->abs.x = g->ref.orig->abs.x + nux;
  g->ref.u->abs.y = g->ref.orig->abs.y + nuy;
  g->ref.u->abs.z = g->ref.orig->abs.z + nuz;

  g->ref.v->abs.x = g->ref.orig->abs.x + nvx;
  g->ref.v->abs.y = g->ref.orig->abs.y + nvy;
  g->ref.v->abs.z = g->ref.orig->abs.z + nvz;
}


inline bool same_point (point_3d p1, point_3d p2)
{
  if (approx_equal (p1.x, p2.x))
    if (approx_equal (p1.y, p2.y))
      if (approx_equal (p1.z, p2.z))
        return true;

  return false;
}


void delete_point (point_type** first, point_type* del)
{
  if (*first != NULL)
    if (del == *first) {
      *first = (*first)->next;
      free(del);
    }
    else {
      point_type* prev;
      for (prev = *first; prev->next != del; prev = prev->next)
        ;
      prev->next = del->next;
      free (del);
    }
}


void delete_group (tri_group** first, tri_group* del)
{
  if (*first != NULL)
    if (del == *first) {
      *first = (*first)->next;
      for (tri_type* t = del->first; t != NULL; ) {
        tri_type* fr = t;
        t = t->next;
        free(fr);
      }
      free(del);
    }
    else {
      tri_group* prev;
      for (prev = *first; prev->next != del; prev = prev->next)
        ;
      prev->next = del->next;
      for (tri_type* t = del->first; t != NULL; ) {
        tri_type* fr = t;
        t = t->next;
        free(fr);
      }
      free (del);
    }
}


bool point_is_used (tri_group* first_g, point_type* p)
{
  for (tri_group* g = first_g; g != NULL; g = g->next) {
    if (g->ref.u == p)
      return true;
    if (g->ref.v == p)
      return true;
    if (g->ref.orig == p)
      return true;
    for (tri_type* t = g->first; t != NULL; t = t->next) {
      if (t->t3d.p1 == p)
        return true;
      if (t->t3d.p2 == p)
        return true;
      if (t->t3d.p3 == p)
        return true;
    }
  }
  
  return false;
}


void show_plane_points (point_type* first, tri_group* g, view_type view)
{
  for (point_type* p = first; p != NULL; p = p->next)
    if (point_relative_to_plane (g->plane, p->abs) == 0)
      p->visible = true;
    else
      p->visible = false;
}


void point_follow_cursor
  (camera_type curs, camera_type prev_curs, bool alt, point_type* p)
{
  double dx = curs.pos.x - prev_curs.pos.x;
  double dy = curs.pos.y - prev_curs.pos.y;
  double dz = curs.pos.z - prev_curs.pos.z;
  double ax = prev_curs.angle_xz;
  double dax = curs.angle_xz - prev_curs.angle_xz;
  double day = curs.angle_yz - prev_curs.angle_yz;
//  double day1 = atan(tan(day) * cos(prev_curs.angle_xz));
//  double day2 = atan(tan(day) * sin(prev_curs.angle_xz));
  
  p->abs.x += dx;
  p->abs.y += dy;
  p->abs.z += dz;
//  rotate (&p->abs.z, &p->abs.y, curs.pos.z, curs.pos.y, -day1);
//  rotate (&p->abs.x, &p->abs.y, curs.pos.x, curs.pos.y, -day2);
//  rotate (&p->abs.x, &p->abs.z, curs.pos.x, curs.pos.z, -dax);

  if (!alt) {
    rotate (&p->abs.x, &p->abs.z, curs.pos.x, curs.pos.z, ax);
    rotate (&p->abs.z, &p->abs.y, curs.pos.z, curs.pos.y, -day);
//  rotate (&result.x, &result.y, curs.pos.x, curs.pos.y, -day2);
    rotate (&p->abs.x, &p->abs.z, curs.pos.x, curs.pos.z, -curs.angle_xz);
  }
}


/*
point_3d point_follow_cursor
  (camera_type curs, camera_type prev_curs, point_3d p)
{
  double dx = curs.pos.x - prev_curs.pos.x;
  double dy = curs.pos.y - prev_curs.pos.y;
  double dz = curs.pos.z - prev_curs.pos.z;
  double ax = prev_curs.angle_xz;
  double dax = curs.angle_xz - prev_curs.angle_xz;
  double day = curs.angle_yz - prev_curs.angle_yz;

//  double day1 = atan2(sin(day) * cos(ax), cos(day));
  double day1 = day;atan2(tan(day), cos(dax));
//  double day2 = atan2(sin(day) * sin(ax), cos(day));
  double day2 = 0;//atan2(tan(day), sin(dax));
  
  point_3d result;
  
  bprint ("1: ", (180 / PI) * day1, 9);
  bprint ("2: ", (180 / PI) * day2, 10);
  
  result.x = p.x + dx;
  result.y = p.y + dy;
  result.z = p.z + dz;
  rotate (&result.x, &result.z, curs.pos.x, curs.pos.z, ax);
  rotate (&result.z, &result.y, curs.pos.z, curs.pos.y, -day);
//  rotate (&result.x, &result.y, curs.pos.x, curs.pos.y, -day2);
  rotate (&result.x, &result.z, curs.pos.x, curs.pos.z, -curs.angle_xz);

  return result;
}
*/

int all_p_sel_group (tri_group* g)
{
  bool some = false;
  bool all = true;
  
  for (tri_type* t = g->first; t != NULL; t = t->next) {
    if (t->t3d.p1->select)
      some = true;
    else
      all = false;
      
    if (t->t3d.p2->select)
      some = true;
    else
      all = false;
      
    if (t->t3d.p3->select)
      some = true;
    else
      all = false;
  }

  if (some)
    if (all)
      return 2;
    else
      return 1;
  else
    return 0;
}


int get_vis_side (tri_group group)
{
  if (group.plane.y_plane)
    if (approx_equal (group.ref.u->abs.x, group.ref.orig->abs.x))
      if (group.ref.u->abs.z > group.ref.orig->abs.z)
        if (group.ref.v->abs.x < group.ref.orig->abs.x)
          return 1;
        else
          return -1;
      else
        if (group.ref.v->abs.x < group.ref.orig->abs.x)
          return -1;
        else
          return 1;
    else if (group.ref.u->abs.x > group.ref.orig->abs.x)
      if (group.ref.v->abs.z > group.ref.orig->abs.z)
        return 1;
      else
        return -1;
    else
      if (group.ref.v->abs.z > group.ref.orig->abs.z)
        return -1;
      else
        return 1;
  else if (group.plane.m1_inf)
    if (approx_equal (group.ref.u->abs.z, group.ref.orig->abs.z))
      if (group.ref.u->abs.y < group.ref.orig->abs.y)
        if (group.ref.v->abs.z < group.ref.orig->abs.z)
          return -1;
        else
          return 1;
      else
        if (group.ref.v->abs.z < group.ref.orig->abs.z)
          return 1;
        else
          return -1;
    else
      return point_grt_line (
        group.ref.orig->abs.z, group.ref.orig->abs.y,
        group.ref.u->abs.z, group.ref.u->abs.y,
        group.ref.v->abs.z, group.ref.v->abs.y);
  else
    if (approx_equal (group.ref.u->abs.x, group.ref.orig->abs.x))
      if (group.ref.u->abs.y < group.ref.orig->abs.y)
        if (group.ref.v->abs.x < group.ref.orig->abs.x)
          return 1;
        else
          return -1;
      else
        if (group.ref.v->abs.x < group.ref.orig->abs.x)
          return -1;
        else
          return 1;
    else
      return point_les_line (
        group.ref.orig->abs.x, group.ref.orig->abs.y,
        group.ref.u->abs.x, group.ref.u->abs.y,
        group.ref.v->abs.x, group.ref.v->abs.y);
}


inline int point_grt_line
  (double x0, double y0, double x1, double y1, double x2, double y2)
{
  double n1 = (x1 - x0) * (y2 - y0);
  double n2 = (y1 - y0) * (x2 - x0);

  if (n1 > n2)
    return 1;
  else
    return -1;
}


inline int point_les_line
  (double x0, double y0, double x1, double y1, double x2, double y2)
{
  double n1 = (x1 - x0) * (y2 - y0);
  double n2 = (y1 - y0) * (x2 - x0);

  if (n1 > n2)
    return -1;
  else
    return 1;
}


     void get_vis_angle (tri_group group, double* ax, double* ay, double* axy)
{
  if (group.plane.y_plane)
    if (group.vis_side == -1) {
      *ax = nan();
      *ay = PI / 2;
      *axy = PI / 2;
    }
    else {
      *ax = nan();
      *ay = -PI / 2;
      *axy = -PI / 2;
    }
  else if (group.plane.m1_inf)
    if (group.vis_side == -1) {
      *ax = PI;
      *ay = nan();
      *axy = atan2(-1, group.plane.m2) + PI / 2;
    }
    else {
      *ax = 0;
      *ay = nan();      
      *axy = -(atan2(-1, group.plane.m2) + PI / 2);
    }
  else
    if (group.vis_side == -1) {
      *ax = atan(group.plane.m1) - PI / 2;
//      *ay = atan(group.plane.m2) - PI / 2;
      *ay = atan2(-1, group.plane.m2) + PI / 2;
      *axy = atan(tan(*ay) * sin(*ax));
    }
    else {
      *ax = atan(group.plane.m1) + PI / 2;
      *ay = (atan2(-1, group.plane.m2) + PI / 2);
      *axy = atan(tan(*ay) * sin(*ax));
//      *ay = atan(group.plane.m2) + PI / 2;
    }
}


void update_refs
  (tri_group prev_group, tri_group* curr_group,
   tri_type curr_tri, point_type fixed)
{
  double ax, ay, axy;
  double pax, pay, paxy;

//  bprint ("pax: ", (180 / PI) * pax, 8);
//  bprint ("pay: ", (180 / PI) * pay, 9);
  
  get_trans_angles
  (prev_group, prev_group.ref.orig->abs,
               prev_group.ref.u->abs,
               prev_group.ref.v->abs,
   *curr_group, curr_tri.t3d.p1->abs,
                curr_tri.t3d.p2->abs,
                curr_tri.t3d.p3->abs,
   &pax, &pay, &paxy, &ax, &ay, &axy);

   rotate (&curr_group->ref.orig->abs.x, &curr_group->ref.orig->abs.z,
     fixed.abs.x, fixed.abs.z, -pax);
   rotate (&curr_group->ref.u->abs.x, &curr_group->ref.u->abs.z,
     fixed.abs.x, fixed.abs.z, -pax);
   rotate (&curr_group->ref.v->abs.x, &curr_group->ref.v->abs.z,
     fixed.abs.x, fixed.abs.z, -pax);

   rotate (&curr_group->ref.orig->abs.z, &curr_group->ref.orig->abs.y,
     fixed.abs.z, fixed.abs.y, -pay);
   rotate (&curr_group->ref.u->abs.z, &curr_group->ref.u->abs.y,
     fixed.abs.z, fixed.abs.y, -pay);
   rotate (&curr_group->ref.v->abs.z, &curr_group->ref.v->abs.y,
     fixed.abs.z, fixed.abs.y, -pay);

   rotate (&curr_group->ref.orig->abs.x, &curr_group->ref.orig->abs.y,
     fixed.abs.x, fixed.abs.y, axy - paxy);
   rotate (&curr_group->ref.u->abs.x, &curr_group->ref.u->abs.y,
     fixed.abs.x, fixed.abs.y, axy - paxy);
   rotate (&curr_group->ref.v->abs.x, &curr_group->ref.v->abs.y,
     fixed.abs.x, fixed.abs.y, axy - paxy);

   rotate (&curr_group->ref.orig->abs.z, &curr_group->ref.orig->abs.y,
     fixed.abs.z, fixed.abs.y, ay);
   rotate (&curr_group->ref.u->abs.z, &curr_group->ref.u->abs.y,
     fixed.abs.z, fixed.abs.y, ay);
   rotate (&curr_group->ref.v->abs.z, &curr_group->ref.v->abs.y,
     fixed.abs.z, fixed.abs.y, ay);

//*
   rotate (&curr_group->ref.orig->abs.x, &curr_group->ref.orig->abs.z,
     fixed.abs.x, fixed.abs.z, ax);
   rotate (&curr_group->ref.u->abs.x, &curr_group->ref.u->abs.z,
     fixed.abs.x, fixed.abs.z, ax);
   rotate (&curr_group->ref.v->abs.x, &curr_group->ref.v->abs.z,
     fixed.abs.x, fixed.abs.z, ax);
}


point_type* get_unsel_point (tri_group& g)
{
  for (tri_type* t = g.first; t != NULL; t = t->next) {
    if (!t->t3d.p1->select)
      return t->t3d.p1;
    else if (!t->t3d.p2->select)
      return t->t3d.p2;
    else if (!t->t3d.p3->select)
      return t->t3d.p3;
    else
      ;
  }

  return NULL;
}


void get_trans_angles
  (tri_group& ga, point_3d p1a, point_3d p2a, point_3d p3a,
   tri_group& gb, point_3d p1b, point_3d p2b, point_3d p3b,
   double* axa, double* aya, double* axya,
   double* axb, double* ayb, double* axyb)
{
  if (ga.plane.y_plane)
    if (gb.plane.y_plane) {
      *axa = *aya = *axya = *axb = *ayb = *axyb = 0;
      return;
    }
    else if (gb.plane.m1_inf) {
      *axa = *aya = *axya =
      *axb = *ayb = *axyb = 0;
      if (ga.vis_side == -1)
        *axya = 0;
      else
        *axya = -PI;
      if (gb.vis_side == -1)
        *axyb = -atan(-1 / gb.plane.m2);
      else
        *axyb = -PI - atan(-1 / gb.plane.m2);
    }
    else {
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      bprint ("ay: ", (180 / PI) * *ayb, 8);
      *axya = -*axb;
      return;
    }
  else if (ga.plane.m1_inf)
    if (gb.plane.y_plane) {
      *axa = *aya = *axya =
      *axb = *ayb = *axyb = 0;
      if (ga.vis_side == -1)
        *axya = -atan(-1 / ga.plane.m2);
      else
        *axya = -PI - atan(-1 / ga.plane.m2);

      if (gb.vis_side == -1)
        *axyb = 0;
      else
        *axyb = -PI;
      return;
    }
    else if (gb.plane.m1_inf) {
//      *axa = *aya = *axya = *axb = *ayb = *axyb = 0;
//      *axya = -atan(-1 / ga.plane.m2) - -atan(-1 / gb.plane.m2);
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      *axya = *axyb = 0;            
      return;      
    }
    else {
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      *axya = *axyb = 0;      
      return;
    }
  else
    if (gb.plane.y_plane) {
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      bprint ("ay: ", (180 / PI) * *ayb, 8);      
      *axyb = -*axa;
      return;
    }
    else if (gb.plane.m1_inf) {
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      *axya = *axyb = 0;
//      *axyb = -*axa;
      return;
    }
    else {
      get_rel_angles (ga, p1a, p2a, p3a, axa, aya, axya);
      get_rel_angles (gb, p1b, p2b, p3b, axb, ayb, axyb);
      bprint ("ay: ", (180 / PI) * *ayb, 8);      
      return;
    }
}


void get_rel_angles
  (tri_group& g, point_3d p1, point_3d p2, point_3d p3,
   double* ax, double* ay, double* axy)
{
  plane_type plane;
  calc_plane (&plane, p1, p2, p3);
  int v = get_vis_side (p1, p2, p3, plane);
  
  if (plane.y_plane) {
    *ax = 0;//atan2(g.ref.u->abs.z - g.ref.orig->abs.z,
    *ay = -PI / 2;
    *axy = 0;    
    if (v == -1)
      *ay += PI;
  }
  else if (plane.m1_inf) {
    *ax = -PI / 2;
    *ay = -atan(plane.m2);
    *axy = -PI / 2;
    if (v == 1) {
      *ax += PI;
      *axy = PI / 2;
      *ay = -*ay;
    }
  }
  else if (approx_zero(plane.m1)) {
    *ax = 0;
    *axy = 0;
    if (approx_equal(p1.z, p2.z))
      if (approx_equal(p1.y, p2.y))
        *ay = -atan((p2.z - p3.z) / (p2.y - p3.y));
      else
        *ay = -atan((p2.z - p1.z) / (p2.y - p1.y));
    else
      *ay = -atan((p2.z - p1.z) / (p2.y - p1.y));
    if (v == 1) {
      *ax += PI;
      *ay = -*ay;
    }
  }
  else
    if (approx_equal(p2.x, p3.x)) {
      *ax = atan(plane.m1);
      *axy = 0;
      double m2 = -1 / plane.m1;
      double b2 = p3.z - m2 * p3.x;
      x_line_3d ref = calc_x_line_3d (p1, p2);
      point_3d pb;
      pb.x = (b2 - ref.bz) / (ref.mz - m2);
      pb.y = ref.my * pb.x + ref.by;
      pb.z = ref.mz * pb.x + ref.bz;

      if (pb.z > p3.z)
        *ay = atan(edge_fall_m (pb, p3));
      else
        *ay = atan(-edge_fall_m (pb, p3));

      if (v == 1) {
        *ax += PI;
        *ay = -*ay;
      }
    }
    else {
      *ax = atan(plane.m1);
      *axy = 0;
      double m2 = -1 / plane.m1;
      double b2 = p1.z - m2 * p1.x;
      x_line_3d ref = calc_x_line_3d (p2, p3);
      point_3d pb;
      pb.x = (b2 - ref.bz) / (ref.mz - m2);
      pb.y = ref.my * pb.x + ref.by;
      pb.z = ref.mz * pb.x + ref.bz;

      if (pb.z > p1.z)
        *ay = atan(edge_fall_m (pb, p1));
      else
        *ay = atan(-edge_fall_m (pb, p1));

      if (v == 1) {
        *ax += PI;
        *ay = -*ay;
      }
    }

/*
  if (plane.y_plane) {
    *ay = PI / 4;
  }
  else if (plane.m1_inf) {
  }
  else {
    *ax = atan(plane.m1);
    *ay = atan2(plane.m2, sin(atan(plane.m1) + PI / 2));
  }
*/
}


bool coplanar (tri_group g)
{
  plane_type plane;
  calc_plane (&plane, g.first->t3d.p1->abs,
                      g.first->t3d.p2->abs,
                      g.first->t3d.p3->abs);
                      
  for (tri_type* tri = g.first; tri != NULL; tri = tri->next) {
    if (point_approx_to_plane (plane, tri->t3d.p1->abs))
      return false;
    if (point_approx_to_plane (plane, tri->t3d.p2->abs))
      return false;    
    if (point_approx_to_plane (plane, tri->t3d.p3->abs))
      return false;
  }
  
  return true;
}


point_3d* copy_sel_points (point_type* first_point)
{
  int pct = 0;
  
  for (point_type* p = first_point; p != NULL; p = p->next)
    if (p->select)
      pct++;

  point_3d* p_backup = (point_3d*) malloc(pct * sizeof(point_3d));

  pct = 0;
  for (point_type* p = first_point; p != NULL; p = p->next)
    if (p->select) {
      p_backup[pct] = p->abs;
      pct++;
    }

  return p_backup;
}


bool is_crap (tri_group* g)
{
  if (point_approx_to_plane (g->plane, g->first->t3d.p1->abs))
    return true;
  if (point_approx_to_plane (g->plane, g->first->t3d.p2->abs))
    return true;
  if (point_approx_to_plane (g->plane, g->first->t3d.p3->abs))
    return true;

  return false;
}


void kill_crap (tri_group** first)
{
  for (tri_group* g = *first; g != NULL; g = g->next)
    if (is_crap (g))
      delete_group (first, g);
}


double alt_fall_a (plane_type plane)
{
//  return atan(1 / sqrt(((plane.m2 * plane.m2) / (plane.m1 * plane.m1)) +
//    plane.m2 * plane.m2));
  return atan(plane.m2 / cos(atan2(-1,plane.m1)));
}


void snap_angle (double* a, double grad)
{
  if (!approx_equal (*a, -PI / 2))
    if (!approx_equal (*a, -3 * PI / 2))
      *a = deg_to_rad (grad * int(rad_to_deg (*a +
      deg_to_rad(grad / 2)) / grad));
}


#endif //!INCLUDE_MAKE
