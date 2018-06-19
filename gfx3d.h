#ifndef INCLUDE_GFX3D
#define INCLUDE_GFX3D

#define Z_CUTOFF .5001
#define DITHER_ZOOM 1.1
#define TEXELS_PER_FOOT 16
#define FEET_PER_UNIT 1
#define TEXELS_PER_UNIT TEXELS_PER_FOOT * FEET_PER_UNIT

#include </programs/zmap/math2d.h>
#include </programs/zmap/math3d.h>
#include </programs/zmap/gfx2d.h>

#ifdef MODE_MAKE
struct point_type {
  point_3d abs, rel;
  bool select;
  bool ref;
  point_type* next;
  int n;
  bool visible;
};
#else
struct point_type {
  point_3d abs, rel;
  point_2d scr;
};
#endif

struct rot_data {
  double m;
  double x_sub, y_sub;
  double inv_dm;
};

struct camera_type {
  point_3d pos;
  double angle_xz;
  double angle_yz;
  double angle_xy;
};

typedef camera_type ray_type;

struct temp_tri {
  point_3d p1, p2, p3;
};

struct view_type {
  camera_type camera;
  point_2d center;
  dbl_pair zoom;
  dbl_pair window[2];
  void* dest;
  rot_data rot_xz, rot_yz, rot_xy;
  double vis_lm, vis_rm, vis_tm, vis_bm;
  double dither_dist;
};

void init_view (view_type* view, void* dest, double field_x, double field_y);
void init_zoom
  (view_type* view, point_2d center, double field_x, double field_y);
inline void calc_rot_data (rot_data* dat, double a);
inline void group_rot (rot_data d, double* x, double* y);
#ifdef MODE_MAKE
inline void calc_rel_3d_points
   (point_3d center, double angle_xz, double angle_yz, double angle_xy,
   point_type* first);
#else
inline void calc_rel_3d_points
   (point_3d center, double angle_xz, double angle_yz, double angle_xy,
   point_type* p, point_type* last);
#endif
inline point_2d map_to_scrn (view_type view, point_type p);
inline point_2d map_to_scrn (view_type view, point_3d p);
inline point_2d map_to_scrn (view_type view, point_type* p);
inline void rel_to_scr (point_3d rel, point_2d* scr, view_type view);
inline line_3d calc_line_3d (point_type* p1, point_type* p2);
inline double edge_fall_m (point_3d p1, point_3d p2);
inline int calc_aim_line (camera_type camera, line_3d* line, point_3d* ref);
inline point_3d calc_view_ref (camera_type camera);
inline bool same_direction (ray_type ray, point_3d p1, point_3d p2);
inline void calc_rel_point
   (point_3d center, double angle_xz, double angle_yz, double angle_xy,
   point_type* p);
inline void deflect_2d
  (double* res_x, double* res_y,
  double l_m, double nt_x, double nt_y, double imag_x, double imag_y);
inline int get_vis_side (point_3d& p1, point_3d& p2, point_3d& p3);
inline int get_vis_side (point_3d& p1, point_3d& p2, point_3d& p3, plane_type& plane);
inline int point_left_of_edge
   (double pnx, double pny,
   double p1x, double p1y,
   double p2x, double p2y);
inline int point_right_of_edge
  (double pnx, double pny,
   double p1x, double p1y,
   double p2x, double p2y);
   

void init_view (view_type* view, void* dest, double field_x, double field_y)
{
  view->camera.pos.x = 0;
  view->camera.pos.y = 0;
  view->camera.pos.z = 0;
  view->camera.angle_xz = 0;
  view->camera.angle_yz = 0;
  view->camera.angle_xy = 0;
  view->center.x = _screen.x_res / 2;
  view->center.y = _screen.y_res / 2;
  view->window[0].x = 0;
  view->window[0].y = 0;
  view->window[1].x = _screen.x_res;
  view->window[1].y = _screen.y_res;
  view->vis_lm = -tan(field_x / 2);
  view->vis_rm = tan(field_x / 2);
  view->vis_tm = -tan(field_y / 2);
  view->vis_bm = tan(field_y / 2);
  view->dest = dest;
  init_zoom (view, view->center, field_x / 2, field_y / 2);
  view->dither_dist = DITHER_ZOOM * view->zoom.x / TEXELS_PER_UNIT;
}


void init_zoom
  (view_type* view, point_2d center, double field_x, double field_y)
{
  view->zoom.x = center.x / tan(field_x);
  view->zoom.y = center.y / tan(field_y);
}


inline void group_rot (rot_data d, double* x, double* y)
{
  double bx = *x * d.x_sub;
  double by = *y * d.y_sub;
  *x = (by - bx) * d.inv_dm;
  *y = d.m * *x + bx;
}


inline void calc_rot_data (rot_data* dat, double a)
{
  dat->m = -1 / tan(a);
  dat->inv_dm = 1 / (dat->m - tan(a));
  dat->x_sub = sin(a) - dat->m * cos(a);
  dat->y_sub = sin(a + PI * .5) - tan(a) * cos(a + PI * .5);
}


#ifdef MODE_MAKE
inline void calc_rel_3d_points
   (point_3d center,
   double angle_xz, double angle_yz, double angle_xy,
   point_type* first)
{
  for (point_type* i = first; i != NULL; i = i->next) {
    i->rel.x = i->abs.x - center.x;
    i->rel.y = i->abs.y - center.y;
    i->rel.z = i->abs.z - center.z;
  }

  if (!approx_equal(angle_xz, 0))
  if (!approx_equal(angle_xz, 2 * PI))
    if (approx_equal(angle_xz, 1 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.x;
        i->rel.x = -i->rel.z;
        i->rel.z = temp;
      }
    else if (approx_equal(angle_xz, 2 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        i->rel.x = -i->rel.x;
        i->rel.z = -i->rel.z;
      }
    else if (approx_equal(angle_xz, 3 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.x;
        i->rel.x = i->rel.z;
        i->rel.z = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_xz);
      for (point_type* i = first; i != NULL; i = i->next)
        group_rot (rot, &i->rel.x, &i->rel.z);
    }

  if (!approx_equal(angle_yz, 0))
  if (!approx_equal(angle_yz, 2 * PI))
    if (approx_equal(angle_yz, 1 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.z;
        i->rel.z = -i->rel.y;
        i->rel.y = temp;
      }
    else if (approx_equal(angle_yz, 2 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        i->rel.z = -i->rel.z;
        i->rel.y = -i->rel.y;
      }
    else if (approx_equal(angle_yz, 3 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.z;
        i->rel.z = i->rel.y;
        i->rel.y = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_yz);
      for (point_type* i = first; i != NULL; i = i->next)
        group_rot (rot, &i->rel.z, &i->rel.y);
    }

  if (!approx_equal(angle_xy, 0))
  if (!approx_equal(angle_xy, 2 * PI))
    if (approx_equal(angle_xy, 1 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.x;
        i->rel.x = -i->rel.y;
        i->rel.y = temp;
      }
    else if (approx_equal(angle_xy, 2 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        i->rel.x = -i->rel.x;
        i->rel.y = -i->rel.y;
      }
    else if (approx_equal(angle_xy, 3 * PI / 2))
      for (point_type* i = first; i != NULL; i = i->next) {
        double temp = i->rel.x;
        i->rel.x = i->rel.y;
        i->rel.y = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_xy);
      for (point_type* i = first; i != NULL; i = i->next)
        group_rot (rot, &i->rel.x, &i->rel.y);
    }
}
#else
inline void calc_rel_3d_points
   (point_3d center,
   double angle_xz, double angle_yz, double angle_xy,
   point_type* p, point_type* last)
{
  for (point_type* i = p; i <= last; i++) {
    i->rel.x = i->abs.x - center.x;
    i->rel.y = i->abs.y - center.y;
    i->rel.z = i->abs.z - center.z;
  }

  if (!approx_equal(angle_xz, 0))
    if (approx_equal(angle_xz, 1 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.x;
        i->rel.x = -i->rel.z;
        i->rel.z = temp;
      }
    else if (approx_equal(angle_xz, 2 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        i->rel.x = -i->rel.x;
        i->rel.z = -i->rel.z;
      }
    else if (approx_equal(angle_xz, 3 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.x;
        i->rel.x = i->rel.z;
        i->rel.z = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_xz);
      for (point_type* i = p; i <= last; i++)
        group_rot (rot, &i->rel.x, &i->rel.z);
    }

  if (!approx_equal(angle_yz, 0))
    if (approx_equal(angle_yz, 1 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.z;
        i->rel.z = -i->rel.y;
        i->rel.y = temp;
      }
    else if (approx_equal(angle_yz, 2 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        i->rel.z = -i->rel.z;
        i->rel.y = -i->rel.y;
      }
    else if (approx_equal(angle_yz, 3 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.z;
        i->rel.z = i->rel.y;
        i->rel.y = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_yz);
      for (point_type* i = p; i <= last; i++)
        group_rot (rot, &i->rel.z, &i->rel.y);
    }

  if (!approx_equal(angle_xy, 0))
    if (approx_equal(angle_xy, 1 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.x;
        i->rel.x = -i->rel.y;
        i->rel.y = temp;
      }
    else if (approx_equal(angle_xy, 2 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        i->rel.x = -i->rel.x;
        i->rel.y = -i->rel.y;
      }
    else if (approx_equal(angle_xy, 3 * PI / 2))
      for (point_type* i = p; i <= last; i++) {
        double temp = i->rel.x;
        i->rel.x = i->rel.y;
        i->rel.y = -temp;
      }
    else {
      rot_data rot;
      calc_rot_data (&rot, angle_xy);
      for (point_type* i = p; i <= last; i++)
        group_rot (rot, &i->rel.x, &i->rel.y);
    }
}
#endif


inline point_2d map_to_scrn (view_type view, point_type p)
{
  point_2d p2d;

  double inv_z = 1 / p.rel.z;
  p2d.x = view.center.x + view.zoom.x * p.rel.x * inv_z;
  p2d.y = view.center.y + view.zoom.y * p.rel.y * inv_z;
  
  return p2d;
}


inline point_2d map_to_scrn (view_type view, point_type* p)
{
  point_2d p2d;

  double inv_z = 1 / p->rel.z;
  p2d.x = view.center.x + view.zoom.x * p->rel.x * inv_z;
  p2d.y = view.center.y + view.zoom.y * p->rel.y * inv_z;
  
  return p2d;
}


inline point_2d map_to_scrn (view_type view, point_3d p)
{
  point_2d p2d;

  double inv_z = 1 / p.z;
  p2d.x = view.center.x + view.zoom.x * p.x * inv_z;
  p2d.y = view.center.y + view.zoom.y * p.y * inv_z;
  
  return p2d;
}


inline void rel_to_scr (point_3d rel, point_2d* scr, view_type view)
{
  double inv_z = 1 / rel.z;
  scr->x = view.center.x + rel.x * view.zoom.x * inv_z;
  scr->y = view.center.y + rel.y * view.zoom.y * inv_z;
}


inline line_3d calc_line_3d (point_type* p1, point_type* p2)
{
  return calc_line_3d (p1->rel, p2->rel);
}


inline double edge_fall_m (point_3d p1, point_3d p2)
{
  return sqrt(sqrd(p2.x - p1.x) + sqrd(p2.z - p1.z)) / (p2.y - p1.y);
}


inline point_3d calc_view_ref (camera_type camera)
{
  point_3d ref;
  
  ref.x = camera.pos.x + sin(camera.angle_xz) * cos(camera.angle_yz);
  ref.y = camera.pos.y - sin(camera.angle_yz);
  ref.z = camera.pos.z + cos(camera.angle_xz) * cos(camera.angle_yz);

  return ref;
}


inline int calc_aim_line (camera_type camera, line_3d* line, point_3d* ref)
{
  *ref = calc_view_ref (camera);
  
  if (fabs(sin(camera.angle_yz)) == 1)
    return 0;
  else
    if (fabs(sin(camera.angle_xz)) == 1)
      if (sin(camera.angle_yz) == 0)
        return -1;
      else {
        line->mx = (ref->y - camera.pos.y) / (ref->x - camera.pos.x);
        line->bx = camera.pos.y - line->mx * camera.pos.x;
        return 1;
      }
    else {
      *line = calc_line_3d (*ref, camera.pos);
      return 2;
    }
}


inline bool same_direction (ray_type ray, point_3d p1, point_3d p2)
{
  if (approx_equal (sin(ray.angle_yz), 1))
    if (p2.y < ray.pos.y)
      return true;
    else
      return false;
  else if (approx_equal (sin(ray.angle_yz), -1))
    if (p2.y > ray.pos.y)
      return true;
    else
      return false;
  else {
    dbl_pair a1, a2;

    if (approx_equal (sin(ray.angle_yz), 0))
      if (approx_equal (sin(ray.angle_xz), -1))
        if (p2.x < ray.pos.x)
          return true;
        else
          return false;
      else if (approx_equal (sin(ray.angle_xz), 1))
        if (p2.x > ray.pos.x)
          return true;
        else
          return false;
      else {
        a1 = angles (ray.pos, p1);
        a2 = angles (ray.pos, p2);
      }
    else {
      a1 = angles (ray.pos, p1);
      a2 = angles (ray.pos, p2);
    }

    if (equivalent_angles (a1.x, a2.x))
      if (equivalent_angles (a1.y, a2.y))
        return true;
      else
        return false;
    else
      return false;
  }
}


inline void calc_rel_point
   (point_3d center, double angle_xz, double angle_yz, double angle_xy,
   point_type* p)
{
  p->rel.x = p->abs.x - center.x;
  p->rel.y = p->abs.y - center.y;
  p->rel.z = p->abs.z - center.z;

  rotate (&p->rel.x, &p->rel.z, 0, 0, angle_xz);
  rotate (&p->rel.z, &p->rel.y, 0, 0, angle_yz);
  rotate (&p->rel.x, &p->rel.y, 0, 0, angle_xy);
}


inline void deflect_2d
  (double* res_x, double* res_y,
  double l_m, double nt_x, double nt_y, double imag_x, double imag_y)
{
  if (is_approx_zero (l_m)) {
    *res_x = imag_x;
    *res_y = nt_y;
  }
  else {
    double m2 = -1 / l_m;
    double b1 = nt_y - l_m * nt_x;
    double b2 = imag_y - m2 * imag_x;
    *res_x = (b2 - b1) / (l_m - m2);
    *res_y = l_m * *res_x + b1;
  }
}


inline int get_vis_side (point_3d& p1, point_3d& p2, point_3d& p3)
{
  plane_type plane;
  calc_plane (&plane, p1, p2, p3);
  return get_vis_side (p1, p2, p3, plane);
}


inline int get_vis_side
  (point_3d& p1, point_3d& p2, point_3d& p3, plane_type& plane)
{
  if (plane.y_plane)
    if (approx_equal (p1.x, p2.x))
      if (approx_equal (p1.x, p3.x)) { // 1-2-3
//        bprint ("1-2-3", 8);
        return 1;
      }
      else if (p3.x > p1.x) { // 1-2,3
//        bprint ("1-2,3", 8);
        if (p1.z < p2.z)
          return -1;
        else
          return 1;
      }
      else { // 3,1-2
//        bprint ("3,1-2", 8);
        if (p1.z > p2.z)
          return -1;
        else
          return 1;
      }
    else if (p1.x > p2.x)
      if (approx_equal (p3.x, p1.x)) { // 2,1-3
//        bprint ("2,1-3", 8);
        if (p1.z < p3.z)
          return -1;
        else
          return 1;
      }
      else if (approx_equal (p3.x, p2.x)) { // 2-3,1
//        bprint ("2-3,1", 8);
        if (p2.z < p3.z)
          return -1;
        else
          return 1;
      }
      else if (p3.x > p1.x) { // 2,1,3
//        bprint ("2,1,3", 8);
        return point_left_of_edge (p1.z, p1.x, p2.z, p2.x, p3.z, p3.x);
      }
      else if (p3.x < p2.x) { // 3,2,1
//        bprint ("3,2,1", 8);
        return point_left_of_edge
          (p2.z, p2.x, p3.z, p3.x, p1.z, p1.x);
      }
      else { // 2,3,1
//        bprint ("2,3,1", 8);
        return point_right_of_edge
          (p3.z, p3.x, p2.z, p2.x, p1.z, p1.x);
      }
    else
      if (approx_equal (p3.x, p1.x)) { // 1-3,2
//        bprint ("1-3,2", 8);
        if (p3.z < p1.z)
          return -1;
        else
          return 1;
      }
      else if (approx_equal (p3.x, p2.x)) { // 1,2-3
//        bprint ("1,2-3", 8);
        if (p3.z < p2.z)
          return -1;
        else
          return 1;
      }
      else if (p3.x > p2.x) { // 1,2,3
//        bprint ("1,2,3", 8);
        return point_right_of_edge
          (p2.z, p2.x, p1.z, p1.x, p3.z, p3.x);
      }
      else if (p3.x < p1.x) { // 3,1,2
//        bprint ("3,1,2", 8);
        return point_right_of_edge
          (p1.z, p1.x, p3.z, p3.x, p2.z, p2.x);
      }
      else { // 1,3,2
//        bprint ("1,3,2", 8);
        return point_left_of_edge
          (p3.z, p3.x, p1.z, p1.x, p2.z, p2.x);
      }
  else if (plane.m1_inf)
    if (approx_equal (p1.y, p2.y))
      if (approx_equal (p1.y, p3.y)) { // 1-2-3
//        bprint ("1-2-3", 8);
        return -1;
      }
      else if (p3.y > p1.y) { // 1-2,3
//        bprint ("1-2,3", 8);
        if (p1.z < p2.z)
          return 1;
        else
          return -1;
      }
      else { // 3,1-2
//        bprint ("3,1-2", 8);
        if (p1.z > p2.z)
          return 1;
        else
          return -1;
      }
    else if (p1.y > p2.y)
      if (approx_equal (p3.y, p1.y)) { // 2,1-3
//        bprint ("2,1-3", 8);
        if (p1.z < p3.z)
          return 1;
        else
          return -1;
      }
      else if (approx_equal (p3.y, p2.y)) { // 2-3,1
//        bprint ("2-3,1", 8);
        if (p2.z < p3.z)
          return 1;
        else
          return -1;
      }
      else if (p3.y > p1.y) { // 2,1,3
//        bprint ("2,1,3", 8);
        return point_right_of_edge (p1.z, p1.y, p2.z, p2.y, p3.z, p3.y);
      }
      else if (p3.y < p2.y) { // 3,2,1
//        bprint ("3,2,1", 8);
        return point_right_of_edge
          (p2.z, p2.y, p3.z, p3.y, p1.z, p1.y);
      }
      else { // 2,3,1
//        bprint ("2,3,1", 8);
        return point_left_of_edge
          (p3.z, p3.y, p2.z, p2.y, p1.z, p1.y);
      }
    else
      if (approx_equal (p3.y, p1.y)) { // 1-3,2
//        bprint ("1-3,2", 8);
        if (p3.z < p1.z)
          return 1;
        else
          return -1;
      }
      else if (approx_equal (p3.y, p2.y)) { // 1,2-3
//        bprint ("1,2-3", 8);
        if (p3.z < p2.z)
          return 1;
        else
          return -1;
      }
      else if (p3.y > p2.y) { // 1,2,3
//        bprint ("1,2,3", 8);
        return point_left_of_edge
          (p2.z, p2.y, p1.z, p1.y, p3.z, p3.y);
      }
      else if (p3.y < p1.y) { // 3,1,2
//        bprint ("3,1,2", 8);
        return point_left_of_edge
          (p1.z, p1.y, p3.z, p3.y, p2.z, p2.y);
      }
      else { // 1,3,2
//        bprint ("1,3,2", 8);
        return point_right_of_edge
          (p3.z, p3.y, p1.z, p1.y, p2.z, p2.y);
      }
  else
    if (approx_equal (p1.y, p2.y))
      if (approx_equal (p1.y, p3.y)) { // 1-2-3
//        bprint ("1-2-3", 8);
        return 1;
      }
      else if (p3.y > p1.y) { // 1-2,3
//        bprint ("1-2,3", 8);
        if (p1.x < p2.x)
          return -1;
        else
          return 1;
      }
      else { // 3,1-2
//        bprint ("3,1-2", 8);
        if (p1.x > p2.x)
          return -1;
        else
          return 1;
      }
    else if (p1.y > p2.y)
      if (approx_equal (p3.y, p1.y)) { // 2,1-3
//        bprint ("2,1-3", 8);
        if (p1.x < p3.x)
          return -1;
        else
          return 1;
      }
      else if (approx_equal (p3.y, p2.y)) { // 2-3,1
//        bprint ("2-3,1", 8);
        if (p2.x < p3.x)
          return -1;
        else
          return 1;
      }
      else if (p3.y > p1.y) { // 2,1,3
//        bprint ("2,1,3", 8);
        return point_left_of_edge (p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
      }
      else if (p3.y < p2.y) { // 3,2,1
//        bprint ("3,2,1", 8);
        return point_left_of_edge
          (p2.x, p2.y, p3.x, p3.y, p1.x, p1.y);
      }
      else { // 2,3,1
//        bprint ("2,3,1", 8);
        return point_right_of_edge
          (p3.x, p3.y, p2.x, p2.y, p1.x, p1.y);
      }
    else
      if (approx_equal (p3.y, p1.y)) { // 1-3,2
//        bprint ("1-3,2", 8);
        if (p3.x < p1.x)
          return -1;
        else
          return 1;
      }
      else if (approx_equal (p3.y, p2.y)) { // 1,2-3
//        bprint ("1,2-3", 8);
        if (p3.x < p2.x)
          return -1;
        else
          return 1;
      }
      else if (p3.y > p2.y) { // 1,2,3
//        bprint ("1,2,3", 8);
        return point_right_of_edge
          (p2.x, p2.y, p1.x, p1.y, p3.x, p3.y);
      }
      else if (p3.y < p1.y) { // 3,1,2
//        bprint ("3,1,2", 8);
        return point_right_of_edge
          (p1.x, p1.y, p3.x, p3.y, p2.x, p2.y);
      }
      else { // 1,3,2
//        bprint ("1,3,2", 8);
        return point_left_of_edge
          (p3.x, p3.y, p1.x, p1.y, p2.x, p2.y);
      }
}


inline int point_left_of_edge
  (double pnx, double pny,
   double p1x, double p1y,
   double p2x, double p2y)
{
  double m = (p2x - p1x) / (p2y - p1y);
  double b = p1x - m * p1y;
  
  if (pnx < m * pny + b)
    return -1;
  else
    return 1;
}


inline int point_right_of_edge
  (double pnx, double pny,
   double p1x, double p1y,
   double p2x, double p2y)
{
  double m = (p2x - p1x) / (p2y - p1y);
  double b = p1x - m * p1y;
  
  if (pnx > m * pny + b)
    return -1;
  else
    return 1;
}


#endif // !INCLUDE_GFX3D
