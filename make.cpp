unsigned char* _font;
bool _collis = true;
#define MODE_MAKE

#include "make.h"
#include "colis.h"
#include "control.h"
#include "gfx3d.h"
#include "math3d.h"
#include "textri.h"
#include "zmap.h"
#include <time.h>


void main()
{
  for (int i = 0; i < 50; i++) {
    _spec[i].abs = p3d(0,0,0);
    _spec[i].select = false;
    _spec[i].ref = false;
    _spec[i].visible = true;
  }
  
  long long curr_t, it;
  tri_group* first_group = NULL;
  point_type* first_point = NULL;
  cylinder_3d cyl;
  view_type view;
  void* scrn_buff;
  _exp_shad = mk_rgb_light(256, 256, 256);
  double b, t;
  tri_group g;
//  g.vis_side = 1;
//  g.plane.y_plane = false;
//  g.plane.m1_inf = false;
//  g.plane.m1 = 1;
//  g.plane.m2 = 1;
  cyl.radius = 1;
  cyl.height = 6.5;
  
  init_keyb();
  if (!init_mouse()) {
    printf ("Can't initialize mouse.");
    getchar();
    return 0;
  }
  set_mouse_play();

  gen_mult_tbl();
  create_shift_tbl();
  create_dither_tbl();
  create_nx88_tbl();
  set_vesa_mode(TRUE_COLOR_320x200);
//  set_vesa_mode(TRUE_COLOR_640x480);
  byte* font = _font = load_text();
  scrn_buff = create_scrn_buff();
  init_view (&view, scrn_buff, deg_to_rad(55), deg_to_rad(55) / 1.5);
  view.dither_dist = 0.0;
  uv_acc_type* acc = create_uv_acc();
  z_buffer_type* zmap = (z_buffer_type*)malloc (sizeof(z_buffer_type)
    * _screen.x_res * _screen.y_res);

  _first_mkb = NULL;
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
  _collis = true;
  
  bool page = false;
  clear_scrn_buff(scrn_buff);
  
  camera_type cursor = view.camera;
  _mk_dat.view_rad = 10.0;
  _mk_dat.plane_mode = false;
  _mk_dat.new_group = true;

  int mx, my;
  load_scene("3dfx.dat", &cursor, &mx, &my, &first_point, &first_group);
  set_mouse_coord (mx, my);
  
  kill_crap (&first_group);
  
  it = ((long long)uclock() / (long long)UCLOCKS_PER_SEC);
  while (true) {
    _pxl_ct = 0;
    wait_vert_retrace();
    show_page (page);
    clear_scrn_buff(scrn_buff);

    calc_view (cursor, &view.camera);
    
    calc_rel_3d_points (
      view.camera.pos,
      view.camera.angle_xz, view.camera.angle_yz, view.camera.angle_xy,
      first_point);

    for (int i = -_screen.x_res - 1;
    i <= _screen.y_res + _screen.x_res + 1; i++)
      acc[i].init = 0;

    int n;
    int ct;
    
    clear_zmap (view, zmap);

    tri_group* g;
    int last_vis = 0;
    for (g = first_group, ct = 1; g != NULL; g = g->next, ct++) {
      zmap_group (*g, zmap, view, ct);
//      if (_pxl_ct == 64000) {
//        last_vis = ct;
//        break;
//      }
    }
    
    for (g = first_group, ct = 1; g != NULL; g = g->next, ct++) {
      draw_tri_group (*g, &view, acc, zmap, ct);
//      if (ct == last_vis)
//        break;
    }

//    point_3d p3dd;
//    tri_type* t = ray_int_tri (view.camera, first_group, &p3dd);

//    if (t == NULL)
//      bprint ("side: NA", 9);
//    else {
//      int boo = is_front_side_of_tri (t);
//      bprint ("y ", p3dd.y, 6);
//      bprint ("side: ", boo, 9);
//    }

//    p3dd.x = view.camera.pos.x;
//    p3dd.z = view.camera.pos.z;
//    p3dd.y = 0;

    bprint (curr_bmp->filename, 1);
    curr_t = (long long)((1000 * uclock()) / UCLOCKS_PER_SEC);
    bprint ("", 1000 / double(curr_t - it), 2);
    it = (long long)((1000 * uclock()) / UCLOCKS_PER_SEC);

    // rad_to_deg
    bprint ("", rad_to_deg(cursor.angle_xz), 3);
    bprint ("", rad_to_deg(cursor.angle_yz), 4);
//    double theta2 = atan(tan(cursor.angle_yz) * cos(cursor.angle_xz));
//    bprint ("", rad_to_deg(theta2), 5);
//    bprint ("", (double)count_tris(first_group), 6);

//    if (first_group != NULL)
//      if (first_group->next != NULL)
//        show_case (tri_rel_tri (*first_group, *(first_group->next), *(first_group->first), *(first_group->next->first)), 6);

    bprint ("vs: ", first_group->vis_side, 9);
//    bprint (": ", first_group->vis_side, 10);
    bprint ("x: ", cursor.pos.x, 11);
    bprint ("y: ", cursor.pos.y, 12);
    bprint ("z: ", cursor.pos.z, 13);
    if (_collis)
      bprint ("Colis: on", 14);
    else
      bprint ("Colis: off", 14);
      
    int pct = 0;
    for (point_type* pt = first_point; pt != NULL; pt = pt->next) {
      show_corner (*pt, view);
      pct++;
    }

//    bprint ("pct: ", pct, 14);
    double aa;

    for (tri_group* g = first_group; g != NULL; g = g->next)
      if (g->select)
        aa = alt_fall_a (g->plane);

    bprint ("aa: ", (180 / PI) * aa, 9);
    show_text(font, view.dest);

    show_cursor (view);

//point_3d deflect_path
//  (plane_type plane, point_3d ntersect, point_3d curr, point_3d tent)

//          for (tri_type* t = g->first; t != NULL; t = t->next)
//            if (check_colis_y_line_tri_3d
//              (*g, *t, &result_head, cyl, curr_head, cursor.pos, curr_feet, tentat_feet))
//              cursor.pos = result_head;
        
    move_player (view, &cursor, cyl, &first_point, &first_group);

    for (int i = 0; i < 50; i++)
      calc_rel_point (view.camera.pos, view.camera.angle_xz,
        view.camera.angle_yz, view.camera.angle_xy, &_spec[i]);

//    pxl (map_to_scrn (view, abs_to_rel(_test_p[0])), 255, 255, 255, scrn_buff);
//    pxl (map_to_scrn (view, abs_to_rel(_test_p[1])), 255, 255, 255, scrn_buff);

//inline point_3d offset_3d
//  (double r, point_3d p, double xz, double yz, double xy)
    
    for (int i = 0; i < 50; i++)
      show_corner_special (_spec[i], view);

    page = !page;
    copy_scrn_buff(scrn_buff, page);
  }
  set_vesa_mode(0x3);
  key_delete();
}

