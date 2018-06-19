#ifndef INCLUDE_CONTROL
#define INCLUDE_CONTROL

#define MOUSE_MAX 2048 * 2

#include <math.h>
#include <go32.h>
#include <dpmi.h>
#include <dos.h>
#include </programs/zmap/gfx2d.h>

struct key_type {
  bool press;
  bool change;
};

struct {
  key_type button[2];
  double x, y;
  double x_sensitiv;
  double y_sensitiv;
} _mouse_status;

key_type _keyb_stat[128];
_go32_dpmi_seginfo _old_key_handler, _new_key_handler;

inline void mouse_button (int raw_dat, int n);
void key_handler();
void log_key(int k);
void init_keyb();
void key_delete();
bool init_mouse();
void show_cursor();
void hide_cursor();
void set_mouse_bound (int x1, int y1, int x2, int y2);
void set_mouse_coord (int x, int y);
inline double units_to_rad (double units, double units_per_rot);
void get_mouse_stat_game();
void get_mouse_stat_screen();
void set_mouse_play();
void set_mouse_screen();
void get_mouse_change (int* x, int* y);
void get_mouse_stat_raw(int* x, int* y, int* b);


inline void mouse_button (int raw_dat, int n)
{
  if (raw_dat & n) {
    if (!_mouse_status.button[n - 1].press)
      _mouse_status.button[n - 1].change = true;
    else
      _mouse_status.button[n - 1].change = false;
    _mouse_status.button[n - 1].press = true;
  }      
  else {
    _mouse_status.button[n - 1].change = false;
    _mouse_status.button[n - 1].press = false;
  }
}


void key_handler()
{
  unsigned char al, ah;
  int raw_key;

  asm("cli; pushal");
  raw_key = inportb(0x60);
  al = inportb(0x61);
  al |= 0x82;
  outportb(0x61, al);
  al &= 0x7f;
  outportb(0x61, al);

  if (raw_key < 128) {
    if (!_keyb_stat[raw_key & 127].press)
      _keyb_stat[raw_key & 127].change = true;
    _keyb_stat[raw_key & 127].press = true;
  }
  else
    _keyb_stat[raw_key & 127].press = false;

  outportb(0x20, 0x20);
  asm("popal; sti");
}


void init_keyb()
{
  // Install new key handler
  _new_key_handler.pm_offset = (int)key_handler;
  _new_key_handler.pm_selector = _go32_my_cs();
  _go32_dpmi_get_protected_mode_interrupt_vector(0x9, &_old_key_handler);
  _go32_dpmi_allocate_iret_wrapper(&_new_key_handler);
  _go32_dpmi_set_protected_mode_interrupt_vector(0x9, &_new_key_handler);

  for (int i = 0; i < 128; i++) {
    _keyb_stat[i].press = false;
    _keyb_stat[i].change = false;
  }
}


void key_delete() // Put standard key handler back
{
  _go32_dpmi_set_protected_mode_interrupt_vector(0x9, &_old_key_handler);
}


bool init_mouse()
{
  union REGS r;
  short result;
  
  r.x.ax = 0x0;
  int86(0x33, &r, &r);
  return r.x.ax;
}


void set_mouse_screen()
{
  set_mouse_bound (0, 0, _screen.x_res, _screen.y_res);
  set_mouse_coord(_screen.x_res / 2, _screen.y_res / 2);

  _mouse_status.x_sensitiv = 1;
  _mouse_status.y_sensitiv = 1;
  
  for (int i = 0; i < 2; i++) {
    _mouse_status.button[0].press = false;
    _mouse_status.button[0].change = false;
  }

  show_cursor();  
}


void set_mouse_play()
{
  hide_cursor();
  set_mouse_bound (0, 0, MOUSE_MAX, MOUSE_MAX);
  set_mouse_coord(MOUSE_MAX / 2, MOUSE_MAX / 2);

  _mouse_status.x_sensitiv = 1;
  _mouse_status.y_sensitiv = 1;

  for (int i = 0; i < 2; i++) {
    _mouse_status.button[0].press = false;
    _mouse_status.button[0].change = false;
  }
}


void show_cursor()
{
  union REGS r;
  
  r.x.ax = 0x1;
  int86(0x33, &r, &r);
}


void hide_cursor()
{
  union REGS r;
  
  r.x.ax = 0x2;
  int86(0x33, &r, &r);
}


void set_mouse_bound (int x1, int y1, int x2, int y2)
{
  union REGS r;

  r.x.ax = 0x7;
  r.x.cx = short(x1);
  r.x.dx = short(x2);
  int86(0x33, &r, &r);

  r.x.ax = 0x8;
  r.x.cx = short(y1);
  r.x.dx = short(y2);
  int86(0x33, &r, &r);
}


void get_mouse_stat_screen()
{
  short mx, my;
  union REGS r;
  short b;
  
  r.x.ax = 0x3;
  int86(0x33, &r, &r);
  mx = r.x.cx;
  my = r.x.dx;
  b = r.x.bx;

  _mouse_status.x = mx * _mouse_status.x_sensitiv;
  _mouse_status.y = my * _mouse_status.y_sensitiv;
  mouse_button (b, 1);
  mouse_button (b, 2);
}


void get_mouse_stat_raw(int* x, int* y, int* b)
{
  union REGS r;
  
  r.x.ax = 0x3;
  int86(0x33, &r, &r);
  *x = r.x.cx;
  *y = r.x.dx;
  *b = r.x.bx;
}


void get_mouse_stat_game()
{
  short mx, my;
  union REGS r;
  short b;
  
  r.x.ax = 0x3;
  int86(0x33, &r, &r);
  mx = r.x.cx;
  my = r.x.dx;
  b = r.x.bx;

  if (mx == MOUSE_MAX) {
    set_mouse_coord (0, my);
    _mouse_status.x = PI;
  }
  else if (mx == 0) {
    set_mouse_coord (MOUSE_MAX, my);
    _mouse_status.x = PI;
  }
  else
    _mouse_status.x = (double(mx - MOUSE_MAX / 2) / double(MOUSE_MAX / 2))
      * PI * _mouse_status.x_sensitiv;

  if (my == MOUSE_MAX)
    _mouse_status.y = PI / 2;
  else if (my == 0)
    _mouse_status.y = (PI * 2) - PI / 2;
  else
    _mouse_status.y = (double(my - MOUSE_MAX / 2) / double(MOUSE_MAX / 2))
      * (PI / 2) * _mouse_status.y_sensitiv;

  mouse_button (b, 1);
  mouse_button (b, 2);
}


void get_mouse_change (int* x, int* y)
{
  short mx, my;
  union REGS r;
  short b;
  
  r.x.ax = 0x3;
  int86(0x33, &r, &r);
  mx = r.x.cx;
  my = r.x.dx;
  b = r.x.bx;

  *x = mx - MOUSE_MAX / 2;
  *y = my - MOUSE_MAX / 2;
  set_mouse_coord(MOUSE_MAX / 2, MOUSE_MAX / 2);
  
  mouse_button (b, 1);
  mouse_button (b, 2);
}


void set_mouse_coord (int x, int y)
{
  union REGS r;

  r.x.ax = 0x4;
  r.x.cx = x;
  r.x.dx = y;
  int86(0x33, &r, &r);
}


inline double units_to_rad (double units, double units_per_rot)
{
  return units * ((PI * 2) / units_per_rot);
}


#endif // !INCLUDE_CONTROL
