#ifndef INCLUDE_GFX2D
#define INCLUDE_GFX2D

#include <stdio.h>
#include <dpmi.h>
#include <io.h>
#include <stdlib.h>
#include <string.h>
#include <dos.h>
#include <go32.h>
#include <sys/nearptr.h>
#include <sys/farptr.h>
#include <c:/programs/zmap/math2d.h>

#define HIGH_COLOR_320x200 0x10e
#define TRUE_COLOR_320x200 0x10f
#define HIGH_COLOR_640x480 0x111
#define TRUE_COLOR_640x480 0x112
#define HIGH_COLOR_800x600 0x114
#define TRUE_COLOR_800x600 0x115
#define HIGH_COLOR_1024x768 0x117
#define TRUE_COLOR_1024x768 0x118
#define BYTES_PER_PXL 3

#ifndef byte
#define byte unsigned char
#endif


struct true_clr {
  byte r, g, b;
};

struct {
  int x_res, y_res;
  point_2d center;
  point_2d window[2];
  int bytes_per_pxl;
  int bank_ct;
} _screen;

struct bmp_type {
  unsigned char* texel; // 0
  unsigned long shift_addr; // 4
  unsigned long bitmask; // 8
  unsigned long row_bytes; // 12
  unsigned long last_x; // 16
  unsigned long last_y; // 20
  unsigned long wrap_y; // 24
  double side; // 28
};

// 14-bytes
struct {
  short bf_type;
  long bf_size;
  long bf_reserved;
  long bf_off_bits; 
} bmp_file_header;

// 40 bytes
struct {
  long bi_size;
  long bi_width;
  long bi_height;
  short bi_planes;
  short bi_bit_count;
  long bi_compression;
  long bi_size_image;
  long bi_x_in_meter;
  long bi_y_in_meter;
  long bi_clr_used;
  long bi_clr_important;
} bmp_info_header;

struct VESA_INFO {
  char legit[4] __attribute__ ((packed));
  unsigned short version __attribute__ ((packed));
  unsigned long oem_string_ptr __attribute__ ((packed));
  unsigned char capabilities[4] __attribute__ ((packed));
  unsigned long video_mode_ptr __attribute__ ((packed));
  unsigned short total_memory __attribute__ ((packed));
  unsigned short oem_software_rev __attribute__ ((packed));
  unsigned long oem_vendor_name_Ptr __attribute__ ((packed));
  unsigned long oem_product_name_ptr __attribute__ ((packed));
  unsigned long oem_product_rev_ptr __attribute__ ((packed));
  unsigned char reserved[222] __attribute__ ((packed));
  unsigned char oem_data[256] __attribute__ ((packed));
} vesa_info;

struct MODE_INFO
{
  unsigned short mode_attributes __attribute__ ((packed));
  unsigned char win_a_attributes __attribute__ ((packed));
  unsigned char win_b_attributes __attribute__ ((packed));
  unsigned short win_granularity __attribute__ ((packed));
  unsigned short win_size __attribute__ ((packed));
  unsigned short win_a_segment __attribute__ ((packed));
  unsigned short win_b_segment __attribute__ ((packed));
  unsigned long win_func_ptr __attribute__ ((packed));
  unsigned short bytes_per_scanline __attribute__ ((packed));
  unsigned short x_res __attribute__ ((packed));
  unsigned short y_res __attribute__ ((packed));
  unsigned char x_char_size __attribute__ ((packed));
  unsigned char y_char_size __attribute__ ((packed));
  unsigned char number_of_planes __attribute__ ((packed));
  unsigned char bits_per_pixel __attribute__ ((packed));
  unsigned char number_of_banks __attribute__ ((packed));
  unsigned char memory_model __attribute__ ((packed));
  unsigned char bank_size __attribute__ ((packed));
  unsigned char number_of_image_pages __attribute__ ((packed));
  unsigned char reserved_page __attribute__ ((packed));
  unsigned char red_mask_size __attribute__ ((packed));
  unsigned char red_mask_pos __attribute__ ((packed));
  unsigned char green_mask_size __attribute__ ((packed));
  unsigned char green_mask_pos __attribute__ ((packed));
  unsigned char blue_mask_size __attribute__ ((packed));
  unsigned char blue_mask_pos __attribute__ ((packed));
  unsigned char reserved_mask_size __attribute__ ((packed));
  unsigned char reserved_mask_pos __attribute__ ((packed));
  unsigned char direct_color_mode_info __attribute__ ((packed));
  unsigned char* phys_base_ptr __attribute__ ((packed));
  unsigned long off_screen_mem_offset __attribute__ ((packed));
  unsigned short off_screen_mem_size __attribute__ ((packed));
  unsigned char reserved[206] __attribute__ ((packed));
} _mode_info;

char _text_buff[20][40];
unsigned long _shift_tbl[10][1024];

void create_shift_tbl();
void show_text (byte* font, void* buff);
void bprint (char* s, double n, int row);
void bprint (char* s, int row);
void bprint (char* s1, char* s2, int row);
inline double calc_dydx_b (point_2d point, double m);
inline double calc_dxdy_b (point_2d point, double m);
void show_page (bool page);
inline void wait_vert_retrace();
inline void write_row_to_scrn (int source, int dest_offset, int len);
void copy_scrn_buff (void* buff, int page);
void pxl (point_2d p, int r, int g, int b, void* buff);
void pxl (int x, int y, int red, int green, int blue, void* buff);
void pxl (int x, int y, true_clr c);
void pxl (int x, int y, int red, int green, int blue);
void pxl (int x, int y, true_clr c, void* buff);
int set_vesa_mode(int mode);
int get_mode_info();
int set_vesa_bank(int bank_number);
int load_bmp (const char* filename, bmp_type* bitmap);
int set_palette_mode (int mode);
int point (int x, int y);
int set_logical_line_len (int len);
int set_first_pxl (int x, int y);
void printg (char* str, int column, int row, byte* font, void* buff);
void gprint (char* s, int column, int row, byte* font, void* buff);
void gprint (char* s, double v, int column, int row, byte* font, void* buff);
byte* load_text();
void* create_scrn_buff();
void clear_scrn_buff(void* buff);
void clear_screen();
long total_video_mem();
int video_mode_pages (int x_res, int y_res);
void draw_line (point_2d p1, point_2d p2, true_clr color);
void draw_line (point_2d p1, point_2d p2, int r, int g, int b);


inline double calc_dxdy_b (point_2d point, double m)
{
  return point.x - point.y * m;
}


inline double calc_dydx_b (point_2d point, double m)
{
  return point.y - point.x * m;
}


void draw_line (line_equat line, true_clr color)
{
  if (is_abs_inf(line.m)) {
    for (int y = 0; y < 199; y++)
      pxl ((int)line.b, y, color);
  }
  else if (fabs(line.m) < 1) {
    double y = line.b;
    for (int x = 0; x < 320; x++, y += line.m)
      pxl (x, int(y), color);
  }
}


void draw_line (point_2d p1, point_2d p2, int r, int g, int b)
{
  true_clr c;
  c.r = r;
  c.g = g;
  c.b = b;

  draw_line (p1, p2, c);
}


void draw_line (point_2d p1, point_2d p2, true_clr color)
{
  line_equat line;
  
  calc_2d_dydx_line (&line, p1, p2);
  
  if (p1.x == p2.x)
    if (p1.y < p2.y)
      for (int y = (int)p1.y; y <= (int)p2.y; y++)
        pxl ((int)line.b, y, color);
    else
      for (int y = (int)p2.y; y <= (int)p1.y; y++)
        pxl ((int)line.b, y, color);
  else if (fabs(line.m) < 1)
    if (p1.x < p2.x) {
      double y = p1.y;
      for (int x = (int)p1.x; x <= (int)p2.x; x++, y += line.m)
        pxl (x, int(y), color);
    }
    else {
      double y = p2.y;
      for (int x = (int)p2.x; x <= (int)p1.x; x++, y += line.m)
        pxl (x, int(y), color);
    }
  else {
    calc_2d_dxdy_line (&line, p1, p2);
    if (p1.y < p2.y) {
      double x = p1.x;
      for (int y = (int)p1.y; y <= (int)p2.y; y++, x += line.m)
        pxl (int(x), y, color);
    }
    else {
      double x = p2.x;
      for (int y = (int)p2.y; y <= (int)p1.y; y++, x += line.m)
        pxl (int(x), y, color);
    }
  }
}


void bprint (char* s, double n, int row)
{
  sprintf (_text_buff[row], "%s%f", s, n);
}


void bprint (char* s1, char* s2, int row)
{
  sprintf (_text_buff[row], "%s%s", s1, s2);
}


void bprint (char* s, int row)
{
  sprintf (_text_buff[row], "%s", s);
}


void show_text (byte* font, void* buff)
{
  for (int n = 1; n < 20; n++)
    printg (_text_buff[n], 1, n, font, buff);
}


void gprint (char* s, int column, int row, byte* font, void* buff)
{
  char str[1000];
  sprintf (str, s);
  printg (str, column, row, font, buff);
}


void gprint (char* s, double v, int column, int row, byte* font, void* buff)
{
  char str[1000];
  sprintf (str, s, v);
  printg (str, column, row, font, buff);
}


void printg (char* str, int column, int row, byte* font, void* buff)
{
  for (int i = 0; str[i] != '\0'; i++) {
    int c = (8 * 8 * _screen.bytes_per_pxl) * ((int)str[i] - 32);
    int pos = i * 8;
    for (int y = 0; y < 8; y++)
      for (int x = 0; x < 8; x++)
        if (font[c + 0 + x * _screen.bytes_per_pxl + y * 24] ||
            font[c + 1 + x * _screen.bytes_per_pxl + y * 24] ||
            font[c + 2 + x * _screen.bytes_per_pxl + y * 24])
          pxl (x + pos, y + row * 8,
            font[c + 0 + x * _screen.bytes_per_pxl + y * 24],
            font[c + 1 + x * _screen.bytes_per_pxl + y * 24],
            font[c + 2 + x * _screen.bytes_per_pxl + y * 24], buff
        );
  }
}


int load_bmp (const char* filename, bmp_type* bitmap)
{
  int handle, len;
 
  if (_dos_open(filename, 0x02, &handle)) {
    printf ("Can't load bitmap \"%s\"'\n", filename);
    getchar();
    exit(5);
  }

  len = filelength(handle);
  bitmap->texel = (unsigned char*) malloc(len);
  if (bitmap->texel == NULL) {
    printf ("Can't init %i bytes for bitmap \"%s\"\n", len, filename);
    getchar();
    exit(5);
  }
  _dos_read(handle, bitmap->texel, len, NULL);
  _dos_close(handle);

  int side = log_2((int)sqrt(len / 3));
  int l = 1 << side;
  bitmap->side = double(l);
  bitmap->shift_addr = (unsigned long)(&_shift_tbl[side - 1]);
  bitmap->last_x = l - 1;
  bitmap->last_y = (l - 1) << side;
  bitmap->row_bytes = l * 3;
  bitmap->wrap_y = l * 3 - (l * l * 3);
  bitmap->bitmask = (l - 1) << 16;

  return 0;
}


void copy_scrn_buff (void* buff, int page)
{
//  __djgpp_nearptr_enable();
  
  int bytes = _screen.x_res * _screen.y_res * _screen.bytes_per_pxl;
  
  if (page) {
    int bank = bytes >> 16;
    int i;
    
    set_vesa_bank(bank);
    write_row_to_scrn
      (int(buff), bytes & 0xFFFF, 65536 - (bytes & 0xFFFF));

    bank++;
    i = 65536 - (bytes & 0xFFFF);
    
    for (; i < bytes - 65536; i += 65536, bank++) {
      set_vesa_bank(bank);
      write_row_to_scrn (int(buff) + i, 0, 65536);
    }

    set_vesa_bank(bank);
    write_row_to_scrn (int(buff) + i, 0, (bytes * 2) & 0xFFFF);
  }
  else {
    int bank = 0;
    int i = 0;

    for (; i < bytes - 65536; i += 65536, bank++) {
      set_vesa_bank(bank);
      write_row_to_scrn (int(buff) + i, 0, 65536);
    }
    
    set_vesa_bank(bank);
    write_row_to_scrn (int(buff) + i, 0, bytes & 65535);
  }
//  __djgpp_nearptr_disable();
}


/*
void show_gfx2 (view_type* view)
{
  int page, rem, p;
  int source;
  int dest_offset;
  int y;
  unsigned char* temp = (unsigned char*)_scrn_buff;
 
//  if (view->vis_page)
//    page = 0;
//  else
    page = 0;
  
//  while (inp(0x3DA) & 0x08);
//  while (!(inportb(0x3DA) & 0x08));

  __djgpp_nearptr_enable();
  source = (int)temp;
  
  y = 0;
  set_vesa_bank(0);
  for (dest_offset = 0; dest_offset < 65536 - 640 * 3; dest_offset += 640 * 3) {
    write_row_to_scrn (source, 320 * 3 * page + dest_offset, 320 * 3);
    source += 320 * 3;
    y++;
  }

  for (p = 1; p <= 4; p++) {
    int rem = (p * 65536) - (3 * (2 * 320 * y));
    write_row_to_scrn (source, 320 * 3 * page + dest_offset, rem);
    set_vesa_bank(p);
    write_row_to_scrn (rem + int(temp) + 3 * (y * 320), 320 * 3 * page, 320 * 3 - rem);
    y++;
    source = int(temp) + 3 * (y * 320);
    dest_offset = (3 * (2 * 320 * y)) & 65535;
    for (; dest_offset < 65536 - 640 * 3; ) {
      write_row_to_scrn (source, 320 * 3 * page + dest_offset, 320 * 3);
      source += 320 * 3;
      dest_offset += 640 * 3;
      y++;
    }
  }
/*
  rem = 2 * 65536 - (3 * (2 * 320 * y));
  write_row_to_scrn (source, dest_offset, rem);

  set_vesa_bank(2);
  write_row_to_scrn (rem + int(temp) + 3 * (y * 320), 0, 320 * 3 - rem);
  y++;  
  source = int(temp) + 3 * (y * 320);
  dest_offset = (3 * (2 * 320 * y)) & 65535;
  for (; dest_offset < 65536 - 640 * 3; ) {
    write_row_to_scrn (source, dest_offset, 320 * 3);
    source += 320 * 3;
    dest_offset += 640 * 3;
    y++;
  }
  rem = 3 * 65536 - (3 * (2 * 320 * y));
  write_row_to_scrn (source, dest_offset, rem);

  set_vesa_bank(3);
  write_row_to_scrn (rem + int(temp) + 3 * (y * 320), 0, 320 * 3 - rem);
  y++;  
  source = int(temp) + 3 * (y * 320);
  dest_offset = (3 * (2 * 320 * y)) & 65535;
  for (; dest_offset < 65536 - 640 * 3; ) {
    write_row_to_scrn (source, dest_offset, 320 * 3);
    source += 320 * 3;
    dest_offset += 640 * 3;
    y++;
  }

  y++;
  set_vesa_bank(4);
  source = int(temp) + 3 * (y * 320);
  dest_offset = (3 * (2 * 320 * y)) & 65535;
  for (; dest_offset < 65536 - 640 * 3; ) {
    write_row_to_scrn (source, dest_offset, 320 * 3);
    source += 320 * 3;
    dest_offset += 640 * 3;
    y++;
  }

  rem = (p * 65536) - (3 * (2 * 320 * y));
  write_row_to_scrn (source, 320 * 3 * page + dest_offset, rem);
  set_vesa_bank(p);
  write_row_to_scrn (rem + int(temp) + 3 * (y * 320), 320 * 3 * page, 320 * 3 - rem);
  y++;
  source = int(temp) + 3 * (y * 320);
  dest_offset = (3 * (2 * 320 * y)) & 65535;

  for (; y < 200; ) {
    write_row_to_scrn (source, 320 * 3 * page + dest_offset, 320 * 3);
    source += 320 * 3;
    dest_offset += 640 * 3;
    y++;
  }

  __djgpp_nearptr_disable();
}
*/


inline void write_row_to_scrn (int source, int dest_offset, int len)
{
  unsigned short old_sel = _fargetsel();
  _farsetsel(_dos_ds);

  asm volatile (
    "movl %0, %%esi \n"
    "movl $0xA0000, %%edi\n"
    "addl %1, %%edi \n"
//    "subl ___djgpp_base_address, %%edi \n"
    "movl %2, %%ecx \n"
    "subl $4, %%ecx \n"
    ".align 4 \n"
    
    "write_row_to_scrn_l1: \n"
    "  movl (%%esi, %%ecx), %%eax \n"
    "  .byte 0x64 \n"
    "  movl %%eax, (%%edi, %%ecx) \n"
    "sub $4, %%ecx \n"
    "jnl write_row_to_scrn_l1 \n"
  :
  :"g" (source), "g" (dest_offset), "g" (len)
  :"eax", "ecx", "edi", "esi", "memory"
  );
  _farsetsel(old_sel);
}


void pxl (int x, int y, true_clr c, void* buff)
{
  pxl (x, y, c.r, c.g, c.b, buff);
}


void pxl (int x, int y, true_clr c)
{
  pxl (x, y, c.r, c.g, c.b);
}


void pxl (point_2d p, int r, int g, int b, void* buff)
{
  pxl (p.x, p.y, r, g, b, buff);
}


void pxl (int x, int y, int red, int green, int blue, void* buff)
{
 unsigned char* temp = (unsigned char*)buff;
 
 if (x >= 0)
   if (y >= 0)
     if (x < _screen.x_res)
       if (y < _screen.y_res) {
         int o = _screen.bytes_per_pxl * (x + y * _screen.x_res);
         temp[0 + o] = blue;
         temp[1 + o] = green;
         temp[2 + o] = red;
       }
}

 
int get_vesa_info()
{
  __dpmi_regs r;
  long dosbuf = __tb & 0xFFFFF;

  for (int c = 0; c < sizeof(VESA_INFO); c++)
    _farpokeb (_dos_ds, dosbuf + c, 0);

  dosmemput ("VBE2", 4, dosbuf);
  r.x.ax = 0x4F00;
  r.x.di = dosbuf & 0xF;
  r.x.es = (dosbuf >> 4) & 0xFFFF;
  __dpmi_int (0x10, &r);

  if (r.h.ah)
    return -1;
  else {
    dosmemget (dosbuf, sizeof(VESA_INFO), &vesa_info);
    if (strncmp(vesa_info.legit, "VESA", 4))
      return -1;
    else
      return 0;
  }
}

 
int set_vesa_mode (int mode)
{
  __dpmi_regs r;
  r.x.ax = 0x4F02;
  r.x.bx = mode;
  __dpmi_int(0x10, &r);

  if (r.h.ah)
    return -1;
  else {
    if (mode == TRUE_COLOR_320x200) {
      _screen.x_res = 320;
      _screen.y_res = 200;
      _screen.bytes_per_pxl = 3;
    }
    else if (mode == TRUE_COLOR_640x480) {
      _screen.x_res = 640;
      _screen.y_res = 480;
      _screen.bytes_per_pxl = 3;
    }
    else if (mode == TRUE_COLOR_800x600) {
      _screen.x_res = 800;
      _screen.y_res = 600;
      _screen.bytes_per_pxl = 3;
    }
    else if (mode == HIGH_COLOR_320x200) {
      _screen.x_res = 320;
      _screen.y_res = 200;
      _screen.bytes_per_pxl = 2;
    }
    else if (mode == HIGH_COLOR_640x480) {
      _screen.x_res = 640;
      _screen.y_res = 480;
      _screen.bytes_per_pxl = 2;
    }
    else if (mode == HIGH_COLOR_800x600) {
      _screen.x_res = 800;
      _screen.y_res = 600;
      _screen.bytes_per_pxl = 2;
    }
    _screen.center.x = _screen.x_res / 2;
    _screen.center.y = _screen.y_res / 2;

    return 0;
  }
}


int get_mode_info()
{
  __dpmi_regs r;
  long dosbuf = __tb & 0xFFFFF;

  for (int c = 0; c < sizeof(MODE_INFO); c++)
    _farpokeb(_dos_ds, dosbuf + c, 0);

  r.x.ax = 0x4F01;
  r.x.di = dosbuf & 0xF;
  r.x.es = (dosbuf >> 4) & 0xFFFF;
  r.x.cx = 257;
  __dpmi_int(0x10, &r);
  
  if (r.h.ah)
    return 0;
  else {
    dosmemget(dosbuf, sizeof(MODE_INFO), &_mode_info);
    return 1;
  }
}


int set_vesa_bank(int bank_number)
{
  __dpmi_regs r;

  r.x.ax = 0x4F05;
  r.x.bx = 0;
  r.x.dx = bank_number;
  __dpmi_int(0x10, &r);

  return r.h.ah;
}


int set_palette_mode (int mode)
{
  __dpmi_regs r;

  r.x.ax = 0x4F08;
  r.x.bx = mode << 8;
  __dpmi_int(0x10, &r);

  return r.h.ah;  
}


int set_logical_line_len (int len)
{
  __dpmi_regs r;

  r.x.ax = 0x4F06;
  r.h.bl = 0;
  r.x.cx = len;
  __dpmi_int(0x10, &r);

  return r.h.ah;  
}


int set_first_pxl (int x, int y)
{
  __dpmi_regs r;

  r.x.ax = 0x4F07;
  r.h.bl = 0x80;
  r.h.bh = 0x00;
  r.x.cx = x;
  r.x.dx = y;
  __dpmi_int(0x10, &r);

  return r.h.ah;  
}


long total_video_mem ()
{
  return vesa_info.total_memory * 64 * 1000;
}


int video_mode_pages (int x_res, int y_res)
{
  long one_page = x_res * y_res * 3;
  
  if (total_video_mem() < one_page)
    return 0; // nope, not gonna happen
  else if (total_video_mem() < one_page * 2)
    return 1; // yes, double buffer
  else
    return 2; // yes, page flipping
}


char* load_bmp (const char* filename, int width, int height)
{
  int handle, len;
  char* bitmap;

  if (_dos_open(filename, 0x02, &handle)) {
    printf ("Can't load bitmap \"%s\"'\n", filename);
    getchar();
    return NULL;
  }
  len = filelength(handle);
  if ((bitmap = (char*) calloc(len, sizeof(char*))) ==
  NULL) {
    printf ("Can't init %i bytes for bitmap \"%s\"\n", len, filename);
    getchar();
    return NULL;
  }
  _dos_read(handle, bitmap, len, NULL);
  _dos_close(handle);

  return bitmap;
}


void pxl (int x, int y, int red, int green, int blue)
{
  long l = _screen.bytes_per_pxl * (x + y * _screen.x_res);
  
  set_vesa_bank(l >> 16);
  _farpokeb (_dos_ds, 0xa0000 + (l & 65535), blue);
  set_vesa_bank((l + 1) >> 16);
  _farpokeb (_dos_ds, 0xa0000 + ((l + 1) & 65535), green);
  set_vesa_bank((l + 2) >> 16);
  _farpokeb (_dos_ds, 0xa0000 + ((l + 2) & 65535), red);
}


int point (int x, int y)
{
  long l = _screen.bytes_per_pxl * (x + y * _screen.x_res);
  
  set_vesa_bank(l >> 16);
  return _farpeekb (_dos_ds, 0xa0000 + (l & 65535));
}


void create_shift_tbl()
{
  for (int bits = 0; bits < 10; bits++)
    for (int n = 0; n < 1024; n++)
      _shift_tbl[bits][n] = n << (bits + 1);
}


byte* load_text()
{
  int handle, len;
  unsigned char* p;
  
  if (_dos_open("text.dat", 0x02, &handle)) {
    printf ("Can't load text");
    getchar();
    exit(1);
  }

  len = filelength(handle);
  p = (unsigned char*) malloc(len);
  if (p == NULL) {
    printf ("crap");
    getchar();
    exit(1);
  }
  _dos_read(handle, p, len, NULL);
  _dos_close(handle);

  return p;
}


void* create_scrn_buff()
{
  return malloc(_screen.x_res * _screen.y_res * _screen.bytes_per_pxl);
}


void clear_scrn_buff(void* buff)
{
  for (int y = 0; y < _screen.y_res; y++)
    for (int x = 0; x < _screen.x_res; x++)
      pxl (x, y, 0, 0, 0, buff);
}


void clear_screen()
{
  for (int y = 0; y < _screen.y_res; y++)
    for (int x = 0; x < _screen.x_res; x++)
      pxl (x, y, 0, 0, 0);
}


void show_page (bool page)
{
  if (page)
    set_first_pxl (0, _screen.y_res);
  else
    set_first_pxl (0, 0);
}


inline void wait_vert_retrace()
{
  while (inp(0x3DA) & 0x08);
  while (!(inportb(0x3DA) & 0x08));
}


#endif // !INCLUDE_GFX2D
