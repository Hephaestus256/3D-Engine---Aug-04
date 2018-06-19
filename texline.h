#ifndef INCLUDE_TEXLINE
#define INCLUDE_TEXLINE

#include "/programs/include/fixed.h"
#include "/programs/zmap/texclass.h"
#include "/programs/zmap/gfx2d.h"

unsigned long _pxl_ct;
long _temp[50];
static long _eax, _ebx, _ecx, _esp, _ebp, _edx, _esi, _edi;
long _ret;
long _ret2;
long _temp1, _temp2, _temp3, _temp4, _temp5,
     _temp6, _temp7, _temp8, _temp9, _temp10;
long _temp11, _temp12, _temp13, _temp14, _temp15,
     _temp16, _temp17, _temp18, _temp19, _temp20;
long _temp21, _temp22, _temp23, _temp24, _temp25;

long _x2_buff;
long _texdat;
long _scrn;
long _f;

bmp_type _bmp;
int _b;
unsigned long _bk_light[3];

unsigned char _mult_tbl[8 * 256][256];
unsigned char _dither[256][256][4];
int _nx88[1024];
double _dtemp[50];

void update_accumulator();
void update_vert_acc();
inline void gen_mult_tbl();
void create_dither_tbl();
inline int sect_top (double y);
inline int sect_bot (double y);
inline dbl init_lef_edge (lin_relat edge, int ind);
inline dbl init_rig_edge (lin_relat edge, int ind);
void create_nx88_tbl();
void create_bmap_end_tbl();
inline unsigned long bit_not (unsigned long n);
void print_bin (int n);
inline long mk_rgb_light (int red, int green, int blue);

long _exp_shad;


inline void tri_group_class::draw_horz_sect
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, zm, istart, iend;
  static int y, p, x1, x2, bm_yc, ilen;
  dbl lef, rig, d_lef, d_rig;
  static int lx, rx;
  
  y = sect_top(t);
  p = group_num;
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  sng m = cv_sng(scan_m);
  int yc = y * _screen.x_res;
  _dtemp[1] = double(y) + .5;
  sng scan_base = cv_sng((double(y) + .5) * scan_m - .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      long scan = (lef.whole << 16) - scan_base.valu;
      f = scan & 0x0000FFFF;
//      acc_start = (int)tex_init.acc + (scan >> 16) * 88;
      acc_start = (int)tex_init.acc - (scan_base.valu >> 16) * 88;
      line_start = (int)tex_init.dest + yc * 3;
//      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      _temp3 = y + 1;
      zm = (int)zmap + _screen.x_res * y * 8;
      istart = lef.whole;
      iend = rig.whole;
//      if (_collis)
//        if ((y & 31) == 0)
//          for (int x = lef.whole; x <= rig.whole + 1; x++)
//            tex_init.acc[x].last_ind = -1;
            
      asm volatile (
        "cli; pushal \n"
        "movl %%edi, __edi \n"
        "movl %%ebp, __ebp \n"
        "movl %%esi, __esi \n"
        "movl %8, %%edi \n"
        "movl %5, %%esi \n"
        "movl %7, %%ebx \n"

        ".align 4 \n"
        "rep_h: \n"
          "movl 4(%%ebx, %%edi, 8), %%eax \n"
          "movl %6, %%edx \n"
          "cmpl %%eax, %%edx \n"
          "jne skip_zmap_h \n"
            "call tex_pxl \n"
          "skip_zmap_h: \n"
        "incl %%edi \n"
        "movl %9, %%eax \n"
        "cmpl %%edi, %%eax \n"
        "jge rep_h \n"

        "movl __ebp, %%ebp \n"
        "movl __edi, %%edi \n"
        "movl __esi, %%esi \n"
        "popal; sti \n"
      :
      :"m" (ilen), "m" (rx),
       "m" (line_start), "m" (bm_yc), "m" (y), "m" (f), "m" (p), "m" (zm),
       "m" (istart), "m" (iend)
      :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base.valu += m.valu;
    _dtemp[1]++;
  }


        asm volatile (
          "jmp skip_tl \n"
          "tex_pxl: \n"
          "movl %%esp, __esp \n"          
          "movl %%ebx, __temp8 \n"
          "movl %4, %%ecx \n"
          
          "movl __nx88(, %%edi, 4), %%eax \n"
          "movl %1, %%ebp \n"
          "addl %%eax, %%ebp \n"
          
          "movl 0 + 16(%%ebp), %%eax \n"
          "cmpl %%ecx, %%eax \n"
          "je s_gfx_hseg12 \n"
          "  leal s_gfx_hseg12, %%eax \n"
          "  movl %%eax, __ret \n"
          "  jmp update_acc \n"
          "s_gfx_hseg12: \n"

          "movl 88 + 16(%%ebp), %%eax \n"
          "cmpl %%ecx, %%eax \n"
          "je s_gfx_hseg13 \n"
          "  leal s_gfx_hseg13, %%eax \n"
          "  movl %%eax, __ret \n"
          "  jmp update_next_acc \n"
          "s_gfx_hseg13: \n"

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
          "addl %%edx, %%ecx \n"

          "movl %%ebx, %%eax \n"
          "imull %%esi \n"
          "shll $16, %%edx \n"
          "shrl $16, %%eax \n"
          "addl %%eax, %%edx \n"
          "addl %%esp, %%edx \n"

          "movl __bmp + 8, %%esp \n"
          "andl %%esp, %%ecx \n"
          "andl %%esp, %%edx \n"
          "shrl $16, %%ecx \n"
          "shrl $16, %%edx \n"
          "movl __bmp + 4, %%eax \n"
          "movl (%%eax, %%edx, 4), %%edx \n"
          
          "addl %%ecx, %%edx \n"
          "movl __bmp + 0, %%ebx \n"
          "leal (%%edx, %%edx, 2), %%edx \n"
          "movl (%%ebx, %%edx), %%edx \n"
          "movl __exp_shad, %%eax \n"

          "movl %%eax, %%ebx \n"
          "shll $8, %%ebx \n"
          "movb %%dl, %%bl \n"
          "andl $0x0003FFFF, %%ebx \n"
          "movb __mult_tbl(%%ebx), %%cl \n"

          "movl %%eax, %%ebx \n"
          "shrl $2, %%ebx \n"
          "movb %%dh, %%bl \n"
          "andl $0x0007FFFF, %%ebx \n"
          "movb __mult_tbl(%%ebx), %%ch \n"
    
          "shrl $13, %%eax \n"
          "shrl $16, %%edx \n"
          "movb %%dl, %%al \n"
          "movb __mult_tbl(%%eax), %%dl \n"

          "movl %2, %%ebx \n"
          "leal (%%edi, %%edi, 2), %%eax \n"

          "movw %%cx, 0(%%eax, %%ebx) \n"
          "movb %%dl, 2(%%eax, %%ebx) \n"
/*
          "movl %2, %%esp \n"
          "leal (%%edi, %%edi, 2), %%eax \n"
          "addl %%esp, %%eax \n"

          "movl __bk_light + 0, %%ecx \n"
          "movb 0(%%ebx, %%edx), %%cl \n"
          "movb __mult_tbl(%%ecx), %%cl \n"
          "movb %%cl, 0(%%eax) \n"

          "movl __bk_light + 4, %%ecx \n"
          "movb 1(%%ebx, %%edx), %%cl \n"
          "movb __mult_tbl(%%ecx), %%cl \n"
          "movb %%cl, 1(%%eax) \n"

          "movl __bk_light + 8, %%ecx \n"
          "movb 2(%%ebx, %%edx), %%cl \n"
          "movb __mult_tbl(%%ecx), %%cl \n"
          "movb %%cl, 2(%%eax) \n"
*/

          "movl __temp3, %%eax \n"
          "movl %%eax, 16(%%ebp) \n"
          "movl __temp8, %%ebx \n"

//          "movl __pxl_ct, %%eax \n"
//          "incl %%eax \n"
//          "movl %%eax, __pxl_ct \n"
          
          "movl __esp, %%esp \n"                    
//          "skip_zmap_h: \n"
//          "movl __temp9, %%eax \n"
          "ret \n"
          "skip_tl: \n"
      :
      :"m" (tex_init.dest), "m" (acc_start),
       "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (p), "m" (zm),
       "m" (x1), "m" (x2)
      :"memory"
     );
//struct bmp_type {
//  unsigned char* texel; // 0
//  unsigned long w_shift; // 4
//  unsigned long x_bitmask, y_bitmask; // 8, 12
//  int width, height;
//};
}


inline void tri_group_class::draw_vert_pos_sect
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, zm;
  static int x, y, p;
  static sng m;  
  static dbl lef, rig, d_lef, d_rig;
  
  p = group_num;
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  m = cv_sng(scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(y + .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      sng scan = scan_base - cv_sng(scan_m * (lef.whole + .5));
      f = scan.valu & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan.valu >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      x = lef.whole;
      zm = (int)zmap + (lef.whole + _screen.x_res * y) * 8 + 4;
      
      asm (
        "cli \n"
        "pushal \n"
        "push %%ebp \n"
        "push %%esi \n"
        "movl %%esp, __esp \n"
        "movl %5, %%esi \n"
        "movl %2, %%edi \n"
        "movl %1, %%ebp \n"
        "movl %9, %%ebx \n"
        
        ".align 4 \n"
        "gfx_v_line_pos_rep1: \n"
        "movl (%%ebx), %%eax \n"
        "movl %8, %%ecx \n"
        "cmpl %%eax, %%ecx \n"
        "jne skip_zmap_vp \n"
        "movl %%ebx, __temp8 \n"
      
        // update u & v
        "movl 16(%%ebp), %%eax \n"

        "movl %7, %%edx \n"
        "cmp %%eax, %%edx \n"
        "je gfx_v_line_pos_skip1 \n"
        "  leal gfx_v_line_pos_skip1, %%eax \n"
        "  movl %%eax, __ret \n"
        "  leal %7, %%ebx \n"
        "  jmp update_v_acc \n"
        "gfx_v_line_pos_skip1: \n"

        // calc u
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
        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%ebx \n"
        "movl %%ebx, __temp1 \n"

        // calc v
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
        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%ebx \n"
        "movl __bmp + 4, %%ecx \n"
        "movl (%%ecx, %%ebx, 4), %%ebx \n" // lookup in shift table

        // calc pointer to bitmap
        "movl __temp1, %%eax \n"
        "addl %%ebx, %%eax \n"
        "leal (%%eax, %%eax, 2), %%eax \n" // eax *= 3
        "movl __bmp, %%ebx \n"

        "movl __bk_light + 0, %%ecx \n"
        "movb 0(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 0(%%edi) \n"
    
        "movl __bk_light + 4, %%ecx \n"
        "movb 1(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 1(%%edi) \n"
    
        "movl __bk_light + 8, %%ecx \n"
        "movb 2(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 2(%%edi) \n"

        // update last x
        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, 16(%%ebp) \n"

        "movl __temp8, %%ebx \n"
        "skip_zmap_vp: \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, %7 \n"
        
        "addl $8, %%ebx \n"
        "addl $3, %%edi \n"
        "subw %6, %%si \n"
        "jnc gfx_v_line_pos_skip2 \n"
        "  subl $88, %%ebp \n"
        "gfx_v_line_pos_skip2: \n"

        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax\n"
        "jge gfx_v_line_pos_rep1 \n"

        "movl __esp, %%esp \n"
        "pop %%esi \n"
        "pop %%ebp \n"
        "popal \n"
        "sti \n"
        :
        :"m" (tex_init.dest), "m" (acc_start),
         "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (m),
         "m" (x), "m" (p), "m" (zm)
        :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base++;
  }
/*
struct bmp_type {
  unsigned char* texel; // 0
  unsigned long w_shift; // 4
  unsigned long x_bitmask, y_bitmask; // 8, 12
  int width, height; // 16, 24
};
*/
}


inline void tri_group_class::draw_vert_neg_sect
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, zm;
  static int x, y, p;
  static sng m;  
  static dbl lef, rig, d_lef, d_rig;
  
  p = group_num;
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  m = cv_sng(-scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(y + .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      sng scan = scan_base - cv_sng(scan_m * (lef.whole + .5));
      f = scan.valu & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan.valu >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      x = lef.whole;
      zm = (int)zmap + (lef.whole + _screen.x_res * y) * 8 + 4;
      
      asm (
        "cli \n"
        "pushal \n"
        "push %%ebp \n"
        "push %%esi \n"
        "movl %%esp, __esp \n"
        "movl %5, %%esi \n"
        "movl %2, %%edi \n"
        "movl %1, %%ebp \n"
        "movl %9, %%ebx \n"
        ".align 4 \n"
        "gfx_v_line_neg_rep1: \n"
        "movl (%%ebx), %%eax \n"
        "movl %8, %%ecx \n"
        "cmpl %%eax, %%ecx \n"
        "jne skip_zmap_vn \n"
        "movl %%ebx, __temp8 \n"
      
        // update u & v
        "movl 16(%%ebp), %%eax \n"

        "movl %7, %%edx \n"
        "cmp %%eax, %%edx \n"
        "je gfx_v_line_neg_skip1 \n"
        "  leal gfx_v_line_neg_skip1, %%eax \n"
        "  movl %%eax, __ret \n"
        "  leal %7, %%ebx \n"
        "  jmp update_v_acc \n"
        "gfx_v_line_neg_skip1: \n"

        // calc u
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
        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%ebx \n"
        "movl %%ebx, __temp1 \n"

        // calc v
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

        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%ebx \n"
        "movl __bmp + 4, %%ecx \n"
        "movl (%%ecx, %%ebx, 4), %%ebx \n"

        "movl __temp1, %%eax \n"
        "addl %%ebx, %%eax \n"
        "leal (%%eax, %%eax, 2), %%eax \n"
        "movl __bmp, %%ebx \n"

        "movl __bk_light + 0, %%ecx \n"
        "movb 0(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 0(%%edi) \n"
    
        "movl __bk_light + 4, %%ecx \n"
        "movb 1(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 1(%%edi) \n"
    
        "movl __bk_light + 8, %%ecx \n"
        "movb 2(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 2(%%edi) \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, 16(%%ebp) \n"

        "movl __temp8, %%ebx \n"
        "skip_zmap_vn: \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, %7 \n"
        
        "addl $8, %%ebx \n"
        "addl $3, %%edi \n"
        "addw %6, %%si \n"
        "jnc gfx_v_line_neg_skip2 \n"
        "  addl $88, %%ebp \n"
        "gfx_v_line_neg_skip2: \n"

        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax \n"
        "jge gfx_v_line_neg_rep1 \n"

        "movl __esp, %%esp \n"
        "pop %%esi \n"
        "pop %%ebp \n"
        "popal \n"
        "sti \n"
        :
        :"m" (tex_init.dest), "m" (acc_start),
         "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (m),
         "m" (x), "m" (p), "m" (zm)
        :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base++;
  }
/*
struct bmp_type {
  unsigned char* texel; // 0
  unsigned long w_shift; // 4
  unsigned long x_bitmask, y_bitmask; // 8, 12
  int width, height; // 16, 24
};
*/
}


inline void tri_group_class::draw_horz_sect_dither
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, zm, p;
  static int y;
  dbl lef, rig, d_lef, d_rig;

  p = group_num;
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  sng m = cv_sng(scan_m);
  int yc = y * _screen.x_res;
  _dtemp[1] = double(y) + .5;
  sng scan_base = cv_sng((double(y) + .5) * scan_m - .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
//      sng s = cv_sng((lef.whole + .5) - scan_m * (y + .5));
      long scan = (lef.whole << 16) - scan_base.valu;
      f = scan & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      _temp3 = y + 1;
      zm = (int)zmap + (lef.whole + _screen.x_res * y) * 8 + 4;
      
      asm (
        "cli; pushal \n"
        "movl %%edi, __edi \n"
        "movl %%ebp, __ebp \n"
        "movl %%esp, __esp \n"
        "movl %%esi, __esi \n"
        "movl %2, %%edi \n"
//        "movl %5, %%esi \n"
        "movl %1, %%ebp \n"
        "movl %6, %%ebx \n"
        
        "l_gfx_h_seg11d: \n"
          "movl %4, %%ecx \n"
          "movl %5, %%esi \n"

          "movl (%%ebx), %%eax \n"
          "movl %7, %%edx \n"
          "cmpl %%eax, %%edx \n"
          "jne skip_zmap_hd \n"

          "movl %%ebx, __ebx \n"
          "movl 0 + 16(%%ebp), %%eax \n"
          "cmpl %%ecx, %%eax \n"
          "je s_gfx_hseg12d \n"
          "  leal s_gfx_hseg12d, %%eax \n"
          "  movl %%eax, __ret \n"
          "  jmp update_acc \n"
          "s_gfx_hseg12d: \n"

          "movl 88 + 16(%%ebp), %%eax \n"
          "cmpl %%ecx, %%eax \n"
          "je s_gfx_hseg13d \n"
          "  leal s_gfx_hseg13d, %%eax \n"
          "  movl %%eax, __ret \n"
          "  jmp update_next_acc \n"
          "s_gfx_hseg13d: \n"

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
          "addl %%edx, %%ecx \n"

          "movl %%ebx, %%eax \n"
          "imull %%esi \n"
          "shll $16, %%edx \n"
          "shrl $16, %%eax \n"
          "addl %%eax, %%edx \n"
          "addl %%esp, %%edx \n"

    "subl $32768, %%ecx \n"
    "subl $32768, %%edx \n"
    "xorl %%esi, %%esi \n"
    "movb %%dh, %%cl \n"
    "movw %%cx, %%si \n"

    "movl __bmp + 8, %%esp \n"
    "andl %%esp, %%ecx \n"
    "andl %%esp, %%edx \n"
    "shrl $16, %%ecx \n"
    "shrl $16, %%edx \n"
    "movl __bmp + 4, %%eax \n"
    "movl (%%eax, %%edx, 4), %%edx \n"

    "movl __bmp + 12, %%esp \n"
    "cmpl __bmp + 20, %%edx \n"
    "jne h_last_y \n"
    "  movl __bmp + 24, %%esp \n"
    "h_last_y: \n"

    "cmpl __bmp + 16, %%ecx \n"
    "je h_last_x \n"
    "  addl %%ecx, %%edx \n"
    "  leal (%%edx, %%edx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "jmp h_endif1 \n"

    "h_last_x: \n"
    "  addl %%ecx, %%edx \n"
    "  leal (%%edx, %%edx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    "  movl __bmp + 12, %%ecx \n"
    "  subl %%ecx, %%ebx \n"
//    "  xorl %%ecx, %%ecx \n"
    
    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "h_endif1: \n"

    "movl __bk_light + 0, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 0(%%edi) \n"

    "movl __bk_light + 4, %%ecx \n"
    "movb %%ah, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 1(%%edi) \n"

    "shrl $16, %%eax \n"
    
    "movl __bk_light + 8, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 2(%%edi) \n"

//    "movb %%al, 0(%%edi) \n"
//    "movb %%ah, 1(%%edi) \n"
//    "shrl $16, %%eax \n"
//    "movb %%al, 2(%%edi) \n"
        "movl __ebx, %%ebx \n"
        "movl __temp3, %%eax \n"
        "movl %%eax, 16(%%ebp) \n"

        "skip_zmap_hd: \n"
          "addl $3, %%edi \n"
          "addl $88, %%ebp \n"
          "addl $8, %%ebx \n"
        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax\n"
        "jge l_gfx_h_seg11d \n"

        "movl __esp, %%esp \n"
        "movl __ebp, %%ebp \n"
        "movl __edi, %%edi \n"
        "movl __esi, %%esi \n"
        "popal; sti \n"
      :
      :"m" (tex_init.dest), "m" (acc_start),
       "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (zm), "m" (p)
      :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base.valu += m.valu;
    _dtemp[1]++;
  }
//struct bmp_type {
//  unsigned char* texel; // 0
//  unsigned long w_shift; // 4
//  unsigned long x_bitmask, y_bitmask; // 8, 12
//  int width, height;
//};
}


inline void tri_group_class::draw_vert_pos_sect_dither
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, p, zm;
  static int x, y;
  static sng m;  
  static dbl lef, rig, d_lef, d_rig;
  
  p = group_num;
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  m = cv_sng(scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(y + .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      sng scan = scan_base - cv_sng(scan_m * (lef.whole + .5));
      f = scan.valu & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan.valu >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      x = lef.whole;
      zm = (int)zmap + (lef.whole + _screen.x_res * y) * 8 + 4;
      
      asm (
        "pushal \n"
        "push %%ebp \n"
        "push %%esi \n"
        "movl %%esp, __esp \n"
        "movl %5, %%esi \n"
        "movl %2, %%edi \n"
        "movl %1, %%ebp \n"
        "movl %8, %%ebx \n"        
        "gfx_v_line_pos_rep1d: \n"
      
        "movl (%%ebx), %%eax \n"
        "movl %9, %%edx \n"        
        "cmpl %%eax, %%edx \n"
        "jne skip_zmap_vpd \n"
        "movl %%ebx, __ebx \n"
        // update u & v
        "movl 16(%%ebp), %%eax \n"

        "movl %7, %%edx \n"
        "cmp %%eax, %%edx \n"
        "je gfx_v_line_pos_skip1d \n"
        "  leal gfx_v_line_pos_skip1d, %%eax \n"
        "  movl %%eax, __ret \n"
        "  leal %7, %%ebx \n"
        "  jmp update_v_acc \n"
        "gfx_v_line_pos_skip1d: \n"

        // calc u
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
        "movl %%ebx, __temp1 \n"

        // calc v
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
        "addl %%edx, %%ebx \n" // ebx = v 16.16
        
        "movl __temp1, %%eax \n" // eax = u 16.16

        "movl %%esi, __temp2 \n"
        "subl $32768, %%eax \n"
        "subl $32768, %%ebx \n"
        "xorl %%esi, %%esi \n"
        "movb %%bh, %%al \n"
        "movw %%ax, %%si \n"
        
        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%eax \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%eax \n"
        "shrl $16, %%ebx \n"
        "movl __bmp + 4, %%ecx \n"
        "movl (%%ecx, %%ebx, 4), %%ebx \n"

    "movl __bmp + 12, %%esp \n"
    "cmpl __bmp + 20, %%ebx \n"
    "jne last_y \n"
    "  movl __bmp + 24, %%esp \n"
    "last_y: \n"

    "cmpl __bmp + 16, %%eax \n"
    "je last_x \n"
    "  addl %%eax, %%ebx \n"
    "  leal (%%ebx, %%ebx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "jmp endif1 \n"

    "last_x: \n"
    "  addl %%eax, %%ebx \n"
    "  leal (%%ebx, %%ebx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    "  movl __bmp + 12, %%ecx \n"
    "  subl %%ecx, %%ebx \n"
    
    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "endif1: \n"

    "movl __bk_light + 0, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 0(%%edi) \n"

    "movl __bk_light + 4, %%ecx \n"
    "movb %%ah, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 1(%%edi) \n"

    "shrl $16, %%eax \n"
    
    "movl __bk_light + 8, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 2(%%edi) \n"

//    "movb %%al, 0(%%edi) \n"
//    "movb %%ah, 1(%%edi) \n"
//    "shrl $16, %%eax \n"
//    "movb %%al, 2(%%edi) \n"

    "movl __temp2, %%esi \n"

        "movl %7, %%eax \n"
        "incl %%eax \n"
        "movl %%eax, 16(%%ebp) \n"
        
        "movl __ebx, %%ebx \n"
        "skip_zmap_vpd: \n"

        "movl %7, %%eax \n"
        "incl %%eax \n"
        "movl %%eax, %7 \n"
        
        "subw %6, %%si \n"
        "jnc gfx_v_line_pos_skip2d \n"
        "  subl $88, %%ebp \n"
        "gfx_v_line_pos_skip2d: \n"
        "addl $3, %%edi \n"
        "addl $8, %%ebx \n"
        
        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax \n"
        "jge gfx_v_line_pos_rep1d \n"

        "movl __esp, %%esp \n"
        "pop %%esi \n"
        "pop %%ebp \n"
        "popal \n"        
        :
        :"m" (tex_init.dest), "m" (acc_start),
         "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (m),
         "m" (x), "m" (zm), "m" (p)
        :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base++;
  }
/*
struct bmp_type {
  unsigned char* texel; // 0
  unsigned long w_shift; // 4
  unsigned long x_bitmask, y_bitmask; // 8, 12
  int width, height; // 16, 24
};
*/
}


inline void tri_group_class::draw_vert_neg_sect_dither
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f, p, zm;
  static int x, y;
  static sng m;  
  static dbl lef, rig, d_lef, d_rig;
  
  p = group_num;
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  m = cv_sng(-scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(y + .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      sng scan = scan_base - cv_sng(scan_m * (lef.whole + .5));
      f = scan.valu & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan.valu >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      x = lef.whole;
      zm = (int)zmap + (lef.whole + _screen.x_res * y) * 8 + 4;
      
      asm (
        "pushal \n"
        "push %%ebp \n"
        "push %%esi \n"
        "movl %%esp, __esp \n"
        "movl %5, %%esi \n"
        "movl %2, %%edi \n"
        "movl %1, %%ebp \n"
        "movl %8, %%ebx \n"
        "gfx_v_line_neg_rep1d: \n"
      
        "movl (%%ebx), %%eax \n"
        "movl %9, %%edx \n"
        "cmpl %%eax, %%edx \n"
        "jne skip_zmap_vnd \n"
        // update u & v
        "movl %%ebx, __ebx \n"
        "movl 16(%%ebp), %%eax \n"

        "movl %7, %%edx \n"
        "cmp %%eax, %%edx \n"
        "je gfx_v_line_neg_skip1d \n"
        "  leal gfx_v_line_neg_skip1d, %%eax \n"
        "  movl %%eax, __ret \n"
        "  leal %7, %%ebx \n"
        "  jmp update_v_acc \n"
        "gfx_v_line_neg_skip1d: \n"

        // calc u
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
        "movl %%ebx, __temp1 \n"

        // calc v
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
        "addl %%edx, %%ebx \n" // ebx = v 16.16
        
        "movl __temp1, %%eax \n" // eax = u 16.16

        "movl %%esi, __temp2 \n"
        "subl $32768, %%eax \n"
        "subl $32768, %%ebx \n"
        "xorl %%esi, %%esi \n"
        "movb %%bh, %%al \n"
        "movw %%ax, %%si \n"
        
        "movl __bmp + 8, %%esp \n"
        "andl %%esp, %%eax \n"
        "andl %%esp, %%ebx \n"
        "shrl $16, %%eax \n"
        "shrl $16, %%ebx \n"
        "movl __bmp + 4, %%ecx \n"
        "movl (%%ecx, %%ebx, 4), %%ebx \n"

    "movl __bmp + 12, %%esp \n"
    "cmpl __bmp + 20, %%ebx \n"
    "jne neg_last_y \n"
    "  movl __bmp + 24, %%esp \n"
    "neg_last_y: \n"

    "cmpl __bmp + 16, %%eax \n"
    "je neg_last_x \n"
    "  addl %%eax, %%ebx \n"
    "  leal (%%ebx, %%ebx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "jmp neg_endif1 \n"

    "neg_last_x: \n"
    "  addl %%eax, %%ebx \n"
    "  leal (%%ebx, %%ebx, 2), %%edx \n"
    "  movl __bmp + 0, %%ebx \n"
    "  addl %%edx, %%ebx \n"
    
    "  xorl %%eax, %%eax \n"
    "  xorl %%edx, %%edx \n"
    "  xorl %%ecx, %%ecx \n"

    // A
    "  movb 0 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  shll $16, %%eax \n"
    "  movb 0*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%al \n"
    "  movb 0*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%ah \n"

    // C
    "  movb 2 + __dither(,%%esi, 4), %%cl \n"
    "  movb 0*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 0*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 0*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"

    "  movl __bmp + 12, %%ecx \n"
    "  subl %%ecx, %%ebx \n"
    
    // B
    "  movb 1 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    
    // D
    "  movb 3 + __dither(,%%esi, 4), %%cl \n"
    "  movb 1*3 + 2(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  shll $16, %%edx \n"
    "  movb 1*3 + 0(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dl \n"
    "  movb 1*3 + 1(%%esp, %%ebx), %%ch \n"
    "  movb __mult_tbl(%%ecx), %%dh \n"
    "  addl %%edx, %%eax \n"
    "neg_endif1: \n"

//    "movb %%al, 0(%%edi) \n"
//    "movb %%ah, 1(%%edi) \n"
//    "shrl $16, %%eax \n"
//    "movb %%al, 2(%%edi) \n"

    "movl __bk_light + 0, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 0(%%edi) \n"

    "movl __bk_light + 4, %%ecx \n"
    "movb %%ah, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 1(%%edi) \n"

    "shrl $16, %%eax \n"
    
    "movl __bk_light + 8, %%ecx \n"
    "movb %%al, %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 2(%%edi) \n"

    "movl __temp2, %%esi \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, 16(%%ebp) \n"

        "movl __ebx, %%ebx \n"
        "skip_zmap_vnd: \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, %7 \n"
        
        "addw %6, %%si \n"
        "jnc gfx_v_line_neg_skip2d \n"
        "  addl $88, %%ebp \n"
        "gfx_v_line_neg_skip2d: \n"
        "addl $8, %%ebx \n"
        "addl $3, %%edi \n"
        
        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax \n"
        "jge gfx_v_line_neg_rep1d \n"

        "movl __esp, %%esp \n"
        "pop %%esi \n"
        "pop %%ebp \n"
        "popal \n"        
        :
        :"m" (tex_init.dest), "m" (acc_start),
         "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (m),
         "m" (x), "m" (zm), "m" (p)
        :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base++;
  }
/*
struct bmp_type {
  unsigned char* texel; // 0
  unsigned long w_shift; // 4
  unsigned long x_bitmask, y_bitmask; // 8, 12
  int width, height; // 16, 24
};
*/
}


/*
inline void tri_group_class::draw_vert_neg_sect_dither
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f;
  static int x, y;
  static sng m;  
  static dbl lef, rig, d_lef, d_rig;
  
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  m = cv_sng(-scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(y + .5);

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      sng scan = scan_base - cv_sng(scan_m * (lef.whole + .5));
      f = scan.valu & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan.valu >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      x = lef.whole;

      asm (
        "pushal \n"
        "push %%ebp \n"
        "push %%esi \n"
        "movl %%esp, __esp \n"
        "movl %5, %%esi \n"
        "movl %2, %%edi \n"
        "movl %1, %%ebp \n"
        "gfx_v_line_neg_rep1d: \n"
      
        // update u & v
        "movl 16(%%ebp), %%eax \n"

        "movl %7, %%edx \n"
        "cmp %%eax, %%edx \n"
        "je gfx_v_line_neg_skip1d \n"
        "  leal gfx_v_line_neg_skip1d, %%eax \n"
        "  movl %%eax, __ret \n"
        "  leal %7, %%ebx \n"
        "  jmp update_v_acc \n"
        "gfx_v_line_neg_skip1d: \n"

        // calc u
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
        "andl __bmp + 8, %%ebx \n"
        "shrl $16, %%ebx \n"
        "movl %%ebx, __temp1 \n"

        // calc v
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
        "andl __bmp + 12, %%ebx \n"
        "movl __bmp + 4, %%ecx \n"
        "shrl %%cl, %%ebx \n"
        "movl __temp1, %%eax \n"
        "addl %%ebx, %%eax \n"
        "leal (%%eax, %%eax, 2), %%eax \n"
        "movl __bmp, %%ebx \n"

        "movl __bk_light + 0, %%ecx \n"
        "movb 0(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 0(%%edi) \n"
    
        "movl __bk_light + 4, %%ecx \n"
        "movb 1(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 1(%%edi) \n"
    
        "movl __bk_light + 8, %%ecx \n"
        "movb 2(%%ebx, %%eax), %%cl \n"
        "movb __mult_tbl(%%ecx), %%cl \n"
        "movb %%cl, 2(%%edi) \n"

        "addl $3, %%edi \n"

        "movl %7, %%eax \n"
        "inc %%eax \n"
        "movl %%eax, 16(%%ebp) \n"
        "movl %%eax, %7 \n"
      
        "addw %6, %%si \n"
        "jnc gfx_v_line_neg_skip2d \n"
        "  addl $88, %%ebp \n"
        "gfx_v_line_neg_skip2d: \n"

        "movl %3, %%eax \n"
        "cmpl %%edi, %%eax\n"
        "jge gfx_v_line_neg_rep1d \n"

        "movl __esp, %%esp \n"
        "pop %%esi \n"
        "pop %%ebp \n"
        "popal \n"        
        :
        :"m" (tex_init.dest), "m" (acc_start),
         "m" (line_start), "m" (line_end), "m" (y), "m" (f), "m" (m),
         "m" (x)
        :"memory"
      );
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base++;
  }
}
*/

/*
inline void tri_group_class::draw_vert_neg_sect
  (double t, double b, lin_relat l, lin_relat r)
{
  static unsigned long line_start, line_end;
  static unsigned long acc_start, f;
  static int y;
  dbl lef, rig, d_lef, d_rig;
  
  y = sect_top(t);
  lef = init_lef_edge (l, y);
  rig = init_rig_edge (r, y);
  d_lef = cv_dbl(l.m);
  d_rig = cv_dbl(r.m);

  double scan_m = tex_init.scan_m;
  sng m = cv_sng(scan_m);
  int yc = y * _screen.x_res;
  sng scan_base = cv_sng(scan_m * (y + .5) - .5);
  _dtemp[1] = double(y) + .5;

  for (; y <= sect_bot(b); y++, yc += _screen.x_res) {
    if (lef.whole <= rig.whole) {
      long scan = (lef.whole << 16) - scan_base.valu;
      f = scan & 0x0000FFFF;
      acc_start = (int)tex_init.acc + (scan >> 16) * 88;
      line_start = (int)tex_init.dest + (lef.whole + yc) * 3;
      line_end = (int)tex_init.dest + (rig.whole + yc) * 3;
      gfx_v_line_neg (y, lef.whole, rig.whole, scan_m);
    }
    lef.valu += d_lef.valu;
    rig.valu += d_rig.valu;
    scan_base.valu += m.valu;
    _dtemp[1]++;
  }

struct bmp_type {
  unsigned char* texel; // 0
  unsigned long w_shift; // 4
  unsigned long x_bitmask, y_bitmask; // 8, 12
  int width, height; // 16, 24
};

}
*/

void update_accumulator()
{
  asm volatile (
    "update_acc: \n"
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
    "jmp *__ret \n"

    "update_next_acc: \n"
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
    "jmp *__ret \n"
    :
    :
    :"memory"
  );
}


void update_vert_acc()
{
   asm (
     "update_v_acc: \n"
     "  fild (%%ebx) \n" // load x
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

     "  fild (%%ebx) \n" // load x
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
     "  jmp *__ret \n"
     :
     :
     :"memory"
   );
}


/*
void tri_group_class::gfx_v_line_neg (int y, int x1, int x2, double scan_m)
{
  if (x1 > x2)
    return;

  sng m = cv_sng(-scan_m);
  long s = int(65536 * (double(y + .5) - double(x1 + .5) * scan_m));
  int scrn = (int)tex_init.dest + (unsigned int)3 * ((unsigned int)x1 + (unsigned int)y * (unsigned int)320);
  _temp3 = x1;
  _temp4 = x2 + 1;//int(&_texdat[DAT_MID + s_scan.whole]);
  _temp5 = m.valu;
  _temp6 = s & 65535;
  _temp7 = int(tex_init.acc) + (s >> 16) * 88;
  
    asm (
      "pushal \n"
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
      "  leal gfx_v_line_neg_skip1, %%eax \n"
      "  movl %%eax, __ret \n"
      "  leal __temp3, %%ebx \n"
      "  jmp update_v_acc \n"
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
      "movl __bmp, %%ebx \n"

    "  movl __bk_light + 0, %%ecx \n"
    "  movb 0(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 0(%%edi) \n"
    
    "  movl __bk_light + 4, %%ecx \n"
    "  movb 1(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 1(%%edi) \n"
    
    "  movl __bk_light + 8, %%ecx \n"
    "  movb 2(%%ebx, %%eax), %%cl \n"
    "  movb __mult_tbl(%%ecx), %%cl \n"
    "  movb %%cl, 2(%%edi) \n"

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
      "popal \n"
    :
    :"g" (_temp7), "g" (_temp6), "g" (_bmp.texel), "g" (scrn)
    :"eax", "ebx", "ecx", "edx", "edi", "esi", "ebp", "memory"
    );
}
*/

/*
void tri_group_class::gfx_h_seg1 ()
{
  asm volatile (
    "movl %%edi, __edi \n"
    "movl %%ebp, __ebp \n"
    "movl %%esp, __esp \n"
    "movl %%esi, __esi \n"
    "movl %0, %%edi \n"
    "addl __scrn, %%edi \n"
    "movl __f, %%esi \n"
    "movl %1, %%ebp \n"
    "addl __temp4, %%ebp \n"

  "l_gfx_h_seg11: \n"
    "movl __y, %%ecx \n"
    
    "movl 0 + 16(%%ebp), %%eax \n"
    "cmpl %%ecx, %%eax \n"
    "je s_gfx_hseg12 \n"
    "  leal s_gfx_hseg12, %%eax \n"
    "  movl %%eax, __ret \n"
    "  jmp update_acc \n"
    "s_gfx_hseg12: \n"

    "movl 88 + 16(%%ebp), %%eax \n"
    "cmpl %%ecx, %%eax \n"
    "je s_gfx_hseg13 \n"
    "  leal s_gfx_hseg13, %%eax \n"
    "  movl %%eax, __ret \n"
    "  jmp update_next_acc \n"
    "s_gfx_hseg13: \n"

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
    "addl %%edx, %%ecx \n"

    "movl %%ebx, %%eax \n"
    "imull %%esi \n"
    "shll $16, %%edx \n"
    "shrl $16, %%eax \n"
    "addl %%eax, %%edx \n"
    "addl %%esp, %%edx \n"

    "andl $63 << 16, %%ecx \n"
    "andl $63 << 16, %%edx \n"
    "shrl $16, %%ecx \n"
    "shrl $16 - 6, %%edx \n"
    "addl %%ecx, %%edx \n"
    "movl __temp5, %%ebx \n"
    "leal (%%edx, %%edx, 2), %%edx \n"
    "movl __temp3, %%eax \n"

//    "xorl %%ecx, %%ecx \n"

    "movl $128 << 8, %%ecx \n"
    "movb 0(%%ebx, %%edx), %%cl \n"    
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 0(%%edi) \n"
    
    "movl $128 << 8, %%ecx \n"
    "movb 1(%%ebx, %%edx), %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 1(%%edi) \n"
    
    "movl $128 << 8, %%ecx \n"
    "movb 2(%%ebx, %%edx), %%cl \n"
    "movb __mult_tbl(%%ecx), %%cl \n"
    "movb %%cl, 2(%%edi) \n"

    "movl %%eax, 16(%%ebp) \n"
    "movl __temp7, %%ebx \n"

    "s_gfx_h_seg11: \n"

    "addl $8, %%ebx \n"
    "addl $3, %%edi \n"
    "addl $88, %%ebp \n"
  "decl __temp2 \n"
  "jge l_gfx_h_seg11 \n"

  "movl __esp, %%esp \n"
  "movl __ebp, %%ebp \n"
  "movl __edi, %%edi \n"
  "movl __esi, %%esi \n"

  "jmp skip_all \n"

      "update_acc: \n"
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
      "jmp *__ret \n"

      "update_next_acc: \n"
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
      "jmp *__ret \n"
      
      "skip_all: \n"
  :
  :"g" (tex_init.dest), "g" (tex_init.acc)
  :"memory", "eax", "ebx", "ecx", "edx"
  );
}
*/


inline void gen_mult_tbl()
{
  for (int a = 0; a < 256 * 8; a++)
    for (int b = 0; b < 256; b++)
      if (((a * b) >> 8) > 255)
        _mult_tbl[a][b] = 255;
      else
        _mult_tbl[a][b] =
        (unsigned char)floor(((double)a * (double)b) / (double)256);
}


void create_nx88_tbl()
{
  for (int i = 0; i < 1024; i++)
    _nx88[i] = i * 88;
}


void create_dither_tbl()
{
  int x, y;
  double v = 0, h = 0;
  double n1, n2, n3, n4;
  double ifactor = 1.00;//166666666666;
  
  for (y = 0; y < 256; y++, v += 1.0/256.0)
    for (x = 0, h = 0; x < 256; x++, h += 1.0/256.0) {
      n1 = ifactor * 255.0 * (1 - h) * (1 - v);
      n2 = ifactor * 255.0 * h * (1 - v);
      n3 = ifactor * 255.0 * (1 - h) * v;
      n4 = ifactor * 255.0 * h * v;

      if (n1 > 255)
        n1 = 255;
      if (n2 > 255)
        n2 = 255;
      if (n3 > 255)
        n3 = 255;
      if (n4 > 255)
        n4 = 255;
        
      _dither[x][y][0] = (unsigned char)n1;
      _dither[x][y][1] = (unsigned char)n2;
      _dither[x][y][2] = (unsigned char)n3;
      _dither[x][y][3] = (unsigned char)n4;
    }
}


void create_dither_tbl2()
{
/*
  int x, y;
  double a, b;
  for (y = 0, b = 0; b < 1; b += 1 / 256.0, y++)
    for (x = 0, a = 0; a < 1; a += 1 / 256.0, x++) {
      double d1 = ((1 - a) * (1 - b));
      double d2 = (a * (1 - b));
      double d3 = ((1 - a) * b);
      double d4 = (a * b);

      double dt = d1 + d2 + d3 + d4;

      double o1 = 255.0 * (d1 / dt);
      double o2 = 255.0 * (d2 / dt);
      double o3 = 255.0 * (d3 / dt);
      double o4 = 255.0 * (d4 / dt);

      _dither[x][y][0] = (unsigned char)o1;
      _dither[x][y][1] = (unsigned char)o2;
      _dither[x][y][2] = (unsigned char)o3;
      _dither[x][y][3] = (unsigned char)o4;
    }
*/
///*
  for (double y = 0; y < 256; y++)
    for (double x = 0; x < 256; x++) {
      _dither[(int)x][(int)y][0] = (unsigned char) (((255 - x) * (255 - y)) / 256.0);
      _dither[(int)x][(int)y][1] = (unsigned char) ((x * (255 - y)) / 256.0);
      _dither[(int)x][(int)y][2] = (unsigned char) (((255 - x) * y) / 256.0);
      _dither[(int)x][(int)y][3] = (unsigned char) ((x * y) / 256.0);
    }
//*/
}


inline dbl init_lef_edge (lin_relat edge, int ind)
{
  return cv_dbl((edge.m * (ind + .5) + edge.b + .5));
}


inline dbl init_rig_edge (lin_relat edge, int ind)
{
  return cv_dbl((edge.m * (ind + .5) + edge.b - .5));
}


inline int sect_top (double y)
{
  return ifloor(y + .5);
}


inline int sect_bot (double y)
{
  return ifloor(y - .5);
}


inline unsigned long bit_not (unsigned long n)
{
  int result;
  
  asm (
    "movl %1, %%eax \n"
    "notl %%eax \n"
    "movl %%eax, %0 \n"
    :"=g" (result)
    :"g" (n)
    :"eax"
  );

  return result;
}


void print_bin (int n)
{
  for (int i = 31; i >= 0; i--) {
    printf ("%i", bool(n & (1 << i)));
  }
}


inline long mk_rgb_light (int red, int green, int blue)
{
  return (red << 21) + (green << 10) + blue;
}


/*
          "testl $0xFFFF0000, %%eax \n"
          "jz skip_bmap_16_0a \n"
            "testl $0xFF000000, %%eax \n"
            "jz skip_bmap_8_0a \n"
              "testl $0xF0000000, %%eax \n"
              "jz skip_bmap_4_0a \n"
                "testl $(3 << (31 - 1)), %%eax \n"
                "jz skip_bmap_2_0a \n"
                  "testl $(1 << 31), %%eax \n"
                  "jz skip_bmap_1_0a \n"
                    "\n"
                  "skip_bmap_1_0a: \n"
                  "testl $(1 << 30), %%eax \n"
                  "jz skip_bmap_1_1a \n"
                    "\n"
                  "skip_bmap_1_1a: \n"
                "skip_bmap_2_0a: \n"
                "testl $(3 << (31 - 1)), %%eax \n"
                "jz skip_bmap_2_0a \n"
                  "testl $(1 << 31), %%eax \n"
                  "jz skip_bmap_1_0a \n"
                    "\n"
                  "skip_bmap_1_0a: \n"
                  "testl $(1 << 30), %%eax \n"
                  "jz skip_bmap_1_1a \n"
                    "\n"
                  "skip_bmap_1_1a: \n"
                "skip_bmap_2_0a: \n"
                
          "skip_bmap_16_0a: \n"
          "testl $0x0000FFFF, %%eax \n"
          "jz skip_bmap_16_1a \n"

          "skip_bmap_16_1a: \n"
*/
#endif //!INCLUDE_TEXLINE
