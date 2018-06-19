#include "gfx2d.h"

double inv_sqr_light (double x, double y, double z, double k, double i);
double calc_k (double bright, double d, double i);
int dig_nterp (double c, int x, int y);
void show_pic_scaled (int w, int h, long* dat, double scale, void* buff);
double get_inv_squared_rad (double k, double abs_i, double z, double i);
double get_inv_squared_max_z (double k, double abs_i, double i);


int main()
{
  double z = 1024;
  double i = 1;
  double k = calc_k (i, z, .01);

  set_vesa_mode(0x3);
  printf ("%f, ", get_inv_squared_rad (k, i, z, .01));
  printf ("%f \n", get_inv_squared_max_z (k, i, .01));
  getchar();
}
/*
  double a1 = 9.9517894;
  long a2 = int(floor(a1 * 256.0));
  double f = .025121389128;
  int n = int(floor(mul * f + .5));
  double o1 = a1 * f;
  double o2 = double(n * a2) / (mul * 256.0);
  double e = (o1 - o2) * 256.0;

  printf ("%f, %f: %f, %f", o1, o2, e, mul * 256.0 * 4.0);
  getchar();

}
*/
/*
  double fa1 = .5, fa2 = fa1;
  double fa = fa1;
  double df = (fa2 - fa1) / 50;
  double mul = 10;
  
  set_vesa_mode(TRUE_COLOR_800x600);
  
  for (int x = 0; x < 50; x++, fa += df)
    for (int y = 0; y < 50; y++) {
      int c = mul * dig_nterp (fa, x, y);
      pxl (x, y, 100 + c, 100 + c, 100 + c);
    }
  getchar();
  
  fa = fa1;
  int acc = 0;
  
  for (int x = 0; x < 50; x++, fa += df)
    for (int y = 0; y < 50; y++)
      acc += mul * dig_nterp (fa, x, y);

  set_vesa_mode(0x3);
  printf ("%f", acc / 2500.0);
  getchar();
  
  return 0;

  double z = 50;
  double i = 1;
  double k = calc_k (i, z, 1.25);

  set_vesa_mode(TRUE_COLOR_320x200);
  
  void* x2_buff = (void*)malloc (_screen.x_res * _screen.y_res * 3);
  
  double d = _screen.x_res;
  double cx = _screen.x_res / 2, cy = _screen.y_res / 2;

  for (i = 1; i >= 0; i-=0.05) {
//    k = calc_k (1, z, 1);
    double rad = get_inv_squared_rad (k, i, z, .01);
    clear_scrn_buff(x2_buff);
//    show_pic_scaled (_screen.x_res, _screen.y_res, dat, scale, x2_buff);
    for (int y = cy - rad; y < cy + rad; y++)
      for (int x = cx - rad; x < cx + rad; x++) {
        double a = .0 + inv_sqr_light(double(x) - cx, double(y) - cy, z, k, i);
        double c = int(a * 255);//dig_nterp (a * 255, x, y);
        pxl (x, y, c, c, c, x2_buff);
      }
    pxl (cx + rad, _screen.y_res / 2, 255, 255, 255, x2_buff);
//    copy_scrn_buff (x2_buff, 0);
    copy_scrn_buff (x2_buff, 0);
    getchar();
  }
  
  return 0;

  long* dat = (long*)malloc (_screen.x_res * _screen.y_res * 4);

  for (int y = 0; y < _screen.y_res; y++)
    for (int x = 0; x < _screen.x_res; x++) {
      double a = .0 + inv_sqr_light(double(x) - cx, double(y) - cy, z, k, i);
      long c = dig_nterp (a * 255, x, y);
      dat[(x + y * _screen.x_res)] = c;
//      dat[(x + y * _screen.x_res) * 3 + 1] = c;
//      dat[(x + y * _screen.x_res) * 3 + 2] = c;
//      pxl (x * s, y * s, c, c, c);
    }

  for (double scale = 1; scale > .01; scale -= .001) {
    clear_scrn_buff(x2_buff);
    show_pic_scaled (_screen.x_res, _screen.y_res, dat, scale, x2_buff);
    copy_scrn_buff (x2_buff, 0);
//    getchar();
  }
  return 0;
  
}


double inv_sqr_light (double x, double y, double z, double k, double i)
{
  double result = k * i * z * pow(x*x + y*y + z*z, -1.5);

  if (result > 1.0)
    return 1;
  else
    return result;
}


int dig_nterp (double c, int x, int y)
{
  double factor = (c - int(c));
  int a;
  
  if (x % 2)
    a = .5 / factor;
  else
    a = 0;
    
  if (int((y + a + 1) * factor) - int((y + a) * factor) == 1)
//    if (int((y + a + 1) * factor) - int((y + a) * factor) == 1)
    return c + 1;
//    else
//      return c;
  else
    return c;
}


void show_pic_scaled (int w, int h, long* dat, double scale, void* buff)
{
  long u, v;
  long du = 65536.0 / scale, dv = 65536.0 / scale;
  int x, y;
  double mul = 1.0;// / (65536.0 * 127.0);// / (65536.0 * 65536.0);
  
  for (v = 0, y = 0; y < h * scale; y++, v += dv)
    for (u = 0, x = 0; x < w * scale; x++, u += du) {
      int i = ((u >> 16) + (v >> 16) * w);
      int c = int(dat[i] * mul);
      if (c > 255)
        pxl (x, y, 255, 255, 255, buff);
      else
        pxl (x, y, c, c, c, buff);
    }
}
//*/


double get_inv_squared_rad (double k, double abs_i, double z, double i)
{
  return sqrt(pow((k * abs_i * z) / i, 2.0/3.0) - z * z);
}


double get_inv_squared_max_z (double k, double abs_i, double i)
{
  return sqrt((k * abs_i) / i);
}


double calc_k (double factor, double d, double i)
{
  return i / (factor * d * pow(d*d, -1.5));
}
