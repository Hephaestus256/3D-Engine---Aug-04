#ifndef INCLUDE_MATH_2D
#define INCLUDE_MATH_2D

#include <math.h>

#ifndef CMP_PRECISION
#define CMP_PRECISION (1.0 / 1000000.0)
#endif

struct point_2d {
  double x, y;
};

struct lin_relat {
  double m, b;
};

struct line_2d {
  double m, b;
};

struct int_pair {
  int x, y;
};

struct dbl_pair {
  double x, y;
};

struct line_equat {
  double m, b;
  double inv_m, inv_b;
};

struct tri_2d {
  point_2d p1, p2, p3;
};

inline bool equivalent_angles (double a1, double a2);
inline point_2d p2d (double x, double y);
inline bool is_real_num (double n);
inline bool is_pos_inf (double n);
inline bool is_neg_inf (double n);
inline bool is_abs_inf (double n);
inline bool is_nan (double n);
inline double nan ();
inline bool is_approx_zero (double n);
inline bool is_approx_abs_inf (double n);
inline double rad_to_deg (double rad);
inline double deg_to_rad (double deg);
inline bool approx_equal (double a, double b);
inline bool approx_greater (double a, double b);
inline bool approx_lesser (double a, double b);
inline bool approx_lesser_inc (double a, double b);
inline bool approx_greater_inc (double a, double b);
inline double square(double n);
inline void rotate
  (double* x, double* y, double center_x, double center_y, double angle);
inline void calc_line (line_equat* line, point_2d point1, point_2d point2);
inline bool is_real_num (double n);
inline bool quadratic
  (double a, double b, double c, double* r1, double* r2);
inline bool approx_zero (double n);
inline void simplify_angle (double* a);
inline double simplify_angle (double a);
inline point_2d intersect_vert_line (point_2d point, line_equat line);
inline point_2d line_intersect (line_equat line1, line_equat line2);
inline double point_rel_line (point_2d p, line_equat line);
inline void calc_2d_dydx_line
  (line_equat* line, point_2d p1, point_2d p2);
inline void calc_2d_dxdy_line
  (line_equat* line, point_2d p1, point_2d p2);
inline lin_relat calc_2d_dxdy_line (point_2d p1, point_2d p2);
inline double find_dydx_b (point_2d point, double m);
inline double find_dxdy_b (point_2d point, double m);
inline int log_2 (int n);
inline bool approx_between (double n, double a, double b);
inline bool approx_between_sorted (double n, double min, double max);


inline point_2d p2d (double x, double y)
{
  point_2d temp;

  temp.x = x;
  temp.y = y;

  return temp;
}


inline bool approx_lesser (double a, double b)
{
  if ((b - a) > CMP_PRECISION)
    return true;
  else
    return false;
}


inline bool approx_lesser_inc (double a, double b)
{
  if ((b - a) > -CMP_PRECISION)
    return true;
  else
    return false;
}


inline bool approx_greater (double a, double b)
{
  if ((a - b) > CMP_PRECISION)
    return true;
  else
    return false;
}


inline bool approx_greater_inc (double a, double b)
{
  if ((a - b) > -CMP_PRECISION)
    return true;
  else
    return false;
}


inline bool approx_equal (double a, double b)
{
  if (fabs(a - b) <= CMP_PRECISION)
    return true;
  else
    return false;
}


inline double deg_to_rad (double deg)
{
  return deg * (PI / 180.0);
}


inline void rotate
  (double* x, double* y, double center_x, double center_y, double angle)
{
  if (*x == center_x)
    if (*y == center_y)
      return;

  double dist = sqrt(square(*x - center_x) + square(*y - center_y));
  double a = atan2(*y - center_y, *x - center_x);

  *x = center_x + cos(angle + a) * dist;
  *y = center_y + sin(angle + a) * dist;
}


inline double square(double n)
{
  return n * n;
}


inline double rad_to_deg (double rad)
{
  return rad / (PI / 180.0);
}


inline bool is_approx_zero (double n)
{
  if (fabs(n) < CMP_PRECISION)
    return true;
  else
    return false;
}


inline bool is_real_num (double n)
{
  double a = 1, b = 0;
  
  double inf = a / b;
  double nan = b / b;
  
  if (n == inf)
    return false;
  else if (n == -inf)
    return false;
  else if (n == nan)
    return false;
  else
    return true;
}


inline bool is_pos_inf (double n)
{
  double a = 1, b = 0;
  return n == a / b;
}


inline bool is_neg_inf (double n)
{
  double a = -1, b = 0;
  return n == a / b;
}


inline bool is_abs_inf (double n)
{
  double a = 1, b = 0;
  return fabs(n) == a / b;
}


inline bool is_nan (double n)
{
  double a = 0;
  return n == a / a;
}


inline double nan ()
{
  double a = 0;
  return a / a;
}


inline bool is_approx_abs_inf (double n)
{
  if (n > 0)
    return n > 1 / CMP_PRECISION;
  else
    return n < -1 / CMP_PRECISION;
}


inline void calc_line (line_equat* line, point_2d point1, point_2d point2)
{
  line->m = (point2.y - point1.y) / (point2.x - point1.x);
  line->b = point1.y - point1.x * line->m;
  line->inv_m = (point2.x - point1.x) / (point2.y - point1.y);
  line->inv_b = point1.x - point1.y * line->inv_m;
}


inline bool quadratic
  (double a, double b, double c, double* r1, double* r2)
{
  double x = sqrt(b * b - 4 * a * c);
  double den = .5 / a;
  *r1 = (x - b) * den;
  *r2 = -(x + b) * den;
  return b * b - 4 * a * c >= 0;
}


inline bool approx_zero (double n)
{
  return (fabs(n) < CMP_PRECISION);
}


inline bool approx_between_sorted (double n, double min, double max)
{
  if (approx_greater_inc (n, min))
    if (approx_lesser_inc (n, max))
      return true;
    else
      return false;
  else
    return false;
}


inline bool approx_between (double n, double a, double b)
{
  if (approx_equal (a, n))
    return true;
  else if (approx_equal (b, n))
    return true;
  else
    if (a < b)
      if (n < a)
        return false;
      else if (n > b)
        return false;
      else
        return true;
    else
      if (n < b)
        return false;
      else if (n > a)
        return false;
      else
        return true;
        
}


inline void calc_2d_dxdy_line
  (line_equat* line, point_2d p1, point_2d p2)
{
  line->m = (p2.x - p1.x) / (p2.y - p1.y);
  line->b = p1.x - p1.y * line->m;
}


inline lin_relat calc_2d_dxdy_line (point_2d p1, point_2d p2)
{
  lin_relat line;
  
  line.m = (p2.x - p1.x) / (p2.y - p1.y);
  line.b = p1.x - p1.y * line.m;

  return line;
}


inline void calc_2d_dydx_line
  (line_equat* line, point_2d p1, point_2d p2)
{
  line->m = (p2.y - p1.y) / (p2.x - p1.x);

  if (is_abs_inf(line->m))
    line->b = p1.x;
  else
    line->b = p1.y - p1.x * line->m;
}


inline double point_rel_line (point_2d p, line_equat line)
{
  return p.y - (line.m * p.x + line.b);
}


inline double find_dydx_b (point_2d point, double m)
{
  return point.y - point.x * m;
}


inline double find_dxdy_b (point_2d point, double m)
{
  return point.x - point.y * m;
}


inline point_2d line_intersect (line_equat line1, line_equat line2)
{
  point_2d p;

  if (is_pos_inf(line2.m)) {
    p.x = line2.b;
    p.y = line1.m * p.x + line1.b;
  }
  else if (is_neg_inf(line2.m)) {
    p.x = line2.b;
    p.y = line1.m * p.x + line1.b;
  }
  else {
    p.x = (line2.b - line1.b) / (line1.m - line2.m);
    p.y = line1.m * p.x + line1.b;
  }
  
  return p;
}


inline point_2d intersect_vert_line (point_2d point, line_equat line)
{
  point_2d p;

  p.x = point.x;
  p.y = line.m * point.x + line.b;
  
  return p;
}


inline double simplify_angle (double a)
{
  // returns the equivalent angle between 0 and 2pi
  
  double temp = atan2(sin(a), cos(a));
  if (temp < 0)
    return temp + 2 * PI;
  else
    return temp;
    
//  return a - (PI * 2 * floor(a / (PI * 2)));
}


inline void simplify_angle (double* a)
{
  // converts the angle to the equivalent angle between 0 and 2pi
  
  *a -= (PI * 2 * floor(*a / (PI * 2)));
}


inline int log_2 (int n)
{
  for (int i = 0; i < 32; i++)
    if ((1 << i) == n)
      return i;
}



inline bool equivalent_angles (double a1, double a2)
{
  if (approx_equal (sin (a1), sin(a2)))
    if (approx_equal (cos (a1), cos(a2)))
      return true;
    else
      return false;
  else
    return false;
}


#endif // !INCLUDE_MATH_2D
