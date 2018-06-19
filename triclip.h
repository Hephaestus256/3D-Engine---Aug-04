#ifndef INCLUDE_TRICLIP
#define INCLUDE_TRICLIP

#include "math3d.h"
#include "texclass.h"

class view_clip_class {
public:
  int view_clip (tri_3d tri, view_type* v, temp_tri* result);
private:
  inline void edge_clip_left (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_right (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_top (point_3d p1, point_3d p2, point_3d p3);
  inline void edge_clip_bot (point_3d p1, point_3d p2, point_3d p3);
  inline void store_result (point_3d p1, point_3d p2, point_3d p3);
  int tri_ct;
  temp_tri* result;
  view_type* view;
};

inline void z_cut
  (temp_tri, double z, temp_tri* i, temp_tri* o, int* in_ct, int* out_ct);
inline point_3d clip_edge_lr (point_3d p1, point_3d p2, double view_m);
inline point_3d clip_edge_tb (point_3d p1, point_3d p2, double view_m);
inline void store_z_cut_result (temp_tri* result, int* ct,
  point_3d p1, point_3d p2, point_3d p3);


int view_clip_class::view_clip (tri_3d t3d, view_type* v, temp_tri* r)
{
  result = r;
  view = v;
  tri_ct = 0;
  
  if (t3d.p1->rel.z >= Z_CUTOFF)
    if (t3d.p2->rel.z >= Z_CUTOFF)
      if (t3d.p3->rel.z >= Z_CUTOFF)
        edge_clip_left (t3d.p1->rel, t3d.p2->rel, t3d.p3->rel);
      else {
        point_3d a = intrapolate_3d (t3d.p1->rel, t3d.p3->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p2->rel, t3d.p3->rel, Z_CUTOFF);
        edge_clip_left (t3d.p1->rel, b, a);
        edge_clip_left (t3d.p1->rel, t3d.p2->rel, b);
      }
    else
      if (t3d.p3->rel.z >= Z_CUTOFF) {
        point_3d a = intrapolate_3d (t3d.p1->rel, t3d.p2->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p2->rel, t3d.p3->rel, Z_CUTOFF);
        edge_clip_left (t3d.p1->rel, a, t3d.p3->rel);
        edge_clip_left (a, b, t3d.p3->rel);
      }
      else {
        point_3d a = intrapolate_3d (t3d.p1->rel, t3d.p3->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p1->rel, t3d.p2->rel, Z_CUTOFF);
        edge_clip_left (t3d.p1->rel, b, a);
      }
  else
    if (t3d.p2->rel.z >= Z_CUTOFF)
      if (t3d.p3->rel.z >= Z_CUTOFF) {
        point_3d a = intrapolate_3d (t3d.p1->rel, t3d.p2->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p1->rel, t3d.p3->rel, Z_CUTOFF);
        edge_clip_left (a, t3d.p2->rel, t3d.p3->rel);
        edge_clip_left (b, a, t3d.p3->rel);
      }
      else {
        point_3d a = intrapolate_3d (t3d.p1->rel, t3d.p2->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p2->rel, t3d.p3->rel, Z_CUTOFF);
        edge_clip_left (t3d.p2->rel, b, a);
      }
    else
      if (t3d.p3->rel.z >= Z_CUTOFF) {
        point_3d a = intrapolate_3d (t3d.p2->rel, t3d.p3->rel, Z_CUTOFF);
        point_3d b = intrapolate_3d (t3d.p1->rel, t3d.p3->rel, Z_CUTOFF);
        edge_clip_left (a, t3d.p3->rel, b);
      }

  return tri_ct;
}


inline void view_clip_class::edge_clip_left
  (point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.x >= view->vis_lm * p1.z)
    if (p2.x >= view->vis_lm * p2.z)
      if (p3.x >= view->vis_lm * p3.z)
        edge_clip_right (p1, p2, p3);
      else {
        point_3d a = clip_edge_lr (p1, p3, view->vis_lm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_lm);
        edge_clip_right (p1, b, a);
        edge_clip_right (p1, p2, b);
      }
    else
      if (p3.x >= view->vis_lm * p3.z) {
        point_3d a = clip_edge_lr (p1, p2, view->vis_lm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_lm);
        edge_clip_right (p1, a, p3);
        edge_clip_right (a, b, p3);
      }
      else {
        point_3d a = clip_edge_lr (p1, p3, view->vis_lm);
        point_3d b = clip_edge_lr (p1, p2, view->vis_lm);
        edge_clip_right (p1, b, a);
      }
  else
    if (p2.x >= view->vis_lm * p2.z)
      if (p3.x >= view->vis_lm * p3.z) {
        point_3d a = clip_edge_lr (p1, p2, view->vis_lm);
        point_3d b = clip_edge_lr (p1, p3, view->vis_lm);
        edge_clip_right (a, p2, p3);
        edge_clip_right (b, a, p3);
      }
      else {
        point_3d a = clip_edge_lr (p1, p2, view->vis_lm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_lm);
        edge_clip_right (p2, b, a);
      }
    else
      if (p3.x >= view->vis_lm * p3.z) {
        point_3d a = clip_edge_lr (p2, p3, view->vis_lm);
        point_3d b = clip_edge_lr (p1, p3, view->vis_lm);
        edge_clip_right (a, p3, b);
      }
}


inline void view_clip_class::edge_clip_right
  (point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.x <= view->vis_rm * p1.z)
    if (p2.x <= view->vis_rm * p2.z)
      if (p3.x <= view->vis_rm * p3.z)
        edge_clip_top (p1, p2, p3);
      else {
        point_3d a = clip_edge_lr (p1, p3, view->vis_rm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_rm);
        edge_clip_top (p1, b, a);
        edge_clip_top (p1, p2, b);
      }
    else
      if (p3.x <= view->vis_rm * p3.z) {
        point_3d a = clip_edge_lr (p1, p2, view->vis_rm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_rm);
        edge_clip_top (p1, a, p3);
        edge_clip_top (a, b, p3);
      }
      else {
        point_3d a = clip_edge_lr (p1, p3, view->vis_rm);
        point_3d b = clip_edge_lr (p1, p2, view->vis_rm);
        edge_clip_top (p1, b, a);
      }
  else
    if (p2.x <= view->vis_rm * p2.z)
      if (p3.x <= view->vis_rm * p3.z) {
        point_3d a = clip_edge_lr (p1, p2, view->vis_rm);
        point_3d b = clip_edge_lr (p1, p3, view->vis_rm);
        edge_clip_top (a, p2, p3);
        edge_clip_top (b, a, p3);
      }
      else {
        point_3d a = clip_edge_lr (p1, p2, view->vis_rm);
        point_3d b = clip_edge_lr (p2, p3, view->vis_rm);
        edge_clip_top (p2, b, a);
      }
    else
      if (p3.x <= view->vis_rm * p3.z) {
        point_3d a = clip_edge_lr (p2, p3, view->vis_rm);
        point_3d b = clip_edge_lr (p1, p3, view->vis_rm);
        edge_clip_top (a, p3, b);
      }
}


inline void view_clip_class::edge_clip_top
  (point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.y >= view->vis_tm * p1.z)
    if (p2.y >= view->vis_tm * p2.z)
      if (p3.y >= view->vis_tm * p3.z)
        edge_clip_bot (p1, p2, p3);
      else {
        point_3d a = clip_edge_tb (p1, p3, view->vis_tm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_tm);
        edge_clip_bot (p1, b, a);
        edge_clip_bot (p1, p2, b);
      }
    else
      if (p3.y >= view->vis_tm * p3.z) {
        point_3d a = clip_edge_tb (p1, p2, view->vis_tm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_tm);
        edge_clip_bot (p1, a, p3);
        edge_clip_bot (a, b, p3);
      }
      else {
        point_3d a = clip_edge_tb (p1, p3, view->vis_tm);
        point_3d b = clip_edge_tb (p1, p2, view->vis_tm);
        edge_clip_bot (p1, b, a);
      }
  else
    if (p2.y >= view->vis_tm * p2.z)
      if (p3.y >= view->vis_tm * p3.z) {
        point_3d a = clip_edge_tb (p1, p2, view->vis_tm);
        point_3d b = clip_edge_tb (p1, p3, view->vis_tm);
        edge_clip_bot (a, p2, p3);
        edge_clip_bot (b, a, p3);
      }
      else {
        point_3d a = clip_edge_tb (p1, p2, view->vis_tm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_tm);
        edge_clip_bot (p2, b, a);
      }
    else
      if (p3.y >= view->vis_tm * p3.z) {
        point_3d a = clip_edge_tb (p2, p3, view->vis_tm);
        point_3d b = clip_edge_tb (p1, p3, view->vis_tm);
        edge_clip_bot (a, p3, b);
      }
}


inline void view_clip_class::edge_clip_bot
  (point_3d p1, point_3d p2, point_3d p3)
{
  if (p1.y <= view->vis_bm * p1.z)
    if (p2.y <= view->vis_bm * p2.z)
      if (p3.y <= view->vis_bm * p3.z)
        store_result (p1, p2, p3);
      else {
        point_3d a = clip_edge_tb (p1, p3, view->vis_bm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_bm);
        store_result (p1, b, a);
        store_result (p1, p2, b);
      }
    else
      if (p3.y <= view->vis_bm * p3.z) {
        point_3d a = clip_edge_tb (p1, p2, view->vis_bm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_bm);
        store_result (p1, a, p3);
        store_result (a, b, p3);
      }
      else {
        point_3d a = clip_edge_tb (p1, p3, view->vis_bm);
        point_3d b = clip_edge_tb (p1, p2, view->vis_bm);
        store_result (p1, b, a);
      }
  else
    if (p2.y <= view->vis_bm * p2.z)
      if (p3.y <= view->vis_bm * p3.z) {
        point_3d a = clip_edge_tb (p1, p2, view->vis_bm);
        point_3d b = clip_edge_tb (p1, p3, view->vis_bm);
        store_result (a, p2, p3);
        store_result (b, a, p3);
      }
      else {
        point_3d a = clip_edge_tb (p1, p2, view->vis_bm);
        point_3d b = clip_edge_tb (p2, p3, view->vis_bm);
        store_result (p2, b, a);
      }
    else
      if (p3.y <= view->vis_bm * p3.z) {
        point_3d a = clip_edge_tb (p2, p3, view->vis_bm);
        point_3d b = clip_edge_tb (p1, p3, view->vis_bm);
        store_result (a, p3, b);
      }
}


inline void view_clip_class::store_result
  (point_3d p1, point_3d p2, point_3d p3)
{
  result[tri_ct].p1 = p1;
  result[tri_ct].p2 = p2;
  result[tri_ct].p3 = p3;
  tri_ct++;
}


inline point_3d clip_edge_lr (point_3d p1, point_3d p2, double view_m)
{
  point_3d temp;

  if (is_approx_zero (p1.z - p2.z)) {
    temp.x = view_m * p1.z;
    temp.z = p1.z;
    double m = (p2.y - p1.y) / (p2.x - p1.x);
    double b = p1.y - m * p1.x;
    temp.y = m * temp.x + b;
  }
  else {
    line_3d l3d = calc_line_3d (p1, p2);
    temp.z = l3d.bx / (view_m - l3d.mx);
    temp.x = l3d.mx * temp.z + l3d.bx;
    temp.y = l3d.my * temp.z + l3d.by;
  }
  
  return temp;
}


inline point_3d clip_edge_tb (point_3d p1, point_3d p2, double view_m)
{
  point_3d temp;

  if (is_approx_zero (p1.z - p2.z)) {
    temp.y = view_m * p1.z;
    temp.z = p1.z;
    double m = (p2.x - p1.x) / (p2.y - p1.y);
    double b = p1.x - m * p1.y;
    temp.x = m * temp.y + b;
  }
  else {
    line_3d l3d = calc_line_3d (p1, p2);
    temp.z = l3d.by / (view_m - l3d.my);
    temp.x = l3d.mx * temp.z + l3d.bx;
    temp.y = l3d.my * temp.z + l3d.by;
  }
  
  return temp;
}


inline void z_cut
  (temp_tri tri, double z, temp_tri* i, temp_tri* o, int* in_ct, int* out_ct)
{
  *in_ct = 0;
  *out_ct = 0;
  
  if (tri.p1.z >= z)
    if (tri.p2.z >= z)
      if (tri.p3.z >= z)
        store_z_cut_result (o, out_ct, tri.p1, tri.p2, tri.p3);
      else {
        point_3d a = intrapolate_3d (tri.p1, tri.p3, z);
        point_3d b = intrapolate_3d (tri.p2, tri.p3, z);
        store_z_cut_result (o, out_ct, tri.p1, b, a);
        store_z_cut_result (o, out_ct, tri.p1, tri.p2, b);
        store_z_cut_result (i, in_ct, a, b, tri.p3);
      }
    else
      if (tri.p3.z >= z) {
        point_3d a = intrapolate_3d (tri.p1, tri.p2, z);
        point_3d b = intrapolate_3d (tri.p2, tri.p3, z);
        store_z_cut_result (o, out_ct, tri.p1, a, tri.p3);
        store_z_cut_result (o, out_ct, a, b, tri.p3);
        store_z_cut_result (i, in_ct, b, a, tri.p2);
      }
      else {
        point_3d a = intrapolate_3d (tri.p1, tri.p3, z);
        point_3d b = intrapolate_3d (tri.p1, tri.p2, z);
        store_z_cut_result (o, out_ct, tri.p1, b, a);
        store_z_cut_result (i, in_ct, tri.p2, tri.p3, a);
        store_z_cut_result (i, in_ct, b, tri.p2, a);
      }
  else
    if (tri.p2.z >= z)
      if (tri.p3.z >= z) {
        point_3d a = intrapolate_3d (tri.p1, tri.p2, z);
        point_3d b = intrapolate_3d (tri.p1, tri.p3, z);
        store_z_cut_result (o, out_ct, a, tri.p2, tri.p3);
        store_z_cut_result (o, out_ct, b, a, tri.p3);
        store_z_cut_result (i, in_ct, a, b, tri.p1);
      }
      else {
        point_3d a = intrapolate_3d (tri.p1, tri.p2, z);
        point_3d b = intrapolate_3d (tri.p2, tri.p3, z);
        store_z_cut_result (o, out_ct, tri.p2, b, a);
        store_z_cut_result (i, in_ct, a, b, tri.p1);
        store_z_cut_result (i, in_ct, b, tri.p3, tri.p1);
      }
    else
      if (tri.p3.z >= z) {
        point_3d a = intrapolate_3d (tri.p2, tri.p3, z);
        point_3d b = intrapolate_3d (tri.p1, tri.p3, z);
        store_z_cut_result (o, out_ct, a, tri.p3, b);
        store_z_cut_result (i, in_ct, a, b, tri.p2);
        store_z_cut_result (i, in_ct, b, tri.p1, tri.p2);
      }
      else
        store_z_cut_result (i, in_ct, tri.p1, tri.p2, tri.p3);
}


inline void store_z_cut_result (temp_tri* result, int* ct,
  point_3d p1, point_3d p2, point_3d p3)
{
  result[*ct].p1 = p1;
  result[*ct].p2 = p2;
  result[*ct].p3 = p3;
  (*ct)++;
}


#endif // !INCLUDE_TRICLIP
