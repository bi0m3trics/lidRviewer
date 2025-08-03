#ifndef DRAWER_H
#define DRAWER_H

#include <Rcpp.h>
#include <SDL2/SDL.h>

#include "Octree.h"
#include "camera.h"

using namespace Rcpp;

enum Attribute{Z, I, RGB, CLASS};
enum RenderMode{RENDER_POINTS, RENDER_LINES, RENDER_MIXED};

class Drawer
{
public:
  Drawer(SDL_Window*, DataFrame, std::string hnof);
  bool draw();
  void resize();
  void setPointSize(float);
  void setAttribute(Attribute x);
  void display_hide_spatial_index() { draw_index = !draw_index; camera.changed = true; };
  void display_hide_edl() { lightning = !lightning; camera.changed = true; };
  void point_size_plus() { point_size++; camera.changed = true; };
  void point_size_minus() { point_size--; camera.changed = true; };
  void budget_plus() { point_budget += 500000; camera.changed = true; };
  void budget_minus() { if (point_budget > 500000) point_budget -= 500000; camera.changed = true; };
  Camera camera;
  Octree index;

  float point_size;
  bool lightning;

private:
  void edl();
  void draw_points();  // New: separated point rendering
  void draw_lines();   // New: line rendering
  bool is_visible(const Node& octant);
  void compute_cell_visibility();
  void query_rendered_point();
  void traverse_and_collect(const Key& key, std::vector<Node*>& visible_octants);
  void init_viewport();

  bool draw_index;
  uint32_t npoints;
  int point_budget;
  int rgb_norm;

  double minx;
  double miny;
  double minz;
  double maxx;
  double maxy;
  double maxz;
  double xcenter;
  double ycenter;
  double zcenter;
  double xrange;
  double yrange;
  double zrange;
  double range;
  double zqmin;
  double zqmax;
  double minattr;
  double maxattr;
  double attrrange;

  DataFrame df;
  NumericVector x;
  NumericVector y;
  NumericVector z;
  IntegerVector r;
  IntegerVector g;
  IntegerVector b;
  NumericVector line_id;  // New: line segment IDs for line rendering (numeric to handle NA values)

  IntegerVector attri;
  NumericVector attrd;

  Attribute attr;
  RenderMode render_mode;  // New: rendering mode
  bool has_line_data;      // New: whether we have line data
  bool has_point_data;     // New: whether we have point data
  std::vector<int> pp;
  std::vector<Node*> visible_octants;

  SDL_Window *window;
  float zNear;
  float zFar;
  float fov;
  int width;
  int height;
};

#endif //DRAWER_H
