/* TODO: Reduce all incoming rectangular boundaries from
         5->4 points to allow quick bbox test for IsInside */

#include <CXYVals.h>
#include <set>

namespace {

template<typename T>
T abs(T v) {
  return (v > T(0) ? v : -v);
}

template<typename T>
T avg(T v1, T v2) {
  return (v1 + v2)/T(2);
}

template<typename T>
bool realEq(T r1, T r2, T tol=1E-3) {
  if (r1 == r2 || abs(r1 - r2) < tol) return true;

  return (abs((r1 - r2)/(abs(r2) > abs(r1) ? r2 : r1)) <= tol);
}

bool PointInsideEvenOdd(double x, double y, const double *xp, const double *yp, int np) {
  double xinters, x1, y1, x2, y2;

  int counter = 0;

  int i2 = np - 1;

  x2 = xp[i2]; y2 = yp[i2];

  // iterate through all lines of the polygon
  for (int i1 = 0; i1 < (int) np; i2 = i1++, x2 = x1, y2 = y1) {
    x1 = xp[i1]; y1 = yp[i1];

    // intersect current line with horizontal line at inside point
    if (y > std::min(y1, y2)) {
      if (y <= std::max(y1, y2)) {
        if (x <= std::max(x1, x2)) {
          if (y1 != y2) {
            // if we have an intersection, increase count
            xinters = (y - y1)*(x2 - x1)/(y2 - y1) + x1;

            if (x1 == x2 || x <= xinters)
              ++counter;
          }
        }
      }
    }
  }

  // if odd then success
  return ((counter % 2) != 0);
}

}

//------

/*
 * Get unique X and Y Values for rectilinear polygon
 */

CXYVals::
CXYVals(const Polygons &polygons)
{
  init(polygons);
}

CXYVals::
CXYVals(const Polygon &polygon)
{
  init(polygon);
}

CXYVals::
CXYVals(const Polygon &polygon1, const Polygon &polygon2)
{
  init(polygon1, polygon2);
}

CXYVals::
CXYVals(const std::vector<double> &x, std::vector<double> &y)
{
  init(x, y);
}

CXYVals::
CXYVals(const std::vector<double> &x1, std::vector<double> &y1,
        const std::vector<double> &x2, std::vector<double> &y2)
{
  init(x1, y1, x2, y2);
}

CXYVals::
CXYVals(const double *x, const double *y, int num_xy)
{
  init(x, y, num_xy);
}

CXYVals::
CXYVals(const double *x1, const double *y1, int num_xy1,
        const double *x2, const double *y2, int num_xy2)
{
  init(x1, y1, num_xy1, x2, y2, num_xy2);
}

void
CXYVals::
init(const Polygon &polygon)
{
  init(&polygon.x[0], &polygon.y[0], polygon.x.size());
}

void
CXYVals::
init(const Polygon &polygon1, const Polygon &polygon2)
{
  assert(polygon1.size() == polygon2.size());

  init(&polygon1.x[0], &polygon1.y[0], polygon1.size(),
       &polygon2.x[0], &polygon2.y[0], polygon2.size());
}

void
CXYVals::
init(const std::vector<double> &x, std::vector<double> &y)
{
  assert(x.size() == y.size());

  init(&x[0], &y[0], x.size());
}

void
CXYVals::
init(const std::vector<double> &x1, std::vector<double> &y1,
     const std::vector<double> &x2, std::vector<double> &y2)
{
  assert(x1.size() == y1.size());
  assert(x2.size() == y2.size());

  init(&x1[0], &y1[0], x1.size(), &x2[0], &y2[0], x2.size());
}

void
CXYVals::
init(const double *x, const double *y, int num_xy)
{
  // Allocate return value for maximum possible size
  std::vector<double> xvals, yvals;

  xvals.reserve(num_xy);
  yvals.reserve(num_xy);

  //---

  int num_xvals = 0;
  int num_yvals = 0;

  for (int i = 0; i < num_xy; ++i) {
    insertValue(x[i], &xvals[0], num_xvals);
    insertValue(y[i], &yvals[0], num_yvals);
  }

  if (num_xy) {
    xvals_ = std::vector<double>(&xvals[0], &xvals[num_xvals]);
    yvals_ = std::vector<double>(&yvals[0], &yvals[num_yvals]);
  }
  else {
    xvals_.clear();
    yvals_.clear();
  }
}

void
CXYVals::
init(const Polygons &polygons)
{
  int num_xy = 0;

  for (uint i = 0; i < polygons.size(); ++i)
    num_xy += polygons[i].size();

  std::vector<double> x, y;

  x.resize(num_xy);
  y.resize(num_xy);

  num_xy = 0;

  for (uint i = 0; i < polygons.size(); ++i) {
    const Polygon &polygon = polygons[i];

    int n = polygon.size();

    memcpy(&x[num_xy], &polygon.x[0], n*sizeof(double));
    memcpy(&y[num_xy], &polygon.y[0], n*sizeof(double));

    num_xy += n;
  }

  init(&x[0], &y[0], num_xy);
}

void
CXYVals::
init(const double *x1, const double *y1, int num_xy1,
     const double *x2, const double *y2, int num_xy2)
{
  int num_xy = num_xy1 + num_xy2;

  std::vector<double> x, y;

  x.resize(num_xy);
  y.resize(num_xy);

  if (num_xy1) {
    memcpy(&x[0      ], x1, num_xy1*sizeof(double));
    memcpy(&y[0      ], y1, num_xy1*sizeof(double));
  }

  if (num_xy2) {
    memcpy(&x[num_xy1], x2, num_xy2*sizeof(double));
    memcpy(&y[num_xy1], y2, num_xy2*sizeof(double));
  }

  init(&x[0], &y[0], num_xy);
}

void
CXYVals::
insertValue(double v, double *vals, int &num_vals)
{
  // handle empty list
  if (num_vals == 0) {
    vals[0] = v;

    ++num_vals;

    return;
  }

  // Find insertion position
  int l = 0;
  int h = num_vals - 1;

  if (realEq(v, vals[l]) || realEq(v, vals[h]))
    return;

  int pos;

  if      (v <= vals[l]) {
    pos = 0;
  }
  else if (v >= vals[h]) {
    pos = num_vals;
  }
  else {
    int m = (l + h)/2;

    while (l != m) {
      if (realEq(v, vals[m]))
        return;

      if (v < vals[m]) h = m;
      else             l = m;

      m = (l + h)/2;
    }

    pos = h;
  }

  // make room
  for (int i = num_vals - 1; i >= pos; --i)
    vals[i + 1] = vals[i];

  // Add to array
  vals[pos] = v;

  ++num_vals;
}

void
CXYVals::
clear()
{
  xvals_.clear();
  yvals_.clear();
}

//------

CXYValsInside::
CXYValsInside(const Polygons &polygons, bool init) :
 CXYVals(polygons), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(polygons);
}

CXYValsInside::
CXYValsInside(const Polygon &polygon, bool init) :
 CXYVals(polygon), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(polygon);
}

CXYValsInside::
CXYValsInside(const Polygon &polygon1, const Polygon &polygon2, bool init) :
 CXYVals(polygon1, polygon2), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(polygon1, polygon2);
}

CXYValsInside::
CXYValsInside(const std::vector<double> &x, std::vector<double> &y, bool init) :
 CXYVals(x, y), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(x, y);
}

CXYValsInside::
CXYValsInside(const std::vector<double> &x1, std::vector<double> &y1,
              const std::vector<double> &x2, std::vector<double> &y2, bool init) :
 CXYVals(x1, y1, x2, y2), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(x1, y1, x2, y2);
}

CXYValsInside::
CXYValsInside(const double *x, const double *y, int num_xy, bool init) :
 CXYVals(x, y, num_xy), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(x, y, num_xy);
}

CXYValsInside::
CXYValsInside(const double *x1, const double *y1, int num_xy1,
              const double *x2, const double *y2, int num_xy2, bool init) :
 CXYVals(x1, y1, num_xy1, x2, y2, num_xy2), num_xvals_(0), num_yvals_(0)
{
  if (init)
    initValues(x1, y1, num_xy1, x2, y2, num_xy2);
}

void
CXYValsInside::
initValues(const Polygons &polygons)
{
  CXYVals::init(polygons);

  initMem();

  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  if (polygons.size() < 8*sizeof(InsideValue)) {
    for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
      double ym = avg(yvals_[iy], yvals_[iy + 1]);

      for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
        double xm = avg(xvals_[ix], xvals_[ix + 1]);

        for (uint k = 0; k < polygons.size(); ++k) {
          const Polygon &polygon = polygons[k];

          InsideValue ival = (1U << k);

          if (polygon.isInside(xm, ym))
            inside_[ix][iy] |= ival;
        }
      }
    }
  }
  else {
    for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
      double ym = avg(yvals_[iy], yvals_[iy + 1]);

      for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
        double xm = avg(xvals_[ix], xvals_[ix + 1]);

        for (uint k = 0; k < polygons.size(); ++k) {
          const Polygon &polygon = polygons[k];

          if (polygon.isInside(xm, ym))
            inside_[ix][iy] = INSIDE1;
        }
      }
    }
  }
}

void
CXYValsInside::
initValues(const Polygon &polygon)
{
  initValues(&polygon.x[0], &polygon.y[0], polygon.size());
}

void
CXYValsInside::
initValues(const Polygon &polygon1, const Polygon &polygon2)
{
  assert(polygon1.size() == polygon2.size());

  initValues(&polygon1.x[0], &polygon1.y[0], polygon1.size(),
             &polygon2.x[0], &polygon2.y[0], polygon2.size());
}

void
CXYValsInside::
initValues(const std::vector<double> &x, std::vector<double> &y)
{
  assert(x.size() == y.size());

  initValues(&x[0], &y[0], x.size());
}

void
CXYValsInside::
initValues(const std::vector<double> &x1, std::vector<double> &y1,
           const std::vector<double> &x2, std::vector<double> &y2)
{
  assert(x1.size() == y1.size());
  assert(x2.size() == y2.size());

  initValues(&x1[0], &y1[0], x1.size(), &x2[0], &y2[0], x2.size());
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const Polygon &polygon)
{
  initValues(xyvals, &polygon.x[0], &polygon.y[0], polygon.size());
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const double *x, const double *y, int num_xy)
{
  std::vector<double> x1(&x[0], &x[num_xy - 1]);
  std::vector<double> y1(&y[0], &y[num_xy - 1]);

  initValues(xyvals, x1, y1);
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const std::vector<double> &x, std::vector<double> &y)
{
  Polygons polygons;

  if (xyvals.getPolygons(polygons)) {
    polygons.push_back(Polygon(x, y));

    initValues(polygons);
  }
  else
    initValues(x, y);
}

void
CXYValsInside::
initValues(const double *x, const double *y, int num_xy)
{
  CXYVals::init(x, y, num_xy);

  initMem();

  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (PointInsideEvenOdd(xm, ym, x, y, num_xy))
        inside_[ix][iy] |= INSIDE1;
    }
  }
}

void
CXYValsInside::
initValues(const double *x1, const double *y1, int num_xy1,
           const double *x2, const double *y2, int num_xy2)
{
  CXYVals::init(x1, y1, num_xy1, x2, y2, num_xy2);

  initMem();

  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (PointInsideEvenOdd(xm, ym, x1, y1, num_xy1))
        inside_[ix][iy] |= INSIDE1;

      if (PointInsideEvenOdd(xm, ym, x2, y2, num_xy2))
        inside_[ix][iy] |= INSIDE2;
    }
  }
}

void
CXYValsInside::
initMem()
{
  num_xvals_ = xvals_.size();
  num_yvals_ = yvals_.size();

  // create inside flag grid
  if (num_xvals_ > 0 && num_yvals_ > 0) {
    inside_.resize(num_xvals_ - 1);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix)
      inside_[ix].resize(num_yvals_ - 1);
  }
  else
    inside_.clear();
}

#define POLY_IS_INSIDE(ix,iy,v) \
  ((ix) >= 0 && (ix) < num_xvals_ - 1 && \
   (iy) >= 0 && (iy) < num_yvals_ - 1 && \
   (inside_[ix][iy] == (v)))

class CXYEdge {
 private:
  class Point {
   private:
    int i_, j_;

   public:
    Point(int i, int j) :
     i_(i), j_(j) {
    }

    friend bool operator==(const Point &p1, const Point &p2) {
      return (p1.i_ == p2.i_ && p1.j_ == p2.j_);
    }

    friend bool operator<(const Point &p1, const Point &p2) {
      return (p1.i_ < p2.i_ || (p1.i_ == p2.i_ && p1.j_ < p2.j_));
    }
  };

  Point p1_;
  Point p2_;

 public:
  CXYEdge(int i1, int j1, int i2, int j2) :
   p1_(i1, j1), p2_(i2, j2) {
  }

  friend bool operator==(const CXYEdge &e1, const CXYEdge &e2) {
    return (e1.p1_ == e2.p1_ && e1.p2_ == e2.p2_);
  }

  friend bool operator<(const CXYEdge &e1, const CXYEdge &e2) {
    return (e1.p1_ < e2.p1_ || (e1.p1_ == e2.p1_ && e1.p2_ < e2.p2_));
  }
};

class CXYEdgeList {
 private:
  typedef std::set<CXYEdge> EdgeSet;

  EdgeSet edgeSet_;

 public:
  CXYEdgeList() { }

  void addEdge(int i1, int j1, int i2, int j2) {
    edgeSet_.insert(CXYEdge(i1, j1, i2, j2));
  }

  bool hasEdge(int i1, int j1, int i2, int j2) {
    return (edgeSet_.find(CXYEdge(i1, j1, i2, j2)) != edgeSet_.end());
  }
};

void
CXYValsInside::
setInside(InsideValue val)
{
  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      inside_[ix][iy] = val;
    }
  }
}

void
CXYValsInside::
combineInside(InsideValue val)
{
  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      if (inside_[ix][iy])
        inside_[ix][iy] = val;
      else
        inside_[ix][iy] = 0;
    }
  }
}

void
CXYValsInside::
clear()
{
  CXYVals::clear();

  inside_.clear();

  num_xvals_ = 0;
  num_yvals_ = 0;
}

bool
CXYValsInside::
getPolygons(Polygons &polygons, bool check_consistent) const
{
  // save inside
  InsideArray inside = inside_;

  CXYValsInside *th = const_cast<CXYValsInside *>(this);

  //---

  Polygon polygon;

  while (getPolygon(polygon, check_consistent)) {
    polygons.push_back(polygon);

    th->removePolygon(polygon);
  }

  //---

  // restore inside
  th->inside_ = inside;

  return ! polygons.empty();
}

bool
CXYValsInside::
getPolygon(Polygon &polygon, bool check_consistent) const
{
  return getPolygon(INSIDE1, polygon, check_consistent);
}

bool
CXYValsInside::
getPolygon(InsideValue inside_val, Polygon &polygon, bool check_consistent) const
{
  return getPolygon(inside_val, polygon.x, polygon.y, check_consistent);
}

bool
CXYValsInside::
getPolygon(std::vector<double> &x, std::vector<double> &y, bool check_consistent) const
{
  return getPolygon(INSIDE1, x, y, check_consistent);
}

bool
CXYValsInside::
getPolygon(double **xo, double **yo, int *num_xyo, bool check_consistent) const
{
  return getPolygon(INSIDE1, xo, yo, num_xyo, check_consistent);
}

bool
CXYValsInside::
getPolygon(InsideValue inside_val, double **xo, double **yo, int *num_xyo,
           bool check_consistent) const
{
  static std::vector<double> rect_temp_x, rect_temp_y;

  bool rc = getPolygon(inside_val, rect_temp_x, rect_temp_y, check_consistent);

  if (rc) {
    if (xo     ) *xo      = &rect_temp_x[0];
    if (yo     ) *yo      = &rect_temp_y[0];
    if (num_xyo) *num_xyo = rect_temp_x.size();
  }
  else {
    if (xo     ) *xo      = 0;
    if (yo     ) *yo      = 0;
    if (num_xyo) *num_xyo = 0;
  }

  return rc;
}

bool
CXYValsInside::
getPolygon(InsideValue inside_val, std::vector<double> &x, std::vector<double> &y,
           bool check_consistent) const
{
  typedef std::pair<int,int> Coord;
  typedef std::vector<Coord> Coords;

  class PolyCoords {
   public:
    PolyCoords(const CXYValsInside *inside) :
     inside_(inside), num_xy_(0) {
      max_xy_ = inside_->numXVals()*inside_->numYVals();

      x_.reserve(max_xy_);
      y_.reserve(max_xy_);
    }

    void addPoint(int i, int j) {
      x_[num_xy_] = inside_->xval(i);
      y_[num_xy_] = inside_->yval(j);

      ++num_xy_;

      coords_.push_back(Coord(i, j));
    }

    int numXY() const { return num_xy_; }

    bool isFull() const { return (num_xy_ >= max_xy_); }

    int isInside(double xm, double ym) const {
      return PointInsideEvenOdd(xm, ym, &x_[0], &y_[0], num_xy_);
    }

    void getXY(std::vector<double> &x, std::vector<double> &y) {
      x = std::vector<double>(&x_[0], &x_[num_xy_]);
      y = std::vector<double>(&y_[0], &y_[num_xy_]);
    }

   private:
    const CXYValsInside *inside_;
    std::vector<double>  x_, y_;
    int                  max_xy_;
    int                  num_xy_;
    Coords               coords_;
  };

  //---

  PolyCoords polyCoords(this);

  CXYEdgeList edgeList;

  //---

  // Find Bottom Left of inside coords
  bool found = false;

  int ix, iy;

  for (iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (ix = 0; ix < num_xvals_ - 1; ++ix) {
      if (inside_[ix][iy] == inside_val) {
        found = true;
        break;
      }
    }

    if (found)
      break;
  }

  if (! found)
    return false;

  int ix1 = ix;
  int iy1 = iy;

  //---

  // Walk alternate X and Y directions to trace boundary
  int ix2 = ix1;
  int iy2 = iy1;

  polyCoords.addPoint(ix2, iy2);

  while (! polyCoords.isFull()) {
    int save_num_xy = polyCoords.numXY();

    //---

    bool r_ok = (ix2 < num_xvals_ - 1 && ! edgeList.hasEdge(ix2, iy2, ix2 + 1, iy2));
    bool l_ok = (ix2 >              0 && ! edgeList.hasEdge(ix2 - 1, iy2, ix2, iy2));

    if (r_ok) {
      // Walk left to right
      int ix = ix2;

      for ( ; ix < num_xvals_ - 1; ++ix)
        if (  POLY_IS_INSIDE(ix, iy2 - 1, inside_val) ||
            ! POLY_IS_INSIDE(ix, iy2    , inside_val))
          break;

      if (ix > ix2) {
        for ( ; ix2 < ix; ++ix2)
          edgeList.addEdge(ix2, iy2, ix2 + 1, iy2);

        if (ix1 == ix2 && iy1 == iy2)
          break;

        polyCoords.addPoint(ix2, iy2);

        if (polyCoords.isFull())
          break;

        l_ok = false;
      }
    }

    if (l_ok) {
      // Walk right to left
      int ix = ix2;

      for ( ; ix >= 0; --ix)
        if (! POLY_IS_INSIDE(ix - 1, iy2 - 1, inside_val) ||
              POLY_IS_INSIDE(ix - 1, iy2    , inside_val))
          break;

      if (ix < ix2) {
        for ( ; ix2 > ix; --ix2)
          edgeList.addEdge(ix2 - 1, iy2, ix2, iy2);

        if (ix1 == ix2 && iy1 == iy2)
          break;

        polyCoords.addPoint(ix2, iy2);

        if (polyCoords.isFull())
          break;
      }
    }

    //---

    bool u_ok = (iy2 < num_yvals_ - 1 && ! edgeList.hasEdge(ix2, iy2, ix2, iy2 + 1));
    bool d_ok = (iy2 >              0 && ! edgeList.hasEdge(ix2, iy2 - 1, ix2, iy2));

    if (u_ok) {
      // Walk bottom to top
      int iy = iy2;

      for ( ; iy < num_yvals_ - 1; ++iy)
        if (! POLY_IS_INSIDE(ix2 - 1, iy, inside_val) ||
              POLY_IS_INSIDE(ix2    , iy, inside_val))
          break;

      if (iy > iy2) {
        for ( ; iy2 < iy; ++iy2)
          edgeList.addEdge(ix2, iy2, ix2, iy2 + 1);

        if (ix1 == ix2 && iy1 == iy2)
          break;

        polyCoords.addPoint(ix2, iy2);

        if (polyCoords.isFull())
          break;

        d_ok = false;
      }
    }

    if (d_ok) {
      // Walk top to bottom
      int iy = iy2;

      for ( ; iy >= 0; --iy)
        if (  POLY_IS_INSIDE(ix2 - 1, iy - 1, inside_val) ||
            ! POLY_IS_INSIDE(ix2    , iy - 1, inside_val))
          break;

      if (iy < iy2) {
        for ( ; iy2 > iy; --iy2)
          edgeList.addEdge(ix2, iy2 - 1, ix2, iy2);

        if (ix1 == ix2 && iy1 == iy2)
          break;

        polyCoords.addPoint(ix2, iy2);

        if (polyCoords.isFull())
          break;
      }
    }

    //---

    // Avoid endless loop
    if (polyCoords.numXY() == save_num_xy)
      break;
  }

  //---

  // Ensure new boundary is consistent with the inside flags
  if (check_consistent) {
    for (iy = 0; iy < num_yvals_ - 1; ++iy) {
      double ym = avg(yvals_[iy], yvals_[iy + 1]);

      for (ix = 0; ix < num_xvals_ - 1; ++ix) {
        double xm = avg(xvals_[ix], xvals_[ix + 1]);

        bool inside1 = (inside_[ix][iy] == inside_val);
        bool inside2 = polyCoords.isInside(xm, ym);

        if (inside1 != inside2)
          return false;
      }
    }
  }

  //---

  polyCoords.getXY(x, y);

  return true;
}

void
CXYValsInside::
removePolygon(const Polygon &polygon)
{
  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (polygon.isInside(xm, ym))
        inside_[ix][iy] = 0;
    }
  }
}

void
CXYValsInside::
fill(bool disconnected)
{
  // if less than 2 polygons we are done
  Polygons polygons;

  if (! getPolygons(polygons) || polygons.size() < 2)
    return;

print(std::cerr); std::cerr << std::endl;
  // set unique number for each polygon
  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      inside_[ix][iy] = 0;

      for (uint k = 0; k < polygons.size(); ++k) {
        const Polygon &polygon = polygons[k];

        if (polygon.isInside(xm, ym)) {
          inside_[ix][iy] = k + 1;
          break;
        }
      }
    }
  }
print(std::cerr); std::cerr << std::endl;

  // connected different regions horizontally and vertically
  int numFilled = 0;

  while (fillH() || fillV()) {
print(std::cerr); std::cerr << std::endl;
    ++numFilled;
  }

  // reset inside to single value
  combineInside(1);
print(std::cerr); std::cerr << std::endl;

  // if still disconnected then force connection on empty rows and refill
  if (numFilled == 0 && disconnected) {
    fillDisconnected();

    return fill(disconnected);
  }
}

bool
CXYValsInside::
fillH()
{
  bool found = false;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    int ix = 0;

    // find first non-zero
    while (ix < num_xvals_ - 1 && ! inside_[ix][iy])
      ++ix;

    if (ix >= num_xvals_ - 1)
      continue;

    while (ix < num_xvals_ - 1) {
      InsideValue val = inside_[ix][iy];

      // skip to last non-zero
      while (ix < num_xvals_ - 1 && inside_[ix][iy] == val)
        ++ix;

      if (ix >= num_xvals_ - 1 || inside_[ix][iy])
        continue;

      int ix1 = ix;

      // find next non-zero
      while (ix < num_xvals_ - 1 && ! inside_[ix][iy])
        ++ix;

      if (ix >= num_xvals_ - 1)
        continue;

      int ix2 = ix - 1;

      // must be different value
      if (inside_[ix][iy] == val)
        continue;

      for (int ii = ix1; ii <= ix2; ++ii)
        inside_[ii][iy] = val;

      found  = true;
    }
  }

  return found;
}

bool
CXYValsInside::
fillV()
{
  bool found = false;

  for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
    int iy = 0;

    // find first non-zero
    while (iy < num_yvals_ - 1 && ! inside_[ix][iy])
      ++iy;

    if (iy >= num_yvals_ - 1)
      continue;

    while (iy < num_yvals_ - 1) {
      InsideValue val = inside_[ix][iy];

      // skip to last non-zero
      while (iy < num_yvals_ - 1 && inside_[ix][iy] == val)
        ++iy;

      if (iy >= num_yvals_ - 1 || inside_[ix][iy])
        continue;

      int iy1 = iy;

      // find next non-zero
      while (iy < num_yvals_ - 1 && ! inside_[ix][iy])
        ++iy;

      if (iy >= num_yvals_ - 1)
        continue;

      int iy2 = iy - 1;

      // must be different value
      if (inside_[ix][iy] == val)
        continue;

      for (int jj = iy1; jj <= iy2; ++jj)
        inside_[ix][jj] = val;

      found  = true;
    }
  }

  return found;
}

void
CXYValsInside::
fillDisconnected()
{
  // find empty row (must be one)
  int iy = 0;

  for ( ; iy < num_yvals_ - 1; ++iy)
    if (isEmptyRow(iy))
      break;

  // must be found and not first or last
  assert(iy > 0 && iy < num_yvals_ - 2);

  // find end of first block on previous row
  int ix1 = 0;

  for ( ; ix1 < num_xvals_ - 1; ++ix1)
    if (inside_[ix1][iy - 1])
      break;

  while (ix1 < num_xvals_ - 2 && inside_[ix1 + 1][iy - 1])
    ++ix1;

  // find start of first block on next row
  int ix2 = 0;

  for ( ; ix2 < num_xvals_ - 1; ++ix2)
    if (inside_[ix2][iy + 1])
      break;

  assert(ix1 != ix2);

  if (ix1 < ix2) {
    for (int ix = ix1; ix <= ix2; ++ix)
      inside_[ix][iy] = 1;
  }
  else {
    for (int ix = ix2; ix <= ix1; ++ix)
      inside_[ix][iy] = 1;
  }

  print(std::cerr); std::cerr << std::endl;
}

bool
CXYValsInside::
isEmptyRow(int iy) const
{
  int num_x = getNumXVals();

  for (int ix = 0; ix < num_x - 1; ++ix)
    if (inside_[ix][iy])
      return false;

  return true;
}

void
CXYValsInside::
print(std::ostream &os) const
{
  int num_x = getNumXVals();
  int num_y = getNumYVals();

  for (int iy = 0; iy < num_y - 1; ++iy) {
    for (int ix = 0; ix < num_x - 1; ++ix)
      os << inside_[ix][iy];

    os << std::endl;
  }
}

bool
CXYValsInside::Polygon::
isInside(double xm, double ym) const
{
  return PointInsideEvenOdd(xm, ym, &x[0], &y[0], x.size());
}

//---------

bool
CXYValsSelfTest(const char *)
{
  static double x[5] = { -0.1, 24.7, 50.3, 75.2, 100.5 };
  static double y[5] = {  0.1, 45.8, 45.8, 77.3,  99.9 };

  CXYValsInside xyvals(x, y, 5);

  CXYValsInside::InsideArray &inside = xyvals.getInside();

  int num_x = xyvals.getNumXVals();
  int num_y = xyvals.getNumYVals();

  assert(num_x == 5 /* "num_x" */);
  assert(num_y == 4 /* "num_y" */);

  for (int iy = 0; iy < num_y - 1; ++iy) {
    for (int ix = 0; ix < num_x - 1; ++ix) {
      inside[ix][iy] = (ix + iy) & 1;
    }
  }

  double *xo, *yo;
  int     num_xyo;

  bool rc = xyvals.getPolygon(CXYValsInside::INSIDE1, &xo, &yo, &num_xyo, false);

  assert(rc /* "getPolygon" */);

  assert(num_xyo == 20 /* "num_xyo" */);

  return true;
}
