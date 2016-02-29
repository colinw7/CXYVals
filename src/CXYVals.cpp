/* TODO: Reduce all incoming rectangular boundaries from
         5->4 points to allow quick bbox test for IsInside */

#include <CXYVals.h>
#include <map>
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
init(const Polygons &polygons)
{
  init1(polygons);
}

void
CXYVals::
init(const Polygon &polygon)
{
  Polygons polygons;

  polygons.push_back(polygon);

  init1(polygons);
}

void
CXYVals::
init(const Polygon &polygon1, const Polygon &polygon2)
{
  assert(polygon1.size() == polygon2.size());

  Polygons polygons;

  polygons.push_back(polygon1);
  polygons.push_back(polygon2);

  init1(polygons);
}

void
CXYVals::
init(const std::vector<double> &x, std::vector<double> &y)
{
  init(Polygon(x, y));
}

void
CXYVals::
init(const std::vector<double> &x1, std::vector<double> &y1,
     const std::vector<double> &x2, std::vector<double> &y2)
{
  init(Polygon(x1, y1), Polygon(x2, y2));
}

void
CXYVals::
init(const double *x, const double *y, int num_xy)
{
  init(Polygon(x, y, num_xy));
}

void
CXYVals::
init(const double *x1, const double *y1, int num_xy1,
     const double *x2, const double *y2, int num_xy2)
{
  init(Polygon(x1, y1, num_xy1), Polygon(x2, y2, num_xy2));
}

void
CXYVals::
init1(const Polygons &polygons)
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

  //---

  // Allocate return value for maximum possible size
  std::vector<double> xvals, yvals;

  xvals.reserve(num_xy);
  yvals.reserve(num_xy);

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
 CXYVals(polygons), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(polygons);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const Polygon &polygon, bool init) :
 CXYVals(polygon), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(polygon);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const Polygon &polygon1, const Polygon &polygon2, bool init) :
 CXYVals(polygon1, polygon2), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(polygon1, polygon2);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const std::vector<double> &x, std::vector<double> &y, bool init) :
 CXYVals(x, y), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(x, y);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const std::vector<double> &x1, std::vector<double> &y1,
              const std::vector<double> &x2, std::vector<double> &y2, bool init) :
 CXYVals(x1, y1, x2, y2), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(x1, y1, x2, y2);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const double *x, const double *y, int num_xy, bool init) :
 CXYVals(x, y, num_xy), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(x, y, num_xy);
  else
    initMem();
}

CXYValsInside::
CXYValsInside(const double *x1, const double *y1, int num_xy1,
              const double *x2, const double *y2, int num_xy2, bool init) :
 CXYVals(x1, y1, num_xy1, x2, y2, num_xy2), num_xvals_(0), num_yvals_(0), orValues_(false)
{
  if (init)
    initValues(x1, y1, num_xy1, x2, y2, num_xy2);
  else
    initMem();
}

void
CXYValsInside::
initValues(const Polygons &polygons)
{
  initValues1(polygons);
}

void
CXYValsInside::
initValues(const Polygon &polygon)
{
  Polygons polygons;

  polygons.push_back(polygon);

  initValues1(polygons);
}

void
CXYValsInside::
initValues(const Polygon &polygon1, const Polygon &polygon2)
{
  assert(polygon1.size() == polygon2.size());

  Polygons polygons;

  polygons.push_back(polygon1);
  polygons.push_back(polygon2);

  initValues1(polygons);
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const Polygon &polygon)
{
  initValues1(xyvals, polygon);
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const double *x, const double *y, int num_xy)
{
  initValues1(xyvals, Polygon(x, y, num_xy));
}

void
CXYValsInside::
initValues(const CXYValsInside &xyvals, const std::vector<double> &x, std::vector<double> &y)
{
  initValues1(xyvals, Polygon(x, y));
}

void
CXYValsInside::
initValues(const std::vector<double> &x, std::vector<double> &y)
{
  initValues(Polygon(x, y));
}

void
CXYValsInside::
initValues(const std::vector<double> &x1, std::vector<double> &y1,
           const std::vector<double> &x2, std::vector<double> &y2)
{
  initValues(Polygon(x1, y1), Polygon(x2, y2));
}

void
CXYValsInside::
initValues(const double *x, const double *y, int num_xy)
{
  initValues(Polygon(x, y, num_xy));
}

void
CXYValsInside::
initValues(const double *x1, const double *y1, int num_xy1,
           const double *x2, const double *y2, int num_xy2)
{
  initValues(Polygon(x1, y1, num_xy1), Polygon(x2, y2, num_xy2));
}

void
CXYValsInside::
initValues1(const CXYValsInside &xyvals, const Polygon &polygon)
{
  Polygons polygons;

  (void) xyvals.getPolygons(polygons);

  polygons.push_back(polygon);

  initValues1(polygons);
}

void
CXYValsInside::
initValues1(const Polygons &polygons)
{
  CXYVals::init1(polygons);

  initMem();

  if (num_xvals_ <= 0 || num_yvals_ <= 0)
    return;

  setPolygonInsideValues(polygons);
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

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      inside_[ix].resize(num_yvals_ - 1);

      for (int iy = 0; iy < num_yvals_ - 1; ++iy)
        inside_[ix][iy] = 0;
    }
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
  // set all inside values to specified value
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
  // set all inside values to specified value if non-zero
  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      if (inside_[ix][iy])
        inside_[ix][iy] = val;
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

  // set all non-zero inside values to INSIDE1
  th->combineInside(INSIDE1);

  // extract each polygon using INSIDE1 value
  Polygon polygon;

  while (getPolygon(INSIDE1, polygon, check_consistent)) {
    polygons.push_back(polygon);

    // remove polygon
    th->setPolygonValue(polygon, 0);
  }

  //---

  // restore inside
  th->inside_ = inside;

  return ! polygons.empty();
}

bool
CXYValsInside::
getPolygons(InsideValue val, Polygons &polygons, bool check_consistent) const
{
  // save inside
  InsideArray inside = inside_;

  CXYValsInside *th = const_cast<CXYValsInside *>(this);

  //---

  // extract each polygon using specified value
  Polygon polygon;

  while (getPolygon(val, polygon, check_consistent)) {
    polygons.push_back(polygon);

    // remove polygon
    th->setPolygonValue(polygon, val ? 0 : 1);
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
  class PolyCoords {
   public:
    PolyCoords(const CXYValsInside *inside) :
     inside_(inside), num_xy_(0) {
      max_xy_ = inside_->numXVals()*inside_->numYVals();

      x_.reserve(max_xy_);
      y_.reserve(max_xy_);
    }

    bool addPoint(int i, int j) {
      Coord c(i, j);

      Coords::const_iterator p = coords_.find(c);

      if (p != coords_.end()) {
        int ii = (*p).second;

        //std::cerr << "Loop from " << ii << " to " << num_xy_ << std::endl;

        for (int i = ii; i < num_xy_; ++i) {
          x_[i - ii] = x_[i];
          y_[i - ii] = y_[i];
        }

        num_xy_ -= ii;

        return false;
      }

      x_[num_xy_] = inside_->xval(i);
      y_[num_xy_] = inside_->yval(j);

      coords_[c] = num_xy_;

      ++num_xy_;

      return true;
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
    typedef std::map<Coord,int> Coords;

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

  (void) polyCoords.addPoint(ix2, iy2);

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

        if (! polyCoords.addPoint(ix2, iy2))
          break;

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

        if (! polyCoords.addPoint(ix2, iy2))
          break;

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

        if (! polyCoords.addPoint(ix2, iy2))
          break;

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

        if (! polyCoords.addPoint(ix2, iy2))
          break;

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
setPolygonValue(const Polygon &polygon, InsideValue val)
{
  // grid cells inside polygon to value
  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (polygon.isInside(xm, ym))
        inside_[ix][iy] = val;
    }
  }
}

double
CXYValsInside::
polygonArea(const Polygon &polygon) const
{
  // get area of specified polygon
  double area = 0;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (polygon.isInside(xm, ym))
        area += CXYVals::area(ix, iy);
    }
  }

  return area;
}

double
CXYValsInside::
valueArea(InsideValue val) const
{
  // get area of cells of specified value
  double area = 0;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      if (inside_[ix][iy] == val)
        area += CXYVals::area(ix, iy);
    }
  }

  return area;
}

CXYValsInside::InsideValue
CXYValsInside::
polygonValue(const Polygon &polygon) const
{
  // get value of cells of specified value
  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    double ym = avg(yvals_[iy], yvals_[iy + 1]);

    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      double xm = avg(xvals_[ix], xvals_[ix + 1]);

      if (polygon.isInside(xm, ym))
        return inside_[ix][iy];
    }
  }

  return 0;
}

int
CXYValsInside::
getNumPolygons() const
{
  // get number of disconnected polygons
  Polygons polygons;

  if (! getPolygons(polygons))
    return 0;

  return polygons.size();
}

int
CXYValsInside::
setPolygonInsideValues()
{
  // set unique values for all disconnected polygons
  Polygons polygons;

  if (! getPolygons(polygons))
    return 0;

  setPolygonInsideValues(polygons);

  return polygons.size();
}

void
CXYValsInside::
setPolygonInsideValues(const Polygons &polygons)
{
  // set unique number for each polygon
  if (isOrValues() && polygons.size() < 8*sizeof(InsideValue)) {
    for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
      double ym = avg(yvals_[iy], yvals_[iy + 1]);

      for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
        double xm = avg(xvals_[ix], xvals_[ix + 1]);

        inside_[ix][iy] = 0;

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
  }
}

bool
CXYValsInside::
fill(bool step)
{
  if (fill1(step, FILL_OUT))
    return true;

  if (fill1(step, FILL_ALL))
    return true;

  return false;
}

bool
CXYValsInside::
fill1(bool step, FillMode mode)
{
  bool rc = true;

  while (getNumPolygons() >= 2) {
    // find polygon with smallest area
    InsideValue s = findSmallest();
    if (s == 0) break;

    // connect smallest horizontally
    if (fillHorizontal(s, mode)) {
      if (setPolygonInsideValues() == 1)
        break;

      if (step)
        break;

      continue;
    }

    // connect smallest vertically
    if (fillVertical(s, mode)) {
      if (setPolygonInsideValues() == 1)
        break;

      if (step)
        break;

      continue;
    }

    // connect polygons across empty rows
    if (fillDisconnectedRows()) {
      if (setPolygonInsideValues() == 1)
        break;

      if (step)
        break;

      continue;
    }

    // connect polygons across empty columns
    if (fillDisconnectedColumns()) {
      if (setPolygonInsideValues() == 1)
        break;

      if (step)
        break;

      continue;
    }

    rc = false;

    break;
  }

  setPolygonInsideValues();

  return rc;
}

CXYValsInside::InsideValue
CXYValsInside::
findSmallest() const
{
  typedef std::map<InsideValue,double> IArea;

  // get area of each unique cell value (non-zero)
  IArea iarea;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      InsideValue ival = inside_[ix][iy];

      if (ival == 0)
        continue;

      iarea[ival] += CXYVals::area(ix, iy);
    }
  }

  // find cell value with miniumum area
  InsideValue minValue = 0;
  double      minArea  = 0;

  for (IArea::const_iterator p = iarea.begin(); p != iarea.end(); ++p) {
    double area = (*p).second;

    if (minValue == 0 || area < minArea) {
      minValue = (*p).first;
      minArea  = area;
    }
  }

  return minValue;
}

bool
CXYValsInside::
fillHorizontal(InsideValue val, FillMode mode)
{
  int x1, y1, x2, y2;

  polygonBounds(val, x1, y1, x2, y2);

  //---

  bool found = false;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    int ix = 0;

    // find first non-zero
    while (ix < num_xvals_ - 1 && ! inside_[ix][iy])
      ++ix;

    if (ix >= num_xvals_ - 1)
      continue;

    while (ix < num_xvals_ - 1) {
      InsideValue val1 = inside_[ix][iy];

      // skip to last non-zero
      while (ix < num_xvals_ - 1 && inside_[ix][iy] == val1)
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
      InsideValue val2 = inside_[ix][iy];

      if (val1 == val2)
        continue;

      // if specified value then one needs to match
      if (val1 != val && val2 != val)
        continue;

      // must fill from bounding edge to bounding edge
      if (mode == FILL_OUT) {
        if (val1 == val) {
          if (ix1 != x2)
            continue;

          int x3, y3, x4, y4;

          polygonBounds(val2, x3, y3, x4, y4);

          if (ix2 + 1 != x3)
            continue;
        }
        else {
          if (ix2 + 1 != x1)
            continue;

          int x3, y3, x4, y4;

          polygonBounds(val1, x3, y3, x4, y4);

          if (ix1 != x4)
            continue;
        }
      }

      // fill zero values
      for (int ii = ix1; ii <= ix2; ++ii)
        inside_[ii][iy] = val;

      found = true;
    }
  }

  return found;
}

bool
CXYValsInside::
fillVertical(InsideValue val, FillMode mode)
{
  int x1, y1, x2, y2;

  polygonBounds(val, x1, y1, x2, y2);

  //---

  bool found = false;

  for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
    int iy = 0;

    // find first non-zero
    while (iy < num_yvals_ - 1 && ! inside_[ix][iy])
      ++iy;

    if (iy >= num_yvals_ - 1)
      continue;

    while (iy < num_yvals_ - 1) {
      InsideValue val1 = inside_[ix][iy];

      // skip to last non-zero
      while (iy < num_yvals_ - 1 && inside_[ix][iy] == val1)
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
      InsideValue val2 = inside_[ix][iy];

      if (val1 == val2)
        continue;

      // if specified value then one needs to match
      if (val1 != val && val2 != val)
        continue;

      // must fill from bounding edge to bounding edge
      if (mode == FILL_OUT) {
        if (val1 == val) {
          if (iy1 != y2)
            continue;

          int x3, y3, x4, y4;

          polygonBounds(val2, x3, y3, x4, y4);

          if (iy2 + 1 != y3)
            continue;
        }
        else {
          if (iy2 + 1 != y1)
            continue;

          int x3, y3, x4, y4;

          polygonBounds(val1, x3, y3, x4, y4);

          if (iy1 != y4)
            continue;
        }
      }

      // fill zero values
      for (int jj = iy1; jj <= iy2; ++jj)
        inside_[ix][jj] = val;

      found  = true;
    }
  }

  return found;
}

void
CXYValsInside::
polygonBounds(InsideValue val, int &x1, int &y1, int &x2, int &y2) const
{
  // find indices of edges of polygon of specified value
  x1 = -1; y1 = -1;
  x2 = -1; y2 = -1;

  for (int iy = 0; iy < num_yvals_ - 1; ++iy) {
    for (int ix = 0; ix < num_xvals_ - 1; ++ix) {
      if (inside_[ix][iy] == val) {
        if (x1 == -1 || ix     < x1) x1 = ix;
        if (y1 == -1 || iy     < y1) y1 = iy;
        if (x2 == -1 || ix + 1 > x2) x2 = ix + 1;
        if (y2 == -1 || iy + 1 > y2) y2 = iy + 1;
      }
    }
  }
}

bool
CXYValsInside::
fillDisconnectedRows()
{
  // find empty row (must be one)
  int iy = 0;

  for ( ; iy < num_yvals_ - 1; ++iy)
    if (isEmptyRow(iy))
      break;

  // must be found and not first or last
  if (iy <= 0 || iy >= num_yvals_ - 2)
    return false;

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

  if (ix1 < ix2) {
    for (int ix = ix1; ix <= ix2; ++ix)
      inside_[ix][iy] = 1;
  }
  else {
    for (int ix = ix2; ix <= ix1; ++ix)
      inside_[ix][iy] = 1;
  }

  return true;
}

bool
CXYValsInside::
fillDisconnectedColumns()
{
  // find empty column (must be one)
  int ix = 0;

  for ( ; ix < num_xvals_ - 1; ++ix)
    if (isEmptyColumn(ix))
      break;

  // must be found and not first or last
  if (ix <= 0 || ix >= num_xvals_ - 2)
    return false;

  // find end of first block on previous column
  int iy1 = 0;

  for ( ; iy1 < num_yvals_ - 1; ++iy1)
    if (inside_[ix - 1][iy1])
      break;

  while (iy1 < num_yvals_ - 2 && inside_[ix - 1][iy1 + 1])
    ++iy1;

  // find start of first block on next column
  int iy2 = 0;

  for ( ; iy2 < num_yvals_ - 1; ++iy2)
    if (inside_[ix + 1][iy2])
      break;

  if (iy1 < iy2) {
    for (int iy = iy1; iy <= iy2; ++iy)
      inside_[ix][iy] = 1;
  }
  else {
    for (int iy = iy2; iy <= iy1; ++iy)
      inside_[ix][iy] = 1;
  }

  return true;
}

bool
CXYValsInside::
isEmptyRow(int iy) const
{
  // check if row is empty
  int num_x = getNumXVals();

  for (int ix = 0; ix < num_x - 1; ++ix)
    if (inside_[ix][iy])
      return false;

  return true;
}

bool
CXYValsInside::
isEmptyColumn(int ix) const
{
  // check if column is empty
  int num_y = getNumYVals();

  for (int iy = 0; iy < num_y - 1; ++iy)
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
CXYValsInside::
unitTest(const std::string &)
{
  static double x[5] = { -0.1, 24.7, 50.3, 75.2, 100.5 };
  static double y[5] = {  0.1, 45.8, 45.8, 77.3,  99.9 };

  CXYValsInside xyvals(x, y, 5);

  std::cerr << "Init" << std::endl;
  xyvals.print(std::cerr);
  std::cerr << std::endl;

  //---

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

  std::cerr << "Set Values" << std::endl;
  xyvals.print(std::cerr);
  std::cerr << std::endl;

  //---

  CXYValsInside::Polygons polygons;

  bool rc = xyvals.getPolygons(polygons, false);
  assert(rc);

  assert(polygons.size() == 6);

  std::cerr << polygons.size() << " Polygons" << std::endl;
  for (uint i = 0; i < polygons.size(); ++i) {
    std::cerr << "Polygon: " << i << std::endl;

    const CXYValsInside::Polygon &polygon = polygons[i];

    for (int j = 0; j < polygon.size(); ++j)
      std::cerr << " {" << polygon.x[j] << "," << polygon.y[j] << "}" << std::endl;
  }
  std::cerr << std::endl;

  //---

  return true;
}
