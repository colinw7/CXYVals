#ifndef CXY_VALS_H
#define CXY_VALS_H

#include <vector>
#include <iostream>
#include <cstring>
#include <cassert>
#include <sys/types.h>

/*
 * Store x, y values in ordered grid of unique values
 */
class CXYVals {
 public:
  enum class FillMode {
    FILL_OUT,
    FILL_ALL
  };

  struct Polygon {
    Polygon() { }

    Polygon(const std::vector<double> &x1, std::vector<double> &y1) :
     x(x1), y(y1) {
      assert(x.size() == y.size());
    }

    Polygon(const double *x1, const double *y1, int n1) {
      x.resize(n1);
      y.resize(n1);

      memcpy(&x[0], x1, n1*sizeof(double));
      memcpy(&y[0], y1, n1*sizeof(double));
    }

    int size() const { return int(x.size()); }

    bool isInside(double xm, double ym) const;

    void centroid(double &xc, double &yc) const {
      xc = 0; yc = 0;

      for (int i = 0; i < size(); ++i) {
        xc += x[i];
        yc += y[i];
      }

      xc /= size();
      yc /= size();
    }

    void print(std::ostream &os) const {
      for (int i = 0; i < size(); ++i)
        os << " (" << x[i] << "," << y[i] << ")";
      os << std::endl;
    }

    std::vector<double> x;
    std::vector<double> y;
  };

  using Polygons = std::vector<Polygon>;

  using Coord = std::pair<int,int>;

 public:
  CXYVals(const Polygons &polygons);

  CXYVals(const Polygon &polygon);
  CXYVals(const Polygon &polygon1, const Polygon &polygon2);

  CXYVals(const std::vector<double> &x, std::vector<double> &y);
  CXYVals(const std::vector<double> &x1, std::vector<double> &y1,
          const std::vector<double> &x2, std::vector<double> &y2);

  CXYVals(const double *x=0, const double *y=0, int num_xy=0);
  CXYVals(const double *x1, const double *y1, int num_xy1,
          const double *x2, const double *y2, int num_xy2);

 ~CXYVals() { }

  const std::vector<double> &xvals() const { return xvals_; }
  const std::vector<double> &yvals() const { return yvals_; }

  double xval(int ix) const { return xvals_[ix]; }
  double yval(int iy) const { return yvals_[iy]; }

  double xsize(int ix) const { return xvals_[ix + 1] - xvals_[ix]; }
  double ysize(int iy) const { return yvals_[iy + 1] - yvals_[iy]; }

  double area(int ix, int iy) const { return xsize(ix)*ysize(iy); }

  const double *getXVals   () const { return &xvals_[0]; }
  int           getNumXVals() const { return int(xvals_.size()); }

  const double *getYVals   () const { return &yvals_[0]; }
  int           getNumYVals() const { return int(yvals_.size()); }

  void init(const Polygons &polygons);

  void init(const Polygon &polygon);
  void init(const Polygon &polygon1, const Polygon &polygon2);

  void init(const std::vector<double> &x, std::vector<double> &y);
  void init(const std::vector<double> &x1, std::vector<double> &y1,
            const std::vector<double> &x2, std::vector<double> &y2);

  void init(const double *x, const double *y, int num_xy);
  void init(const double *x1, const double *y1, int num_xy1,
            const double *x2, const double *y2, int num_xy2);

  void clear();

 protected:
  void init1(const Polygons &polygons);

  void insertValue(double v, double *vals, int &num_vals);

 protected:
  std::vector<double> xvals_;
  std::vector<double> yvals_;
};

//------

/*
 * track inside state of each cell of x, y grid to perform functions
 * on rectilinear poylgons.
 */
class CXYValsInside : public CXYVals {
 public:
  enum {
    INSIDE1 = (1U<<0),
    INSIDE2 = (1U<<1)
  };

  typedef uint                                    InsideValue;
  typedef std::vector< std::vector<InsideValue> > InsideArray;

 public:
  CXYValsInside(const Polygons &polygons, bool init=false);

  CXYValsInside(const Polygon &polygon, bool init=false);
  CXYValsInside(const Polygon &polygon1, const Polygon &polygon2, bool init=false);

  CXYValsInside(const std::vector<double> &x, std::vector<double> &y,
                bool init=false);
  CXYValsInside(const std::vector<double> &x1, std::vector<double> &y1,
                const std::vector<double> &x2, std::vector<double> &y2,
                bool init=false);

  CXYValsInside(const double *x=0, const double *y=0, int num_xy=0,
                bool init=false);
  CXYValsInside(const double *x1, const double *y1, int num_xy1,
                const double *x2, const double *y2, int num_xy2,
                bool init=false);

 ~CXYValsInside() { }

  bool isOrValues() const { return orValues_; }
  void setOrValues(bool b) { orValues_ = b; }

  void initValues(const Polygons &polygons);

  void initValues(const Polygon &polygon);
  void initValues(const Polygon &polygon1, const Polygon &polygon2);

  void initValues(const CXYValsInside &xyvals, const Polygon &polygon);
  void initValues(const CXYValsInside &xyvals, const std::vector<double> &x,
                  std::vector<double> &y);
  void initValues(const CXYValsInside &xyvals, const double *x, const double *y, int num_xy);

  void initValues(const std::vector<double> &x, std::vector<double> &y);
  void initValues(const std::vector<double> &x1, std::vector<double> &y1,
                  const std::vector<double> &x2, std::vector<double> &y2);

  void initValues(const double *x, const double *y, int num_xy);
  void initValues(const double *x1, const double *y1, int num_xy1,
                  const double *x2, const double *y2, int num_xy2);

  const InsideArray &getInside() const { return inside_; }

  InsideArray &getInside() { return inside_; }

  int numXVals() const { return num_xvals_; }
  int numYVals() const { return num_yvals_; }

  bool isInside1(int ix, int iy) const { return (inside_[ix][iy] & CXYValsInside::INSIDE1); }
  bool isInside2(int ix, int iy) const { return (inside_[ix][iy] & CXYValsInside::INSIDE2); }

  bool isInside(int ix, int iy, InsideValue val=INSIDE1) const {
    return (inside_[ix][iy] == val);
  }

  void setInside(InsideValue val=INSIDE1);

  InsideValue insideVal(int ix, int iy) const { return inside_[ix][iy]; }
  void setInsideVal(int ix, int iy, InsideValue val=INSIDE1);

  void combineInside(InsideValue val=INSIDE1);

  void clear();

  bool getPolygons(Polygons &polygons, bool check_consistent=false) const;
  bool getPolygons(InsideValue ival, Polygons &polygons, bool check_consistent=false) const;

  bool getPolygon(Polygon &polygon, bool check_consistent=false) const;
  bool getPolygon(InsideValue ival, Polygon &polygon, bool check_consistent=false) const;

  bool getPolygon(std::vector<double> &x, std::vector<double> &y,
                  bool check_consistent=false) const;
  bool getPolygon(InsideValue ival, std::vector<double> &x, std::vector<double> &y,
                  bool check_consistent=false) const;

  bool getPolygon(double **xo, double **yo, int *num_xyo, bool check_consistent=false) const;
  bool getPolygon(InsideValue ival, double **xo, double **yo, int *num_xyo,
                  bool check_consistent=false) const;

  void setPolygonValue(const Polygon &poly, InsideValue val);

  double polygonArea(const Polygon &polygon) const;

  double valueArea(InsideValue ival) const;

  InsideValue polygonValue(const Polygon &polygon) const;

  int getNumPolygons() const;

  int  setPolygonInsideValues();
  void setPolygonInsideValues(const Polygons &polygons);

  bool fill(bool step=false);

  InsideValue findSmallest() const;

  bool fillDisconnectedRows();
  bool fillDisconnectedColumns();

  void print(std::ostream &os) const;

  static bool unitTest(const std::string &args="");

 private:
  void initValues1(const CXYValsInside &xyvals, const Polygon &polygon);
  void initValues1(const Polygons &polygons);

  void initMem();

  bool fill1(bool step, FillMode mode);

  bool fillHorizontal(InsideValue val, FillMode mode);
  bool fillVertical  (InsideValue val, FillMode mode);

  bool isOnEdge(const Polygon &poly) const;

  void polygonBounds(InsideValue val, int &x1, int &y1, int &x2, int &y2) const;

  bool isEmptyRow   (int iy) const;
  bool isEmptyColumn(int ix) const;

 private:
  InsideArray inside_;
  int         num_xvals_ { 0 };
  int         num_yvals_ { 0 };
  bool        orValues_ { false };
};

#endif
