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
  struct Polygon {
    Polygon() { }

    Polygon(const std::vector<double> &x1, std::vector<double> &y1) :
     x(x1), y(y1) {
      assert(x.size() == y.size());
    }

    int size() const { return x.size(); }

    bool isInside(double xm, double ym) const;

    std::vector<double> x;
    std::vector<double> y;
  };

  typedef std::vector<Polygon> Polygons;

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

  double xval(int i) const { return xvals_[i]; }
  double yval(int i) const { return yvals_[i]; }

  const double *getXVals   () const { return &xvals_[0]; }
  int           getNumXVals() const { return xvals_.size(); }

  const double *getYVals   () const { return &yvals_[0]; }
  int           getNumYVals() const { return yvals_.size(); }

  void init(const Polygon &polygon);
  void init(const Polygon &polygon1, const Polygon &polygon2);

  void init(const std::vector<double> &x, std::vector<double> &y);
  void init(const std::vector<double> &x1, std::vector<double> &y1,
            const std::vector<double> &x2, std::vector<double> &y2);

  void init(const double *x, const double *y, int num_xy);
  void init(const double *x1, const double *y1, int num_xy1,
            const double *x2, const double *y2, int num_xy2);

  void init(const Polygons &polygons);

  void clear();

 private:
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

  bool isInside1(int i, int j) const { return (inside_[i][j] & CXYValsInside::INSIDE1); }
  bool isInside2(int i, int j) const { return (inside_[i][j] & CXYValsInside::INSIDE2); }

  bool isInside(int i, int j, InsideValue val=INSIDE1) const { return (inside_[i][j] == val); }

  void setInside(InsideValue val=INSIDE1);

  void combineInside(InsideValue val=INSIDE1);

  void clear();

  bool getPolygons(Polygons &polygons, bool check_consistent=false) const;
  bool getPolygons(InsideValue inside_val, Polygons &polygons, bool check_consistent=false) const;

  bool getPolygon(Polygon &polygon, bool check_consistent=false) const;
  bool getPolygon(InsideValue inside_val, Polygon &polygon, bool check_consistent=false) const;

  bool getPolygon(std::vector<double> &x, std::vector<double> &y,
                  bool check_consistent=false) const;
  bool getPolygon(InsideValue inside_val, std::vector<double> &x, std::vector<double> &y,
                  bool check_consistent=false) const;

  bool getPolygon(double **xo, double **yo, int *num_xyo, bool check_consistent=false) const;
  bool getPolygon(InsideValue inside_val, double **xo, double **yo, int *num_xyo,
                  bool check_consistent=false) const;

  void removePolygon(const Polygon &poly);

  void fill(bool disconnected=false);

  void print(std::ostream &os) const;

  static bool unitTest(const std::string &args="");

 private:
  void initMem();

  bool fillH();
  bool fillV();

  bool fillHoles();

  void fillDisconnected();

  bool isEmptyRow(int iy) const;

 private:
  InsideArray inside_;
  int         num_xvals_ { 0 };
  int         num_yvals_ { 0 };
};

#endif
