#include "Array.hh"
#include <stdio.h>
#include <stdlib.h>

// Structural methods.

Array::Array(int cols, int rows, float missing) :
  _ncols(cols), _nrows(rows), _missing(missing), _size(_ncols * _nrows)
{
  _data = new float[_size];
  _rows = new float*[_nrows];
  for (int idx = 0; idx < _nrows; ++idx) _rows[idx] = _data + idx * _ncols;
  for (int idx = 0; idx < _size; ++idx)  _data[idx] = _missing;
}

Array::~Array()
{
  delete [] _data;
  delete [] _rows;
}

void Array::reset(void)
{
  for (int idx = 0; idx < _size; ++idx) _data[idx] = _missing;
}

Array &Array::operator=(float val)
{
  for (int idx = 0; idx < _size; ++idx) _data[idx] = val;
  return *this;
}


// Rotate the entries in an array in the x-direction by a given
// offset.  The effect is basically (for all valid r and c):
//
//    _rows[r][(c + shift) % _ncols] <- _rows[r][c]

void Array::rotate(int shift)
{
//  float tmp[shift > 0 ? shift : -shift];  //Doug 11/10: changed to lines below in order to make compatatble with fortran compiler
  float *tmp;
  if (shift>0) {
    tmp=(float*)malloc(sizeof(float)*shift);
  } else  {
    tmp=(float*)malloc(sizeof(float)*(-shift));
  }
  for (int r = 0; r < _nrows; ++r) {
    float *row = _rows[r];
    if (shift < 0) {
      for (int idx = 0; idx < -shift; ++idx)
        tmp[idx] = row[idx];
      for (int idx = -shift; idx < _ncols; ++idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < -shift; ++idx)
        row[_ncols + shift + idx] = tmp[idx];
    } else {
      for (int idx = 0; idx < shift; ++idx)
        tmp[idx] = row[_ncols - shift + idx];
      for (int idx = _ncols - shift - 1; idx >= 0; --idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < shift; ++idx)
        row[idx] = tmp[idx];
    }
  }
}


// Reverse order of rows.

void Array::flip_rows(void)
{
  for (int r = 0; r < _nrows / 2; ++r)
    for (int c = 0; c < _ncols; ++c) {
      float tmp = _rows[r][c];
      _rows[r][c] = _rows[_nrows - r - 1][c];
      _rows[_nrows - r - 1][c] = tmp;
    }
}


// Fix up any missing values in an array: rows with some missing and
// some non-missing values have the non-missing values adjacent to
// missing regions copied into the missing regions, then rows that are
// all missing values have values copied into them from adjacent rows
// that have no missing values.

void Array::fix_missing(void)
{
  // Fill in missing values for rows with only some missing values.

  int empty_rows = 0;
  for (int row = 0; row < _nrows; ++row) {
    if (all_missing(row)) { ++empty_rows; continue; }

    // Repeat until all the missing values in this row are filled in.

    int missing_count = count_missing(row);
    while (missing_count > 0) {

      // Find start and end points with missing values.

      float *vals = _rows[row];
      int stidx = _ncols - 1;
      while (stidx > 0) {
        if (vals[stidx] == _missing && vals[stidx - 1] != _missing) break;
        --stidx;
      }
      int enidx = (stidx + 1) % _ncols;
      while (vals[enidx] == _missing)
        enidx = (enidx + 1) % _ncols;


      // Fill in the missing values.
          
      float fill_val = vals[(stidx - 1 + _ncols) % _ncols];
      int idx = stidx;
      int swidx = (enidx > stidx) ?
        (enidx + stidx) / 2 :
        ((enidx + stidx + _ncols) / 2) % _ncols;
      while (idx != enidx) {
        if (idx == swidx) fill_val = vals[enidx];
        vals[idx] = fill_val;
        --missing_count;
        idx = (idx + 1) % _ncols;
      }
    }
  }


  // Fill in missing values for empty rows.

  while (empty_rows > 0) {

    // Find start and end points with missing values: note that at
    // this point, all rows are either all missing values, or have had
    // their missing values fixed.  This means we can detect missing
    // value rows just by checking a single value.

    int stidx = _nrows - 1;
    while (stidx > 0) {
      if (_rows[stidx][0] == _missing && _rows[stidx - 1][0] != _missing)
        break;
      --stidx;
    }
    int enidx = (stidx + 1) % _nrows;
    while (_rows[enidx][0] == _missing) enidx = (enidx + 1) % _nrows;


    // Fill in the missing values.
       
    int idx = stidx;
    int fill_row = (stidx - 1 + _nrows) % _nrows;
    int swidx = (enidx > stidx) ?
      (enidx + stidx) / 2 :
      ((enidx + stidx + _nrows) / 2) % _nrows;
    while (idx != enidx) {
      if (idx == swidx) fill_row = enidx;
      copy_row(fill_row, idx);
      --empty_rows;
      idx = (idx + 1) % _nrows;
    }
  }
}


// Mask an array according to a boolean mask -- entries are set to the
// missing value wherever the mask is false.

void Array::mask(const BoolArray &mask)
{
  if (mask.rows() != rows() || mask.cols() != cols())
    return;
  for (int r = 0; r < rows(); ++r)
    for (int c = 0; c < cols(); ++c)
      if (!mask(c, r)) _rows[r][c] = _missing;
}


// Regrid data from one array to another using bilinear interpolation.

void Array::regrid(const Array &in, const float *in_x, const float *in_y,
                   Array &out, const float *out_x, const float *out_y)
{
//   // We don't support extrapolation!
//   if (out_x[0] < in_x[0] || out_x[out.cols() - 1] > in_x[in.cols() - 1] ||
//       out_y[0] < in_y[0] || out_y[out.rows() - 1] > in_y[in.rows() - 1])
//     throw ArrayException("Extrapolation not supported in Array::regrid!");
      

  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point, i.e. the points that should
  // be used for interpolation.  Similarly for the y values.

//  int xbrack[out.cols()], ybrack[out.rows()]; //Doug 11/10: replaced with lines below the next to make compatatible with intel compiler
  int inidx = 0;
  int *xbrack;
  int *ybrack;
  xbrack=(int*)malloc(sizeof(int)*out.cols());
  ybrack=(int*)malloc(sizeof(int)*out.rows());
  for (int outidx = 0; outidx < out.cols(); ++outidx) {
    while (inidx <= in.cols() - 3 && in_x[inidx + 1] < out_x[outidx]) ++inidx;
    xbrack[outidx] = inidx;
  }
  inidx = 0;
  for (int outidx = 0; outidx < out.rows(); ++outidx) {
    while (inidx <= in.rows() - 3 && in_y[inidx + 1] < out_y[outidx]) ++inidx;
    ybrack[outidx] = inidx;
  }


  // Do bilinear interpolation, by taking the input x and y values
  // bracketing the output point (i.e. the corner points of the
  // smallest rectangle defined by the input axes that contains the
  // output point), and using the fractional x and y distances of the
  // output point across the bracketing intervals as interpolation
  // fractions.

  for (int row = 0; row < out.rows(); ++row) {
    int iym = ybrack[row], iyp = ybrack[row] + 1;
    float yfrac = (out_y[row] - in_y[iym]) / (in_y[iyp] - in_y[iym]);
    for (int col = 0; col < out.cols(); ++col) {
      int ixm = xbrack[col], ixp = xbrack[col] + 1;
      float xfrac = (out_x[col] - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
      out(col, row) =
        in(ixm, iym) * (1 - xfrac) * (1 - yfrac) +
        in(ixp, iym) * xfrac       * (1 - yfrac) +
        in(ixm, iyp) * (1 - xfrac) * yfrac       +
        in(ixp, iyp) * xfrac       * yfrac;
    }
  }
}


// Regrid data from one array to another using bilinear interpolation,
// assuming that the x-coordinate is circular (very useful for
// geophysical data...).

void Array::regrid_circular_x(const Array &in,
                              const float *in_x, const float *in_y,
                              Array &out,
                              const float *out_x, const float *out_y,
                              float min_x, float max_x)
{
  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point, i.e. the points that should
  // be used for interpolation.  Similarly for the y values.

  int *xbrack = new int[out.cols()], *ybrack = new int[out.rows()];
  int inidx = 0;
  for (int outidx = 0; outidx < out.cols(); ++outidx) {
    while (inidx <= in.cols() - 2 && in_x[inidx + 1] < out_x[outidx])
      ++inidx;
    xbrack[outidx] = inidx;
  }
  inidx = 0;
  for (int outidx = 0; outidx < out.rows(); ++outidx) {
    while (inidx <= in.rows() - 3 && in_y[inidx + 1] < out_y[outidx]) ++inidx;
    ybrack[outidx] = inidx;
  }


  // Do bilinear interpolation, by taking the input x and y values
  // bracketing the output point (i.e. the corner points of the
  // smallest rectangle defined by the input axes that contains the
  // output point), and using the fractional x and y distances of the
  // output point across the bracketing intervals as interpolation
  // fractions.

  for (int row = 0; row < out.rows(); ++row) {
    int iym = ybrack[row], iyp = ybrack[row] + 1;
    float yfrac = (out_y[row] - in_y[iym]) / (in_y[iyp] - in_y[iym]);
    for (int col = 0; col < out.cols(); ++col) {
      int ixm = xbrack[col], ixp = xbrack[col] + 1;
      if (ixm == in.cols() - 1) {
        float xfrac = (out_x[col] - in_x[ixm]) / (max_x - in_x[ixm]);
        out(col, row) =
          in(ixm, iym) * (1 - xfrac) * (1 - yfrac) +
          in(0, iym)   * xfrac       * (1 - yfrac) +
          in(ixm, iyp) * (1 - xfrac) * yfrac       +
          in(0, iyp)   * xfrac       * yfrac;
      } else {
        float xfrac = (out_x[col] - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
        out(col, row) =
          in(ixm, iym) * (1 - xfrac) * (1 - yfrac) +
          in(ixp, iym) * xfrac       * (1 - yfrac) +
          in(ixm, iyp) * (1 - xfrac) * yfrac       +
          in(ixp, iyp) * xfrac       * yfrac;
      }
    }
  }

  delete [] xbrack;
  delete [] ybrack;
}


// Regrid data from one array to another using area-weighted
// averaging.  (This is designed to be used when regridding from a
// finer to a coarser grid.)

void Array::regrid_average(const Array &in,
                           const float *in_x, const float *in_y,
                           Array &out,
                           const float *out_x1, const float *out_y1)
{
//  float out_x[out.cols() + 1], out_y[out.rows() + 1]; //Doug 11/10: replaced with lines below to make compatatble with intel compiler
  float *out_x;
  float *out_y;
  out_x=(float*)malloc(sizeof(float)*out.cols()+1);
  out_y=(float*)malloc(sizeof(float)*out.rows()+1);

  for (int idx = 1; idx < out.cols() - 1; ++idx)
    out_x[idx] = (out_x1[idx - 1] + out_x1[idx]) / 2.0;
  for (int idx = 1; idx < out.rows() - 1; ++idx)
    out_y[idx] = (out_y1[idx - 1] + out_y1[idx]) / 2.0;
  out_x[0] = out_x1[0] - (out_x[1] - out_x1[0]);
  out_y[0] = out_y1[0] - (out_y[1] - out_y1[0]);
  int tc = out.cols() - 1, tr = out.rows() - 1;
  out_x[tc] = out_x1[tc] + (out_x[tc] - out_x1[tc]);
  out_y[tr] = out_y1[tr] + (out_y[tr] - out_y1[tr]);


  if (out_x[0] < in_x[0]) out_x[0] = in_x[0];
  if (out_y[0] < in_y[0]) out_y[0] = in_y[0];
  if (out_x[out.cols()] > in_x[in.cols() - 1])
    out_x[out.cols()] = in_x[in.cols() - 1];
  if (out_y[out.rows()] > in_y[in.rows() - 1])
    out_y[out.rows()] = in_y[in.rows() - 1];


  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point.  Similarly for the y values.

//  int xbrack[out.cols() + 1], ybrack[out.rows() + 1]; //Doug 11/10: replaced with lines below to make compatatble with intel compiler
  int idx = 0;
  int *xbrack;
  int *ybrack;
  xbrack=(int*)malloc(sizeof(int)*out.cols()+1);
  ybrack=(int*)malloc(sizeof(int)*out.rows()+1);
  for (int outidx = 0; outidx < out.cols() + 1; ++outidx) {
    while (idx <= in.cols() - 3 && in_x[idx + 1] < out_x[outidx]) ++idx;
    xbrack[outidx] = idx;
  }
  idx = 0;
  for (int outidx = 0; outidx < out.rows() + 1; ++outidx) {
    while (idx <= in.rows() - 3 && in_y[idx + 1] < out_y[outidx]) ++idx;
    ybrack[outidx] = idx;
  }


  // Do averaging.

  for (int row = 0; row < out.rows(); ++row) {
    int iyfbm = ybrack[row], iyfbp = ybrack[row] + 1;
    int iyftm = ybrack[row + 1], iyftp = ybrack[row + 1] + 1;
    for (int col = 0; col < out.cols(); ++col) {
      int ixflm = xbrack[col], ixflp = xbrack[col] + 1;
      int ixfrm = xbrack[col + 1], ixfrp = xbrack[col + 1] + 1;

      double sum = 0.0, area = 0.0, darea;
      double ps = 0.0, pa = 0.0, ns = 0.0, na = 0.0;

      // Fully covered cells.

      for (int iy = iyfbp; iy < iyftm; ++iy)
        for (int ix = ixflp; ix < ixfrm; ++ix) {
          darea = (in_y[iy + 1] - in_y[iy]) * (in_x[ix + 1] - in_x[ix]);
          sum += in(ix, iy) * darea;
          area += darea;
          if (in(ix, iy) > 0) { ps += in(ix, iy) * darea; pa += darea; }
          else                { ns += in(ix, iy) * darea; na += darea; }
        }

      // Edges.

      double frac1 =
        1.0 - (out_x[col] - in_x[ixflm]) / (in_x[ixflp] - in_x[ixflm]);
      double frac2 =
        (out_x[col + 1] - in_x[ixfrm]) / (in_x[ixfrp] - in_x[ixfrm]);
      for (int iy = iyfbp; iy < iyftm; ++iy) {
        double darea1 = frac1 * (in_y[iy + 1] - in_y[iy]) *
          (in_x[ixflp] - in_x[ixflm]);
        double darea2 = frac2 * (in_y[iy + 1] - in_y[iy]) *
          (in_x[ixfrp] - in_x[ixfrm]);
        sum += in(ixflm, iy) * darea1 + in(ixfrm, iy) * darea2;
        area += darea1 + darea2;
        if (in(ixflm, iy) > 0) { ps += in(ixflm, iy) * darea; pa += darea; }
        else                   { ns += in(ixflm, iy) * darea; na += darea; }
        if (in(ixfrm, iy) > 0) { ps += in(ixfrm, iy) * darea; pa += darea; }
        else                   { ns += in(ixfrm, iy) * darea; na += darea; }
      }

      frac1 = 1.0 - (out_y[row] - in_y[iyfbm]) / (in_y[iyfbp] - in_y[iyfbm]);
      frac2 = (out_y[row + 1] - in_y[iyftm]) / (in_y[iyftp] - in_y[iyftm]);
      for (int ix = ixflp; ix < ixfrm; ++ix) {
        double darea1 = frac1 * (in_x[ix + 1] - in_x[ix]) *
          (in_y[iyfbp] - in_y[iyfbm]);
        double darea2 = frac2 * (in_x[ix + 1] - in_x[ix]) *
          (in_y[iyftp] - in_y[iyftm]);
        sum += in(ix, iyfbm) * darea1 + in(ix, iyftm) * darea2;
        area += darea1 + darea2;
        if (in(ix, iyfbm) > 0) { ps += in(ix, iyfbm) * darea; pa += darea; }
        else                   { ns += in(ix, iyfbm) * darea; na += darea; }
        if (in(ix, iyftm) > 0) { ps += in(ix, iyftm) * darea; pa += darea; }
        else                   { ns += in(ix, iyftm) * darea; na += darea; }
      }


      // Corners.

      double fracbl =
        (1.0 - (out_x[col] - in_x[ixflm]) / (in_x[ixflp] - in_x[ixflm])) *
        (1.0 - (out_y[row] - in_y[iyfbm]) / (in_y[iyfbp] - in_y[iyfbm]));
      double fracbr =
        (out_x[col + 1] - in_x[ixfrm]) / (in_x[ixfrp] - in_x[ixfrm]) *
        (1.0 - (out_y[row] - in_y[iyfbm]) / (in_y[iyfbp] - in_y[iyfbm]));
      double fractl =
        (1.0 - (out_x[col] - in_x[ixflm]) / (in_x[ixflp] - in_x[ixflm])) *
        (out_y[row + 1] - in_y[iyftm]) / (in_y[iyftp] - in_y[iyftm]);
      double fractr =
        (out_x[col + 1] - in_x[ixfrm]) / (in_x[ixfrp] - in_x[ixfrm]) *
        (out_y[row + 1] - in_y[iyftm]) / (in_y[iyftp] - in_y[iyftm]);
      double dareabl = fracbl * (in_y[iyfbp] - in_y[iyfbm]) *
        (in_x[ixflp] - in_x[ixflm]);
      double dareabr = fracbr * (in_y[iyfbp] - in_y[iyfbm]) *
        (in_x[ixfrp] - in_x[ixfrm]);
      double dareatl = fractl * (in_y[iyftp] - in_y[iyftm]) *
        (in_x[ixflp] - in_x[ixflm]);
      double dareatr = fractr * (in_y[iyftp] - in_y[iyftm]) *
        (in_x[ixfrp] - in_x[ixfrm]);
      sum += in(ixflm, iyfbm) * dareabl + in(ixfrm, iyfbm) * dareabr
           + in(ixflm, iyftm) * dareatl + in(ixfrm, iyftm) * dareatr;
      area += dareabl + dareabr + dareatl + dareatr;
      if (in(ixflm, iyfbm) > 0) {
        ps += in(ixflm, iyfbm) * dareabl; pa += dareabl;
      } else {
        ns += in(ixflm, iyfbm) * dareabl; na += dareabl; }
      if (in(ixfrm, iyfbm) > 0) {
        ps += in(ixfrm, iyfbm) * dareabr; pa += dareabr;
      } else {
        ns += in(ixfrm, iyfbm) * dareabr; na += dareabr; }
      if (in(ixflm, iyftm) > 0) {
        ps += in(ixflm, iyftm) * dareatl; pa += dareatl;
      } else {
        ns += in(ixflm, iyftm) * dareatl; na += dareatl; }
      if (in(ixfrm, iyftm) > 0) {
        ps += in(ixfrm, iyftm) * dareatr; pa += dareatr;
      } else {
        ns += in(ixfrm, iyftm) * dareatr; na += dareatr; }

      if (pa > 0.0 && na > 0.0) {
        if (pa > na)
          out(col, row) = ps / pa;
        else
          out(col, row) = ns / na;
      } else
        out(col, row) = sum / area;
    }
  }
}


// Regrid data from one array to another using a discrete area-based
// interpolation method (used for things like soil types).

void Array::regrid_discrete(const Array &in,
                            const float *in_x, const float *in_y,
                            Array &out,
                            const float *out_x, const float *out_y)
{
//   // We don't support extrapolation!
//   if (out_x[0] < in_x[0] || out_x[out.cols() - 1] > in_x[in.cols() - 1] ||
//       out_y[0] < in_y[0] || out_y[out.rows() - 1] > in_y[in.rows() - 1])
//     throw ArrayException("Extrapolation not supported "
//                          "in Array::regrid_discrete!");
      

  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point, i.e. the points that should
  // be used for interpolation.  Similarly for the y values.

//  int xbrack[out.cols()], ybrack[out.rows()]; //Doug 11/10: rep;aced with lines below to make compatable with intel compiler
  int *xbrack;
  int *ybrack;
  xbrack=(int*)malloc(sizeof(int)*out.cols());  
  ybrack=(int*)malloc(sizeof(int)*out.rows());
  int inidx = 0;
  for (int outidx = 0; outidx < out.cols(); ++outidx) {
    while (inidx <= in.cols() && in_x[inidx + 1] < out_x[outidx]) ++inidx;
    xbrack[outidx] = inidx;
  }
  inidx = 0;
  for (int outidx = 0; outidx < out.rows(); ++outidx) {
    while (inidx <= in.rows() && in_y[inidx + 1] < out_y[outidx]) ++inidx;
    ybrack[outidx] = inidx;
  }


  // Do discrete interpolation, by finding the maximum overlap between
  // the output grid square and the four possible input grid squares
  // that it may cover.  This is done using the input x and y values
  // bracketing the output point (i.e. the corner points of the
  // smallest rectangle defined by the input axes that contains the
  // output point), and using the fractional x and y distances of the
  // output point across the bracketing intervals to calculate
  // proportional area overlaps.

  for (int row = 0; row < out.rows(); ++row) {
    int iym = ybrack[row], iyp = ybrack[row] + 1;
    float yfrac = (out_y[row] - in_y[iym]) / (in_y[iyp] - in_y[iym]);
    for (int col = 0; col < out.cols(); ++col) {
      int ixm = xbrack[col], ixp = xbrack[col] + 1;
      float xfrac = (out_x[col] - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
      int idx = 0;
      float max_area = (1 - xfrac) * (1 - yfrac);
      float area = xfrac * (1 - yfrac);
      if (area > max_area) { idx = 1; max_area = area; }
      area = (1 - xfrac) * yfrac;
      if (area > max_area) { idx = 2; max_area = area; }
      area = xfrac * yfrac;
      if (area > max_area) { idx = 3; max_area = area; }
      switch (idx) {
      case 0: out(col, row) = in(ixm, iym);  break;
      case 1: out(col, row) = in(ixp, iym);  break;
      case 2: out(col, row) = in(ixm, iyp);  break;
      case 3: out(col, row) = in(ixp, iyp);  break;
      }
    }
  }
}


// Is a row all missing values?

bool Array::all_missing(int row)
{
  float *vals = _rows[row];
  for (int col = 0; col < _ncols; ++col)
    if (vals[col] != _missing)
      return false;
  return true;
}


// Count the missing values in a row.

int Array::count_missing(int row)
{
  int count = 0;
  float *vals = _rows[row];
  for (int col = 0; col < _ncols; ++col)
    if (vals[col] == _missing) ++count;
  return count;
}


// Copy one row to another.

void Array::copy_row(int from_row, int to_row)
{
  float *from_vals = _rows[from_row], *to_vals = _rows[to_row];
  for (int col = 0; col < _ncols; ++col)
    to_vals[col] = from_vals[col];
}


Array3::Array3(int cols, int rows, int blocks, float missing) :
  _ncols(cols), _nrows(rows), _nblocks(blocks), _missing(missing),
  _size(_ncols * _nrows * _nblocks)
{
  _blocks = new Array* [_nblocks];
  for (int bl = 0; bl < _nblocks; ++bl)
    _blocks[bl] = new Array(cols, rows, missing);
}

Array3::~Array3()
{
  for (int bl = 0; bl < _nblocks; ++bl) delete _blocks[bl];
  delete [] _blocks;
}



// Structural methods.

IntArray::IntArray(int cols, int rows, int missing) :
  _ncols(cols), _nrows(rows), _missing(missing), _size(_ncols * _nrows)
{
  _data = new int[_size];
  _rows = new int*[_nrows];
  for (int idx = 0; idx < _nrows; ++idx) _rows[idx] = _data + idx * _ncols;
  for (int idx = 0; idx < _size; ++idx)  _data[idx] = _missing;
}

IntArray::~IntArray()
{
  delete [] _data;
  delete [] _rows;
}

void IntArray::reset(void)
{
  for (int idx = 0; idx < _size; ++idx) _data[idx] = _missing;
}

IntArray &IntArray::operator=(int val)
{
  for (int idx = 0; idx < _size; ++idx) _data[idx] = val;
  return *this;
}


// Rotate the entries in an array in the x-direction by a given
// offset.  The effect is basically (for all valid r and c):
//
//    _rows[r][(c + shift) % _ncols] <- _rows[r][c]

void IntArray::rotate(int shift)
{
//  int tmp[shift > 0 ? shift : -shift];//Doug 11/10: changed to lines below in order to make compatatble with fortran compiler
  int *tmp;
  if (shift>0) {
    tmp=(int*)malloc(sizeof(int)*shift);
  } else  {
    tmp=(int*)malloc(sizeof(int)*(-shift));
  }

  for (int r = 0; r < _nrows; ++r) {
    int *row = _rows[r];
    if (shift < 0) {
      for (int idx = 0; idx < -shift; ++idx)
        tmp[idx] = row[idx];
      for (int idx = -shift; idx < _ncols; ++idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < -shift; ++idx)
        row[_ncols + shift + idx] = tmp[idx];
    } else {
      for (int idx = 0; idx < shift; ++idx)
        tmp[idx] = row[_ncols - shift + idx];
      for (int idx = _ncols - shift - 1; idx >= 0; --idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < shift; ++idx)
        row[idx] = tmp[idx];
    }
  }
}


// Fix up any missing values in an array: rows with some missing and
// some non-missing values have the non-missing values adjacent to
// missing regions copied into the missing regions, then rows that are
// all missing values have values copied into them from adjacent rows
// that have no missing values.

void IntArray::fix_missing(void)
{
  // Fill in missing values for rows with only some missing values.

  int empty_rows = 0;
  for (int row = 0; row < _nrows; ++row) {
    if (all_missing(row)) { ++empty_rows; continue; }

    // Repeat until all the missing values in this row are filled in.

    int missing_count = count_missing(row);
    while (missing_count > 0) {

      // Find start and end points with missing values.

      int *vals = _rows[row];
      int stidx = _ncols - 1;
      while (stidx > 0) {
        if (vals[stidx] == _missing && vals[stidx - 1] != _missing) break;
        --stidx;
      }
      int enidx = (stidx + 1) % _ncols;
      while (vals[enidx] == _missing)
        enidx = (enidx + 1) % _ncols;


      // Fill in the missing values.
          
      int fill_val = vals[(stidx - 1 + _ncols) % _ncols];
      int idx = stidx;
      int swidx = (enidx > stidx) ?
        (enidx + stidx) / 2 :
        ((enidx + stidx + _ncols) / 2) % _ncols;
      while (idx != enidx) {
        if (idx == swidx) fill_val = vals[enidx];
        vals[idx] = fill_val;
        --missing_count;
        idx = (idx + 1) % _ncols;
      }
    }
  }


  // Fill in missing values for empty rows.

  while (empty_rows > 0) {

    // Find start and end points with missing values: note that at
    // this point, all rows are either all missing values, or have had
    // their missing values fixed.  This means we can detect missing
    // value rows just by checking a single value.

    int stidx = _nrows - 1;
    while (stidx > 0) {
      if (_rows[stidx][0] == _missing && _rows[stidx - 1][0] != _missing)
        break;
      --stidx;
    }
    int enidx = (stidx + 1) % _nrows;
    while (_rows[enidx][0] == _missing) enidx = (enidx + 1) % _nrows;


    // Fill in the missing values.
       
    int idx = stidx;
    int fill_row = (stidx - 1 + _nrows) % _nrows;
    int swidx = (enidx > stidx) ?
      (enidx + stidx) / 2 :
      ((enidx + stidx + _nrows) / 2) % _nrows;
    while (idx != enidx) {
      if (idx == swidx) fill_row = enidx;
      copy_row(fill_row, idx);
      --empty_rows;
      idx = (idx + 1) % _nrows;
    }
  }
}


// Mask an array according to a boolean mask -- entries are set to the
// missing value wherever the mask is false.

void IntArray::mask(const BoolArray &mask)
{
  if (mask.rows() != rows() || mask.cols() != cols())
    return;
  for (int r = 0; r < rows(); ++r)
    for (int c = 0; c < cols(); ++c)
      if (!mask(c, r)) _rows[r][c] = _missing;
}


int IntArray::interpolate(float x, float y, float *in_x, float *in_y)
{
  int xbrack, ybrack;
  int inidx = 0;
  while (inidx <= cols() - 3 && in_x[inidx + 1] < x) ++inidx;
  xbrack = inidx;
  inidx = 0;
  while (inidx <= rows() - 3 && in_y[inidx + 1] < y) ++inidx;
  ybrack = inidx;

  int iym = ybrack, iyp = ybrack + 1;
  float yfrac = (y - in_y[iym]) / (in_y[iyp] - in_y[iym]);
  int ixm = xbrack, ixp = xbrack + 1;
  float xfrac = (x - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
  float out =
    (*this)(ixm, iym) * (1 - xfrac) * (1 - yfrac) +
    (*this)(ixp, iym) * xfrac       * (1 - yfrac) +
    (*this)(ixm, iyp) * (1 - xfrac) * yfrac       +
    (*this)(ixp, iyp) * xfrac       * yfrac;

  return static_cast<int>(out);
}


// Regrid data from one array to another using a discrete area-based
// interpolation method (used for things like soil types).

void IntArray::regrid_discrete(const IntArray &in,
                               const float *in_x, const float *in_y,
                               IntArray &out,
                               const float *out_x, const float *out_y)
{
//   // We don't support extrapolation!
//   if (out_x[0] < in_x[0] || out_x[out.cols() - 1] > in_x[in.cols() - 1] ||
//       out_y[0] < in_y[0] || out_y[out.rows() - 1] > in_y[in.rows() - 1])
//     throw ArrayException("Extrapolation not supported "
//                          "in IntArray::regrid_discrete!");
      

  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point, i.e. the points that should
  // be used for interpolation.  Similarly for the y values.

//  int xbrack[out.cols()], ybrack[out.rows()]; //Doug 11/10: relaced with line below

  int *xbrack;
  int *ybrack;
  xbrack=(int*)malloc(sizeof(int)*out.cols());
  ybrack=(int*)malloc(sizeof(int)*out.rows());

  int inidx = 0;
  for (int outidx = 0; outidx < out.cols(); ++outidx) {
    while (inidx <= in.cols() && in_x[inidx + 1] < out_x[outidx]) ++inidx;
    xbrack[outidx] = inidx;
  }
  inidx = 0;
  for (int outidx = 0; outidx < out.rows(); ++outidx) {
    while (inidx <= in.rows() && in_y[inidx + 1] < out_y[outidx]) ++inidx;
    ybrack[outidx] = inidx;
  }


  // Do discrete interpolation, by finding the maximum overlap between
  // the output grid square and the four possible input grid squares
  // that it may cover.  This is done using the input x and y values
  // bracketing the output point (i.e. the corner points of the
  // smallest rectangle defined by the input axes that contains the
  // output point), and using the fractional x and y distances of the
  // output point across the bracketing intervals to calculate
  // proportional area overlaps.

  for (int row = 0; row < out.rows(); ++row) {
    int iym = ybrack[row], iyp = ybrack[row] + 1;
    float yfrac = (out_y[row] - in_y[iym]) / (in_y[iyp] - in_y[iym]);
    for (int col = 0; col < out.cols(); ++col) {
      int ixm = xbrack[col], ixp = xbrack[col] + 1;
      float xfrac = (out_x[col] - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
      int idx = 0;
      float max_area = (1 - xfrac) * (1 - yfrac);
      float area = xfrac * (1 - yfrac);
      if (area > max_area) { idx = 1; max_area = area; }
      area = (1 - xfrac) * yfrac;
      if (area > max_area) { idx = 2; max_area = area; }
      area = xfrac * yfrac;
      if (area > max_area) { idx = 3; max_area = area; }
      if (ixm >= in.cols()) ixm = in.cols() - 1;
      if (ixp >= in.cols()) ixp = in.cols() - 1;
      if (iym >= in.rows()) iym = in.rows() - 1;
      if (iyp >= in.rows()) iyp = in.rows() - 1;
      switch (idx) {
      case 0: out(col, row) = in(ixm, iym);  break;
      case 1: out(col, row) = in(ixp, iym);  break;
      case 2: out(col, row) = in(ixm, iyp);  break;
      case 3: out(col, row) = in(ixp, iyp);  break;
      }
    }
  }
}


// Is a row all missing values?

bool IntArray::all_missing(int row)
{
  int *vals = _rows[row];
  for (int col = 0; col < _ncols; ++col)
    if (vals[col] != _missing)
      return false;
  return true;
}


// Count the missing values in a row.

int IntArray::count_missing(int row)
{
  int count = 0;
  int *vals = _rows[row];
  for (int col = 0; col < _ncols; ++col)
    if (vals[col] == _missing) ++count;
  return count;
}


// Copy one row to another.

void IntArray::copy_row(int from_row, int to_row)
{
  int *from_vals = _rows[from_row], *to_vals = _rows[to_row];
  for (int col = 0; col < _ncols; ++col)
    to_vals[col] = from_vals[col];
}



BoolArray::BoolArray(int cols, int rows) :
  _ncols(cols), _nrows(rows), _size(_ncols * _nrows)
{
  _data = new bool[_size];
  _rows = new bool*[_nrows];
  for (int idx = 0; idx < _nrows; ++idx) _rows[idx] = _data + idx * _ncols;
  for (int idx = 0; idx < _size; ++idx)  _data[idx] = false;
}

BoolArray::~BoolArray()
{
  delete [] _data;
  delete [] _rows;
}


// Rotate the entries in an array in the x-direction by a given
// offset.  The effect is basically (for all valid r and c):
//
//    _rows[r][(c + shift) % _ncols] <- _rows[r][c]

void BoolArray::rotate(int shift)
{
//  bool tmp[shift > 0 ? shift : -shift];  //Doug 11/10: changed to lines below in order to make compatatble with fortran compiler
  bool *tmp;
  if (shift>0) {
    tmp=(bool*)malloc(sizeof(bool)*shift);
  } else  {
    tmp=(bool*)malloc(sizeof(bool)*(-shift));
  }
  for (int r = 0; r < _nrows; ++r) {
    bool *row = _rows[r];
    if (shift < 0) {
      for (int idx = 0; idx < -shift; ++idx)
        tmp[idx] = row[idx];
      for (int idx = -shift; idx < _ncols; ++idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < -shift; ++idx)
        row[_ncols + shift + idx] = tmp[idx];
    } else {
      for (int idx = 0; idx < shift; ++idx)
        tmp[idx] = row[_ncols - shift + idx];
      for (int idx = _ncols - shift - 1; idx >= 0; --idx)
        row[idx + shift] = row[idx];
      for (int idx = 0; idx < shift; ++idx)
        row[idx] = tmp[idx];
    }
  }
}


// Regrid data from one array to another using a discrete area-based
// interpolation method (used for things like soil types).

void BoolArray::regrid_discrete(const BoolArray &in,
                                const float *in_x, const float *in_y,
                                BoolArray &out,
                                const float *out_x, const float *out_y)
{
//   // We don't support extrapolation!
//   if (out_x[0] < in_x[0] || out_x[out.cols() - 1] > in_x[in.cols() - 1] ||
//       out_y[0] < in_y[0] || out_y[out.rows() - 1] > in_y[in.rows() - 1])
//     throw ArrayException("Extrapolation not supported "
//                          "in IntArray::regrid_discrete!");
      

  // Find bracketing index information: brack_x[col] gives the index
  // into the in_x array of the greatest value smaller than
  // corresponding value in the out_x array.  brack_x[col] and
  // brack_x[col] + 1 thus give the indices into the in_x array of the
  // points bracketing the output point, i.e. the points that should
  // be used for interpolation.  Similarly for the y values.

//  int xbrack[out.cols()], ybrack[out.rows()];//Doug 11/10: replaced with lines below the next to make compatatible with intel compiler
  int inidx = 0;
  int *xbrack;
  int *ybrack;
  xbrack=(int*)malloc(sizeof(int)*out.cols());
  ybrack=(int*)malloc(sizeof(int)*out.rows());
  for (int outidx = 0; outidx < out.cols(); ++outidx) {
    while (inidx <= in.cols() && in_x[inidx + 1] < out_x[outidx]) ++inidx;
    xbrack[outidx] = inidx;
  }
  inidx = 0;
  for (int outidx = 0; outidx < out.rows(); ++outidx) {
    while (inidx <= in.rows() && in_y[inidx + 1] < out_y[outidx]) ++inidx;
    ybrack[outidx] = inidx;
  }


  // Do discrete interpolation, by finding the maximum overlap between
  // the output grid square and the four possible input grid squares
  // that it may cover.  This is done using the input x and y values
  // bracketing the output point (i.e. the corner points of the
  // smallest rectangle defined by the input axes that contains the
  // output point), and using the fractional x and y distances of the
  // output point across the bracketing intervals to calculate
  // proportional area overlaps.

  for (int row = 0; row < out.rows(); ++row) {
    int iym = ybrack[row], iyp = ybrack[row] + 1;
    float yfrac = (out_y[row] - in_y[iym]) / (in_y[iyp] - in_y[iym]);
    for (int col = 0; col < out.cols(); ++col) {
      int ixm = xbrack[col], ixp = xbrack[col] + 1;
      float xfrac = (out_x[col] - in_x[ixm]) / (in_x[ixp] - in_x[ixm]);
      int idx = 0;
      float max_area = (1 - xfrac) * (1 - yfrac);
      float area = xfrac * (1 - yfrac);
      if (area > max_area) { idx = 1; max_area = area; }
      area = (1 - xfrac) * yfrac;
      if (area > max_area) { idx = 2; max_area = area; }
      area = xfrac * yfrac;
      if (area > max_area) { idx = 3; max_area = area; }
      switch (idx) {
      case 0: out(col, row) = in(ixm, iym);  break;
      case 1: out(col, row) = in(ixp, iym);  break;
      case 2: out(col, row) = in(ixm, iyp);  break;
      case 3: out(col, row) = in(ixp, iyp);  break;
      }
    }
  }
}