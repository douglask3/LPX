
#ifndef _ARRAY_HH_
#define _ARRAY_HH_

class BoolArray;

class ArrayException {
public:

  ArrayException(const char *msg) : _msg(msg) { }
  const char *msg(void) const { return _msg; }

private:

  const char *_msg;
};


class Array {
public:

  Array(int cols, int rows, float missing = 0.0);

  ~Array();

  int rows(void) const { return _nrows; }
  int cols(void) const { return _ncols; }
  float missing(void) const { return _missing; }
  long size(void) const { return _size; }

  float *data(void) { return _data; }
  float *operator()(int row) { return _rows[row]; }
  float &operator()(int col, int row) { return _rows[row][col]; }
  const float &operator()(int col, int row) const { return _rows[row][col]; }

  Array &operator=(float val);

  void reset(void);

  void rotate(int shift);
  void flip_rows(void);
  void fix_missing(void);
  void mask(const BoolArray &mask);

  static void regrid(const Array &in, const float *in_x, const float *in_y,
                     Array &out, const float *out_x, const float *out_y);

  static void regrid_circular_x(const Array &in,
                                const float *in_x, const float *in_y,
                                Array &out,
                                const float *out_x, const float *out_y,
                                float min_x, float max_x);

  static void regrid_average(const Array &in,
                             const float *in_x, const float *in_y,
                             Array &out,
                             const float *out_x, const float *out_y);

  static void regrid_discrete(const Array &in,
                              const float *in_x, const float *in_y,
                              Array &out,
                              const float *out_x, const float *out_y);

private:

  bool all_missing(int row);
  int count_missing(int row);
  void copy_row(int from_row, int to_row);

  int _ncols;
  int _nrows;
  float _missing;

  long _size;
  float *_data;
  float **_rows;
};


class Array3 {
public:

  Array3(int cols, int rows, int blocks, float missing = 0.0);
  ~Array3();

  int rows(void) const { return _nrows; }
  int cols(void) const { return _ncols; }
  int blocks(void) const { return _nblocks; }
  float missing(void) const { return _missing; }
  long size(void) const { return _size; }

  Array &operator()(int block) { return *_blocks[block]; }
  float *operator()(int row, int block) { return (*_blocks[block])(row); }
  float &operator()(int col, int row, int block)
  { return (*_blocks[block])(col, row); }
  const float &operator()(int col, int row, int block) const
  { return (*_blocks[block])(col, row); }

private:

  int _ncols;
  int _nrows;
  int _nblocks;
  float _missing;

  long _size;
  Array **_blocks;
};


class IntArray {
public:

  IntArray(int cols, int rows, int missing = 0);

  ~IntArray();

  int rows(void) const { return _nrows; }
  int cols(void) const { return _ncols; }
  int missing(void) const { return _missing; }
  long size(void) const { return _size; }

  int *data(void) { return _data; }
  int *operator()(int row) { return _rows[row]; }
  int &operator()(int col, int row) { return _rows[row][col]; }
  const int &operator()(int col, int row) const { return _rows[row][col]; }

  IntArray &operator=(int val);

  void reset(void);

  void rotate(int shift);

  void fix_missing(void);
  void mask(const BoolArray &mask);

  int interpolate(float x, float y, float *in_x, float *in_y);

  static void regrid_discrete(const IntArray &in,
                              const float *in_x, const float *in_y,
                              IntArray &out,
                              const float *out_x, const float *out_y);

private:

  bool all_missing(int row);
  int count_missing(int row);
  void copy_row(int from_row, int to_row);

  int _ncols;
  int _nrows;
  int _missing;

  long _size;
  int *_data;
  int **_rows;
};


class BoolArray {
public:

  BoolArray(int cols, int rows);
  BoolArray(const BoolArray &other);

  ~BoolArray();

  int rows(void) const { return _nrows; }
  int cols(void) const { return _ncols; }
  long size(void) const { return _size; }

  bool &operator()(int col, int row) { return _rows[row][col]; }
  const bool &operator()(int col, int row) const { return _rows[row][col]; }

  void rotate(int shift);

  static void regrid_discrete(const BoolArray &in,
                              const float *in_x, const float *in_y,
                              BoolArray &out,
                              const float *out_x, const float *out_y);

private:

  int _ncols;
  int _nrows;

  long _size;
  bool *_data;
  bool **_rows;
};


#endif
