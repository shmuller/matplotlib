/* -*- mode: c++; c-basic-offset: 4 -*- */

#include "agg_py_path_iterator.h"
#include "agg_py_transforms.h"
#include "path_converters.h"

#include <limits>
#include <math.h>

#include "CXX/Extensions.hxx"

#include "agg_conv_contour.h"
#include "agg_conv_curve.h"
#include "agg_conv_stroke.h"
#include "agg_conv_transform.h"
#include "agg_path_storage.h"
#include "agg_trans_affine.h"

struct XY
{
    double x;
    double y;

    XY(double x_, double y_) : x(x_), y(y_) {}
};

//
// The following function was found in the Agg 2.3 examples (interactive_polygon.cpp).
// It has been generalized to work on (possibly curved) polylines, rather than
// just polygons.  The original comments have been kept intact.
//  -- Michael Droettboom 2007-10-02
//
//======= Crossings Multiply algorithm of InsideTest ========================
//
// By Eric Haines, 3D/Eye Inc, erich@eye.com
//
// This version is usually somewhat faster than the original published in
// Graphics Gems IV; by turning the division for testing the X axis crossing
// into a tricky multiplication test this part of the test became faster,
// which had the additional effect of making the test for "both to left or
// both to right" a bit slower for triangles than simply computing the
// intersection each time.  The main increase is in triangle testing speed,
// which was about 15% faster; all other polygon complexities were pretty much
// the same as before.  On machines where division is very expensive (not the
// case on the HP 9000 series on which I tested) this test should be much
// faster overall than the old code.  Your mileage may (in fact, will) vary,
// depending on the machine and the test data, but in general I believe this
// code is both shorter and faster.  This test was inspired by unpublished
// Graphics Gems submitted by Joseph Samosky and Mark Haigh-Hutchinson.
// Related work by Samosky is in:
//
// Samosky, Joseph, "SectionView: A system for interactively specifying and
// visualizing sections through three-dimensional medical image data",
// M.S. Thesis, Department of Electrical Engineering and Computer Science,
// Massachusetts Institute of Technology, 1993.
//
// Shoot a test ray along +X axis.  The strategy is to compare vertex Y values
// to the testing point's Y and quickly discard edges which are entirely to one
// side of the test ray.  Note that CONVEX and WINDING code can be added as
// for the CrossingsTest() code; it is left out here for clarity.
//
// Input 2D polygon _pgon_ with _numverts_ number of vertices and test point
// _point_, returns 1 if inside, 0 if outside.
template<class T>
static void
point_in_path_impl(const void* const points_, const size_t s0,
                   const size_t s1, const size_t n, T& path,
                   npy_bool* const inside_flag)
{
    int *yflag0;
    int *subpath_flag;
    int yflag1;
    double vtx0, vty0, vtx1, vty1;
    double tx, ty;
    double sx, sy;
    double x, y;
    size_t i;
    int all_done;
    const char *const points = (const char * const)points_;

    yflag0 = (int *)malloc(n * sizeof(int));
    subpath_flag = (int *)malloc(n * sizeof(int));

    path.rewind(0);

    for (i = 0; i < n; ++i) {
        inside_flag[i] = 0;
    }

    unsigned code = 0;
    do
    {
        if (code != agg::path_cmd_move_to)
        {
            code = path.vertex(&x, &y);
            if (code == agg::path_cmd_stop ||
                (code & agg::path_cmd_end_poly) == agg::path_cmd_end_poly) {
                continue;
            }
        }

        sx = vtx0 = vtx1 = x;
        sy = vty0 = vty1 = y;

        for (i = 0; i < n; ++i) {
            ty = *(double *)(points + s0 * i + s1);

            // get test bit for above/below X axis
            yflag0[i] = (vty0 >= ty);

            subpath_flag[i] = 0;
        }

        do
        {
            code = path.vertex(&x, &y);

            // The following cases denote the beginning on a new subpath
            if (code == agg::path_cmd_stop ||
                (code & agg::path_cmd_end_poly) == agg::path_cmd_end_poly)
            {
                x = sx;
                y = sy;
            }
            else if (code == agg::path_cmd_move_to)
            {
                break;
            }

            for (i = 0; i < n; ++i) {
                tx = *(double *)(points + s0 * i);
                ty = *(double *)(points + s0 * i + s1);

                yflag1 = (vty1 >= ty);
                // Check if endpoints straddle (are on opposite sides) of
                // X axis (i.e. the Y's differ); if so, +X ray could
                // intersect this edge.  The old test also checked whether
                // the endpoints are both to the right or to the left of
                // the test point.  However, given the faster intersection
                // point computation used below, this test was found to be
                // a break-even proposition for most polygons and a loser
                // for triangles (where 50% or more of the edges which
                // survive this test will cross quadrants and so have to
                // have the X intersection computed anyway).  I credit
                // Joseph Samosky with inspiring me to try dropping the
                // "both left or both right" part of my code.
                if (yflag0[i] != yflag1) {
                    // Check intersection of pgon segment with +X ray.
                    // Note if >= point's X; if so, the ray hits it.  The
                    // division operation is avoided for the ">=" test by
                    // checking the sign of the first vertex wrto the test
                    // point; idea inspired by Joseph Samosky's and Mark
                    // Haigh-Hutchinson's different polygon inclusion
                    // tests.
                    if (((vty1 - ty) * (vtx0 - vtx1) >=
                         (vtx1 - tx) * (vty0 - vty1)) == yflag1) {
                        subpath_flag[i] ^= 1;
                    }
                }

                // Move to the next pair of vertices, retaining info as
                // possible.
                yflag0[i] = yflag1;
            }

            vtx0 = vtx1;
            vty0 = vty1;

            vtx1 = x;
            vty1 = y;
        }
        while (code != agg::path_cmd_stop &&
               (code & agg::path_cmd_end_poly) != agg::path_cmd_end_poly);

        all_done = 1;
        for (i = 0; i < n; ++i) {
            tx = *(double *)(points + s0 * i);
            ty = *(double *)(points + s0 * i + s1);

            yflag1 = (vty1 >= ty);
            if (yflag0[i] != yflag1) {
                if (((vty1 - ty) * (vtx0 - vtx1) >=
                     (vtx1 - tx) * (vty0 - vty1)) == yflag1) {
                    subpath_flag[i] ^= 1;
                }
            }
            inside_flag[i] |= subpath_flag[i];
            if (inside_flag[i] == 0) {
                all_done = 0;
            }
        }

        if (all_done) {
            goto exit;
        }
    }
    while (code != agg::path_cmd_stop);

 exit:

    free(yflag0);
    free(subpath_flag);
}

inline void
points_in_path(const void* const points, const size_t s0,
               const size_t s1, const size_t n,
               const double r, PathIterator& path,
               const agg::trans_affine& trans,
               npy_bool* result)
{
    typedef agg::conv_transform<PathIterator> transformed_path_t;
    typedef PathNanRemover<transformed_path_t> no_nans_t;
    typedef agg::conv_curve<no_nans_t> curve_t;
    typedef agg::conv_contour<curve_t> contour_t;

    size_t i;
    for (i = 0; i < n; ++i) {
        result[i] = 0;
    }

    if (path.total_vertices() < 3)
    {
        return;
    }

    transformed_path_t trans_path(path, trans);
    no_nans_t no_nans_path(trans_path, true, path.has_curves());
    curve_t curved_path(no_nans_path);
    contour_t contoured_path(curved_path);
    contoured_path.width(r);
    point_in_path_impl(points, s0, s1, n, contoured_path, result);
}

inline bool
point_in_path(const double x, const double y, const double r,
              PathIterator& path, const agg::trans_affine& trans)
{
    double points[2];
    npy_bool result;

    points[0] = x;
    points[1] = y;

    points_in_path(points, 0, sizeof(double), 1, r, path, trans, &result);
    return result;
}

inline void
points_on_path(const void* const points, const size_t s0,
               const size_t s1, const size_t n,
               const double r, PathIterator& path,
               const agg::trans_affine& trans,
               npy_bool* result)
{
    typedef agg::conv_transform<PathIterator> transformed_path_t;
    typedef PathNanRemover<transformed_path_t> no_nans_t;
    typedef agg::conv_curve<no_nans_t> curve_t;
    typedef agg::conv_stroke<curve_t> stroke_t;

    transformed_path_t trans_path(path, trans);
    no_nans_t nan_removed_path(trans_path, true, path.has_curves());
    curve_t curved_path(nan_removed_path);
    stroke_t stroked_path(curved_path);
    stroked_path.width(r * 2.0);
    point_in_path_impl(points, s0, s1, n, stroked_path, result);
}

inline bool
point_on_path(const double x, const double y, const double r,
              PathIterator& path, const agg::trans_affine& trans)
{
    double points[2];
    npy_bool result;

    points[0] = x;
    points[1] = y;

    points_on_path(points, 0, sizeof(double), 1, r, path, trans, &result);
    return result;
}

PyObject *_point_in_path(PyObject *self, PyObject *_args)
{
    double x, y, r;
    PyObject *_path, *_trans;
    if (!PyArg_ParseTuple(_args, "dddOO", &x, &y, &r, &_path, &_trans)) {
        return NULL;
    }
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);
    
    if (::point_in_path(x, y, r, path, trans)) {
        return PyInt_FromLong(1);
    }
    return PyInt_FromLong(0);
}

PyObject *_points_in_path(PyObject *self, PyObject *_args)
{
    double r;
    PyObject *_points_array, *_path, *_trans;
    if (!PyArg_ParseTuple(_args, "OdOO", &_points_array, &r, &_path, &_trans)) {
        return NULL;
    }

    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);

    PyArrayObject* points_array;
    points_array = (PyArrayObject*)PyArray_FromObject(_points_array, PyArray_DOUBLE, 2, 2);
    if (points_array == NULL || PyArray_DIM(points_array, 1) != 2) {
        PyErr_SetString(PyExc_TypeError, 
                "Argument 0 to points_in_path must be an Nx2 numpy array");
        Py_XDECREF(points_array);
        return NULL;
    }
        
    npy_intp n = PyArray_DIM(points_array, 0);
    PyObject* result = PyArray_ZEROS(1, &n, PyArray_BOOL, 0);
    if (result == NULL) {
        PyErr_SetString(PyExc_MemoryError,
                "Could not allocate memory for result");
        Py_DECREF(points_array);
        return NULL;
    }

    ::points_in_path(PyArray_DATA(points_array),
                     PyArray_STRIDE(points_array, 0),
                     PyArray_STRIDE(points_array, 1),
                     n, r, path, trans,
                     (npy_bool *)PyArray_DATA(result));
    Py_DECREF(points_array);

    return result;
}

PyObject *_point_on_path(PyObject *self, PyObject *_args)
{
    double x, y, r;
    PyObject *_path, *_trans;
    if (!PyArg_ParseTuple(_args, "dddOO", &x, &y, &r, &_path, &_trans)) {
        return NULL;
    }
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);

    if (::point_on_path(x, y, r, path, trans))
    {
        return PyInt_FromLong(1);
    }
    return PyInt_FromLong(0);
}

void
update_limits(double x, double y,
              double* x0, double* y0, double* x1, double* y1,
              double* xm, double* ym)
{
    if (x < *x0) *x0 = x;
    if (y < *y0) *y0 = y;
    if (x > *x1) *x1 = x;
    if (y > *y1) *y1 = y;
    /* xm and ym are the minimum positive values in the data, used
       by log scaling */
    if (x > 0.0 && x < *xm) *xm = x;
    if (y > 0.0 && y < *ym) *ym = y;
}

void
get_path_extents(PathIterator& path, const agg::trans_affine& trans,
                 double* x0, double* y0, double* x1, double* y1,
                 double* xm, double* ym)
{
    typedef agg::conv_transform<PathIterator> transformed_path_t;
    typedef PathNanRemover<transformed_path_t> nan_removed_t;
    typedef agg::conv_curve<nan_removed_t> curve_t;
    double x, y;
    unsigned code;

    transformed_path_t tpath(path, trans);
    nan_removed_t nan_removed(tpath, true, path.has_curves());

    nan_removed.rewind(0);

    while ((code = nan_removed.vertex(&x, &y)) != agg::path_cmd_stop)
    {
        if ((code & agg::path_cmd_end_poly) == agg::path_cmd_end_poly)
        {
            continue;
        }
        update_limits(x, y, x0, y0, x1, y1, xm, ym);
    }
}

PyObject *_get_path_extents(PyObject *self, PyObject *_args)
{
    PyObject *_path, *_trans;
    if (!PyArg_ParseTuple(_args, "OO", &_path, &_trans)) {
        return NULL;
    }
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);

    npy_intp extent_dims[] = { 2, 2, 0 };
    double* extents_data = NULL;
    double xm, ym;
    PyArrayObject* extents = NULL;

    extents = (PyArrayObject*)PyArray_SimpleNew(2, extent_dims, PyArray_DOUBLE);
    if (extents == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
        return NULL;
    }

    extents_data = (double*)PyArray_DATA(extents);

    extents_data[0] = std::numeric_limits<double>::infinity();
    extents_data[1] = std::numeric_limits<double>::infinity();
    extents_data[2] = -std::numeric_limits<double>::infinity();
    extents_data[3] = -std::numeric_limits<double>::infinity();
    /* xm and ym are the minimum positive values in the data, used by log scaling */
    xm = std::numeric_limits<double>::infinity();
    ym = std::numeric_limits<double>::infinity();

    // shmuller: never throws, so remove try-catch block
    ::get_path_extents(path, trans, &extents_data[0], &extents_data[1],
                       &extents_data[2], &extents_data[3], &xm, &ym);
    
    return (PyObject*)extents;
}

PyObject *_update_path_extents(PyObject *self, PyObject *_args)
{
    PyObject *_path, *_trans, *_bbox, *_minpos, *_ignore;
    if (!PyArg_ParseTuple(_args, "OOOOO", 
                &_path, &_trans, &_bbox, &_minpos, &_ignore)) {
        return NULL;
    }
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);
    
    double x0, y0, x1, y1;
    if (!py_convert_bbox(_bbox, x0, y0, x1, y1))
    {
        PyErr_SetString(PyExc_ValueError,
                "Must pass Bbox object as arg 3 of update_path_extents");
        return NULL;
    }
    bool ignore = PyObject_IsTrue(_ignore) != 0;

    double xm, ym;
    PyArrayObject* input_minpos = NULL;
    input_minpos = (PyArrayObject*)PyArray_FromObject(_minpos, PyArray_DOUBLE, 1, 1);
    if (!input_minpos || PyArray_DIM(input_minpos, 0) != 2)
    {
        PyErr_SetString(PyExc_TypeError,
                "Argument 4 to update_path_extents must be a length-2 numpy array.");
        Py_XDECREF(input_minpos);
        return NULL;
    }
    xm = *(double*)PyArray_GETPTR1(input_minpos, 0);
    ym = *(double*)PyArray_GETPTR1(input_minpos, 1);
    Py_DECREF(input_minpos);

    npy_intp extent_dims[] = { 2, 2, 0 };
    double* extents_data = NULL;
    npy_intp minpos_dims[] = { 2, 0 };
    double* minpos_data = NULL;
    PyArrayObject* extents = NULL;
    PyArrayObject* minpos = NULL;
    bool changed = false;

        extents = (PyArrayObject*)PyArray_SimpleNew
                  (2, extent_dims, PyArray_DOUBLE);
        if (extents == NULL)
        {   
            PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
            return NULL;
        }
        minpos = (PyArrayObject*)PyArray_SimpleNew
                 (1, minpos_dims, PyArray_DOUBLE);
        if (minpos == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
            Py_DECREF(extents);
            return NULL;
        }

        extents_data = (double*)PyArray_DATA(extents);
        minpos_data = (double*)PyArray_DATA(minpos);

        if (ignore)
        {
            extents_data[0] = std::numeric_limits<double>::infinity();
            extents_data[1] = std::numeric_limits<double>::infinity();
            extents_data[2] = -std::numeric_limits<double>::infinity();
            extents_data[3] = -std::numeric_limits<double>::infinity();
            minpos_data[0] = std::numeric_limits<double>::infinity();
            minpos_data[1] = std::numeric_limits<double>::infinity();
        }
        else
        {
            if (x0 > x1)
            {
                extents_data[0] = std::numeric_limits<double>::infinity();
                extents_data[2] = -std::numeric_limits<double>::infinity();
            }
            else
            {
                extents_data[0] = x0;
                extents_data[2] = x1;
            }
            if (y0 > y1)
            {
                extents_data[1] = std::numeric_limits<double>::infinity();
                extents_data[3] = -std::numeric_limits<double>::infinity();
            }
            else
            {
                extents_data[1] = y0;
                extents_data[3] = y1;
            }
            minpos_data[0] = xm;
            minpos_data[1] = ym;
        }

        ::get_path_extents(path, trans, &extents_data[0], &extents_data[1],
                           &extents_data[2], &extents_data[3], &minpos_data[0],
                           &minpos_data[1]);

        changed = (extents_data[0] != x0 ||
                   extents_data[1] != y0 ||
                   extents_data[2] != x1 ||
                   extents_data[3] != y1 ||
                   minpos_data[0]  != xm ||
                   minpos_data[1]  != ym);
    
    return Py_BuildValue("(NNi)", extents, minpos, changed ? 1 : 0);
}

PyObject *_get_path_collection_extents(PyObject *self, PyObject *_args)
{
    PyObject *_master_transform, *_paths, *_transforms, *_offsets, *_offset_trans;
    if (!PyArg_ParseTuple(_args, "OOOOO", 
            &_master_transform, &_paths, &_transforms, &_offsets, &_offset_trans)) {
        return NULL;
    }
    agg::trans_affine master_transform = py_to_agg_transformation_matrix
        (_master_transform, false);
    agg::trans_affine offset_trans     = py_to_agg_transformation_matrix
        (_offset_trans, false);
    
    double x0, y0, x1, y1, xm, ym;

        PyObject* __paths = PySequence_Fast(_paths, "paths must be a sequence");
        if (__paths == NULL) 
        {
            return NULL;
        }

        PyObject* __transforms = PySequence_Fast(_transforms, "transforms must be a sequence");
        if (__transforms == NULL) 
        {
            Py_DECREF(__paths);
            return NULL;
        }
        
        PyArrayObject* offsets = (PyArrayObject*)PyArray_FromObject(
            _offsets, PyArray_DOUBLE, 0, 2);
        if (!offsets ||
            (PyArray_NDIM(offsets) == 2 && PyArray_DIM(offsets, 1) != 2) ||
            (PyArray_NDIM(offsets) == 1 && PyArray_DIM(offsets, 0) != 0))
        {
            PyErr_SetString(PyExc_ValueError, "Offsets array must be Nx2");
            Py_DECREF(__paths);
            Py_DECREF(__transforms);
            Py_XDECREF(offsets);
            return NULL;
        }

        PyObject** paths_arr = PySequence_Fast_ITEMS(__paths);
        PyObject** transforms_arr = PySequence_Fast_ITEMS(__transforms);

        size_t Npaths      = PySequence_Fast_GET_SIZE(__paths);
        size_t Noffsets    = PyArray_DIM(offsets, 0);
        size_t N           = std::max(Npaths, Noffsets);
        size_t Ntransforms = std::min<size_t>(PySequence_Fast_GET_SIZE(__transforms), N);
        size_t i;

        // Convert all of the transforms up front
        typedef std::vector<agg::trans_affine> transforms_t;
        transforms_t transforms;
        transforms.reserve(Ntransforms);
        for (i = 0; i < Ntransforms; ++i)
        {
            agg::trans_affine trans = py_to_agg_transformation_matrix
                (transforms_arr[i], false);
            trans *= master_transform;
            transforms.push_back(trans);
        }

        // The offset each of those and collect the mins/maxs
        x0 = std::numeric_limits<double>::infinity();
        y0 = std::numeric_limits<double>::infinity();
        x1 = -std::numeric_limits<double>::infinity();
        y1 = -std::numeric_limits<double>::infinity();
        xm = std::numeric_limits<double>::infinity();
        ym = std::numeric_limits<double>::infinity();
        agg::trans_affine trans;

        if (transforms.size() <= 1 && Npaths == 1)
        {
            PathIterator path(paths_arr[0]);
            if (Ntransforms)
            {
                trans = transforms[0];
            }
            else
            {
                trans = master_transform;
            }

            double bx0 = std::numeric_limits<double>::infinity();
            double by0 = std::numeric_limits<double>::infinity();
            double bx1 = -std::numeric_limits<double>::infinity();
            double by1 = -std::numeric_limits<double>::infinity();
            double bxm = std::numeric_limits<double>::infinity();
            double bym = std::numeric_limits<double>::infinity();

            ::get_path_extents(path, trans, &bx0, &by0, &bx1, &by1, &bxm, &bym);

            for (i = 0; i < Noffsets; ++i)
            {
                double xo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 0);
                double yo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 1);
                offset_trans.transform(&xo, &yo);
                update_limits(xo + bx0, yo + by0, &x0, &y0, &x1, &y1, &xm, &ym);
                update_limits(xo + bx1, yo + by1, &x0, &y0, &x1, &y1, &xm, &ym);
            }
        } else {
            for (i = 0; i < N; ++i)
            {
                PathIterator path(paths_arr[i % Npaths]);
                if (Ntransforms)
                {
                    trans = transforms[i % Ntransforms];
                }
                else
                {
                    trans = master_transform;
                }

                if (Noffsets)
                {
                    double xo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 0);
                    double yo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 1);
                    offset_trans.transform(&xo, &yo);
                    trans *= agg::trans_affine_translation(xo, yo);
                }

                ::get_path_extents(path, trans, &x0, &y0, &x1, &y1, &xm, &ym);
            }
        }
    
    Py_DECREF(__paths);
    Py_DECREF(__transforms);
    Py_DECREF(offsets);

    return Py_BuildValue("(dddd)", x0, y0, x1, y1);
}

PyObject *_point_in_path_collection(PyObject *self, PyObject *_args)
{
    double x, y, radius;
    PyObject *_master_transform, *_paths, *_transforms, *_offsets, *_offset_trans, *_filled;
    const char *_offset_position;
    if (!PyArg_ParseTuple(_args, "dddOOOOOOs", &x, &y, &radius,
             &_master_transform, &_paths, &_transforms, &_offsets, &_offset_trans, &_filled, 
             &_offset_position)) {
        return NULL;
    }
    agg::trans_affine master_transform = py_to_agg_transformation_matrix
        (_master_transform, false);
    agg::trans_affine offset_trans     = py_to_agg_transformation_matrix
        (_offset_trans, false);

    PyObject* __paths = PySequence_Fast(_paths, "paths must be a sequence");
    if (__paths == NULL) 
    {
        return NULL;
    }
    PyObject* __transforms = PySequence_Fast(_transforms, "transforms must be a sequence");
    if (__transforms == NULL) 
    {
        Py_DECREF(__paths);
        return NULL;
    }

    PyArrayObject* offsets = (PyArrayObject*)PyArray_FromObject(
        _offsets, PyArray_DOUBLE, 0, 2);
    if (!offsets ||
        (PyArray_NDIM(offsets) == 2 && PyArray_DIM(offsets, 1) != 2) ||
        (PyArray_NDIM(offsets) == 1 && PyArray_DIM(offsets, 0) != 0))
    {
        PyErr_SetString(PyExc_ValueError, "Offsets array must be Nx2");
        Py_DECREF(__paths);
        Py_DECREF(__transforms);
        Py_XDECREF(offsets);
        return NULL;
    }

    bool filled = PyObject_IsTrue(_filled) != 0;
    std::string offset_position(_offset_position);

    bool data_offsets = (offset_position == "data");

    PyObject** paths_arr = PySequence_Fast_ITEMS(__paths);
    PyObject** transforms_arr = PySequence_Fast_ITEMS(__transforms);

    size_t Npaths      = PySequence_Fast_GET_SIZE(__paths);
    size_t Noffsets    = PyArray_DIM(offsets, 0);
    size_t N           = std::max(Npaths, Noffsets);
    size_t Ntransforms = std::min<size_t>(PySequence_Fast_GET_SIZE(__transforms), N);
    size_t i;

    // Convert all of the transforms up front
    typedef std::vector<agg::trans_affine> transforms_t;
    transforms_t transforms;
    transforms.reserve(Ntransforms);
    for (i = 0; i < Ntransforms; ++i)
    {
        agg::trans_affine trans = py_to_agg_transformation_matrix
                                  (transforms_arr[i], false);
        trans *= master_transform;
        transforms.push_back(trans);
    }

    agg::trans_affine trans;

    PyObject* result = PyList_New(0);
    if (result == NULL) goto Done;
    
    for (i = 0; i < N; ++i)
    {
        PathIterator path(paths_arr[i % Npaths]);

        if (Ntransforms)
        {
            trans = transforms[i % Ntransforms];
        }
        else
        {
            trans = master_transform;
        }

        if (Noffsets)
        {
            double xo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 0);
            double yo = *(double*)PyArray_GETPTR2(offsets, i % Noffsets, 1);
            offset_trans.transform(&xo, &yo);
            if (data_offsets) {
                trans = agg::trans_affine_translation(xo, yo) * trans;
            } else {
                trans *= agg::trans_affine_translation(xo, yo);
            }
        }

        if (filled)
        {
            if (::point_in_path(x, y, radius, path, trans)) {
                PyObject* _i = PyInt_FromLong((int)i);
                int status = PyList_Append(result, _i);
                Py_DECREF(_i);
                if (status == -1) goto Fail;
            }
        }
        else
        {
            if (::point_on_path(x, y, radius, path, trans)) {
                PyObject* _i = PyInt_FromLong((int)i);
                int status = PyList_Append(result, _i);
                Py_DECREF(_i);
                if (status == -1) goto Fail;
            }
        }
    }

Fail:
    Py_DECREF(result);
    result = NULL;
Done:
    Py_DECREF(__paths);
    Py_DECREF(__transforms);
    Py_DECREF(offsets);

    return result;
}

bool
path_in_path(PathIterator& a, const agg::trans_affine& atrans,
             PathIterator& b, const agg::trans_affine& btrans)
{
    typedef agg::conv_transform<PathIterator> transformed_path_t;
    typedef PathNanRemover<transformed_path_t> no_nans_t;
    typedef agg::conv_curve<no_nans_t> curve_t;

    if (a.total_vertices() < 3)
        return false;

    transformed_path_t b_path_trans(b, btrans);
    no_nans_t b_no_nans(b_path_trans, true, b.has_curves());
    curve_t b_curved(b_no_nans);

    double x, y;
    b_curved.rewind(0);
    while (b_curved.vertex(&x, &y) != agg::path_cmd_stop)
    {
        if (!::point_in_path(x, y, 0.0, a, atrans))
            return false;
    }

    return true;
}

PyObject *_path_in_path(PyObject *self, PyObject *_args)
{
    PyObject *_a, *_atrans, *_b, *_btrans;
    if (!PyArg_ParseTuple(_args, "OOOO", &_a, &_atrans, &_b, &_btrans)) {
        return NULL;
    }
    PathIterator a(_a);
    agg::trans_affine atrans = py_to_agg_transformation_matrix(_atrans, false);

    PathIterator b(_b);
    agg::trans_affine btrans = py_to_agg_transformation_matrix(_btrans, false);

    return PyInt_FromLong(::path_in_path(a, atrans, b, btrans));
}

/** The clip_path_to_rect code here is a clean-room implementation of
    the Sutherland-Hodgman clipping algorithm described here:

  http://en.wikipedia.org/wiki/Sutherland-Hodgman_clipping_algorithm
*/

typedef std::vector<XY> Polygon;

namespace clip_to_rect_filters
{
    /* There are four different passes needed to create/remove
       vertices (one for each side of the rectangle).  The differences
       between those passes are encapsulated in these functor classes.
    */
    struct bisectx
    {
        double m_x;

        bisectx(double x) : m_x(x) {}

        inline void
        bisect(double sx, double sy, double px, double py, double* bx,
               double* by) const
        {
            *bx = m_x;
            double dx = px - sx;
            double dy = py - sy;
            *by = sy + dy * ((m_x - sx) / dx);
        }
    };

    struct xlt : public bisectx
    {
        xlt(double x) : bisectx(x) {}

        inline bool
        is_inside(double x, double y) const
        {
            return x <= m_x;
        }
    };

    struct xgt : public bisectx
    {
        xgt(double x) : bisectx(x) {}

        inline bool
        is_inside(double x, double y) const
        {
            return x >= m_x;
        }
    };

    struct bisecty
    {
        double m_y;

        bisecty(double y) : m_y(y) {}

        inline void
        bisect(double sx, double sy, double px, double py, double* bx,
               double* by) const
        {
            *by = m_y;
            double dx = px - sx;
            double dy = py - sy;
            *bx = sx + dx * ((m_y - sy) / dy);
        }
    };

    struct ylt : public bisecty
    {
        ylt(double y) : bisecty(y) {}

        inline bool
        is_inside(double x, double y) const
        {
            return y <= m_y;
        }
    };

    struct ygt : public bisecty
    {
        ygt(double y) : bisecty(y) {}

        inline bool
        is_inside(double x, double y) const
        {
            return y >= m_y;
        }
    };
}

template<class Filter>
inline void
clip_to_rect_one_step(const Polygon& polygon, Polygon& result, const Filter& filter)
{
    double sx, sy, px, py, bx, by;
    bool sinside, pinside;
    result.clear();

    if (polygon.size() == 0)
    {
        return;
    }

    sx = polygon.back().x;
    sy = polygon.back().y;
    for (Polygon::const_iterator i = polygon.begin(); i != polygon.end(); ++i)
    {
        px = i->x;
        py = i->y;

        sinside = filter.is_inside(sx, sy);
        pinside = filter.is_inside(px, py);

        if (sinside ^ pinside)
        {
            filter.bisect(sx, sy, px, py, &bx, &by);
            result.push_back(XY(bx, by));
        }

        if (pinside)
        {
            result.push_back(XY(px, py));
        }

        sx = px;
        sy = py;
    }
}

template<class Path>
void
clip_to_rect(Path& path,
             double x0, double y0, double x1, double y1,
             bool inside, std::vector<Polygon>& results)
{
    double xmin, ymin, xmax, ymax;
    if (x0 < x1)
    {
        xmin = x0;
        xmax = x1;
    }
    else
    {
        xmin = x1;
        xmax = x0;
    }

    if (y0 < y1)
    {
        ymin = y0;
        ymax = y1;
    }
    else
    {
        ymin = y1;
        ymax = y0;
    }

    if (!inside)
    {
        std::swap(xmin, xmax);
        std::swap(ymin, ymax);
    }

    Polygon polygon1, polygon2;
    double x = 0, y = 0;
    unsigned code = 0;
    path.rewind(0);

    do
    {
        // Grab the next subpath and store it in polygon1
        polygon1.clear();
        do
        {
            if (code == agg::path_cmd_move_to)
            {
                polygon1.push_back(XY(x, y));
            }

            code = path.vertex(&x, &y);

            if (code == agg::path_cmd_stop)
            {
                break;
            }

            if (code != agg::path_cmd_move_to)
            {
                polygon1.push_back(XY(x, y));
            }
        }
        while ((code & agg::path_cmd_end_poly) != agg::path_cmd_end_poly);

        // The result of each step is fed into the next (note the
        // swapping of polygon1 and polygon2 at each step).
        clip_to_rect_one_step(polygon1, polygon2, clip_to_rect_filters::xlt(xmax));
        clip_to_rect_one_step(polygon2, polygon1, clip_to_rect_filters::xgt(xmin));
        clip_to_rect_one_step(polygon1, polygon2, clip_to_rect_filters::ylt(ymax));
        clip_to_rect_one_step(polygon2, polygon1, clip_to_rect_filters::ygt(ymin));

        // Empty polygons aren't very useful, so skip them
        if (polygon1.size())
        {
            results.push_back(polygon1);
        }
    }
    while (code != agg::path_cmd_stop);
}

PyObject *_clip_path_to_rect(PyObject *self, PyObject *_args)
{
    PyObject *_path, *_bbox, *_inside;
    if (!PyArg_ParseTuple(_args, "OOO", &_path, &_bbox, &_inside)) {
        return NULL;
    }
    PathIterator path(_path);
    bool inside = PyObject_IsTrue(_inside) != 0;

    double x0, y0, x1, y1;
    if (!py_convert_bbox(_bbox, x0, y0, x1, y1))
    {
        PyErr_SetString(PyExc_TypeError,
                "Argument 2 to clip_to_rect must be a Bbox object.");
        return NULL;
    }

    std::vector<Polygon> results;
    typedef agg::conv_curve<PathIterator> curve_t;
    curve_t curve(path);

    ::clip_to_rect(curve, x0, y0, x1, y1, inside, results);

    npy_intp dims[2];
    dims[1] = 2;
    PyObject* py_results = PyList_New(results.size());
    if (!py_results)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error creating results list");
        return NULL;
    }
        
        for (std::vector<Polygon>::const_iterator p = results.begin(); p != results.end(); ++p)
        {
            size_t size = p->size();
            dims[0] = (npy_intp)size + 1;
            PyArrayObject* pyarray = (PyArrayObject*)PyArray_SimpleNew(2, dims, PyArray_DOUBLE);
            if (pyarray == NULL)
            {
                PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
                Py_DECREF(py_results);
                return NULL;
            }

            double *data = (double *) PyArray_DATA(pyarray);

            for (size_t i = 0; i < size; ++i)
            {
                data[2*i]   = (*p)[i].x;
                data[2*i+1] = (*p)[i].y;
            }
            data[2*size]   = (*p)[0].x;
            data[2*size+1] = (*p)[0].y;

            // cannot fail
            PyList_SET_ITEM(py_results, p - results.begin(), (PyObject *)pyarray);
        }
    
    return py_results;
}

PyObject *_affine_transform(PyObject *self, PyObject *_args)
{
    PyObject *_vertices, *_transform;
    if (!PyArg_ParseTuple(_args, "OO", &_vertices, &_transform)) {
        return NULL;
    }

    PyArrayObject* vertices = NULL;
    PyArrayObject* transform = NULL;
    PyArrayObject* result = NULL;

    int nd;
    size_t n;
    npy_intp dims[2];
        
        vertices = (PyArrayObject*)PyArray_FromObject(_vertices, PyArray_DOUBLE, 1, 2);
        if (!vertices ||
            (PyArray_NDIM(vertices) == 2 && PyArray_DIM(vertices, 0) != 0 &&
             PyArray_DIM(vertices, 1) != 2) ||
            (PyArray_NDIM(vertices) == 1 &&
             PyArray_DIM(vertices, 0) != 2 && PyArray_DIM(vertices, 0) != 0))
        {
            PyErr_SetString(PyExc_ValueError, "Invalid vertices array.");
            goto Fail;
        }

        transform = (PyArrayObject*)PyArray_FromObject(_transform, PyArray_DOUBLE, 2, 2);
        if (!transform ||
            PyArray_DIM(transform, 0) != 3 ||
            PyArray_DIM(transform, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "Invalid transform.");
            goto Fail;
        }

        double a, b, c, d, e, f;
        {
            size_t stride0 = PyArray_STRIDE(transform, 0);
            size_t stride1 = PyArray_STRIDE(transform, 1);
            char* row0 = PyArray_BYTES(transform);
            char* row1 = row0 + stride0;

            a = *(double*)(row0);
            row0 += stride1;
            c = *(double*)(row0);
            row0 += stride1;
            e = *(double*)(row0);

            b = *(double*)(row1);
            row1 += stride1;
            d = *(double*)(row1);
            row1 += stride1;
            f = *(double*)(row1);
        }

        // PyPy's PyArray_DIMS() is inefficient, avoid where possible
        nd = PyArray_NDIM(vertices);
        n = dims[0] = PyArray_DIM(vertices, 0);
        if (nd == 2) dims[1] = PyArray_DIM(vertices, 1);
        result = (PyArrayObject*)PyArray_SimpleNew(nd, dims, PyArray_DOUBLE);
                 
        if (result == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for path");
            goto Fail;
        }
        if (nd == 2)
        {
            char* vertex_in = PyArray_BYTES(vertices);
            double* vertex_out = (double*)PyArray_DATA(result);
            size_t stride0 = PyArray_STRIDE(vertices, 0);
            size_t stride1 = PyArray_STRIDE(vertices, 1);
            double x, y;
            volatile double t0, t1, t;

            for (size_t i = 0; i < n; ++i)
            {
                x = *(double*)(vertex_in);
                y = *(double*)(vertex_in + stride1);

                t0 = a * x;
                t1 = c * y;
                t = t0 + t1 + e;
                *(vertex_out++) = t;

                t0 = b * x;
                t1 = d * y;
                t = t0 + t1 + f;
                *(vertex_out++) = t;

                vertex_in += stride0;
            }
        }
        else if (n != 0)
        {
            char* vertex_in = PyArray_BYTES(vertices);
            double* vertex_out = (double*)PyArray_DATA(result);
            size_t stride0 = PyArray_STRIDE(vertices, 0);
            double x;
            double y;
            x = *(double*)(vertex_in);
            y = *(double*)(vertex_in + stride0);
            *vertex_out++ = a * x + c * y + e;
            *vertex_out++ = b * x + d * y + f;
        }

Fail:
    Py_XDECREF(vertices);
    Py_XDECREF(transform);

    return (PyObject*)result;
}

PyObject *_count_bboxes_overlapping_bbox(PyObject *self, PyObject *_args)
{
    PyObject *_bbox, *_bboxes;
    if (!PyArg_ParseTuple(_args, "OO", &_bbox, &_bboxes)) {
        return NULL;
    }

    double ax0, ay0, ax1, ay1;
    if (!py_convert_bbox(_bbox, ax0, ay0, ax1, ay1))
    {
        PyErr_SetString(PyExc_ValueError,
                "First argument to count_bboxes_overlapping_bbox must be a Bbox object.");
        return NULL;
    }

    PyObject* __bboxes = PySequence_Fast(_bboxes, "bboxes must be a sequence");
    if (__bboxes == NULL)
    {
        return NULL;
    }
    PyObject** bboxes = PySequence_Fast_ITEMS(__bboxes);
    size_t num_bboxes = PySequence_Fast_GET_SIZE(__bboxes);

    double bx0, by0, bx1, by1;
    long count = 0;

        if (ax1 < ax0)
        {
            std::swap(ax0, ax1);
        }
        if (ay1 < ay0)
        {
            std::swap(ay0, ay1);
        }

        for (size_t i = 0; i < num_bboxes; ++i)
        {
            if (py_convert_bbox(bboxes[i], bx0, by0, bx1, by1))
            {
                if (bx1 < bx0)
                {
                    std::swap(bx0, bx1);
                }
                if (by1 < by0)
                {
                    std::swap(by0, by1);
                }
                if (!((bx1 <= ax0) ||
                      (by1 <= ay0) ||
                      (bx0 >= ax1) ||
                      (by0 >= ay1)))
                {
                    ++count;
                }
            }
            else
            {
                PyErr_SetString(PyExc_ValueError, "Non-bbox object in bboxes list");
                Py_DECREF(__bboxes);
                return NULL;
            }
        }
    
    Py_DECREF(__bboxes);

    return PyInt_FromLong(count);
}

inline bool
segments_intersect(const double& x1, const double& y1,
                   const double& x2, const double& y2,
                   const double& x3, const double& y3,
                   const double& x4, const double& y4)
{
    double den = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
    if (den == 0.0)
    {
        return false;
    }

    double n1 = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
    double n2 = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));

    double u1 = n1 / den;
    double u2 = n2 / den;

    return (u1 >= 0.0 && u1 <= 1.0 &&
            u2 >= 0.0 && u2 <= 1.0);
}

bool
path_intersects_path(PathIterator& p1, PathIterator& p2)
{
    typedef PathNanRemover<PathIterator> no_nans_t;
    typedef agg::conv_curve<no_nans_t> curve_t;

    if (p1.total_vertices() < 2 || p2.total_vertices() < 2)
    {
        return false;
    }

    no_nans_t n1(p1, true, p1.has_curves());
    no_nans_t n2(p2, true, p2.has_curves());

    curve_t c1(n1);
    curve_t c2(n2);

    double x11, y11, x12, y12;
    double x21, y21, x22, y22;

    c1.vertex(&x11, &y11);
    while (c1.vertex(&x12, &y12) != agg::path_cmd_stop)
    {
        c2.rewind(0);
        c2.vertex(&x21, &y21);
        while (c2.vertex(&x22, &y22) != agg::path_cmd_stop)
        {
            if (segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22))
            {
                return true;
            }
            x21 = x22;
            y21 = y22;
        }
        x11 = x12;
        y11 = y12;
    }

    return false;
}

PyObject *_path_intersects_path(PyObject *self, PyObject *_args)
{
    PyObject *_p1, *_p2, *_filled = Py_False;
    if (!PyArg_ParseTuple(_args, "OO|O", &_p1, &_p2, &_filled)) {
        return NULL;
    }

    PathIterator p1(_p1);
    PathIterator p2(_p2);
    bool filled = PyObject_IsTrue(_filled) != 0;
    
    bool isect;
    if (!filled)
    {
        isect = ::path_intersects_path(p1, p2);
    }
    else
    {
        isect = ::path_intersects_path(p1, p2)
             || ::path_in_path(p1, agg::trans_affine(), p2, agg::trans_affine())
             || ::path_in_path(p2, agg::trans_affine(), p1, agg::trans_affine());
    }
    return PyInt_FromLong(isect);
}

int
_add_polygon(PyObject *_polygons, const std::vector<double>& polygon)
{
    if (polygon.size() == 0)
    {
        return 0;
    }
    npy_intp polygon_dims[] = { static_cast<npy_intp>(polygon.size() / 2), 2, 0 };
    PyArrayObject* polygon_array = NULL;
    polygon_array = (PyArrayObject*)PyArray_SimpleNew
                    (2, polygon_dims, PyArray_DOUBLE);
    if (!polygon_array)
    {
        PyErr_SetString(PyExc_MemoryError, "Error creating polygon array");
        return -1;
    }
    double* polygon_data = (double*)PyArray_DATA(polygon_array);
    memcpy(polygon_data, &polygon[0], polygon.size() * sizeof(double));
    int status = PyList_Append(_polygons, (PyObject*)polygon_array);
    Py_DECREF(polygon_array);
    return status;
}

PyObject *_convert_path_to_polygons(PyObject *self, PyObject *_args)
{
    typedef agg::conv_transform<PathIterator>  transformed_path_t;
    typedef PathNanRemover<transformed_path_t> nan_removal_t;
    typedef PathClipper<nan_removal_t>         clipped_t;
    typedef PathSimplifier<clipped_t>          simplify_t;
    typedef agg::conv_curve<simplify_t>        curve_t;

    typedef std::vector<double> vertices_t;

    PyObject *_path, *_trans;
    double width, height;
    if (!PyArg_ParseTuple(_args, "OOdd", &_path, &_trans, &width, &height)) {
        return NULL;
    }
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);

    bool do_clip = width != 0.0 && height != 0.0;

    bool simplify = path.should_simplify();

    transformed_path_t tpath(path, trans);
    nan_removal_t      nan_removed(tpath, true, path.has_curves());
    clipped_t          clipped(nan_removed, do_clip, width, height);
    simplify_t         simplified(clipped, simplify, path.simplify_threshold());
    curve_t            curve(simplified);

    PyObject* polygons = PyList_New(0);
    if (polygons == NULL) return NULL;
    vertices_t polygon;
    double x, y;
    unsigned code;

    polygon.reserve(path.total_vertices() * 2);

    while ((code = curve.vertex(&x, &y)) != agg::path_cmd_stop)
    {
        if ((code & agg::path_cmd_end_poly) == agg::path_cmd_end_poly)
        {
            if (polygon.size() >= 2)
            {
                polygon.push_back(polygon[0]);
                polygon.push_back(polygon[1]);
                if (_add_polygon(polygons, polygon) == -1) goto Fail;
            }
            polygon.clear();
        }
        else
        {
            if (code == agg::path_cmd_move_to)
            {
                if (_add_polygon(polygons, polygon) == -1) goto Fail;
                polygon.clear();
            }
            polygon.push_back(x);
            polygon.push_back(y);
        }
    }

    if (_add_polygon(polygons, polygon) == -1) goto Fail;

    return polygons;

Fail:
    Py_DECREF(polygons);
    return NULL;
}

template<class VertexSource>
void
__cleanup_path(VertexSource& source,
               std::vector<double>& vertices,
               std::vector<npy_uint8>& codes)
{
    unsigned code;
    double x, y;
    do
    {
        code = source.vertex(&x, &y);
        vertices.push_back(x);
        vertices.push_back(y);
        codes.push_back((npy_uint8)code);
    }
    while (code != agg::path_cmd_stop);
}

void
_cleanup_path(PathIterator& path, const agg::trans_affine& trans,
              bool remove_nans, bool do_clip,
              const agg::rect_base<double>& rect,
              e_snap_mode snap_mode, double stroke_width,
              bool do_simplify, bool return_curves,
              double sketch_scale, double sketch_length,
              double sketch_randomness,
              std::vector<double>& vertices,
              std::vector<npy_uint8>& codes)
{
    typedef agg::conv_transform<PathIterator>  transformed_path_t;
    typedef PathNanRemover<transformed_path_t> nan_removal_t;
    typedef PathClipper<nan_removal_t>         clipped_t;
    typedef PathSnapper<clipped_t>             snapped_t;
    typedef PathSimplifier<snapped_t>          simplify_t;
    typedef agg::conv_curve<simplify_t>        curve_t;
    typedef Sketch<curve_t>                    sketch_t;

    transformed_path_t tpath(path, trans);
    nan_removal_t      nan_removed(tpath, remove_nans, path.has_curves());
    clipped_t          clipped(nan_removed, do_clip, rect);
    snapped_t          snapped(clipped, snap_mode, path.total_vertices(), stroke_width);
    simplify_t         simplified(snapped, do_simplify, path.simplify_threshold());

    vertices.reserve(path.total_vertices() * 2);
    codes.reserve(path.total_vertices());

    if (return_curves && sketch_scale == 0.0)
    {
        __cleanup_path(simplified, vertices, codes);
    }
    else
    {
        curve_t curve(simplified);
        sketch_t sketch(curve, sketch_scale, sketch_length, sketch_randomness);
        __cleanup_path(sketch, vertices, codes);
    }
}

PyObject *_cleanup_path(PyObject *self, PyObject *_args)
{
    PyObject *_path, *_trans, *_nans, *_clip, *_snap;
    double stroke_width;
    PyObject *_simplify, *_return_curves, *_sketch_params;
    if (!PyArg_ParseTuple(_args, "OOOOOdOOO", &_path, &_trans, &_nans, &_clip, &_snap,
             &stroke_width, &_simplify, &_return_curves, &_sketch_params)) {
        return NULL;
    }

    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);
    bool remove_nans = PyObject_IsTrue(_nans) != 0;
    bool do_clip = _clip != Py_None;

    agg::rect_base<double> clip_rect;
    if (do_clip) 
    {
        double x1, y1, x2, y2;
        if (!PyArg_ParseTuple(_clip, "dddd", &x1, &y1, &x2, &y2)) {
            return NULL;
        }
        clip_rect.init(x1, y1, x2, y2);
    }

    e_snap_mode snap_mode;
    if (_snap == Py_None)
    {
        snap_mode = SNAP_AUTO;
    }
    else if (PyObject_IsTrue(_snap))
    {
        snap_mode = SNAP_TRUE;
    }
    else
    {
        snap_mode = SNAP_FALSE;
    }

    bool simplify;
    if (_simplify == Py_None)
    {
        simplify = path.should_simplify();
    }
    else
    {
        simplify = PyObject_IsTrue(_simplify) != 0;
    }

    bool return_curves = PyObject_IsTrue(_return_curves) != 0;

    double sketch_scale = 0.0;
    double sketch_length = 0.0;
    double sketch_randomness = 0.0;
    if (_sketch_params != Py_None && 
        !PyArg_ParseTuple(_sketch_params, "ddd", 
            &sketch_scale, &sketch_length, &sketch_randomness)) 
    {
        return NULL;
    }

    std::vector<double> vertices;
    std::vector<npy_uint8> codes;

    _cleanup_path(path, trans, remove_nans, do_clip, clip_rect, snap_mode,
                  stroke_width, simplify, return_curves, sketch_scale,
                  sketch_length, sketch_randomness, vertices, codes);

    npy_intp length = codes.size();
    npy_intp dims[] = { length, 2, 0 };

    PyArrayObject* vertices_obj = NULL;
    PyArrayObject* codes_obj = NULL;
        
        vertices_obj = (PyArrayObject*)PyArray_SimpleNew
                       (2, dims, PyArray_DOUBLE);
        if (vertices_obj == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
            return NULL;
        }

        codes_obj = (PyArrayObject*)PyArray_SimpleNew
                    (1, dims, PyArray_UINT8);
        if (codes_obj == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate result array");
            Py_DECREF(vertices_obj);
            return NULL;
        }

        memcpy(PyArray_DATA(vertices_obj), &vertices[0], sizeof(double) * 2 * length);
        memcpy(PyArray_DATA(codes_obj), &codes[0], sizeof(npy_uint8) * length);

    return Py_BuildValue("(NN)", vertices_obj, codes_obj);
}

PyObject *_convert_to_svg(PyObject *self, PyObject *_args)
{
    PyObject *_path, *_trans, *_clip, *_simplify;
    int precision;
    if (!PyArg_ParseTuple(_args, "OOOOi", 
            &_path, &_trans, &_clip, &_simplify, &precision)) {
        return NULL;
    }
    
    PathIterator path(_path);
    agg::trans_affine trans = py_to_agg_transformation_matrix(_trans, false);
    bool do_clip = PyObject_IsTrue(_clip) != 0;

    agg::rect_base<double> clip_rect;
    if (do_clip) 
    {
        double x1, y1, x2, y2;
        if (!PyArg_ParseTuple(_clip, "dddd", &x1, &y1, &x2, &y2)) {
            return NULL;
        }
        clip_rect.init(x1, y1, x2, y2);
    }

    bool simplify;
    if (_simplify == Py_None)
    {
        simplify = path.should_simplify();
    }
    else
    {
        simplify = PyObject_IsTrue(_simplify) != 0;
    }

    char format[64];
    snprintf(format, 64, "%%.%dg %%.%dg", precision, precision);
   
    typedef agg::conv_transform<PathIterator>  transformed_path_t;
    typedef PathNanRemover<transformed_path_t> nan_removal_t;
    typedef PathClipper<nan_removal_t>         clipped_t;
    typedef PathSimplifier<clipped_t>          simplify_t;

    transformed_path_t tpath(path, trans);
    nan_removal_t      nan_removed(tpath, true, path.has_curves());
    clipped_t          clipped(nan_removed, do_clip, clip_rect);
    simplify_t         simplified(clipped, simplify, path.simplify_threshold());

    size_t buffersize = path.total_vertices() * (precision + 5) * 4;
    char* buffer = (char *)malloc(buffersize);
    char* p = buffer;

    const char codes[] = {'M', 'L', 'Q', 'C'};
    const int  waits[] = {  1,   1,   2,   3};

    int wait = 0;
    unsigned code;
    double x = 0, y = 0;
    while ((code = simplified.vertex(&x, &y)) != agg::path_cmd_stop)
    {
        if (wait == 0)
        {
            *p++ = '\n';

            if (code == 0x4f)
            {
                *p++ = 'z';
                *p++ = '\n';
                continue;
            }

            *p++ = codes[code-1];
            wait = waits[code-1];
        }
        else
        {
            *p++ = ' ';
        }

        p += snprintf(p, buffersize - (p - buffer), format, x, y);
        
        --wait;
    }

    #if PY3K
    PyObject* result = PyUnicode_FromStringAndSize(buffer, p - buffer);
    #else
    PyObject* result = PyString_FromStringAndSize(buffer, p - buffer);
    #endif
    free(buffer);

    return result;
}


CXX_WRAPPED(_point_in_path)
CXX_WRAPPED(_points_in_path)
CXX_WRAPPED(_point_on_path)
CXX_WRAPPED(_get_path_extents)
CXX_WRAPPED(_update_path_extents)
CXX_WRAPPED(_get_path_collection_extents)
CXX_WRAPPED(_point_in_path_collection)
CXX_WRAPPED(_path_in_path)
CXX_WRAPPED(_clip_path_to_rect)
CXX_WRAPPED(_affine_transform)
CXX_WRAPPED(_count_bboxes_overlapping_bbox)
CXX_WRAPPED(_path_intersects_path)
CXX_WRAPPED(_convert_path_to_polygons)
CXX_WRAPPED(_cleanup_path)
CXX_WRAPPED(_convert_to_svg)

static PyMethodDef methods[] = {
    {"point_in_path", &CXX_point_in_path, METH_VARARGS,
         "point_in_path(x, y, path, trans)"},
    {"points_in_path", &CXX_points_in_path, METH_VARARGS,
         "points_in_path(points, path, trans)"},
    {"point_on_path", &CXX_point_on_path, METH_VARARGS,
         "point_on_path(x, y, r, path, trans)"},
    {"get_path_extents", &CXX_get_path_extents, METH_VARARGS,
         "get_path_extents(path, trans)"},
    {"update_path_extents", &CXX_update_path_extents, METH_VARARGS,
         "update_path_extents(path, trans, bbox, minpos)"},
    {"get_path_collection_extents", &CXX_get_path_collection_extents, METH_VARARGS,
         "get_path_collection_extents(trans, paths, transforms, offsets, offsetTrans)"},
    {"point_in_path_collection", &CXX_point_in_path_collection, METH_VARARGS,
         "point_in_path_collection(x, y, r, trans, paths, transforms, offsets, offsetTrans, filled)"},
    {"path_in_path", &CXX_path_in_path, METH_VARARGS,
         "path_in_path(a, atrans, b, btrans)"},
    {"clip_path_to_rect", &CXX_clip_path_to_rect, METH_VARARGS,
         "clip_path_to_rect(path, bbox, inside)"},
    {"affine_transform", &CXX_affine_transform, METH_VARARGS,
         "affine_transform(vertices, transform)"},
    {"count_bboxes_overlapping_bbox", &CXX_count_bboxes_overlapping_bbox, METH_VARARGS,
         "count_bboxes_overlapping_bbox(bbox, bboxes)"},
    {"path_intersects_path", &CXX_path_intersects_path, METH_VARARGS,
         "path_intersects_path(p1, p2)"},
    {"convert_path_to_polygons", &CXX_convert_path_to_polygons, METH_VARARGS,
         "convert_path_to_polygons(path, trans, width, height)"},
    {"cleanup_path", &CXX_cleanup_path, METH_VARARGS,
         "cleanup_path(path, trans, remove_nans, clip, snap, simplify, curves, sketch_params)"},
    {"convert_to_svg", &CXX_convert_to_svg, METH_VARARGS,
         "convert_to_svg(path, trans, clip, simplify, precision)"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
#if PY3K
PyInit__path(void)
#else
init_path(void)
#endif
{
    import_array();

    //static _path_module* _path = NULL;
    //_path = new _path_module;

#if PY3K
    //return _path->module().ptr();
    return
#endif
    Py_InitModule3("_path", methods, "Helper functions for paths");
}


