/* -*- mode: c++; c-basic-offset: 4 -*- */

/* mplutils.h   --
 *
 * $Header$
 * $Log$
 * Revision 1.2  2004/11/24 15:26:12  jdh2358
 * added Printf
 *
 * Revision 1.1  2004/06/24 20:11:17  jdh2358
 * added mpl src
 *
 */

#ifndef _MPLUTILS_H
#define _MPLUTILS_H

#include <Python.h>

#include <string>
#include <iostream>
#include <sstream>

#if PY_MAJOR_VERSION >= 3
#define PY3K 1
#else
#define PY3K 0
#endif

#if PY3K
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, methods, doc) \
          static PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef)
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
  #define MOD_DEF(ob, name, methods, doc) \
          ob = Py_InitModule3(name, methods, doc)
#endif

#if PY3K
#define PyInt_FromLong PyLong_FromLong
#define PyString_FromString PyUnicode_FromString
#endif

inline int setIfL(PyObject *d, const char *name, long value) {
    PyObject* o = PyInt_FromLong(value);
    if (!o) return -1;
    int status = PyDict_SetItemString(d, name, o);
    Py_DECREF(o);
    return status;
}

inline int setLfL(PyObject *d, const char *name, long value) {
    PyObject* o = PyLong_FromLong(value);
    if (!o) return -1;
    int status = PyDict_SetItemString(d, name, o);
    Py_DECREF(o);
    return status;
}

inline int setLfUL(PyObject *d, const char *name, unsigned long value) {
    PyObject* o = PyLong_FromUnsignedLong(value);
    if (!o) return -1;
    int status = PyDict_SetItemString(d, name, o);
    Py_DECREF(o);
    return status;
}

void _VERBOSE(const std::string&);


#undef  CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#undef  MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

inline double mpl_round(double v)
{
    return (double)(int)(v + ((v >= 0.0) ? 0.5 : -0.5));
}

class Printf
{
private :
    char *buffer;
    // prevent copying
    Printf(const Printf&);
    Printf& operator=(const Printf&);
public :
    Printf(const char *, ...);
    ~Printf();
    std::string str()
    {
        return buffer;
    }
    friend std::ostream &operator <<(std::ostream &, const Printf &);
};

#if defined(_MSC_VER) && (_MSC_VER == 1400)

/* Required by libpng and zlib */
#pragma comment(lib, "bufferoverflowU")

/* std::max and std::min are missing in Windows Server 2003 R2
   Platform SDK compiler.  See matplotlib bug #3067191 */
namespace std {

    template <class T> inline T max(const T& a, const T& b)
    {
        return (a > b) ? a : b;
    }

    template <class T> inline T min(const T& a, const T& b)
    {
        return (a < b) ? a : b;
    }

}

#endif

#endif
