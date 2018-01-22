/*
  LFA Lab - Library to simplify local Fourier analysis.
  Copyright (C) 2018  Hannah Rittich

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 
*/

// vim: set filetype=cpp:

// =========================================================
// convert to and from vector
// =========================================================

%{

template <typename T>
T pythonToC(PyObject* obj)
{
    throw std::logic_error("pythonToC: Not implemented for this type.");
}

template <>
double pythonToC<double>(PyObject* obj)
{
    double value = PyFloat_AsDouble(obj);
    // PyFloat_AsDouble returns -1.0 on error
    if (value == -1.0 && PyErr_Occurred()) {
        throw std::runtime_error("Cannot convert to float");
    }
    return value;
}

template <>
long pythonToC<long>(PyObject* obj)
{
    int overflow;

    // only available in >= 3.2
    long value = PyLong_AsLongLongAndOverflow(obj, &overflow);
    if (overflow) {
        throw std::runtime_error("Integer overflow.");
    }

    if(PyErr_Occurred()) {
        throw std::runtime_error("Something is wrong. I mean, really wrong.");
    }

    return value;
}

template <>
int pythonToC<int>(PyObject* obj)
{
    int overflow;

    int value = PyLong_AsLongAndOverflow(obj, &overflow);
    if (overflow) {
        throw std::runtime_error("Integer overflow.");
    }

    if(PyErr_Occurred()) {
        throw std::runtime_error("Something is wrong. I mean, really wrong.");
    }

    return value;
}

template <typename T>
PyObject* CToPython(const T& input) {
    throw std::logic_error("CToPython: Not implemented for this type.");
}

template <>
PyObject* CToPython<double>(const double& input) {
    return PyFloat_FromDouble(input);
}

template <>
PyObject* CToPython<int>(const int& input) {
    return PyLong_FromLong(input);
}

template <typename T>
bool convertibleToC(PyObject* obj) {
    throw std::logic_error("Typecheck not implemented for this type.");
}

template <>
bool convertibleToC<double>(PyObject* obj) {
    return PyFloat_Check(obj);
}
template <>
bool convertibleToC<long>(PyObject* obj) {
    return PyLong_Check(obj);
}

template <typename T>
T pythonToArray(PyObject* input)
{
    using namespace std;

    typedef typename T::Scalar Scalar;

    int l = PySequence_Length(input);
    if (l < 0) {
        throw std::runtime_error("Expected a sequence.");
    }

    if (l >= T::MaxRowsAtCompileTime) {
        stringstream msg;
        msg << "Too many elements in sequence. "
            << "Got " << l
            << " and expected at most " << T::MaxRowsAtCompileTime << ".";
        throw std::runtime_error(msg.str());
    }

    T result(l);
    for (int i = 0; i < l; ++i) {
        PyObject* o = PySequence_GetItem(input, i);
        result[i] = pythonToC<Scalar>(o);
    }
    return result;
}

template <typename T>
PyObject* arrayToPython(const T& vec)
{
    typedef typename T::Scalar Scalar;

    PyObject* tuple = PyTuple_New(vec.rows());
    for (int i = 0; i < vec.rows(); ++i) {
        PyObject* v = CToPython<Scalar>( vec(i) );
        PyTuple_SetItem(tuple, i, v);
    }
    return tuple;
}
%}

%define ARRAY_TYPEMAP(type)
%typemap(typecheck) type {
    $1 = 1;

    if (!PySequence_Check($input))
        $1 = 0;
}
%typemap(in) type {
    try {
        $1 = pythonToArray<type>($input);
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}
%typemap(out) type {
    $result = arrayToPython<type>($1);
}
%enddef

ARRAY_TYPEMAP(ArrayXd)
ARRAY_TYPEMAP(ArrayFd)
ARRAY_TYPEMAP(ArrayXi)
ARRAY_TYPEMAP(ArrayFi)

// =========================================================
// NumPy Typemap

// Convert C++ ==> Python
%typemap(out) MatrixXcd {

    const npy_intp nd = 2; // number of dimensions
    npy_intp dims[nd];
    dims[0] = $1.rows();
    dims[1] = $1.cols();
    PyObject* obj = PyArray_SimpleNew(nd, dims, NPY_COMPLEX128);
    PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj);

    // copy data
    for (int i = 0; i < $1.rows(); ++i) {
        for (int j = 0; j < $1.cols(); ++j) {
            *reinterpret_cast<std::complex<double>* > (PyArray_GETPTR2(arr, i, j)) = $1(i, j);
        }
    }

    $result = obj;
}

%typemap(out) VectorXcd {

    const npy_intp nd = 1; // number of dimensions
    npy_intp dims[nd];
    dims[0] = $1.rows();
    PyObject* obj = PyArray_SimpleNew(nd, dims, NPY_COMPLEX128);
    PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj);

    // copy data
    for (int i = 0; i < $1.rows(); ++i) {
        *reinterpret_cast<std::complex<double>* > (PyArray_GETPTR1(arr, i)) = $1(i);
    }

    $result = obj;
}


