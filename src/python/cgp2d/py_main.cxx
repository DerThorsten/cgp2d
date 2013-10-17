#define PY_ARRAY_UNIQUE_SYMBOL superimg_PyArray_API

//#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array_converters.hxx>


#include "py_cell0.hxx"
#include "py_cell1.hxx"
#include "py_cell2.hxx"

#include "py_cell0vec.hxx"
#include "py_cell1vec.hxx"
#include "py_cell2vec.hxx"

#include "py_cgp2d.hxx"

BOOST_PYTHON_MODULE_INIT(_cgp2d) {
    //using namespace boost::python;
    //using namespace vigra;

    vigra::import_vigranumpy();
    cgp2d::export_cgp2d();
    cgp2d::export_cell0();
    cgp2d::export_cell1();
    cgp2d::export_cell2();

    cgp2d::export_cell0vec();
    cgp2d::export_cell1vec();
    cgp2d::export_cell2vec();
}
