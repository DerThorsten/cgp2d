#define PY_ARRAY_UNIQUE_SYMBOL superimg_PyArray_API
#define NO_IMPORT_ARRAY


#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "cgp2d.hxx"
#include "py_cell_visitor.hxx"

namespace cgp2d {

namespace python = boost::python;

void export_cell2()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);


    ////////////////////////////////////////
    // Region Graph
    ////////////////////////////////////////
    // basic types
    // tgrid and input image type
    typedef TopologicalGrid<vigra::UInt32> TopologicalGridType;
	typedef Cgp<vigra::UInt32,vigra::UInt32> CgpType;

    typedef  vigra::NumpyArray<2 ,vigra::Singleband < vigra::UInt32 > > InputLabelImageType;
    // cgp type and cell types
    typedef CgpType::PointType PointType;
    // bound vector
    typedef std::vector<float> FloatVectorType;
    typedef std::vector<vigra::UInt32> LabelVectorType;
    // point vector
    typedef std::vector<PointType> PointVectorType;
    // geo cells 
    typedef CgpType::GeoCell2 GeoCell2Type;
    typedef CgpType::GeoCells2 GeoCell2VectorType;

    // cells
    python::class_<GeoCell2Type>("Cell2",python::init<>())
        .def(CellTypeSuite<GeoCell2Type>())
    ;
}

} // namespace vigra

