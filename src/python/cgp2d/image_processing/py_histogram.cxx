#define PY_ARRAY_UNIQUE_SYMBOL superimg_PyArray_API
#define NO_IMPORT_ARRAY

#include <string>
#include <cmath>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


#include <boost/array.hpp>

#include <boost/accumulators/accumulators.hpp>

#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "cgp2d.hxx"
#include "cgp2d_python.hxx"


namespace python = boost::python;

namespace cgp2d {



	vigra::NumpyAnyArray  colorHistogram(
		vigra::NumpyArray<2, vigra::TinyVector< float  ,3 > > 	colorImage,
		vigra::NumpyArray<1, vigra::TinyVector< float  ,3 > > 	binning,
		vigra::TinyVector< float  ,3 > 							sigma
	){
		typedef vigra::MultiArray< 5, float  >  Image3dHistType;
		typedef Image3dHistType::difference_type ShapeType;


		// allocate an explicit multidimensional hist
		vigra::MultiArray< 5, float  >  image3dHist(
			ShapeType(colorImage.shape(0),colorImage.shape(1),binning(0)[0],binning(1)[0],binning(2)[0]));


	}	



	void export_histogram(){

		python::def("colorHistogram", vigra::registerConverters( &colorHistogram) );

	}

}