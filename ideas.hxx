



	//modular cgp : 


	//- we have a tgrid an provide a set of 
	//" extract methods"
	
	
	// get the geometry
	Geometry<0,cgp2d::TopolocigalPoint> 				topoGeo0;
	Geometry<0,cgp2d::CartesianPoint>   				cartGeo0;
	Geometry<1,cgp2d::TopolocigalPoint,cgp2d::FILL_GAB> topoGeo1;
	Geometry<1,cgp2d::CartesianPoint>   				cartGeo1;
	Geometry<2,cgp2d::TopolocigalPoint> 				topoGeo2;
	Geometry<2,cgp2d::CartesianPoint,cgp2d::FILL_GAB>  	cartGeo2;

	tgrid.extract_geometry(topoGeo0);
	tgrid.extract_geometry(cartGeo0);
	tgrid.extract_geometry(topoGeo1);
	tgrid.extract_geometry(cartGeo1);
	tgrid.extract_geometry(topoGeo2);
	tgrid.extract_geometry(cartGeo2);

		
	// extract bounds
	Bounds<0> bounds0;
	Bounds<1> bounds1;

	tgrid.extract_bounds(bounds0);
	tgrid.extract_bounds(bounds1);


	// extract bounded by
	BoundedBy<1> boundedBy1;
	BoundedBy<2> boundedBy2;

	tgrid.extract_bounded_by(boundedBy1);
	tgrid.extract_bounded_by(boundedBy2);


	// extract non-hyper edge graph
	AdjacencyGraph<2,1,cgp2d::MERGE_PARALLEL_EDGES> adjacencyGraph21;   // nodes are cell2, edges are cell1 "RAG"
	AdjacencyGraph<0,1> adjacencyGraph01;  								// nodes are junction , edges cell1  exclude "dead edges" ??



	// extract hyper edge graph
	AdjacencyGraph<1,0> adjacencyGraph10;   // nodes are cell1 , hyperedges are cell0



	// algorithms on geometry

	// accumulate from geometry
	Features  features0  = cgp2d.accumulate(0,topoGeo0,"some_acc_settings");
	Features  features1  = cgp2d.accumulate(1,topoGeo1,"some_acc_settings");
	Features  features2  = cgp2d.accumulate(2,topoGeo2,"some_acc_settings");








