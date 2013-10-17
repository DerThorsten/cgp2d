#ifndef CGP2D_HXX
#define CGP2D_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>

/* opengm */
//#include "opengm/config.hxx"
//#include "opengm/utilities/metaprogramming.hxx"

/* this project */
#include "partition.hxx"
#include "macros.hxx"



/****************************************************************************/
/* C e l l T y p e                                                          */
/****************************************************************************/
    
struct CellType {
    enum Values {
        Junction =0,
        Boundary =1,
        Region   =2
    };
};

namespace cgp2d {

/****************************************************************************/
/* L i n e 2 d                                                              */
/****************************************************************************/

template<class T>
class Line2d{
public:
    typedef T ValueType;
    typedef vigra::TinyVector<ValueType,2> PointType;
    Line2d(){}
    template<class P>
    Line2d(const P & pA,const P & pB)
    :   pA_(pA),
        pB_(pB){

    }

    template<class ITER>
    ValueType accDistance( ITER begin,ITER end ){
        ValueType dist = 0.0 ;
        while(begin!=end){
            dist+=distance(*begin);
            ++begin;
        }
        return dist;
    }

    template<class V>
    ValueType distance(const vigra::TinyVector<V,2> & p){
        PointType n = pA_ - pB_;
        n/=vigra::norm(n);
        PointType ap = pA_-p;
        const float apn =  vigra::dot(ap,n);
        return vigra::norm(ap - (apn)*n);
    }


    vigra::TinyVector<float,2> angle()const{
        vigra::TinyVector<float,2> r=(pA_- pB_);
        r/=vigra::norm(r);
        return r;
    }
private:
    PointType pA_;
    PointType pB_;
};  

/****************************************************************************/
/* L i n e P a t h 2 d                                                      */
/****************************************************************************/

template<class T>
class LinePath2d{
public:
    typedef vigra::TinyVector<T,2>         PointType;
    typedef std::vector<PointType>  PointVector;
    typedef Line2d<float>           FLine2d;
    typedef std::vector<FLine2d>    FLine2dVec;

    LinePath2d(const PointVector & inPoints,std::vector<vigra::TinyVector<float,2> >  & angles)
    :   inPoints_(inPoints),
        lineVec_(),
        angles_(angles)
    {
        const size_t numPoints=inPoints_.size();

        if (numPoints==1){
            const PointType & p = inPoints_.front();
            // A|B
            // vertical  boundary point
            if(p[0]%2==1){
                FLine2d line(PointType(0,0),PointType(0,1));
                angles_[0]=line.angle();
            }
            else{
                FLine2d line(PointType(0,0),PointType(1,0));
                angles_[0]=line.angle();
            }
        }
        else if (numPoints<=3){  
            FLine2d line(inPoints_.front(),inPoints_.back());
            for(size_t i=0;i<numPoints;++i){
                angles_[i]=line.angle();
            }
            const float dist = line.accDistance(inPoints_.begin(),inPoints_.end());
            //std::cout<<"tiny ! dist "<<dist<<"\n";
        }
        else{
            //std::cout<<"to big\n";
            size_t end = 0;
            while(end!=inPoints_.size()){
                //std::cout<<"-----\n";
                size_t start = end;
                end = this->getAngle(start);
            }
            
        }



    }

    size_t getAngle(const size_t start){
        //std::cout<<"start "<<start<<"\n";
        size_t end = inPoints_.size();
        const float  threshold = 2.0;
        float distance = 9999999.0f;
        size_t nP  = end -start;
        while(true) {

            //std::cout<<" * end "<<end<<"\n";
            CGP_ASSERT_OP(nP,>=,2);
            CGP_ASSERT_OP(start,!=,end-1);
            FLine2d line(inPoints_[start],inPoints_[end-1]);
            distance = line.accDistance(inPoints_.begin()+start,inPoints_.begin()+end);
            if(nP<=3 || distance < threshold){
                for(size_t i=start;i<end;++i){
                    angles_[i]=line.angle();
                }
                break;
            }
            else{
                end = start+nP/2;
                nP  = end -start;
            }
        }
        //std::cout<<"at ending where end is "<<end<<"\n";
        return end; 
    }
private:
    const PointVector &        inPoints_;
    FLine2dVec                 lineVec_;
    std::vector<vigra::TinyVector<float,2> >  &       angles_;

};

/****************************************************************************/
/* T o p o l o g i c a l G r i d                                            */
/****************************************************************************/

template<class LABEL_TYPE>
class TopologicalGrid{
public:
    typedef partition::Partition<size_t>    UfdType;
    typedef LABEL_TYPE LabelType;
    typedef vigra::MultiArray<2,LabelType>  LabelImageType;
    typedef typename LabelImageType::difference_type ShapeType;

    TopologicalGrid(){

    }
    // constructor 
    template<class INPUT_IMG>
    TopologicalGrid(const INPUT_IMG & seg);

    // query
    const LabelImageType & tgrid()const;
    size_t numCells(const size_t i)const;
    size_t shape(const size_t d)const;

    const ShapeType & shapeTopologicalGrid()const{
        return tShape_;
    }
    const ShapeType & shapeLabeling()const{
        return lShape_;
    }

    LabelType operator()(const size_t tx,const size_t ty)     {return tgrid_(tx,ty);}
    LabelType operator()(const size_t tx,const size_t ty)const{return tgrid_(tx,ty);}

private:

    size_t numCells_[3];
    LabelImageType tgrid_;
    ShapeType tShape_;
    ShapeType lShape_;

};

template<class COORDINATE_TYPE,class LABEL_TYPE>
class Geometry;

template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cgp;

/****************************************************************************/
/* C e l l B a s e                                                          */
/****************************************************************************/

template<class COORDINATE_TYPE,class LABEL_TYPE,int CELLTYPE>
class CellBase {

    // friend classes
    friend class Cgp<COORDINATE_TYPE,LABEL_TYPE>;
public:

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef vigra::TinyVector<float,2> FloatPointType;
    //typedef BoundingBox<CoordinateType> BoundingBoxType;
    typedef LABEL_TYPE LabelType;
    size_t size()const{
        return points_.size();
    }
    const PointType & operator[](const size_t i) const{
        CGP_ASSERT_OP( i , < , points_.size() );
        return points_[i];
    }
    const PointType & operator[](const size_t i) {
        CGP_ASSERT_OP(i,<,points_.size());
        return points_[i];
    }
    LabelType label()const{
        return label_;
    }


    bool operator == (const CellBase & other){
        return label_==other.label_;
    }
    bool operator != (const CellBase & other){
        return label_!=other.label_;
    }



    const CgpType & cgp()const{
        return *cgp_;
    }
    
    size_t cellType()const{
        return static_cast<size_t>(CELLTYPE);
    }


    std::pair<PointType,PointType> boundingBox()const{
        PointType ul=points_[0];
        PointType lr=points_[0];
        for(size_t p=0;p<size();++p){
            ul[0] = points_[p][0] < ul[0] ? points_[p][0] : ul[0];
            ul[1] = points_[p][1] < ul[1] ? points_[p][1] : ul[1];
            lr[0] = points_[p][0] > lr[0] ? points_[p][0] : lr[0];
            lr[1] = points_[p][1] > lr[1] ? points_[p][1] : lr[1];
        }
        return std::pair<PointType,PointType>(ul,lr);
    }
    
    
    
    FloatPointType centerCoordinate()const{
        FloatPointType cp(0.0f,0.0f);
        for(size_t p=0;p<size();++p){
            cp+=points_[p];
        }
        cp/=size();
        return cp;
    }
    
    const std::vector<LabelType> & bounds()const{
        return bounds_;
    }
    const std::vector<LabelType> & boundedBy()const{
        return boundedBy_;
    }

    const std::vector<PointType> & points()const{
        return points_;
    }
protected:
    void push_back_point(const PointType & p){
        points_.push_back(p);
    }

    void push_back_bound(const LabelType l){
        bounds_.push_back(l);
    }

    void push_back_bounded_by(const LabelType l){
        boundedBy_.push_back(l);
    }

    void sortAdjaceny(){
        std::sort(bounds_.begin(),bounds_.end());
        std::sort(boundedBy_.begin(),boundedBy_.end());
    }

    void setLabel(const LabelType l){
        label_=l;
    }




    LabelType label_;
    // coordinates
    std::vector<PointType> points_;

    // bounds
    std::vector<LabelType> bounds_;
    std::vector<LabelType> boundedBy_;
    std::vector<LabelType> adjaceny_;

    // cgp pointer
    CgpType * cgp_;
};


template<class COORDINATE_TYPE,class LABEL_TYPE,int CELLTYPE>
class Cell;

template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cell<COORDINATE_TYPE,LABEL_TYPE,0> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,0>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;


    void getAngles(const TopologicalGridType & tgrid,const size_t radius,std::vector<float> & angles)const{
        const size_t numBoundaries = this->bounds_.size();
        const int r=static_cast<int>(radius);
        const CoordinateType tx=this->points_[0][0],ty=this->points_[0][1];
        const CoordinateType xmin=    static_cast<int>(tx)-r < 0 ? 0 : static_cast<CoordinateType>(static_cast<int>(tx)-r );
        const CoordinateType ymin=    static_cast<int>(ty)-r < 0 ? 0 : static_cast<CoordinateType>(static_cast<int>(ty)-r );
        const CoordinateType xmax=    static_cast<int>(tx)+r+1 > tgrid.shape(0) ?   tgrid.shape(0) : static_cast<CoordinateType>(static_cast<int>(tx)+r+1 );
        const CoordinateType ymax=    static_cast<int>(ty)+r+1 > tgrid.shape(1) ?   tgrid.shape(1) : static_cast<CoordinateType>(static_cast<int>(ty)+r+1 );

        //std::cout<<"min "<<xmin<<" , "<<ymin<<"\n";
        //std::cout<<"max "<<xmax<<" , "<<ymax<<"\n";

        typedef std::pair<PointType,LabelType>  MapItem;
        typedef std::map<LabelType,MapItem > AverageMapType;
        typedef typename AverageMapType::const_iterator AverageMapConstIter;
        typedef typename AverageMapType::iterator AverageMapIter;

        AverageMapType averageMap;
        for(size_t b=0;b<this->bounds_.size();++b){
            MapItem initItem= MapItem(PointType(0,0),0);
            averageMap[this->bounds_[b]]=initItem;
        }

        // collect
        for(CoordinateType tyy=ymin;tyy<ymax;++tyy)
        for(CoordinateType txx=xmin;txx<xmax;++txx){
            // if boundary
            if(  (txx%2==1 && tyy%2==0) || (txx%2==0 && tyy%2==1) ){

                LabelType cell1Label=tgrid(txx,tyy);
                if(cell1Label!=0){
                    AverageMapIter iter=averageMap.find(cell1Label);
                    if(iter!=averageMap.end()){
                        MapItem item=iter->second;
                        item.first+=PointType(txx,tyy);
                        ++item.second;
                        averageMap[cell1Label]=item;
                    }
                }
            }
        }

        angles.resize(numBoundaries);
        size_t index=0;
        for(AverageMapConstIter iter=averageMap.begin();iter!=averageMap.end();++iter,++index){\
            MapItem item=iter->second;
            PointType averagePoint = item.first;


            averagePoint/=item.second;

            const float x=static_cast<float>(tx);
            const float y=static_cast<float>(ty);
            const float ax=static_cast<float>(averagePoint[0]);
            const float ay=static_cast<float>(averagePoint[1]);


            const float rx=ax-x;
            const float ry=ay-y;


            //std::cout<<" num P "<<item.second<<"\n";
            //std::cout<<"point          "<< x << " , "<< y<<"\n";
            //std::cout<<"averge point   "<< ax<< " , "<<ay<<"\n";
            //std::cout<<"relative point "<< rx<< " , "<<ry<<"\n";
    
            const float result = std::atan2 (ry,rx) * 180.0 / M_PI;
            angles[index]=result;
        }
        //std::cout<<"\n";
    }
};

template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cell<COORDINATE_TYPE,LABEL_TYPE,1> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,1>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef vigra::TinyVector<float,2> FloatPointType;

    typedef LinePath2d<CoordinateType> LinePath2dType;

    void getAngles(){
        angles_.resize(this->size());
        LinePath2dType linePath(this->points_,angles_);            
    }

    float getAngles(const size_t start,const size_t center,const size_t end){
        //std::cout<<"s "<<start<<" c "<<center<<" e "<<end<<"\n";

        FloatPointType   sp   = this->points_[start];
        FloatPointType   akkP(0.0,0.0);
        size_t c=0;
        for(size_t i=start+1;i<end;++i){
            ++c;
            akkP+=this->points_[i];
        }
        CGP_ASSERT_OP(c,==,(end-start-1));
        akkP/=(end-start-1);
        akkP-=sp;

        const float result = std::atan2 (akkP[1],akkP[0]) * 180.0 / M_PI;
        return result;
    }


    void sortCells(){
        size_t finishedPoints=0;
        std::vector<bool>  finishedPoint(this->size(),false);
        std::deque<LabelType> sorted;
        PointType pFront = this->points_[0];
        PointType pBack  = this->points_[0];

        // insert first point 
        sorted.push_back(0);
        finishedPoint[0]=true;
        finishedPoints=1;

        while(finishedPoints < this->size()){
            const size_t oldSize=finishedPoints;
            for(size_t i=0;i<this->size();++i){
                if(finishedPoint[i]==false){
                    const PointType & point = this->points_[i];
                    if(adjacent2Cells(pFront,point)){
                        sorted.push_front(i);
                        pFront=point;
                        finishedPoint[i]=true;
                        ++finishedPoints;
                        break;
                    }
                    else if(adjacent2Cells(pBack,point)){
                        sorted.push_back(i);
                        pBack=point;
                        finishedPoint[i]=true;
                        ++finishedPoints;
                        break;
                    }
                }
            }
            if(oldSize+1!=finishedPoints){
                std::cout<<"\n\n\n\n size "<<this->size()<<"\n";
            }
            CGP_ASSERT_OP(oldSize+1,==,finishedPoints);
        }

        std::vector<PointType> points;
        points.reserve(this->size());
        while(sorted.empty()==false){
            points.push_back(this->points_[sorted.front()]);
            sorted.pop_front();
        }
        this->points_=points;

        this->getAngles();
    }
private:
    bool adjacent2Cells(const PointType & pa,const PointType & pb ){
        
        // A|B
        // vertical  boundary point
        if(pa[0]%2==1){
            // six possible neighbours:
            //
            //    |      case  1
            //  --*--    case 2,3
            //    |               <-self
            //  --*--    case 4,5
            //    |      case  6

            //    
            //    1
            //  2 * 3
            //    |
            //  4 * 5
            //    6

            //case 1 (with border check)
            if     (pa[1]!=0 && ( pa[0]  ==pb[0] && pa[1]-2==pb[1] ) ){return true;}
            //case 2 
            else if(pa[1]!=0 && ( pa[0]-1==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 3 
            else if(pa[1]!=0 && ( pa[0]+1==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 4 
            else if(            ( pa[0]-1==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 5 
            else if(            ( pa[0]+1==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 6 
            else if(            ( pa[0]  ==pb[0] && pa[1]+2==pb[1] ) ){return true;}
            
            return false;
        }
        // horizontal boundary
        else{
            // six possible neighbours:
            //
            //      |  |        
            //    --*--*--      
            //      |  | 
            //          
            //      2  4
            //    1 *--* 6
            //      3  5 

            //case 1 (with border check)
            if     (pa[0]!=0  && ( pa[0] -2 ==pb[0] && pa[1]  ==pb[1] ) ){return true;}
            //case 2 
            else if(pa[0]!=0  && ( pa[0] -1 ==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 3 
            else if(pa[0]!=0  && ( pa[0] -1 ==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 4 
            else if(             ( pa[0] +1 ==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 5 
            else if(             ( pa[0] +1 ==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            else if(             ( pa[0] +2 ==pb[0] && pa[1]  ==pb[1] ) ){return true;}

            return false;
        }
    }
public:
    std::vector<vigra::TinyVector<float,2> > angles_;
};

template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cell<COORDINATE_TYPE,LABEL_TYPE,2> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,2>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
private:
};

/****************************************************************************/
/* C g p                                                                    */
/****************************************************************************/


template<class COORDINATE_TYPE, class LABEL_TYPE>
class Cgp;



template<int CELL_TYPE,class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper{
};


template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<0,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,0>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry0(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry0(); }
};




template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<1,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,1>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry1(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry1(); }
};



template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<2,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,2>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry2(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry2(); }
};





template<class COORDINATE_TYPE, class LABEL_TYPE>
class Cgp {
    
    public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef typename TopologicalGridType::ShapeType ShapeType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,0> GeoCell0;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,1> GeoCell1;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,2> GeoCell2;
    typedef std::vector< GeoCell0 > GeoCells0;
    typedef std::vector< GeoCell1 > GeoCells1;
    typedef std::vector< GeoCell2 > GeoCells2;

    typedef std::vector<std::vector< LABEL_TYPE> > CellAdjacencyGraphVectorType;
    typedef std::vector<std::set< LABEL_TYPE> > CellAdjacencyGraphSetType;

    // Constructor
    Cgp(const TopologicalGridType & tgrid );
    // Query
    const GeoCells0  & geometry0()const;
    const GeoCells1  & geometry1()const;
    const GeoCells2  & geometry2()const;

    const TopologicalGridType & tgrid()const;

    size_t numCells(const size_t cellType)const{
        return tgrid_.numCells(cellType);
    }

    size_t shape(const size_t d)const{
        return tgrid_.shape(d);
    }

    const ShapeType & shapeTopologicalGrid()const{
        return tgrid_.shapeTopologicalGrid();
    }

    const ShapeType & shapeLabeling()const{
        return tgrid_.shapeLabeling();
    }


    LabelType operator()(const size_t x,const size_t y) const{
        return tgrid_(x,y);
    }

    LabelType operator()(const size_t x,const size_t y)     {
        return tgrid_(x,y);
    }
    


    template<int CELL_TYPE>
    std::vector< Cell<COORDINATE_TYPE,LABEL_TYPE,CELL_TYPE> > & cells()const{
        return CgpHelper<CELL_TYPE,CoordinateType,LabelType>::getCells(*this);
    }




    template<int CELL_TYPE>
    size_t cellSizeT(const LabelType cellIndex)const{
        return this-> template cells<CELL_TYPE> ()[cellIndex].size();
    }


    size_t cellSize(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                return this-> template cells<0>()[cellIndex].size();
            case 1:
                return this-> template cells<1>()[cellIndex].size();
            case 2:
                return this-> template cells<2>()[cellIndex].size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }



    size_t nBounds(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                return this-> template cells<0>().bounds().size();
            case 1:
                return this-> template cells<1>().bounds().size();
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }
    
    template<int CELL_TYPE>
    size_t nBoundsT(const LabelType cellIndex)const{
        if(CELL_TYPE==2){
            CGP_ASSERT_OP(false,==,true);
        }
        return this-> template cells<CELL_TYPE>().bounds().size();

    }


    size_t nBoundedBy(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return geoCells1_[cellIndex].boundedBy().size();
            case 2:
                return geoCells2_[cellIndex].boundedBy().size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }
    template<int CELL_TYPE>
    size_t nBoundedByT(const LabelType cellIndex)const{
        switch(CELL_TYPE){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return geoCells1_[cellIndex].boundedBy().size();
            case 2:
                return geoCells2_[cellIndex].boundedBy().size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    LabelType bound(const size_t cellType ,const LabelType cellIndex,const size_t boundNr)const{
        switch(cellType){
            case 0:
                return geoCells0_[cellIndex].bounds()[boundNr];
            case 1:
                return geoCells1_[cellIndex].bounds()[boundNr];
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    template<int CELL_TYPE>
    LabelType bound(const LabelType cellIndex,const size_t boundNr)const{
        switch(CELL_TYPE){
            case 0:
                return geoCells0_[cellIndex].bounds()[boundNr];
            case 1:
                return geoCells1_[cellIndex].bounds()[boundNr];
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }


    LabelType boundedBy(const size_t cellType ,const LabelType cellIndex,const size_t boundedByNr)const{
        switch(cellType){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return geoCells1_[cellIndex].boundedBy()[boundedByNr];
            case 2:
                return geoCells2_[cellIndex].boundedBy()[boundedByNr];
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    template<int CELL_TYPE>
    LabelType boundedBy(const LabelType cellIndex,const size_t boundedByNr)const{
        switch(CELL_TYPE){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return geoCells1_[cellIndex].boundedBy()[boundedByNr];
            case 2:
                return geoCells2_[cellIndex].boundedBy()[boundedByNr];
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }



    void cellAdjacencyGraphVector(int cellType,CellAdjacencyGraphVectorType & graph)const{
        CGP_ASSERT_OP(cellType,>=,1);
        CGP_ASSERT_OP(cellType,<=,2);
        const size_t numCells=this->numCells(cellType);
        graph.clear();
        graph.resize(numCells);

        CellAdjacencyGraphSetType  graphSet;
        this->cellAdjacencySetGraph(cellType,graphSet);
        for(size_t ci=0;ci<numCells;++ci){
            graph[ci].assign(graphSet[ci].begin(),graphSet[ci].end());
        }
    }


    void cellAdjacencySetGraph(int cellType,CellAdjacencyGraphSetType & graph)const{
        CGP_ASSERT_OP(cellType,>=,1);
        CGP_ASSERT_OP(cellType,<=,2);

        const size_t numCells=this->numCells(cellType);
        graph.clear();
        graph.resize(numCells);

        if(cellType==1){

            const size_t nBoundaries = numCells;
            // iterate over all boundaries
            for(size_t bi=0;bi<nBoundaries;++bi){

                // iterate over all junctions which bound the boudary
                const size_t nBoundedBy = geoCells1_[bi].boundedBy().size();
                CGP_ASSERT_OP(nBoundedBy,<=,2);

                for(size_t j=0;j<nBoundedBy;++j){
                    const LabelType  ji = geoCells1_[bi].boundedBy()[j]-1.0;

                    // get the number boundaries for the junction
                    const size_t nBounds = geoCells0_[ji].bounds().size();
                    CGP_ASSERT_OP(nBounds,>=,3);
                    CGP_ASSERT_OP(nBounds,<=,4);

                    for(size_t b=0;b<nBounds;++b){

                        // insert other region to to adj.
                        const LabelType biOther = geoCells0_[ji].bounds()[b];
                        if(biOther!=bi){
                            graph[bi].insert(biOther);
                        }
                    }
                }
            }
        }

        if(cellType==2){
            const size_t nRegion = numCells;

            // iterate over all regions
            for(size_t ri=0;ri<nRegion;++ri){

                // iterate over all faces which bound the region
                const size_t nBoundedBy = geoCells2_[ri].boundedBy().size();
                const bool assertTrue  = nBoundedBy>0 || nRegion==1;
                CGP_ASSERT_OP(assertTrue,==,true);

                for(size_t b=0;b<nBoundedBy;++b){
                    const LabelType  bi = geoCells2_[ri].boundedBy()[b]-1.0;

                    // get the 2 regions of the boundary
                    const size_t nBounds = geoCells1_[bi].bounds().size();
                    CGP_ASSERT_OP(nBounds,==,2);
                    for(size_t r=0;r<2;++r){

                        // insert other region to to adj.
                        const LabelType riOther = geoCells1_[bi].bounds()[r];
                        if(riOther!=ri){
                            graph[ri].insert(riOther);
                        }
                    }
                }
            }
        }
    }


    /** given a consecutive list of topological points
     *  (that define a line in 2D),
     *  return a list of cartesian coordinates which can be used
     *  to draw the line.
     */
    std::vector<PointType> cartesianLine(const GeoCell1& c1) const {
        std::vector<PointType> res;
        for(int x=0; x<c1.size(); ++x) {
            const typename GeoCell1::PointType& pt = c1[x];

            PointType p1; //start
            PointType p2; //end
            if(pt[0] % 2 == 0) { //normal in y direction
                p1[0] = (pt[0])/2;
                p2[0] = (pt[0])/2+1;
                p1[1] = (pt[1]+1)/2;
                p2[1] = (pt[1]+1)/2;
            }
            else { //normal in x direction
                p1[0] = (pt[0]+1)/2;
                p2[0] = (pt[0]+1)/2;
                p1[1] = (pt[1])/2;
                p2[1] = (pt[1])/2+1;
            }
            res.push_back(p1);
            res.push_back(p2);
        }
        
        if(res.size() > 2) {
            //sort points such that they form a consecutive line.
            //to do this, look at consecutive pairs of points and reorder them.
            for(int x=0; x<res.size()-2; x+=2) {
                if(res[x+0] == res[x+3]) {
                    std::swap(res[x+0], res[x+1]);
                    std::swap(res[x+2], res[x+3]);
                }
                else if(res[x+1] == res[x+2]) {
                    continue;
                }
                else if(res[x+0] == res[x+2]) {
                    std::swap(res[x+0], res[x+1]);
                }
                else if(res[x+1] == res[x+3]) {
                    std::swap(res[x+2], res[x+3]);
                }
                else {
                    throw std::runtime_error("path is not connected");
                }
                if(res[x+1] != res[x+2]) {
                    throw std::runtime_error("err");
                }
            }
        }
        
        std::vector<PointType> res2;
        for(int x=0; x<res.size(); x+=2) {
            if(x==res.size()-2) {
                res2.push_back(res[x]);
                res2.push_back(res[x+1]);
            }
            else {
                res2.push_back(res[x]); 
            }
        }
        
        return res2;
    }
    
    std::vector<unsigned int> serialize() const {
        //do not wonder that this code looks odd,
        //it was copied and adapted from the 3D code
        
        typedef std::vector<const GeoCell1*> Lines;
        typedef std::map<LABEL_TYPE, Lines> Label2Lines;
        typedef unsigned int data_t;
        typedef std::vector<unsigned int> SerializeType; 
        
        std::map<LABEL_TYPE, std::vector<const GeoCell1*> > label2lines;
        for(typename GeoCells1::const_iterator it=geoCells1_.begin(); it!=geoCells1_.end(); ++it) {
            label2lines[it->label()].push_back(&(*it));
        }
        
        //determine size of serialized data
        int numPoints = 0;
        for(typename Label2Lines::const_iterator it=label2lines.begin(); it!=label2lines.end(); ++it) {
            numPoints += 2;               //<label><size>
            for(typename Lines::const_iterator l= it->second.begin();
                l!=it->second.end(); ++l)
            {
                numPoints += 2*((*l)->size() +1);
            }
            numPoints += 2*it->second.size()-2; //extra room for seperator tokens
        }

        //array in which to serialize into
        SerializeType M(numPoints); 

        data_t* m    = &M[0];
        data_t* mEnd = m + numPoints;

        const data_t sepToken = std::numeric_limits<data_t>::max();

        //now serialize the whole slice into the array
        data_t* k;
        for(typename Label2Lines::const_iterator it = label2lines.begin(); it!=label2lines.end(); ++it) {
            const LABEL_TYPE& label = it->first;
            const std::vector<const GeoCell1*>& lines = it->second;

            //write label
            *m = label; ++m;
            
            //indicate following length
            k = m; //save position
            ++m;

            for(typename Lines::const_iterator l=lines.begin(); l!=lines.end(); ++l) {
                const GeoCell1& c1 = **l;
                
                std::vector<PointType> line = cartesianLine(c1);
                
                for(int x=0; x<line.size(); ++x) {
                    const typename GeoCell1::PointType& pt = line[x];
                    
                    *m = pt[0]; ++m;
                    *m = pt[1]; ++m;
                }
                if(l!=(--lines.end())) {
                    //write seperator token
                    *m = sepToken; ++m;
                    *m = sepToken; ++m;
                }
            }
            if( (m-k-1) % 2 != 0) {throw std::runtime_error("not good"); }
            *(k)   = m-k-1;

        }
        if(m != mEnd) throw std::runtime_error("Oh no!");
        
        return M;
    }

    template<class CELL1_STATE_ITER>
    void merge2Cells (
        CELL1_STATE_ITER statesBegin,
        CELL1_STATE_ITER statesEnd,
        TopologicalGridType & newTGrid
    ) const {
        typedef partition::Partition<LabelType>    UfdType;

        const size_t numCell1=numCells(1);
        const size_t numCell2=numCells(2);

        UfdType cell2Ufd(numCell2);
        // merging
        for(size_t  cell1Index=0;cell1Index<numCell1;++cell1Index,++statesBegin){
            const size_t state = static_cast<size_t>(*statesBegin);
            const bool mergeCell2= state==0;
            if(mergeCell2){
                const LabelType cell2IndexA=geoCells1_[cell1Index].bounds()[0]-1;
                const LabelType cell2IndexB=geoCells1_[cell1Index].bounds()[1]-1;
                cell2Ufd.merge(cell2IndexA,cell2IndexB);
            }
        }
        std::map<LabelType,LabelType> denseRelabeling;
        cell2Ufd.representativeLabeling(denseRelabeling);
        //std::cout<<"number of sets "<<cell2Ufd.numberOfSets()<<"\n";
        typedef vigra::MultiArray<2,LabelType>  LabelImageType;
        typedef typename LabelImageType::difference_type ShapeType;

        ShapeType labelingShape( (shape(0)+1)/2,(shape(1)+1)/2 );
        LabelImageType newLabeling(labelingShape);

        for(size_t y=0;y<labelingShape[1];++y)
        for(size_t x=0;x<labelingShape[0];++x){
            // get old label from tgrid;
            const LabelType oldCell2Label=tgrid_(x*2,y*2);
            const LabelType oldCell1Index=oldCell2Label-1;
            // find old index in udf
            const LabelType foundIndex=cell2Ufd.find(oldCell1Index);
            // relabel to dense new INDEX
            const LabelType newCell2Index=denseRelabeling[foundIndex];
            // finaly the new label
            const LabelType newCell2Label=newCell2Index+1;
            newLabeling(x,y)=newCell2Label;
        }
        newTGrid=TopologicalGridType(newLabeling);
    }

    private:
    const TopologicalGridType & tgrid_;
    std::vector< GeoCell0 >  geoCells0_;
    std::vector< GeoCell1 >  geoCells1_;
    std::vector< GeoCell2 >  geoCells2_;
};


// Implementation ToopogicalGrid
template<class LABEL_TYPE>
template<class INPUT_IMG>
TopologicalGrid<LABEL_TYPE>::TopologicalGrid(const INPUT_IMG & seg)
: tgrid_(ShapeType(seg.shape(0)*2-1,seg.shape(1)*2-1)),
    tShape_(seg.shape(0)*2-1,seg.shape(1)*2-1),
    lShape_(seg.shape(0),seg.shape(1))
{
    const size_t dx=seg.shape(0);
    const size_t dy=seg.shape(1);
    const size_t tdx=seg.shape(0)*2-1;
    const size_t tdy=seg.shape(1)*2-1;
    size_t shape[] = { tdx,tdy};
    // counters
    size_t maxRegionLabel=0; //TODO TODO TODO
    size_t junctionIndex=0;
    size_t boundaryElementLabel=1;
    ////////////////
    // 1. PASS //
    ////////////////
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        //std::cout<<" tx "<<tx<<" ty "<<ty<<"\n";
        // if region
        if(tx%2==0 && ty%2==0){
            size_t label=seg(tx/2,ty/2);
            tgrid_(tx,ty)=label;
            maxRegionLabel=label>maxRegionLabel ? label : maxRegionLabel;
        }
        // if junction
        else if(tx%2!=0 && ty%2!=0){
            //  A|B
            //  _ _
            //  C|D
            // check labels of A,B,C and D
            std::set<LABEL_TYPE> lset;
            lset.insert( seg((tx-1)/2,(ty-1)/2));  // A
            lset.insert( seg((tx+1)/2,(ty-1)/2));  // B
            lset.insert( seg((tx-1)/2,(ty+1)/2));  // A
            lset.insert( seg((tx+1)/2,(ty+1)/2));  // A
            if(lset.size()>=3){
                tgrid_(tx,ty)=junctionIndex+1;
                ++junctionIndex;
            }
            else{
                tgrid_(tx,ty)=0;
            }
        }
        // boundary
        else{
            size_t l0,l1;
            // A|B
            // vertical  boundary 
            if(tx%2==1){
                l0=seg( (tx-1)/2, ty/2 );
                l1=seg( (tx+1)/2, ty/2 );
            }
            // horizontal boundary
            else{
                l0=seg( tx/2, (ty-1)/2);
                l1=seg( tx/2, (ty+1)/2);
            }
            // active boundary ?
            if(l0!=l1){
                //std::cout<<l0<<"/"<<l1<<"\n";
                tgrid_(tx,ty)=boundaryElementLabel;
                ++boundaryElementLabel;
            }
            else
                tgrid_(tx,ty)=0;
        }
    }
    /////////////////
    // 2. PASS //
    /////////////////
    UfdType boundaryUdf(boundaryElementLabel-1);
    const size_t num1Elements=boundaryElementLabel-1;
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        // boundary
        if((tx%2!=0 && ty%2!=1 ) || (tx%2!=1 && ty%2!=0 )) {
            if ( tgrid_(tx,ty)!=0){
                size_t ownIndex=tgrid_(tx,ty);
                // vertical boundary
                if(tx%2==1){
                    // each horizontal boundary has 6 candidate neighbours:
                    //  _|_
                    //  _|_  <= this is the boundary we are looking right now
                    //   |
                    // if junction down is inctive
                    LabelType other;
                    if (ty+1 < tdy  && tgrid_(tx,ty+1)==0){
                        // boundary up is active?
                        other=tgrid_(tx,ty+2);
                        if( other!=0){
                            CGP_ASSERT_OP(other-1, <, num1Elements);
                            boundaryUdf.merge(ownIndex-1,other-1);
                        }
                        // boundary left is active?
                        other=tgrid_(tx-1,ty+1);
                        if( other!=0){
                            CGP_ASSERT_OP(other-1 ,<, num1Elements);
                            boundaryUdf.merge(ownIndex-1,other-1 );
                        }
                        // boundary right is active?
                        if( tgrid_(tx+1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty+1)-1 );
                        }
                    }
                    // if junction up is inctive
                    if(ty > 0 && tgrid_(tx,ty-1)==0){
                        // boundary up is active?
                        if( tgrid_(tx,ty-2)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx,ty-2) -1);
                        }
                        // boundary left is active?
                        if( tgrid_(tx-1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty-1)-1 );
                        }
                        // boundary right is active?
                        if( tgrid_(tx+1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty-1)-1 );
                        }
                    }
                }
                // horizontal boundary 
                else{
                    //   each horizontal boundary has 6 candidate neighbours:
                    //   _|_|_
                    //    | |
                    // 
                    // if left junction inactive?     
                    if(tx >0 && tgrid_( tx-1,ty)==0){
                        // boundary left is active?
                        if( tgrid_(tx-2,ty)!=0){
                            //std::cout<<"merge left \n";
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-2,ty)-1 );
                        }
                        // boundary up is active?
                        if( tgrid_(tx-1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty-1)-1 );
                        }
                        // boundary down is active?
                        if( tgrid_(tx-1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty+1)-1 );
                        }
                    }
                    // if right junction inactive?     
                    if(tx+1<tdx &&tgrid_( tx+1,ty)==0){
                        // boundary right is active?
                        if( tgrid_(tx+2,ty)!=0){
                            //std::cout<<"merge right \n";
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+2,ty)-1 );
                        }
                        // boundary up is active?
                        if( tgrid_(tx+1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty-1)-1 );
                        }
                        // boundary down is active?
                        if( tgrid_(tx+1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty+1)-1 );
                        }
                    }
                }
            }
        }
    }

    // dense relabeling
    std::map<size_t,size_t> relabel;
    boundaryUdf.representativeLabeling(relabel);
    /////////////////
    // 3. PASS //
    /////////////////
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        // boundary
        if((tx%2!=0 && ty%2!=1 ) || (tx%2!=1 && ty%2!=0 )) {
            if(tgrid_(tx,ty)!=0){
                // relabel
                size_t notDenseIndex=boundaryUdf.find( tgrid_(tx,ty)-1 );
                size_t denseIndex=relabel[notDenseIndex];
                tgrid_(tx,ty)=denseIndex+1;
            }
        }
    }

    // update cell counters
    numCells_[CellType::Region]=maxRegionLabel;
    CGP_ASSERT_OP(boundaryUdf.numberOfSets(),==,relabel.size());
    numCells_[CellType::Boundary]=relabel.size();
    numCells_[CellType::Junction]=junctionIndex;
}

template<class LABEL_TYPE>
const typename TopologicalGrid<LABEL_TYPE>::LabelImageType & TopologicalGrid<LABEL_TYPE>::tgrid()const{
    return tgrid_;
}

template<class LABEL_TYPE>
size_t TopologicalGrid<LABEL_TYPE>::numCells(const size_t i)const{
    return numCells_[i];
}

template<class LABEL_TYPE>
size_t TopologicalGrid<LABEL_TYPE>::shape(const size_t d)const{
    return tgrid_.shape(d);
}

// Implementation Cgp
template<class COORDINATE_TYPE,class LABEL_TYPE>
Cgp<COORDINATE_TYPE,LABEL_TYPE>::Cgp(const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::TopologicalGridType  & tgrid )
:   tgrid_(tgrid),
    geoCells0_(tgrid.numCells(0)),
    geoCells1_(tgrid.numCells(1)),
    geoCells2_(tgrid.numCells(2))
{
    // set up geometry
    const typename TopologicalGridType::LabelImageType & grid=tgrid.tgrid();
    for(size_t ty=0;ty<tgrid.shape(1);++ty)
    for(size_t tx=0;tx<tgrid.shape(0);++tx){
        // Cell 2 (=Region)
        if(tx%2==0 && ty%2==0){
            int label = grid(tx,ty);
            CGP_ASSERT_OP(label,>,0)
            CGP_ASSERT_OP(label,<=,tgrid.numCells(2));
            geoCells2_[label-1].push_back_point(PointType(tx,ty));
        }
        // Cell 0 (== Junction)
        else if(tx%2!=0 && ty%2!=0){
            int label = grid(tx,ty);
            if(label!=0){
                CGP_ASSERT_OP(label,>,0)
                CGP_ASSERT_OP(label,<=,tgrid.numCells(0));
                geoCells0_[label-1].push_back_point(PointType(tx,ty));
            }
                
        }
        // Cell 1 (== Boundary)
        else{
            int label = grid(tx,ty);
            if(label!=0){
                CGP_ASSERT_OP(label,>,0)
                CGP_ASSERT_OP(label,<=,tgrid.numCells(1));
                geoCells1_[label-1].push_back_point(PointType(tx,ty));
            }
        }
    }
    // check size of geometry
    CGP_ASSERT_OP(geoCells0_.size(),==,tgrid.numCells(0));
    CGP_ASSERT_OP(geoCells1_.size(),==,tgrid.numCells(1));
    CGP_ASSERT_OP(geoCells2_.size(),==,tgrid.numCells(2));
    // set up bounds and bounded by

    // iterate over all 0-cells / junctions
    for(size_t cell0Index=0;cell0Index<tgrid.numCells(0);++cell0Index){
        const LabelType cell0Label=cell0Index+1;
        // set up label
        geoCells0_[cell0Index].setLabel(cell0Label);
        // get coordinates
        const size_t tx=geoCells0_[cell0Index][0][0];
        const size_t ty=geoCells0_[cell0Index][0][1];
        // Loop over all possible Cell1's / boundaries of the Cell0 / Junction
        const int px[]={ 1, -1, 0, 0};
        const int py[]={ 0,  0, 1,-1};
        for(size_t b=0;b<4;++b){
            LabelType cell1Label=grid(int(tx)+px[b],int(ty)+py[b]);
            // check if Cell1 / boundary is active
            if(cell1Label!=0){

                CGP_ASSERT_OP(cell1Label,>,0)
                CGP_ASSERT_OP(cell1Label,<=,tgrid.numCells(1));

                LabelType cell1Index=cell1Label-1;
                // bounds ( boundaries of a juction)
                geoCells0_[cell0Index].push_back_bound(cell1Label);
                // junctions of a boundaty
                geoCells1_[cell1Index].push_back_bounded_by(cell0Label);
            }
        }
        CGP_ASSERT_OP(geoCells0_[cell0Index].bounds().size(),>=,3);
        CGP_ASSERT_OP(geoCells0_[cell0Index].bounds().size(),<=,4);
    }

    // iterate over all 1-cells / boundaries
    for(size_t cell1Index=0;cell1Index<tgrid.numCells(1);++cell1Index){
        const LabelType cell1Label=cell1Index+1;
    
        // set up label
        geoCells1_[cell1Index].setLabel(cell1Label);
        // get tx and ty of SOME element of the boundary (the first in this case) 
        const size_t tx=geoCells1_[cell1Index][0][0];
        const size_t ty=geoCells1_[cell1Index][0][1];
        // bounds (region labels)
        LabelType cell2LabelA,cell2LabelB;
        // vertical boundary
        if(tx%2==1){
            cell2LabelA=static_cast<LabelType>(grid(tx-1,ty));
            cell2LabelB=static_cast<LabelType>(grid(tx+1,ty));

        }
        else{
            cell2LabelA=static_cast<LabelType>(grid(tx,ty-1));
            cell2LabelB=static_cast<LabelType>(grid(tx,ty+1));
        }
        CGP_ASSERT_OP(cell2LabelA,>,0)
        CGP_ASSERT_OP(cell2LabelA,<=,tgrid.numCells(2));
        CGP_ASSERT_OP(cell2LabelB,>,0)
        CGP_ASSERT_OP(cell2LabelB,<=,tgrid.numCells(2));
        const LabelType cell2IndexA=cell2LabelA-1;
        const LabelType cell2IndexB=cell2LabelB-1;

        // set up bounds (the 2 adj. regions to this boundary)
        geoCells1_[cell1Index].push_back_bound(cell2LabelA);
        geoCells1_[cell1Index].push_back_bound(cell2LabelB);
        // set up bounded by ( n adj. boundaries of a region)
        geoCells2_[cell2IndexA].push_back_bounded_by(cell1Label);
        geoCells2_[cell2IndexB].push_back_bounded_by(cell1Label);

        geoCells1_[cell1Index].sortCells();
    }
    // sortAdjaceny

    // iterate over all 2-cells / regions 
    for(size_t cell2Index=0;cell2Index<tgrid.numCells(2);++cell2Index){
        // set up ptr
        geoCells2_[cell2Index].cgp_=this;
        // set up label
        geoCells2_[cell2Index].setLabel(cell2Index+1);
        // sortAdjaceny
        geoCells2_[cell2Index].sortAdjaceny();
    }
    // iterate over all 1-cells / boundaries
    for(size_t cell1Index=0;cell1Index<tgrid.numCells(1);++cell1Index){
        // set up ptr
        geoCells1_[cell1Index].cgp_=this;
        // sortAdjaceny// sortAdjaceny
        geoCells1_[cell1Index].setLabel(cell1Index+1);
        geoCells1_[cell1Index].sortAdjaceny();
    }
    // iterate over all 0-cells / junctions
    for(size_t cell0Index=0;cell0Index<tgrid.numCells(0);++cell0Index){
        // set up ptr
        geoCells0_[cell0Index].cgp_=this;
        // sortAdjaceny
        geoCells0_[cell0Index].setLabel(cell0Index+1);
        // sortAdjaceny
        geoCells0_[cell0Index].sortAdjaceny();
    }

}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::TopologicalGridType & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::tgrid()const{
    return tgrid_;
}
template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::GeoCells0  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry0()const{
    return geoCells0_;
}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::GeoCells1  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry1()const{
    return geoCells1_;
}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::GeoCells2  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry2()const{
    return geoCells2_;
}

} /* namespace vigra */

#endif /* CGP2D_HXX */