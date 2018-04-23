#ifndef THoughTransform_hxx_seen
#define THoughTransform_hxx_seen

#include <cmath>
#include <iostream>
#include <TVector3.h>
#include <list>
#include <functional>
#include <algorithm>
#include <vector>
#include "TH2F.h"
#include "TPad.h"

#define PI 3.14159265

/// The implementation of Hough Transformation which takes vector of std::pairs.
/// First member in pair correspond to X coodr and second to Y
/// Hough Transform X,Y coord in to Ther,Rad polar space and find point with max contribution, that correspond to the line

namespace CP {
  template <class DataPointer> class THoughCell;
  template <class DataPointer> class THoughTrans;
};

/// Hough Cell class is disabled right now, It can be implementerd if
/// more parameters for each point need to be stored 

template<class DataPointer>
class CP::THoughCell{
  public:
  typedef std::vector<DataPointer> Points;

  THoughCell(){}
  
  virtual ~THoughCell(){}

  THoughCell(int number1, int number2): fN1(number1),fN2(number2){}
  
  int GetIndex1() const {return fN1;}

  int GetIndex2() const {return fN2;}

  const Points& GetPoints() const {return fPoints;}

  std::size_t GetPointsCount(){return fPoints.size();}

  void AddPoint(const DataPointer& point){
    fPoints.push_back(point);
  }

  
  bool operator ==(const THoughCell<DataPointer>& cell) const {
    return cell.GetIndex1() == fN1 && cell.GetIndex2() == fN2;
      }

  
private:
  
  int fN1;
  
  int fN2;
  
  Points fPoints;
 
  
};


/// Actual Hough Transformtation Template 
template <class DataPointer>
class CP::THoughTrans{
public:

  typedef std::vector<CP::THoughCell<DataPointer>&> Cells;
  
    THoughTrans();
  
  virtual ~THoughTrans(){
    if(fHoughHist) delete fHoughHist;
  }

  ///COnstruct Hough Object. Maximum range of theta and radius can be set manually here
  /// The number of bins for each polar axis is provided as an input for the constructtor in main() 
  
  THoughTrans(double angleBin, double radBin) : fAngleBin(angleBin),fRadBin(radBin),fMaxThet(360.0),fMaxRad(5000.0){}

  /// Main Algorithm

  template <class DataIterator>
  void HoughTransform(DataIterator begin, DataIterator end);

  /// Take parameters of found Hough Line pair.fist is a and pair.second is b, ther y=a*x+b

  std::pair<double,double> GetLineParam(int limit);

  ///Get Filled histogram for the Hough transform

  TH2F* GetHist(){return fHoughHist;}
  

private:

  double fAngleBin;

  double fRadBin;

  double fMaxRad;

  double fMaxThet;

  TH2F* fHoughHist;

};



///*********************************************
///Implenmentation of the Algorithm
///*********************************************



template<class DataPointer>
template<class DataIterator>
void CP::THoughTrans<DataPointer>::HoughTransform(DataIterator begin, DataIterator end){

fHoughHist = new TH2F("HH","HH",fAngleBin,-90,90,fRadBin,-fMaxRad,fMaxRad);
  
  for(DataIterator p=begin; p != end;++p){

    //Metric definition for input points
    
    double x = (*p).first;
    double y = (*p).second;

 
    for(int i=1;i<fAngleBin+1;++i){
      double thet=fHoughHist->GetXaxis()->GetBinCenter(i);
      double rad = x*sin(thet*PI/180.0)+y*cos(thet*PI/180.0);
      //  std::cout<<rad<<std::endl;
      int ibin = fHoughHist->FindFixBin(thet,rad);
      int j=fHoughHist->GetYaxis()->FindFixBin(rad);

      //Make Wide Lines
       if(i==1){
	if(j==1){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j>1 && j<fRadBin){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j==fRadBin){
	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);
	}	
      }
        if(i>1 && i<fAngleBin){
	if(j==1){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j>1 && j<fRadBin){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);
	  
	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j==fRadBin){
	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i+1,j);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);
	}
      }
	 if(i==fAngleBin){
	if(j==1){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j>1 && j<fRadBin){
	  ibin = fHoughHist->GetBin(i,j+1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);
	}
	if(j==fRadBin){
	  ibin = fHoughHist->GetBin(i,j-1);
	  fHoughHist->AddBinContent(ibin);

	  ibin = fHoughHist->GetBin(i-1,j);
	  fHoughHist->AddBinContent(ibin);
	}	
	}
    }
  }

}

template<class DataPointer>
std::pair<double,double> CP::THoughTrans<DataPointer>::GetLineParam(int limit){

int maxi=-9999;
  int maxj=-9999;
  int maxW=-9999;
  std::pair<double,double> line_bad = std::make_pair(-9999,-9999);
  for(int i=1;i<fAngleBin+1;++i){
    for(int j=1;j<fRadBin+1;++j){
      int ibin = fHoughHist->GetBin(i,j);
      int weight = fHoughHist->GetBinContent(ibin);
      if(weight>maxW){maxW=weight;maxi=i;maxj=j;}
    }
  }
  if(maxW>limit){
   
  //y=a*x+b <- y = rad/cos(thet)-x*sin(thet)/cos(thet)
    double thet = fHoughHist->GetXaxis()->GetBinCenter(maxi);
  double rad = fHoughHist->GetYaxis()->GetBinCenter(maxj);

  double a = -sin(thet*PI/180.0)/cos(thet*PI/180.0);
  double b = rad/cos(thet*PI/180.0);
  std::pair<double,double> line = std::make_pair(a,b);
  return line;}
  return line_bad;
}




#endif
