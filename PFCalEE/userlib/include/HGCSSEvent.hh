#ifndef _hgcssevent_hh_
#define _hgcssevent_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

class HGCSSEvent{


public:
  HGCSSEvent():
    event_(0)
  {
    
  };

  ~HGCSSEvent(){};

  inline unsigned eventNumber() const{
    return event_;
  };

  inline void eventNumber(const unsigned aNum){
    event_ = aNum;
  };

private:

  unsigned event_;

  ClassDef(HGCSSEvent,1);



};


#endif
