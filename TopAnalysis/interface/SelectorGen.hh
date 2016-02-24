#ifndef _SelectorGen_hh_
#define _SelectorGen_hh__
#include "TLorentzVector.h"

class ColourFlowAnalysisTool;
class MiniEvent_t;

using namespace std;
class SelectorGen
{
  
public:
  
  const MiniEvent_t            * event_ptr_;
  ColourFlowAnalysisTool       * CFAT_;
  void Work();
}
