#include "evidence.h"
#include "graph/node.h"

#include <sstream>
#include <string>
#include <cstdlib>

#include <dai/util.h>
#include <boost/lexical_cast.hpp>

namespace pdbntk {

void Evidence::addEvidenceTabFile( std::istream &is, FactorGraph &fg ) {
  std::map<std::string, Node*> varMap;
  for (std::vector<Node*>::const_iterator v = fg.nodes().begin(); v != fg.nodes().end(); ++v ) {
    std::stringstream s;
    s << (*v)->index();
    varMap[s.str()] = *v;
  }

  addEvidenceTabFile(is, varMap);
}


void Evidence::addEvidenceTabFile(std::istream &is, std::map<std::string, Node*> &varMap) {
  using std::vector;
  std::string line;
  getline( is, line );
  size_t line_number = 2;

  // Parse header
  std::vector<std::string> header_fields;
  header_fields = dai::tokenizeString(line, true);
  std::vector<std::string>::const_iterator p_field = header_fields.begin();
  if( p_field == header_fields.end() )
    DAI_THROWE(INVALID_EVIDENCE_FILE,"Empty header line");

  std::vector<Node*> nodes;
  for( ; p_field != header_fields.end(); ++p_field ) {
    std::map<std::string, Node*>::iterator elem = varMap.find( *p_field );
    if( elem == varMap.end() )
      DAI_THROWE(INVALID_EVIDENCE_FILE,"Variable " + *p_field + " not known");
    nodes.push_back(elem->second);
  }

  getline(is,line);
  if( is.fail() || line.size() > 0 )
    DAI_THROWE(INVALID_EVIDENCE_FILE,"Expecting empty line");

  // Read samples
  while( getline(is, line) ) {
    line_number++;

    std::vector<std::string> fields;
    fields = dai::tokenizeString( line, true, "\t" );
    if( fields.size() != nodes.size() )
      DAI_THROWE(INVALID_EVIDENCE_FILE,"Invalid number of fields in line " + boost::lexical_cast<std::string>(line_number));

    Observation sample;
    for( size_t i = 0; i < nodes.size(); ++i ) {
      if( fields[i].size() > 0 ) { // skip if missing observation
        if( fields[i].find_first_not_of("0123456789") != std::string::npos )
          DAI_THROWE(INVALID_EVIDENCE_FILE,"Invalid state " + fields[i] + " in line " + boost::lexical_cast<std::string>(line_number));
      }
    }
    _samples.push_back( sample );
  } // finished sample line
}


} // end of namespace dai
