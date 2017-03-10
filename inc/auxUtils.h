#include "inc/RooStatsHead.h"
#include "inc/RooFitHead.h"
#include "inc/CommonHead.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

namespace auxUtils {

  std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters) {
    std::vector<std::string> tokens;

    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
  }

}
