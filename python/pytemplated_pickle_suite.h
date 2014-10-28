#ifndef PYTEMPLATED_PICKLE_SUITE_H
#define PYTEMPLATED_PICKLE_SUITE_H

#include <boost/serialization/serialization.hpp>
#include <boost/python.hpp>

// Templated pickle suite
template<class T>
struct TemplatedPickleSuite : boost::python::pickle_suite {
  static std::string getstate( const T& g ) {
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa & g;
    return ss.str();
  }

  static void setstate( T& g, const std::string& state ) {
    std::stringstream ss(state);
    boost::archive::text_iarchive ia(ss);
    ia & g;
  }
};

#endif // PYTEMPLATED_PICKLE_SUITE_H
