#ifndef HASH_H
#define HASH_H

#if defined (__GNUC__) && ( __GNUC__ <= 2 )
#include <hash_map>
using std::hash_map;
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif

#endif //HASH_H
