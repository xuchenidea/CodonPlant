#pragma once

//#define __DEBUG__CODON
#ifdef  __DEBUG__CODON
#define log(a,...) printf((a), __VA_ARGS__)
#else
#define log(a,...)
#endif