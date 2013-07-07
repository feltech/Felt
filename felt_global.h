#ifndef FELT_GLOBAL_H
#define FELT_GLOBAL_H

#include <QtCore/qglobal.h>
#include <omp.h>

#if defined(FELT_LIBRARY)
#  define FELTSHARED_EXPORT Q_DECL_EXPORT
#else
#  define FELTSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // FELT_GLOBAL_H
