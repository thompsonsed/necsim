// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
/**
 * @author Sam Thompson
 * @file CPLCustomHandler.h
 * @brief Contains a custom CPL error handler.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
#ifndef SPECIATIONCOUNTER_CPLCUSTOMHANDLER_H
#define SPECIATIONCOUNTER_CPLCUSTOMHANDLER_H
#ifdef with_gdal
#include <cpl_error.h>
/**
 * @brief Contains a custom CPLErrorHandler for reporting errors from GDAL imports.
 * @param eErrClass the error class (CE_None, CE_Debug, CE_Warning, CE_Failure or CE_Fatal)
 * @param err_no the error number to report
 * @param msg the message to report
 */
void cplCustomErrorHandler(CPLErr eErrClass, int err_no, const char *msg);

#endif // with_gdal


#endif //SPECIATIONCOUNTER_CPLCUSTOMHANDLER_H
