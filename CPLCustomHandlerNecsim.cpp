// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file CPLCustomHandlerNecsim.cpp
 * @brief Contains a custom CPL error handler.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include <sstream>
#include "CPLCustomHandlerNecsim.h"
#include "Logging.h"

#ifdef with_gdal

void cplNecsimCustomErrorHandler(CPLErr eErrClass, int err_no, const char *msg)
{
	stringstream error_msg;
	error_msg << "Gdal error: " << err_no << ". " << msg << endl;
	if(eErrClass == CE_Fatal)
	{
		writeCritical(error_msg.str());
	}
	else if(eErrClass == CE_Failure)
	{
		writeError(error_msg.str());
	}
	else if(eErrClass == CE_Warning)
	{
		writeWarning(error_msg.str());
	}
#ifdef DEBUG
	else
	{
		writeLog(10, error_msg.str());
	}
#endif // DEBUG
}

#endif //with_gdal