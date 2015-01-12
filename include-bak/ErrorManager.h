/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-13 by the Author.                                       *
 * Authors: Samuel Flores, Alex Tek                                           *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to use the Software for non-profit academic research only. Any commercial, *
 * for-profit, or industrial use requires the payment of usage fees.        *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#ifndef ERRORMANAGER_H_
#define ERRORMANAGER_H_

#include <string>
#include <iostream>
#include <sstream>
#include <exception>
#include <streambuf>

#include "ExportMacros.h" 

/**
 * \brief Simple MMB exception. 
 * 
 * Store error messages in a string
 * 
 */
class MMB_EXPORT MMBException : public std::exception
{
public:
    MMBException( const std::string & message)
    {
        std::stringstream ss;
        ss << "MMB_ERROR: " << message;
        this->_message = ss.str();
    }
 
    virtual ~MMBException() throw()
    {
 
    }
 
    virtual const char * what() const throw()
    {
        return this->_message.c_str();
    }
 
private:
    std::string _message;
};


/**
 * \brief Write the error messages in an output stream (std::cerr by default). Treating the error throws a MMBexception. 
 *
 * This was implmentated to keep the program alive and propagate errors when used as a library,
 * while it still exits for the command line version of MMB.
 * Errors can be redirected to a file too if needed.
 *
 * 
 * How-to use:
 * A global instance is availaible as ErrorManager::instance
 * (No singleton, to keep it simple and have potentially different Managers (e.g. for logs))
 * 
 * ErrorManager::instance is used like a stream :
 * ErrorManager::instance << "error message" << endl;
 * By default, it redirects to std::cerr
 * ErrorManager::instance.treatError() will throw a MMBException containing the last error messages
 * ErrorManager::setInstanceOutputStream(std:ostream&) can be used to change the output stream
 */
class MMB_EXPORT ErrorManager {

/**
* << operator overload to redirect val to the member output stream
* and fill a buffer to reconstruct the entire message afterward
*/
template<class T> 
friend ErrorManager& operator<< (ErrorManager& em, T val){
    *(em._output) << val;
    em._buffer << val;
    return em;
}

protected:

    /**
    * output stream
    */
    std::ostream * _output;

    /**
    * stream buffer backup
    */
    std::streambuf * _backup;

    /**
    * message buffer
    */
    std::stringstream _buffer;
    
public:
    /**
     * Initialize the instance with the output stream provided.
     * std::cerr is used by default.
     *
     * \param a pointer to the stream to use
     */
    ErrorManager(std::ostream & output = std::cerr);

    /**
     * Copy constructor. Necessary to properly copy the stringstream member
     *
     * \param a reference to the ErrorManager to copy
     */
    ErrorManager(const ErrorManager & em);

    /**
     * = operator overload. Necessary to properly copy the stringstream member
     *
     * \param a reference to the ErrorManager to copy
     */
    ErrorManager& operator=(const ErrorManager & em);

    void errorLog(const std::string & message)
    {
        *this << message;
    }

    /**
     * Destructor
     * Set the output back to the original output
    */
    ~ErrorManager();

    /**
    * Propagate the error message through a MMBException
    */
    void treatError();

    /**
    * << overload for stream manipulators (like std::endl)
    */
    ErrorManager& operator<<(std::ostream& (*pfun)(std::ostream&));

    /**
     * Recreate the global instance with a redirection to dest
     * @param dest the new output stream
     */
    static void setInstanceOutputStream(std::ostream& dest);

    /**
     * Global ErrorManager
     */
    static ErrorManager instance;

};


#endif
