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

#include <cstdlib>


#include "ErrorManager.h"


using namespace std;

ErrorManager ErrorManager::instance;

ErrorManager::ErrorManager(std::ostream & output) : _output(&output), _backup(output.rdbuf()){

}

ErrorManager::ErrorManager(const ErrorManager & em) : _output(em._output), _backup(em._backup){
    _buffer << em._buffer.rdbuf();
}

ErrorManager& ErrorManager::operator=(const ErrorManager & em)
{
    _output = em._output; 
    _backup = em._backup;
    _buffer << em._buffer.rdbuf();
}

ErrorManager::~ErrorManager(){
    _output->rdbuf(_backup);
}

void ErrorManager::treatError(){
    string str = _buffer.str();
    _buffer.flush();
    _buffer.str("");
    throw MMBException(str);    
}

ErrorManager& ErrorManager::operator<<(std::ostream& (*pfun)(std::ostream&)){
    pfun(*_output);
    return *this;
}


void ErrorManager::setInstanceOutputStream(std::ostream& dest){
    instance = ErrorManager(dest);
}

