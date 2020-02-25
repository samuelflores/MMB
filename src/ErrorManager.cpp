
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

