#include "MMBLogger.h"
#include "ProgressWriter.h"

#include <iostream>

static std::mutex initMutex;

inline
std::string msgPrefix(const MMBLogger::Severity severity) {
    switch (severity) {
    case MMBLogger::Severity::DEBUG:
        return "DEBUG: ";
    case MMBLogger::Severity::INFO:
        return "INFO: ";
    case MMBLogger::Severity::WARNING:
        return "WARNING: ";
    case MMBLogger::Severity::ALWAYS:
        return ""; // No prefix
    case MMBLogger::Severity::CRITICAL:
        return "CRITICAL: ";
    }

    assert("Invalid Severity level");
}

MMBLogger * MMBLogger::s_me(nullptr);

MMBLogger::MMBLogger() :
    _loggingSeverity(Severity::INFO),
    _output(nullptr)
{
}

MMBLogger & MMBLogger::instance() {
    std::lock_guard<std::mutex> lk(initMutex);
    if (s_me == nullptr) {
        s_me = new MMBLogger();
        s_me->setOutput(&std::cout);
    }

    return *s_me;
}

void MMBLogger::log(const Severity severity, const std::ostringstream& oss) {
    assert(_output != nullptr);

    const std::string msg = msgPrefix(severity) + oss.str();

    if (severity >= _loggingSeverity) {
        std::lock_guard<std::mutex> lk{_writeMutex};
        (*_output) << msg;
    }
}

#ifndef MMBLOG_DONT_THROW_ON_CRITICAL
void MMBLogger::logCritical [[noreturn]] (const std::ostringstream& oss) {
    assert(_output != nullptr);

    GlobalProgressWriter::get().update(ProgressWriter::State::FAILED);

    const std::string msg = msgPrefix(Severity::CRITICAL) + oss.str();

    std::lock_guard<std::mutex> lk{_writeMutex};
    (*_output) << msg;
    _output->flush();
    throw MMBException(msg);
}
#endif // MMBLOG_DONT_THROW_ON_CRITICAL

void MMBLogger::setLoggingSeverity(const Severity severity)
{
    _loggingSeverity = severity < Severity::ALWAYS ? severity : Severity::ALWAYS;
}

void MMBLogger::setOutput(std::ostream *output)
{
    std::lock_guard<std::mutex> lk(_writeMutex);

    if (_output != nullptr)
        _output->flush();
    _output = output;
}
