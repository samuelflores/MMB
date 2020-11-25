#include "MMBLogger.h"
#include "Impossible.h"
#include "ProgressWriter.h"

#include <iostream>

static std::mutex initMutex;

inline
constexpr
const char * msgPrefix(const MMBLogger::Severity severity) {
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

    __IMPOSSIBLE__;
}

MMBLogger * MMBLogger::s_me(nullptr);

MMBLogger::MMBLogger() :
    _loggingSeverity(Severity::INFO),
    _output(nullptr),
    _newlinesSinceFlush(0) {
}

MMBLogger & MMBLogger::instance() {
    std::lock_guard<std::mutex> lk(initMutex);
    if (s_me == nullptr) {
        s_me = new MMBLogger();
        s_me->setOutput(&std::cout);
    }

    return *s_me;
}

void MMBLogger::log(const Severity severity, const std::ostringstream& oss, const bool printSeverity) {
    assert(_output != nullptr);

    const std::string msg = (printSeverity ? msgPrefix(severity) : "") + oss.str();

    if (severity >= _loggingSeverity) {
        std::lock_guard<std::mutex> lk{_writeMutex};
        (*_output) << msg;
        maybeFlush(msg);
    }
}

#ifndef MMBLOG_DONT_THROW_ON_CRITICAL
void MMBLogger::logCritical [[noreturn]] (const std::ostringstream& oss) {
    assert(_output != nullptr);

    if (GlobalProgressWriter::isInitialized())
        GlobalProgressWriter::get().update(ProgressWriter::State::FAILED);

    const std::string msg = msgPrefix(Severity::CRITICAL) + oss.str();

    std::lock_guard<std::mutex> lk{_writeMutex};
    (*_output) << msg;
    _output->flush();
    throw MMBException(msg);
}
#endif // MMBLOG_DONT_THROW_ON_CRITICAL

void MMBLogger::maybeFlush(const std::string &msg) {
    for (const auto &ch : msg)
        _newlinesSinceFlush += size_t(ch == '\n');

    if (_newlinesSinceFlush > 5) {
        _output->flush();
        _newlinesSinceFlush = 0;
    }
}

void MMBLogger::setLoggingSeverity(const Severity severity) {
    _loggingSeverity = severity < Severity::ALWAYS ? severity : Severity::ALWAYS;
}

void MMBLogger::setOutput(std::ostream *output) {
    std::lock_guard<std::mutex> lk(_writeMutex);

    if (_output != nullptr)
        _output->flush();
    _output = output;
}
