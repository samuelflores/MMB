#ifndef MMBLOGGER_H_
#define MMBLOGGER_H_

#include "ExportMacros.h"
#include <cassert>
#include <functional>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

using LogFunc = std::function<std::ostringstream ()>;

class MMB_EXPORT MMBException : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
    using std::runtime_error::what;
};

class MMB_EXPORT MMBLogger {
public:
    enum class Severity : uint8_t {
        DEBUG,
        INFO,
        WARNING,
        ALWAYS = 254,     /*!< Non-critical message that is always displayed */
        CRITICAL = 255    /*!< Unrecoverable error message. Logging the message also terminates the program */
    };

    static MMBLogger & instance();

    void flush();
    void log(const Severity severity, const LogFunc &logFunc, const bool printSeverity = true);

    #ifndef MMBLOG_DONT_THROW_ON_CRITICAL
    void logCritical [[noreturn]] (const LogFunc &oss);
    #endif // MMBLOG_DONT_THROW_ON_CRITICAL

    void setOutput(std::ostream *output);
    void setLoggingSeverity(const Severity severity);

private:
    MMBLogger();
    void maybeFlush(const std::string &msg);

    Severity _loggingSeverity;
    std::ostream * _output;
    std::mutex _writeMutex;
    size_t _newlinesSinceFlush;

    static MMBLogger *s_me;
};


template <MMBLogger::Severity S>
struct MMBLoggerDispatcher {
    static void call(const LogFunc &logFunc, const bool printSeverity = true) {
        MMBLogger::instance().log(S, logFunc, printSeverity);
    }
};

#ifndef MMBLOG_DONT_THROW_ON_CRITICAL
template <>
struct MMBLoggerDispatcher<MMBLogger::Severity::CRITICAL> {
    static void call [[noreturn]] (const LogFunc &logFunc) {
        MMBLogger::instance().logCritical(logFunc);
    }
};
#endif // MMBLOG_DONT_THROW_ON_CRITICAL

#define MMBLOG_PLAIN(sev, msg) \
    do { \
        auto f = [&]() { \
            std::ostringstream oss{}; \
            oss << msg; \
            return oss; \
        }; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(f); \
    } while (false)

#define MMBLOG_PLAIN_NOSEV(sev, msg) \
    do { \
        auto f = [&]() { \
            std::ostringstream oss{}; \
            oss << msg; \
            return oss; \
        }; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(f, false); \
    } while (false)

#define MMBLOG_FILE_LINE(sev, msg) \
    do { \
        auto f = [&]() { \
            std::ostringstream oss{}; \
            oss << __FILE__ << ":" << ":" << __LINE__ << ": " << msg; \
            return oss; \
        }; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(f); \
    } while (false)

#define MMBLOG_FILE_FUNC_LINE(sev, msg) \
    do { \
        const auto __FUNC_STR__ = __FUNCTION__; \
        auto f = [&]() { \
            std::ostringstream oss{}; \
            oss << __FILE__ << ":" << __FUNC_STR__ << ":" << __LINE__ << ": " << msg; \
            return oss; \
        }; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(f); \
    } while (false)

#endif // MMBLOGGER_H_
