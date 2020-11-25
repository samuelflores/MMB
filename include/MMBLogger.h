#ifndef MMBLOGGER_H_
#define MMBLOGGER_H_

#include <cassert>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

class MMBException : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
    using std::runtime_error::what;
};

class MMBLogger {
public:
    enum class Severity : uint8_t {
        DEBUG,
        INFO,
        WARNING,
        ALWAYS = 254,     /*!< Non-critical message that is always displayed */
        CRITICAL = 255    /*!< Unrecoverable error message. Logging the message also terminates the program */
    };

    static MMBLogger & instance();

    void log(const Severity severity, const std::ostringstream& oss, const bool printSeverity = true);

    #ifndef MMBLOG_DONT_THROW_ON_CRITICAL
    void logCritical [[noreturn]] (const std::ostringstream& oss);
    #endif // MMBLOG_DONT_THROW_ON_CRITICAL

    void setOutput(std::ostream *output);
    void setLoggingSeverity(const Severity severity);

private:
    MMBLogger();

    Severity _loggingSeverity;
    std::ostream * _output;
    std::mutex _writeMutex;

    static MMBLogger *s_me;
};

template <MMBLogger::Severity S>
struct MMBLoggerDispatcher {
    static void call(const std::ostringstream &oss, const bool printSeverity = true) {
        MMBLogger::instance().log(S, oss, printSeverity);
    }
};

#ifndef MMBLOG_DONT_THROW_ON_CRITICAL
template <>
struct MMBLoggerDispatcher<MMBLogger::Severity::CRITICAL> {
    static void call [[noreturn]] (const std::ostringstream &oss) {
        MMBLogger::instance().logCritical(oss);
    }
};
#endif // MMBLOG_DONT_THROW_ON_CRITICAL

#define MMBLOG_PLAIN(sev, msg) \
    do { \
        std::ostringstream oss{}; \
        oss << msg; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(oss); \
    } while (false)

#define MMBLOG_PLAIN_NOSEV(sev, msg) \
    do { \
        std::ostringstream oss{}; \
        oss << msg; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(oss, false); \
    } while (false)

#define MMBLOG_FILE_LINE(sev, msg) \
    do { \
        std::ostringstream oss{}; \
        oss << __FILE__ << ":" << ":" << __LINE__ << ": " << msg; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(oss); \
    } while (false)

#define MMBLOG_FILE_FUNC_LINE(sev, msg) \
    do { \
        std::ostringstream oss{}; \
        oss << __FILE__ << ":" << __FUNCTION__ << ":" << __LINE__ << ": " << msg; \
        MMBLoggerDispatcher<MMBLogger::Severity::sev>::call(oss); \
    } while (false)

#endif // MMBLOGGER_H_
