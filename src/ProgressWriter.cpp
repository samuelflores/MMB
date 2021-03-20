#include "ProgressWriter.h"
#include "Constexpr.h"
#include "Impossible.h"
#include "MMBLogger.h"

#ifndef _WINDOWS
#include <unistd.h>
#include <fcntl.h>
#endif // _WINDOWS

static
MMB_CONSTEXPR
const char * stateToStr(const ProgressWriter::State s) noexcept {
    switch (s) {
    case ProgressWriter::State::NOT_STARTED:
        return "NotStarted";
    case ProgressWriter::State::PREPARING:
        return "Preparing";
    case ProgressWriter::State::RUNNING:
        return "Running";
    case ProgressWriter::State::FINISHED:
        return "Finished";
    case ProgressWriter::State::FAILED:
        return "Failed";
    }

    __IMPOSSIBLE__;
}

#ifndef _WINDOWS
class Locker {
public:
    Locker(const int fd, const bool wait) :
        _fd(fd) {

        flock lk{};
        lk.l_type = F_WRLCK;
        lk.l_whence = SEEK_SET;
        lk.l_start = 0;
        lk.l_len = 0;

        if (wait)
            lock(lk);
        else
            tryLock(lk);
    }

    ~Locker() {
        if (!_haveLock)
            return;

        flock lk{};
        lk.l_type = F_UNLCK;
        lk.l_whence = SEEK_SET;
        lk.l_start = 0;
        lk.l_len = 0;

        if (fcntl(_fd, F_SETLK, &lk) != 0)
            MMBLOG_PLAIN(CRITICAL, "Failed to release progress file lock");
    }

    bool haveLock() const { return _haveLock; }
private:
    void tryLock(flock &lk) {
        const int ret = fcntl(_fd, F_SETLK, &lk);
        if (ret == -1) {
            if (errno == EACCES || errno == EAGAIN)
                _haveLock = false;
            else
                MMBLOG_PLAIN(CRITICAL, "Error while trying to acquire progress file lock\n");
        } else
            _haveLock = true;
    }

    void lock(flock &lk) {
        const int ret = fcntl(_fd, F_SETLKW, &lk);
        if (ret == -1) {
            if (errno == EINTR)
                _haveLock = false;
            else
                MMBLOG_PLAIN(CRITICAL, "Error while trying to acquire progress file lock\n");
        } else
            _haveLock = true;
    }

    const int _fd;
    bool _haveLock;
};
#endif // _WINDOWS

AbstractProgressWriter::~AbstractProgressWriter() {}

void DummyProgressWriter::setTotalSteps(const int total) {
    // NOOP
}

void DummyProgressWriter::update(const State s) {
    // NOOP
}

void DummyProgressWriter::update(const State s, const int step) {
    // NOOP
}

ProgressWriter::ProgressWriter(const std::string &path) :
    _totalSteps(0) {
#ifdef _WINDOWS
    _path = path;
    HANDLE hOutput = CreateFileA(_path.c_str(), GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hOutput == INVALID_HANDLE_VALUE)
        MMBLOG_PLAIN(CRITICAL, "Cannot open progress file for writing");
    CloseHandle(hOutput);
#else
    _output = creat(path.c_str(), S_IRUSR | S_IWUSR);
    if (_output < 1)
        MMBLOG_PLAIN(CRITICAL, "Cannot open progress file for writing");
#endif // _WINDOWS

    ProgressWriter::update(State::NOT_STARTED, 0);
}

ProgressWriter::~ProgressWriter() {
    write(true);
}

void ProgressWriter::setTotalSteps(const int total) {
    _totalSteps = total;
}

void ProgressWriter::update(const State s) {
    _state = s;

    write(false);
}

void ProgressWriter::update(const State s, const int step) {
    assert(step <= _totalSteps);

    _step = step;

    ProgressWriter::update(s);
}

void ProgressWriter::write(const bool wait) {
#ifdef _WINDOWS
    int len = snprintf(
        nullptr,
        0,
        "{\"state\":\"%s\",\"step\": %d,\"total_steps\":%d}\n",
        stateToStr(_state),
        _step,
        _totalSteps
    );
    if (len < 1) {
        MMBLOG_PLAIN(WARNING, "Failed to create progress report\n");
        return;
    }

    auto text = std::unique_ptr<char[]>(new char[len + 1]);
    int len2 = snprintf(
        text.get(),
        len,
        "{\"state\":\"%s\",\"step\": %d,\"total_steps\":%d}\n",
        stateToStr(_state),
        _step,
        _totalSteps
    );

    if (len2 < 1 || len2 >= len) {
        MMBLOG_PLAIN(WARNING, "Failed to create progress report\n");
        return;
    }

    HANDLE hOutput = CreateFileA(_path.c_str(), GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hOutput == INVALID_HANDLE_VALUE) {
        MMBLOG_PLAIN(WARNING, "Cannot open progress file for writing");
        return;
    }

    DWORD written;
    if (WriteFile(hOutput, text.get(), len, &written, NULL) == FALSE)
        MMBLOG_PLAIN(WARNING, "Cannot write progress to file");
    CloseHandle(hOutput);
#else
    char *text = nullptr;
    int ret = asprintf(
        &text,
        "{\"state\":\"%s\",\"step\": %d,\"total_steps\":%d}\n",
        stateToStr(_state),
        _step,
        _totalSteps
    );

    if (ret < 1) {
        MMBLOG_PLAIN(WARNING, "Failed to create progress report\n");
        return;
    }

    Locker lk(_output, wait);
    if (lk.haveLock()) {
        ftruncate(_output, 0);
	lseek(_output, 0, SEEK_SET);
        if (::write(_output, text, ret) != ret)
            MMBLOG_PLAIN(WARNING, "Failed to write progress report\n");
    }

    free(text);
#endif // _WINDOWS
}

std::unique_ptr<AbstractProgressWriter> GlobalProgressWriter::_writer{nullptr};

void GlobalProgressWriter::close() {
    _writer = nullptr;
}

void GlobalProgressWriter::initialize(const std::string &path) {
    if (path.length() == 0)
        _writer = std::make_unique<DummyProgressWriter>();
    else
        _writer = std::make_unique<ProgressWriter>(path);
}

AbstractProgressWriter & GlobalProgressWriter::get() {
    assert(_writer != nullptr);
    return *_writer;
}

bool GlobalProgressWriter::isInitialized() {
    return _writer != nullptr;
}
