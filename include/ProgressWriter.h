#ifndef ProgressWriter_H
#define ProgressWriter_H

#ifdef _WINDOWS
#include <windows.h>
#endif // _WINDOWS

#include "ExportMacros.h"
#include <memory>
#include <string>

class MMB_EXPORT AbstractProgressWriter {
public:
    enum class State {
        NOT_STARTED,
        PREPARING,
        RUNNING,
        FINISHED,
        FAILED
    };

    virtual ~AbstractProgressWriter();

    virtual void setTotalSteps(const int total) = 0;
    virtual void update(const State s) = 0;
    virtual void update(const State s, const int step) = 0;
};

class MMB_EXPORT DummyProgressWriter : public AbstractProgressWriter {
public:
    // Dummy writer that does absolutely nothing
    void setTotalSteps(const int total) override;
    void update(const State s) override;
    void update(const State s, const int completed) override;
};

class MMB_EXPORT ProgressWriter : public AbstractProgressWriter {
public:
    explicit ProgressWriter(const std::string &path);
    ~ProgressWriter();
    void setTotalSteps(const int total) override;
    void update(const State s) override;
    void update(const State s, const int step) override;

private:
    void write(const bool wait);

#ifdef _WINDOWS
    std::string _path;
#else
    int _output;
#endif // _WINDOWS

    State _state;
    int _step;
    int _totalSteps;
};

class MMB_EXPORT GlobalProgressWriter {
public:
    static void close();
    static void initialize(const std::string &path);
    static AbstractProgressWriter & get();
    static bool isInitialized();

private:
    static std::unique_ptr<AbstractProgressWriter> _writer;
};

#endif // ProgressWriter_H

