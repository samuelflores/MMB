#ifndef ProgressWriter_H
#define ProgressWriter_H

#include <memory>
#include <string>

class AbstractProgressWriter {
protected:
    enum class State {
        NOT_STARTED,
        PREPARING,
        RUNNING,
        FINISHED,
        FAILED
    };


public:
    virtual ~AbstractProgressWriter();

    virtual void setTotalSteps(const int total) = 0;
    virtual void update(const State s) = 0;
    virtual void update(const State s, const int step) = 0;
};

class DummyProgressWriter : public AbstractProgressWriter {
public:
    // Dummy writer that does absolutely nothing
    void setTotalSteps(const int total) override;
    void update(const State s) override;
    void update(const State s, const int completed) override;
};

class ProgressWriter : public AbstractProgressWriter {
public:
    using AbstractProgressWriter::State;

    explicit ProgressWriter(const std::string &path);
    ~ProgressWriter();
    void setTotalSteps(const int total) override;
    void update(const State s) override;
    void update(const State s, const int step) override;

private:
    void write(const bool wait);

    int _output;

    State _state;
    int _step;
    int _totalSteps;
};

class GlobalProgressWriter {
public:
    static void initialize(const std::string &path);
    static AbstractProgressWriter & get();
    static bool isInitialized();

private:
    static std::unique_ptr<AbstractProgressWriter> _writer;
};

#endif // ProgressWriter_H

