#pragma once

#include <string>

#ifdef MPETRACE
#include <mpi.h>
#include <mpe.h>
#endif

class TraceEvent {
public:
    TraceEvent(const std::string &name) {
#ifdef MPETRACE
        MPE_Log_get_state_eventIDs(&_begin, &_end);
        MPE_Describe_state(_begin, _end, name.c_str(), "red");
#endif
    }

    void begin() {
#ifdef MPETRACE
        MPE_Log_event(_begin, 0, NULL);
#endif
    }

    void end() {
#ifdef MPETRACE
        MPE_Log_event(_end, 0, NULL);
#endif
    }

private:
    int _begin, _end;
};

class TraceEventWrapper {
public:
    TraceEventWrapper(TraceEvent &event):
        event(event) {
        event.begin();
    }

    ~TraceEventWrapper() {
        event.end();
    }

private:
    TraceEvent &event;
};
