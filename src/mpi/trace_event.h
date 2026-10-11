#pragma once

#include <string>

#ifdef MPETRACE
#include <mpi.h>
#include <mpe.h>
#endif

#ifdef ITACTRACE
#include <VT.h>
#endif

class TraceEvent {
public:
    TraceEvent(const std::string &name) {
#ifdef MPETRACE
        MPE_Log_get_state_eventIDs(&_begin, &_end);
        MPE_Describe_state(_begin, _end, name.c_str(), "red");
#endif

#ifdef ITACTRACE
        VT_classdef("User", &class_handle);
        VT_funcdef(name.c_str(), class_handle, &state_handle);
#endif
    }

    void begin() {
#ifdef MPETRACE
        MPE_Log_event(_begin, 0, NULL);
#endif

#ifdef ITACTRACE
        VT_begin(state_handle);
#endif
    }

    void end() {
#ifdef MPETRACE
        MPE_Log_event(_end, 0, NULL);
#endif

#ifdef ITACTRACE
        VT_end(state_handle);
#endif
    }

private:
#ifdef MPETRACE
    int _begin, _end;
#endif

#ifdef ITACTRACE
    int class_handle, state_handle;
#endif
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
