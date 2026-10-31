// function_ref.h
// Non-owning reference to a callable with a fixed signature: one pointer to the
// callable plus one call thunk, no allocation. A driver that takes a
// FunctionRef instead of a template parameter compiles once, while each caller
// keeps its own (capturing) lambda. The referenced callable must outlive the
// FunctionRef; only lvalues bind, so a temporary cannot be captured by mistake.

#ifndef TULPAOBS_FUNCTION_REF_H
#define TULPAOBS_FUNCTION_REF_H

#include <memory>
#include <type_traits>
#include <utility>

namespace tulpaObs {

template <class Sig> class FunctionRef;

template <class R, class... Args>
class FunctionRef<R(Args...)> {
public:
    template <class F,
              class = std::enable_if_t<!std::is_same<std::decay_t<F>, FunctionRef>::value>>
    FunctionRef(F& f) noexcept
        : obj_(const_cast<void*>(static_cast<const void*>(std::addressof(f)))),
          call_([](void* obj, Args... args) -> R {
              return (*static_cast<F*>(obj))(std::forward<Args>(args)...);
          }) {}

    R operator()(Args... args) const {
        return call_(obj_, std::forward<Args>(args)...);
    }

private:
    void* obj_;
    R (*call_)(void*, Args...);
};

}  // namespace tulpaObs

#endif  // TULPAOBS_FUNCTION_REF_H
