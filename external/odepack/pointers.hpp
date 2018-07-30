#pragma once

template <class T>
struct lambda_traits : lambda_traits<decltype(&T::operator())>
{ };

template <class T, class R, class... Args>
struct lambda_traits<R(T::*)(Args...) const> {
    typedef R (*pointer)(Args...);

    static pointer function_ptr(T t) {
        static T fn = t;
        return [](Args... args) {
            return fn(args...);
        };
    }
};

template <class T>
inline typename lambda_traits<T>::pointer function_ptr(T t) {
    return lambda_traits<T>::function_ptr(t);
}
