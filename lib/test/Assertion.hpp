#include <iomanip>
#include <iostream>
#include <string>

class Assertion {
   private:
    int err = 0;
    int last_err = 0;

   public:
    Assertion() {}

    void operator()(bool t) {
        last_err = t ? 0 : 1;
        err = err == 0 && last_err == 0 ? 0 : 1;
    }
    void operator()(bool t, std::string msg) {
        operator()(t);
        if (last_status() != 0) { std::cout << msg << '\n'; }
    }
    template <typename T>
    void test(bool t, const T& msg) {
        operator()(t);
        std::cout << std::boolalpha << '[' << t << "] " << msg << '\n';
    }
    int last_status() const { return last_err; }
    int status() const { return err; }
};
