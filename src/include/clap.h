/*
 * Command Line Argument Parser
 */

#ifndef ZQ_CLAP
#define ZQ_CLAP

#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <vector>

#define TYPE_LIST()        \
    PROCESS_TYPE(bool)     \
    PROCESS_TYPE(int)      \
    PROCESS_TYPE(unsigned) \
    PROCESS_TYPE(double)   \
    PROCESS_TYPE(string)

enum class TYPE_CODE {
#define PROCESS_TYPE(TYPE) TYPE##_t,
    TYPE_LIST()
#undef PROCESS_TYPE
};

template <typename T>
struct CLAP {
    using base = T;
    using string = std::string;

    static auto& get_options() {
        static std::unordered_map<std::string,
                                  std::pair<std::size_t, TYPE_CODE>>
            parameter_map{};
        return parameter_map;
    }

    static auto& get_arguments() {
        static std::vector<std::pair<std::size_t, TYPE_CODE>> arguments;
        return arguments;
    }

    template <typename Arg>
    static TYPE_CODE get_type_code() {
#define PROCESS_TYPE(TYPE) \
    if (std::is_same<Arg, TYPE>::value) { return TYPE_CODE::TYPE##_t; }
        TYPE_LIST()
#undef PROCESS_TYPE
        throw std::invalid_argument("Type not supported yet.");
    }

    static void assign_value(char* adr,
                             TYPE_CODE type_code,
                             const std::string& raw_value) {
        std::istringstream iss(raw_value);
#define PROCESS_TYPE(TYPE)                    \
    case (TYPE_CODE::TYPE##_t): {             \
        iss >> *reinterpret_cast<TYPE*>(adr); \
        break;                                \
    }
        switch (type_code) { TYPE_LIST() }
#undef PROCESS_TYPE
    }

    static void parse_input(T& input, int argc, char** argv) {
        std::size_t argument_idx = 0;
        for (int i = 1; i < argc; ++i) {
            std::string option_or_arg{argv[i]};
            if (option_or_arg[0] == '-') {
                // it's an option
                const auto& opts = get_options();
                auto iter = opts.find(option_or_arg);
                if (iter == opts.end()) {
                    std::ostringstream oss;
                    oss << argv[0] << ": unrecognized option '" << option_or_arg
                        << "'\n";
                    throw std::invalid_argument(oss.str());
                }

                // valid option
                auto [offset, type_code] = iter->second;
                std::string raw_value;
                if (type_code == TYPE_CODE::bool_t) {
                    raw_value = "1";
                } else if (i + 1 == argc) {
                    std::ostringstream oss;
                    oss << argv[0] << ": missing value of option '"
                        << option_or_arg << "'\n";
                    throw std::invalid_argument(oss.str());
                } else {
                    raw_value = argv[++i];
                }
                assign_value(reinterpret_cast<char*>(&input) + offset,
                             type_code, raw_value);
            } else {
                // it's an argument
                if (argument_idx == get_arguments().size()) {
                    std::ostringstream oss;
                    oss << argv[0] << ": too many arguments\n";
                    throw std::invalid_argument(oss.str());
                }
                auto [offset, type_code] = get_arguments()[argument_idx];
                assign_value(reinterpret_cast<char*>(&input) + offset,
                             type_code, argv[i]);
                ++argument_idx;
            }
        }
        if (argument_idx < get_arguments().size()) {
            std::ostringstream oss;
            oss << argv[0] << ": too few arguments\n";
            throw std::invalid_argument(oss.str());
        }
    }
#undef TYPE_LIST
};

#define CLAP_BEGIN(NAME)                        \
    struct _clap_##NAME##_ : NAME, CLAP<NAME> { \
        using base = NAME;                      \
        static void define_parameters() {
#define CLAP_END(NAME) \
    }                  \
    }                  \
    ;                  \
    _clap_##NAME##_::define_parameters();

#define CLAP_REGISTER_DIRECT(NAME)                                        \
    {                                                                     \
        get_options().emplace(                                            \
            "--" #NAME, std::make_pair(offsetof(base, NAME),              \
                                       get_type_code<decltype(NAME)>())); \
    }

#define CLAP_REGISTER_EXTEND(NAME, OPTION_NAME)                            \
    {                                                                      \
        get_options().emplace(                                             \
            OPTION_NAME, std::make_pair(offsetof(base, NAME),              \
                                        get_type_code<decltype(NAME)>())); \
    }

// Macro overloading depending on argument number
// https://stackoverflow.com/a/11763277/7255197
#define CLAP_GET_MACRO(_1, _2, NAME, ...) NAME
#define CLAP_REGISTER_OPTION(...)                                             \
    CLAP_GET_MACRO(__VA_ARGS__, CLAP_REGISTER_EXTEND, CLAP_REGISTER_DIRECT, ) \
    (__VA_ARGS__)

#define CLAP_REGISTER_ARG(NAME)                                        \
    {                                                                  \
        get_arguments().emplace_back(offsetof(base, NAME),             \
                                     get_type_code<decltype(NAME)>()); \
    }

#endif  // ZQ_CLAP
