#ifndef ZQ_CLAP
#define ZQ_CLAP

#include <sstream>
#include <type_traits>
#include <unordered_map>

enum class TYPE_CODE { int_t, double_t, string_t };

template <typename T>
struct CLAP {
    using base = T;
    using string = std::string;

    static inline std::unordered_map<std::string,
                                     std::pair<std::size_t, TYPE_CODE>>
        parameter_map{};

#define TYPE_LIST()      \
    PROCESS_TYPE(int)    \
    PROCESS_TYPE(double) \
    PROCESS_TYPE(string)

    template <typename Arg>
    static TYPE_CODE get_type_code() {
#define PROCESS_TYPE(TYPE) \
    if constexpr (std::is_same_v<Arg, TYPE>) { return TYPE_CODE::TYPE##_t; }
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
        (void)input;

        for (int i = 1; i < argc; ++i) {
            std::string option{argv[i]};
            bool valid_option = false;
            for (auto&& [key, config] : parameter_map) {
                if (option == key) {
                    if (++i == argc) {
                        std::ostringstream oss;
                        oss << argv[0] << ": missing value of option '"
                            << argv[i - 1] << "'\n";
                        throw std::invalid_argument(oss.str());
                    }
                    std::string raw_value{argv[i]};
                    assign_value(reinterpret_cast<char*>(&input) + config.first,
                                 config.second, raw_value);

                    valid_option = true;
                    break;
                }
            }
            if (!valid_option) {
                std::ostringstream oss;
                oss << argv[0] << ": invalid option '" << argv[i] << "'\n";
                throw std::invalid_argument(oss.str());
            }
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
        parameter_map.emplace(                                            \
            "--" #NAME, std::make_pair(offsetof(base, NAME),              \
                                       get_type_code<decltype(NAME)>())); \
    }

#define CLAP_REGISTER_EXTEND(NAME, OPTION_NAME)                            \
    {                                                                      \
        parameter_map.emplace(                                             \
            OPTION_NAME, std::make_pair(offsetof(base, NAME),              \
                                        get_type_code<decltype(NAME)>())); \
    }

// Macro overloading depending on argument number
// https://stackoverflow.com/a/11763277/7255197
#define CLAP_GET_MACRO(_1, _2, NAME, ...) NAME
#define CLAP_REGISTER(...)                                                    \
    CLAP_GET_MACRO(__VA_ARGS__, CLAP_REGISTER_EXTEND, CLAP_REGISTER_DIRECT, ) \
    (__VA_ARGS__)

#endif  // ZQ_CLAP
