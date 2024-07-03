/*
 * Command Line Argument Parser
 */

#ifndef ZQ_CLAP
#define ZQ_CLAP

#include <algorithm>
#include <cstdlib>
#include <iostream>
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

struct OptionConfig {
    std::size_t offset;
    TYPE_CODE type_code;
    std::string short_name;
    std::string description;

    bool has_short_name() const { return !short_name.empty(); }
};

template <typename T>
struct CLAP {
    using base = T;
    using string = std::string;
    using option_container_type = std::unordered_map<std::string, OptionConfig>;

    static auto& get_options() {
        static option_container_type options{};
        return options;
    }

    static auto& get_short_name_map() {
        static std::unordered_map<std::string,
                                  option_container_type::const_iterator>
            short_name_map;
        return short_name_map;
    }

    static auto& get_arguments() {
        static std::vector<std::pair<std::size_t, TYPE_CODE>> arguments;
        return arguments;
    }

    static auto& get_usage() {
        static std::string usage;
        return usage;
    }

    static auto& get_description() {
        static std::string description;
        return description;
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
        const std::string l_prefix = "--";
        for (int i = 1; i < argc; ++i) {
            std::string option_or_arg{argv[i]};
            if (option_or_arg.compare("--help") == 0 ||
                option_or_arg.compare("-h") == 0) {
                print_help(argv[0]);
                std::exit(0);
            }
            if (option_or_arg[0] == '-') {
                const auto& opts = get_options();
                std::size_t offset{};
                TYPE_CODE type_code{};
                auto iter = opts.cbegin();
                std::string raw_value;

                bool invalid_option = false;
                bool has_value = false;

                // parse option name
                if (option_or_arg[1] == '-') {
                    // long option
                    std::size_t equal_sign_idx = 0;
                    if ((equal_sign_idx = option_or_arg.find_first_of('=')) !=
                        std::string::npos) {
                        iter =
                            opts.find(option_or_arg.substr(0, equal_sign_idx));
                        raw_value = option_or_arg.substr(equal_sign_idx + 1,
                                                         std::string::npos);
                        has_value = true;
                    } else {
                        iter = opts.find(option_or_arg);
                    }
                    invalid_option = iter == opts.end();
                } else {
                    // short option
                    auto& short_opts = get_short_name_map();
                    auto iter_short = short_opts.find(option_or_arg);
                    if (!(invalid_option = iter_short == short_opts.end())) {
                        iter = iter_short->second;
                    }
                }
                if (invalid_option) {
                    std::ostringstream oss;
                    oss << argv[0] << ": unrecognized option '" << option_or_arg
                        << "'\n";
                    throw std::invalid_argument(oss.str());
                } else {
                    offset = iter->second.offset;
                    type_code = iter->second.type_code;
                }

                // parse option value
                if (!has_value) {
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

    static void print_help(const char* program_name) {
        std::vector<option_container_type::const_iterator> iters;
        for (auto it = get_options().begin(); it != get_options().end(); ++it) {
            iters.push_back(it);
        }
        std::sort(iters.begin(), iters.end(),
                  [](auto it1, auto it2) { return it1->first < it2->first; });
        auto print_space = [](std::size_t n) {
            for (std::size_t i = 0; i < n; ++i) { std::cout << ' '; }
        };

        std::cout << "Usage: " << program_name << " " << get_usage() << '\n';
        std::cout << get_description() << "\n\n";
        for (auto it : iters) {
            std::size_t col = 0;
            if (it->second.has_short_name()) {
                print_space(2);
                std::cout << it->second.short_name << ", ";
                col = 4 + it->second.short_name.size();
            } else {
                print_space(6);
                col = 6;
            }
            std::cout << it->first;  // long name
            col += it->first.size();
            if (col > 28) {
                std::cout << '\n';
                print_space(30);
            } else {
                print_space(30 - col);
            }
            std::cout << it->second.description << '\n';
        }
    }

    static void report_repeat_definition(std::string option_name) {
#ifdef CLAP_DEBUG
        std::ostringstream oss;
        oss << "You can not define option '" << option_name << "' twice!\n";
        throw std::invalid_argument(oss.str());
#else
        (void)option_name;
#endif
    }
#undef TYPE_LIST
};

#define CLAP_BEGIN(NAME)                        \
    struct _clap_##NAME##_ : NAME, CLAP<NAME> { \
        using base = NAME;                      \
        static void define_parameters() {
#define CLAP_END(NAME)                                                     \
    {                                                                      \
        auto result = get_options().emplace(                               \
            "--help", OptionConfig{0, TYPE_CODE::int_t, "",                \
                                   "display this help message and exit"}); \
        if (get_short_name_map().emplace("-h", result.first).second) {     \
            result.first->second.short_name = "-h";                        \
        }                                                                  \
    }                                                                      \
    }                                                                      \
    }                                                                      \
    ;                                                                      \
    _clap_##NAME##_::define_parameters();

#define CLAP_REGISTER_OPT_MINIMAL(NAME)                                 \
    {                                                                   \
        auto result = get_options().emplace(                            \
            "--" #NAME, OptionConfig{offsetof(base, NAME),              \
                                     get_type_code<decltype(NAME)>()}); \
        if (!result.second) { report_repeat_definition(OPTION_NAME); }  \
    }

#define CLAP_REGISTER_OPT_LONG(NAME, OPTION_NAME)                        \
    {                                                                    \
        auto result = get_options().emplace(                             \
            OPTION_NAME, OptionConfig{offsetof(base, NAME),              \
                                      get_type_code<decltype(NAME)>()}); \
        if (!result.second) { report_repeat_definition(OPTION_NAME); }   \
    }

#define CLAP_REGISTER_OPT_LONG_SHORT(NAME, OPTION_NAME, OPTION_NAME_SHORT)     \
    {                                                                          \
        auto result = get_options().emplace(                                   \
            OPTION_NAME,                                                       \
            OptionConfig{offsetof(base, NAME),                                 \
                         get_type_code<decltype(NAME)>(), OPTION_NAME_SHORT}); \
        if (!result.second) { report_repeat_definition(OPTION_NAME); }         \
        get_short_name_map().emplace(OPTION_NAME_SHORT, result.first);         \
    }

#define CLAP_REGISTER_OPT_DESC(NAME, DESC)                             \
    {                                                                  \
        auto result = get_options().emplace(                           \
            "--" #NAME,                                                \
            OptionConfig{offsetof(base, NAME),                         \
                         get_type_code<decltype(NAME)>(), "", DESC});  \
        if (!result.second) { report_repeat_definition(OPTION_NAME); } \
    }

#define CLAP_REGISTER_OPT_LONG_DESC(NAME, OPTION_NAME, DESC)           \
    {                                                                  \
        auto result = get_options().emplace(                           \
            OPTION_NAME,                                               \
            OptionConfig{offsetof(base, NAME),                         \
                         get_type_code<decltype(NAME)>(), "", DESC});  \
        if (!result.second) { report_repeat_definition(OPTION_NAME); } \
    }

#define CLAP_REGISTER_OPT_FULL(NAME, OPTION_NAME, OPTION_NAME_SHORT, DESC) \
    {                                                                      \
        auto result = get_options().emplace(                               \
            OPTION_NAME, OptionConfig{offsetof(base, NAME),                \
                                      get_type_code<decltype(NAME)>(),     \
                                      OPTION_NAME_SHORT, DESC});           \
        if (!result.second) { report_repeat_definition(OPTION_NAME); }     \
        get_short_name_map().emplace(OPTION_NAME_SHORT, result.first);     \
    }

// Macro overloading depending on argument number
// https://stackoverflow.com/a/11763277/7255197
#define CLAP_GET_MACRO_3(_1, _2, _3, NAME, ...) NAME
#define CLAP_REGISTER_OPTION(...)                                       \
    CLAP_GET_MACRO(__VA_ARGS__, CLAP_REGISTER_OPT_LONG_SHORT,           \
                   CLAP_REGISTER_OPT_LONG, CLAP_REGISTER_OPT_MINIMAL, ) \
    (__VA_ARGS__)

#define CLAP_GET_MACRO_4(_1, _2, _3, _4, NAME, ...) NAME
#define CLAP_REGISTER_OPTION_WITH_DESCRIPTION(...)                          \
    CLAP_GET_MACRO_4(__VA_ARGS__, CLAP_REGISTER_OPT_FULL,                   \
                     CLAP_REGISTER_OPT_LONG_DESC, CLAP_REGISTER_OPT_DESC, ) \
    (__VA_ARGS__)

#define CLAP_REGISTER_ARG(NAME)                                        \
    {                                                                  \
        get_arguments().emplace_back(offsetof(base, NAME),             \
                                     get_type_code<decltype(NAME)>()); \
    }

#define CLAP_ADD_USAGE(USAGE) get_usage() = USAGE;
#define CLAP_ADD_DESCRIPTION(DESC) get_description() = DESC;

#endif  // ZQ_CLAP
