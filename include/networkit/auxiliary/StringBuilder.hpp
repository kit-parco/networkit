// no-networkit-format
#ifndef NETWORKIT_AUXILIARY_STRING_BUILDER_HPP_
#define NETWORKIT_AUXILIARY_STRING_BUILDER_HPP_

#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iterator>
#include <type_traits>
#include <tuple>
#include <utility>

#include <networkit/auxiliary/TemplateUtils.hpp>

namespace Aux {

template<typename...T>
std::string toString(const T&...args);

template<typename...T>
std::ostream& printToStream(std::ostream& stream, const T&...args);

template<typename...T>
std::string toStringF(const std::string& format, const T&...args);

template<typename...T>
std::ostream& printToStreamF(std::ostream& stream, const std::string& format, const T&...args);

// Implementation
/////////////////

namespace Impl {

// Categories of how a type might be printable
enum class PrintableCategory {
    Unprintable,
    Iteratable,
    Pair,
    Tuple,
    Streamable,
};

template<typename T> constexpr bool isStreamable();
template<typename T> constexpr bool isPair();
template<typename T> constexpr bool isTuple();
template<typename T> constexpr bool isIteratable();

template<typename T>
constexpr PrintableCategory getPrintableCategory() {
    return
        isStreamable<T>() ? PrintableCategory::Streamable  :
        isPair<T>()       ? PrintableCategory::Pair        :
        isTuple<T>()      ? PrintableCategory::Tuple       :
        isIteratable<T>() ? PrintableCategory::Iteratable  :
        /* else: */         PrintableCategory::Unprintable ;
}
template<PrintableCategory Tag> struct PrintableCategoryTag{};
using IteratableTag   = PrintableCategoryTag< PrintableCategory::Iteratable  >;
using PairTag         = PrintableCategoryTag< PrintableCategory::Pair        >;
using TupleTag        = PrintableCategoryTag< PrintableCategory::Tuple       >;
using StreamableTag   = PrintableCategoryTag< PrintableCategory::Streamable  >;
using UnprintableTag  = PrintableCategoryTag< PrintableCategory::Unprintable >;

template<typename T, typename...Args> void printToStream(std::ostream& stream, const T&, const Args&...);
inline void printToStream(std::ostream&) {}

template<typename T> void printToStreamTagged(std::ostream& stream, const T&, IteratableTag);
template<typename T> void printToStreamTagged(std::ostream& stream, const T&, PairTag);
template<typename T> void printToStreamTagged(std::ostream& stream, const T&, TupleTag);
template<typename T> void printToStreamTagged(std::ostream& stream, const T&, StreamableTag);

// Claim that this function exists somewhere else to keep the errors clean
// (calling this function is not a problem since the error is already caught earlier
// in printToStream, and calling it directly will result in a linker-error)
template<typename T> extern void printToStreamTagged(std::ostream& , const T&, UnprintableTag);

template<typename T, typename...Args>
void printToStream(std::ostream& stream, const T& arg, const Args&...args) {
    static_assert(getPrintableCategory<T>() != PrintableCategory::Unprintable,
            "printToStream must not be called with an unprintable argument");
    printToStreamTagged(stream, arg, PrintableCategoryTag<getPrintableCategory<T>()>{});
    printToStream(stream, args...);
}

inline std::tuple<std::string::const_iterator, bool> printFormatPartToStream(std::ostream& stream,
        std::string::const_iterator begin, std::string::const_iterator end);

inline void printToStreamF(std::ostream& stream, std::string::const_iterator format_begin,
        std::string::const_iterator format_end) {
    bool printArgument;
    using iterator = std::string::const_iterator;
    iterator it;
    std::tie(it, printArgument) = printFormatPartToStream(stream, format_begin, format_end);
    if (printArgument) {
        throw std::invalid_argument{"formatstring requests more arguments then provided"};
    }
}
template<typename T, typename...Args>
void printToStreamF(std::ostream& stream, std::string::const_iterator format_begin,
        std::string::const_iterator format_end, const T& arg, const Args&...args) {
    bool printArgument;
    using iterator = std::string::const_iterator;
    iterator it;
    std::tie(it, printArgument) = printFormatPartToStream(stream, format_begin, format_end);
    if(printArgument) {
        printToStream(stream, arg);
        printToStreamF(stream, it, format_end, args...);
    } else {
        assert(it == format_end);
        return;
    }
}



/**
 * Write the formatstring until the first occurance of "%s"
 * to stream, '%%' will be replaced by '%'.
 *
 * @returns The iterator that points one after the last consumed
 * character and whether a further argument should be printed.
 */
inline auto printFormatPartToStream(std::ostream& stream, std::string::const_iterator begin,
        std::string::const_iterator end)
    -> std::tuple<std::string::const_iterator, bool>
{
    auto it = begin;
    while (it != end) {
        auto nextPercent = std::find(it, end, '%');
        stream.write(&*it, std::distance(it, nextPercent));
        
        if(nextPercent == end) {
            it = nextPercent;
            break;
        }
        
        it = ++nextPercent;
        
        if(it == end) {
            throw std::invalid_argument{"formatstrings must not end on unmatched '%'"};
        } else if (*it == '%') {
            stream.put('%');
            ++it;
        } else if (*it == 's') {
            ++it;
            return std::make_tuple(it, true);
        } else {
            throw std::invalid_argument{"formatstring contains illegal format-specifier"};
        }
    }
    assert(it == end);
    return std::make_tuple(end, false);
}

// Brace Yourself: Templatemetaprogramming is coming
////////////////////////////////////////////////////

// Iteratable
struct IsIteratableHelper {
    static std::false_type isIteratable(...);
    
    template<typename T,
        class Iterator = decltype(std::begin(std::declval<T>())),
        AUX_REQUIRE(EndIteratorValid, isSame<Iterator, decltype(std::end(std::declval<T>()))>()),
        AUX_REQUIRE(HasInputIterator, isBaseOrSame<std::input_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>())
    > static std::true_type isIteratable(const T&);
};
template<typename T> constexpr bool isIteratable() {
    return decltype(IsIteratableHelper::isIteratable(std::declval<T>()))::value;
}

// Pair
template<typename T              > struct IsPairHelper                    : std::false_type {};
template<typename T1, typename T2> struct IsPairHelper<std::pair<T1, T2>> : std::true_type {};
template<typename T> constexpr bool isPair() {
    return IsPairHelper<T>::value;
}

// Tuple
template<typename   T> struct IsTupleHelper                   : std::false_type {};
template<typename...T> struct IsTupleHelper<std::tuple<T...>> : std::true_type {};
template<typename T> constexpr bool isTuple() {
    return IsTupleHelper<T>::value;
}

// Streamable
struct IsStreamableHelper {
    static std::false_type isStreamable(...);
    
    template<typename T,
        typename _streamT = typename std::decay< decltype(std::declval<std::ostream&>() << std::declval<const T&>())>::type,
        typename std::enable_if< std::is_same<std::ostream, _streamT>::value || std::is_base_of<std::ostream, _streamT>::value >* = nullptr
    >
    static std::true_type isStreamable(const T&);
    
};
template<typename T> constexpr bool isStreamable() {
    return decltype(IsStreamableHelper::isStreamable(std::declval<T>()))::value;
}


// And now: implement the actual printing:
//////////////////////////////////////////

// Streamable
template<typename T> void printToStreamTagged(std::ostream& stream, const T& arg, StreamableTag) {
    stream << arg;
}

// Pair
template<typename T> void printToStreamTagged(std::ostream& stream, const T& arg, PairTag) {
    stream << '(';
    printToStream(stream, arg.first);
    stream << ", ";
    printToStream(stream, arg.second);
    stream << ')';
}

// Tuple
template<typename Tuple, int I, int TupleSize> struct printTupleHelper {
    static void print(std::ostream& stream, const Tuple& arg) {
        printToStream(stream, std::get<I-1>(arg));
        stream << ", ";
        printTupleHelper<Tuple, I+1, TupleSize>::print(stream, arg);
    }
};
template<typename Tuple, int I> struct printTupleHelper<Tuple, I, I>{
    static void print(std::ostream& stream, const Tuple& arg) {
        printToStream(stream, std::get<I-1>(arg));
    }
};
template<typename T> void printToStreamTagged(std::ostream& stream, const T& arg, TupleTag) {
    stream << '(';
    printTupleHelper<T, 1, std::tuple_size<T>::value>::print(stream, arg);
    stream << ')';
}

// Iteratable
template<typename T> void printToStreamTagged(std::ostream& stream, const T& arg, IteratableTag) {
    auto it = std::begin(arg);
    auto end = std::end(arg);
    bool firstpass = true;
    stream << '[';
    while(it != end) {
        if(firstpass) {
            firstpass = false;
        }
        else {
            stream << ", ";
        }
        printToStream(stream, *it);
        ++it;
    }
    stream << ']';
}


} // namespace Impl

// Finally: put together the public interface:
//////////////////////////////////////////////

template<typename...T>
std::string toString(const T&...args) {
    std::stringstream stream;
    printToStream(stream, args...);
    return stream.str();
}

template<typename...T>
std::ostream& printToStream(std::ostream& stream, const T&...args) {
    Impl::printToStream(stream, args...);
    return stream;
}

template<typename...T>
std::string toStringF(const std::string& format, const T&...args) {
    std::stringstream stream;
    printToStreamF(stream, format, args...);
    return stream.str();
}

template<typename...T>
std::ostream& printToStreamF(std::ostream& stream, const std::string& format, const T&...args) {
    Impl::printToStreamF(stream, format.begin(), format.end(), args...);
    return stream;
}

} //namespace Aux


#endif // NETWORKIT_AUXILIARY_STRING_BUILDER_HPP_
