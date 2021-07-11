// no-networkit-format
/*
 * MemoryMappedFile.hpp
 *
 *  Created on: 16.07.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef NETWORKIT_IO_MEMORY_MAPPED_FILE_HPP_
#define NETWORKIT_IO_MEMORY_MAPPED_FILE_HPP_

#include <memory>
#include <string>

namespace NetworKit {
struct MemoryMappedFileState;

/**
 * This class is a wrapper to os-dependend file mapping implementations.
 * It supports to map a read-only file into virtual address space.
 * After successful opening, pointers to the begin/end of the virtual
 * address space are accessible via cbegin() / cend(). All file and mapping
 * handlers are freed automatically.
 *
 * This wrapper is non-copyable but can be cheaply moved and swapped.
 */
class MemoryMappedFile final {
public:
    using value_type = char;
    using const_iterator = const value_type*;

    //! Creates a MemoryMappedFile instance in an unmapped state.
    //! A call to open() is required.
    MemoryMappedFile();

    //! Invokes open(path) automatically
    explicit MemoryMappedFile(const std::string& path);

    //! Invokes close
    ~MemoryMappedFile();

    //! It's non copy-able
    MemoryMappedFile(const MemoryMappedFile&) = delete;
    MemoryMappedFile& operator=(const MemoryMappedFile&) = delete;

    //! Takes over the mapping (if any) from o, leaves o in a "closed"
    //! state ready to be opened again.
    MemoryMappedFile(MemoryMappedFile&& o) noexcept;

    //! If *this currently holds a mapping, it is closed and replaced
    //! the other's state which in turn is left closed.
    MemoryMappedFile& operator=(MemoryMappedFile&& o) noexcept;

    //! Opens the file and maps it to cbegin() ... cend()
    //! Opening an empty file is considered an error.
    void open(const std::string& file);

    //! If a file is mapped, it is closed. Otherwise, operation is carried out.
    //! @note This function is automatically called by the destructor.
    void close() noexcept;

    //! If a file is opened, a random-access iterator to the first byte mapped is returned.
    //! If no file is opened, nullptr is returned.
    const_iterator cbegin() const {
        return beginIt;
    }

    //! Analogously to cbegin()
    const_iterator cend() const {return endIt;}

    //! Number of bytes mapped
    size_t size() const {
        return std::distance(beginIt, endIt);
    }

private:
    const_iterator beginIt {nullptr};
    const_iterator endIt {nullptr};

    std::unique_ptr<MemoryMappedFileState> state; //!< used only by windows implementation to keep handlers around.
};

}
#endif // NETWORKIT_IO_MEMORY_MAPPED_FILE_HPP_
