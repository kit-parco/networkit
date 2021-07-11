// no-networkit-format
/*
* MemoryMappedFileGTest.h
*
*  Created on: 28.07.2018
*      Author: Manuel Penschuck (networkit@manuel.jetzt)
*/

#ifdef NETWORKIT_WINDOWS
#define NOMINMAX // windows.h by default defines the maros min/max, which prevent the usage of numeric_limits
#include <windows.h>
#else
#include <cstdlib>
#endif

#include <algorithm>
#include <fstream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <stdexcept>

#include "MemoryMappedFileGTest.hpp"

#include <networkit/io/MemoryMappedFile.hpp>


namespace NetworKit {
    // This object generates a temporary file and fills it with a request number
    // of random bytes. On destruction of the elemen, the file is deleted automatically.
    // Make sure that no application holds a handle to the file at this point, since
    // on Windows removal of the file fails silently in this case.
    class TemporaryFile {
    public:
        using value_type = char;
        using const_iterator = const value_type*;

        TemporaryFile(size_t bytes, std::default_random_engine& prng) : 
            dataBegin(new value_type[bytes]),
            dataEnd(dataBegin.get() + bytes)
        {
            // fill array with random data
            {
                // we use int (rather than char), as char is not an integer type for MSVC
                std::uniform_int_distribution<int> distr{std::numeric_limits<value_type>::min(),
                                                         std::numeric_limits<value_type>::max()};

                std::generate(dataBegin.get(), dataEnd, [&distr, &prng] {
                    return static_cast<value_type>(distr(prng));
                });
            }

            // obtain temporary file name (C++17 offers std::filesystem, which we currently cannot rely on)
#ifdef NETWORKIT_WINDOWS
            {
                char cpath[MAX_PATH];

                if (!GetTempFileName(".",  "Tst", 0, cpath)) {
                    throw std::runtime_error("Could not create temporary file");
                }

                path = cpath;
            }
#else
            {
                char cpath[32];
                strcpy(cpath, "TestMemoryMappedFile_XXXXXX");
                int fd = mkstemp(cpath);
                if (fd == -1)
                    throw std::runtime_error("Could not create temporary file");

                path = cpath;
                ::close(fd);
            }
#endif

            {
                std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::trunc);
                file.write(reinterpret_cast<const char*>(dataBegin.get()), bytes);
            }

        }

        ~TemporaryFile() {
            if (!path.empty()) {
                std::remove(path.c_str());
            }
        }

        const std::string& filename() const { return path; }

        const_iterator cbegin() const {
            return dataBegin.get();
        }

        const_iterator cend() const {
            return dataEnd;
        }

        size_t size() const { return dataEnd - dataBegin.get(); }

        void verify_mapping(const MemoryMappedFile& mmf) const {
            ASSERT_EQ(size(), mmf.size());

            for (size_t i = 0; i < size(); i++) {
                ASSERT_EQ(*(mmf.cbegin() + i), *(cbegin() + i)) << "i=" << i;
            }
        }

    private:
        std::unique_ptr<value_type[]> dataBegin;
        value_type* dataEnd;
        std::string path;
    };


    TEST_P(MemoryMappedFileIOGTest, testReadFile) {
        const size_t bytes = GetParam();
        std::default_random_engine prng(static_cast<unsigned>(bytes));
        TemporaryFile tmpFile(bytes, prng);
        ASSERT_EQ(bytes, tmpFile.size());

        MemoryMappedFile mmf(tmpFile.filename());

        tmpFile.verify_mapping(mmf);
    }

    INSTANTIATE_TEST_SUITE_P(MemoryMappedFileIOGTest, MemoryMappedFileIOGTest,
        ::testing::Values(0, 1,
            (1 <<  2) - 1, (1 <<  2), (1 <<  2) + 1,
            (1 << 10) - 1, (1 << 10), (1 << 10) + 1,
            (1 << 16) - 1, (1 << 16), (1 << 16) + 1,
            (1 << 22) - 1, (1 << 22), (1 << 22) + 1
        ));

    TEST_F(MemoryMappedFileGTest, testMove) {
        const size_t bytes = 1000;
        std::default_random_engine prng(static_cast<unsigned>(bytes));
        TemporaryFile tmpFile1(bytes, prng);
        TemporaryFile tmpFile2(bytes+1, prng);

        MemoryMappedFile mmf1(tmpFile1.filename());
        MemoryMappedFile mmf2(tmpFile2.filename());

        tmpFile1.verify_mapping(mmf1);
        tmpFile2.verify_mapping(mmf2);

        mmf1 = std::move(mmf2);  // <- TESTED

        ASSERT_FALSE(mmf2.size());
        ASSERT_EQ(mmf2.cbegin(), nullptr);
        ASSERT_EQ(mmf2.cend(), nullptr);

        tmpFile2.verify_mapping(mmf1);
    }

    TEST_F(MemoryMappedFileGTest, testSwap) {
        const size_t bytes = 1100;
        std::default_random_engine prng(static_cast<unsigned>(bytes));
        TemporaryFile tmpFile1(bytes, prng);
        TemporaryFile tmpFile2(bytes+1, prng);

        MemoryMappedFile mmf1(tmpFile1.filename());
        MemoryMappedFile mmf2(tmpFile2.filename());

        tmpFile1.verify_mapping(mmf1);
        tmpFile2.verify_mapping(mmf2);

        std::swap(mmf1, mmf2); // <- TESTED

        tmpFile1.verify_mapping(mmf2);
        tmpFile2.verify_mapping(mmf1);
    }

}
