#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <cassert>

#define _USE_MATH_DEFINES
#include <cmath>

#if !defined(_MSC_VER)
#define HWM_WAVE_TABLE_NOEXCEPT noexcept
#else
#define HWM_WAVE_TABLE_NOEXCEPT
#endif

namespace hwm {

class WaveTable
{
public:
    typedef std::uint32_t entry_value_type;
    typedef std::uint64_t expansion_type;
    typedef std::vector<entry_value_type> container;

    //! a scale factor for mapping a range of [0.0 .. 1.0] to a range of the entry_value_type.
    static entry_value_type const kScaling = 0xFFFFFFFF;

    //! Exception type to notify invalid entry size error.
    struct InvalidEntrySizeError
        : std::runtime_error
    {
        InvalidEntrySizeError()
            : std::runtime_error("invalid entry size was specified")
        {}
    };

private:
    //! Size of a wave table entry.
    std::uint32_t kNumEntry;
    //! A shift amout to get an index for higher bit data.
    std::uint32_t kEntryShift;
    //! A masking value to get an index for lower bit data.
    std::uint32_t kEntryMask;

    //! a scale factor for mapping a range [0.0 .. pi/2] to a range of the entry_value_type.
    double kScaleForXToIndex;
    double kScalingTwice;

    container sin1_table;
    container cos1_table;
    container sin2_table;
    container cos2_table;

    bool check_pow2(std::uint32_t x)
    {
        assert(x != 0);

        for( ; ; ) {
            if(x & 0x1) {
                return (x == 1);
            }
            x >>= 1;
        }
    }

    double calc_sin1(std::uint32_t i) const { return sin(2 * M_PI * i / (double)kNumEntry / 4); };
    double calc_sin2(std::uint32_t i) const { return sin(2 * M_PI * i / (double)(kNumEntry * kNumEntry) / 4); };
    double calc_cos1(std::uint32_t i) const { return cos(2 * M_PI * i / (double)kNumEntry / 4); };
    double calc_cos2(std::uint32_t i) const { return cos(2 * M_PI * i / (double)(kNumEntry * kNumEntry) / 4); };

public:
    //! @pre x = [0 .. kNumEntry^2)
    double fast_sin_with_index(std::uint32_t x) const HWM_WAVE_TABLE_NOEXCEPT
    {
        std::uint32_t const a = x >> kEntryShift;
        std::uint32_t const b = x & kEntryMask;

        expansion_type tmp =
            (expansion_type)sin1_table[a] * (expansion_type)cos2_table[b]
            + (expansion_type)cos1_table[a] * (expansion_type)sin2_table[b];

        return tmp / kScalingTwice;
    }

    //! calculate sin value faster.
    /*!
     * @param x is radian. the value range should be [0.0 .. pi/2)
     */
    double fast_sin(double x) const HWM_WAVE_TABLE_NOEXCEPT
    {
        std::uint32_t const index = (std::uint32_t)(x * kScaleForXToIndex + 0.5);
        return fast_sin_with_index(index);
    }

    //! Get size of allocated memory by this instance.
    /*!
     * a result value of this function does not contain of a size of this class.
     * To get total amout of memory using this class, add sizeof(WaveTable) to the result.
     */
    std::uint32_t GetSizeOfAllocatedMemory() const
    {
        return (kNumEntry + 1) * sizeof(entry_value_type) * 4;
    }

    //! constructor
    /*!
     * @pre num_entry must be a power of 2. ex.) 1, 2, 16, 512, 32768, ...
     * and must be smaller than 65536.
     */
    WaveTable(std::uint32_t num_entry = 1)
        :   kNumEntry(0)
        ,   kEntryShift(0)
        ,   kEntryMask(0)
    {
        if(num_entry == 0 || num_entry > 0xFFFF || check_pow2(num_entry) == false) {
            throw InvalidEntrySizeError();
        }

        kNumEntry = num_entry;

        for(auto tmp = kNumEntry; tmp != 1; tmp >>= 1) {
            kEntryShift += 1;
        }

        kEntryMask = kNumEntry - 1;
        kScaleForXToIndex = 2 / M_PI * (kNumEntry * kNumEntry);
        kScalingTwice = (double)kScaling * kScaling;

        //! allocate entries
        sin1_table.resize(kNumEntry+1);
        sin2_table.resize(kNumEntry+1);
        cos1_table.resize(kNumEntry+1);
        cos2_table.resize(kNumEntry+1);

        //! assign sin, cos values to tables.
        for(std::uint32_t i = 0; i < kNumEntry; ++i) {
            double const sin1_value = calc_sin1(i);
            double const sin2_value = calc_sin2(i);
            double const cos1_value = calc_cos1(i);
            double const cos2_value = calc_cos2(i);

            sin1_table[i] = sin1_value * kScaling + 0.5;
            sin2_table[i] = sin2_value * kScaling + 0.5;
            cos1_table[i] = cos1_value * kScaling + 0.5;
            cos2_table[i] = cos2_value * kScaling + 0.5;
        }

        //! values for x = pi/2
        sin1_table[kNumEntry] = kScaling;
        sin2_table[kNumEntry] = 0;
        cos1_table[kNumEntry] = 0;
        cos2_table[kNumEntry] = 0;
    }

    //! Utility function to check WaveTable's correctness and precision.
    void CheckWaveTableInfo(std::ostream &os) const
    {
        os << "### Calculate memory usage" << std::endl;
        os << "Memory usage: " << sizeof(*this) + GetSizeOfAllocatedMemory() << " bytes" << std::endl;

        os.setf(std::ios_base::scientific, std::ios_base::floatfield);

        int const entry_width = log10(kNumEntry) + 1;

        os << "### Show table entry data" << std::endl;
        for(std::uint32_t i = 0; i < kNumEntry; ++i) {
            os
                << "sin1_table[" << std::setw(entry_width) << i << "]: " << std::setw(10) << sin1_table[i] << ", "
                << "sin2_table[" << std::setw(entry_width) << i << "]: " << std::setw(10) << sin2_table[i] << ", "
                << "cos1_table[" << std::setw(entry_width) << i << "]: " << std::setw(10) << cos1_table[i] << ", "
                << "cos2_table[" << std::setw(entry_width) << i << "]: " << std::setw(10) << cos2_table[i] << std::endl;
        }

        struct statistics {
            statistics() : sum_(0), sum_square_(0), num_(0) {}
            void add_value(double v) {
                sum_ += v;
                sum_square_ += (v * v);
                num_ += 1;
            }
            double average() const { assert(num_ != 0); return sum_ / num_; }
            double variance() const { assert(num_ != 0); return sum_square_ / num_ - (average() * average()); }
            double standard_deviation() const { return sqrt(variance()); }
        private:
            double sum_;
            double sum_square_;
            int num_;
        };

        statistics st_sin1;
        statistics st_sin2;
        statistics st_cos1;
        statistics st_cos2;

        os << "### Restore sin, cos values from each table entries." << std::endl;
        for(std::uint32_t i = 0; i < kNumEntry; ++i) {
            double const sin1_value = sin1_table[i] / (double)kScaling;
            double const sin2_value = sin2_table[i] / (double)kScaling;
            double const cos1_value = cos1_table[i] / (double)kScaling;
            double const cos2_value = cos2_table[i] / (double)kScaling;

            double sin1_diff = sin1_value - calc_sin1(i);
            double sin2_diff = sin2_value - calc_sin2(i);
            double cos1_diff = cos1_value - calc_cos1(i);
            double cos2_diff = cos2_value - calc_cos2(i);

            os  << "[" << std::setw(entry_width) << i << "]" << ", "
                << "difference from calc_sin1(i): " << std::setw(17) << std::setprecision(10) << sin1_diff << ", "
                << "difference from calc_sin2(i): " << std::setw(17) << std::setprecision(10) << sin2_diff << ", "
                << "difference from calc_cos1(i): " << std::setw(17) << std::setprecision(10) << cos1_diff << ", "
                << "difference from calc_cos2(i): " << std::setw(17) << std::setprecision(10) << cos2_diff << std::endl;

            st_sin1.add_value(sin1_diff);
            st_sin2.add_value(sin2_diff);
            st_cos1.add_value(cos1_diff);
            st_cos2.add_value(cos2_diff);
        }

        os << "#### Statistics" << std::endl;
        os  << "|------|------------------------|------------------------|------------------------|\n"
            << "|      | average of differences |        variance        |   standard deviation   |\n"
            << "|------|------------------------|------------------------|------------------------|\n"
            << "| sin1 |"
                << std::setw(24) << std::setprecision(10) << st_sin1.average() << "|"
                << std::setw(24) << std::setprecision(10) << st_sin1.variance() << "|"
                << std::setw(24) << std::setprecision(10) << st_sin1.standard_deviation() << "|\n"
            << "|------|------------------------|------------------------|------------------------|\n"
            << "| sin2 |"
                << std::setw(24) << std::setprecision(10) << st_sin2.average() << "|"
                << std::setw(24) << std::setprecision(10) << st_sin2.variance() << "|"
                << std::setw(24) << std::setprecision(10) << st_sin2.standard_deviation() << "|\n"
            << "|------|------------------------|------------------------|------------------------|\n"
            << "| cos1 |"
                << std::setw(24) << std::setprecision(10) << st_cos1.average() << "|"
                << std::setw(24) << std::setprecision(10) << st_cos1.variance() << "|"
                << std::setw(24) << std::setprecision(10) << st_cos1.standard_deviation() << "|\n"
            << "|------|------------------------|------------------------|------------------------|\n"
            << "| cos2 |"
                << std::setw(24) << std::setprecision(10) << st_cos2.average() << "|"
                << std::setw(24) << std::setprecision(10) << st_cos2.variance() << "|"
                << std::setw(24) << std::setprecision(10) << st_cos2.standard_deviation() << "|\n"
            << "|------|------------------------|------------------------|------------------------|\n"
            << std::endl;

        os << "### Compare fast_sin() and sin() of STL" << std::endl;
        std::uint32_t const kAnglePrecision = 100;
        statistics st_fast_sin;

        for(std::uint32_t i = 0; i < kAnglePrecision; ++i) {
            double theta = 2 * M_PI / 4 * i / (double)kAnglePrecision;
            double diff = fast_sin(theta) - sin(theta);

            st_fast_sin.add_value(diff);

            os
                << "theta = " << std::setw(16) << std::setprecision(10) << theta << ", "
                << "fast_sin(theta): " << std::setw(16) << std::setprecision(10) << fast_sin(theta) << ", "
                << "sin(theta): " << std::setw(16) << std::setprecision(10) << sin(theta) << ", "
                << "difference: " << std::setw(17) << std::setprecision(10) << diff << std::endl;
        }

        os << "#### Statistics" << std::endl;
        os  << "|------------------------|------------------------|------------------------|\n"
            << "| average of differences |        variance        |   standard deviation   |\n"
            << "|------------------------|------------------------|------------------------|\n"
            << "|"
                << std::setw(24) << std::setprecision(10) << st_fast_sin.average() << "|"
                << std::setw(24) << std::setprecision(10) << st_fast_sin.variance() << "|"
                << std::setw(24) << std::setprecision(10) << st_fast_sin.standard_deviation() << "|\n"
            << "|------------------------|------------------------|------------------------|\n"
            << std::endl;
    }
};

} // namespace hwm

