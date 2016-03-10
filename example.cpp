#include <sstream>
#include <boost/timer.hpp>
#include "WaveTable.hpp"

int main(int argc, char **argv)
{
    if(argc != 2) {
        std::cout << "usage: a.out <Table Entry Size (which should be power of 2)>" << std::endl;
        return 1;
    }

    int tmp_num_entry = 0;
    std::stringstream ss;
    ss << argv[1];
    ss >> tmp_num_entry;

    hwm::WaveTable table;
    try {
        table = hwm::WaveTable(tmp_num_entry);
    } catch(std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    //! Output WaveTable Information.
    table.CheckWaveTableInfo(std::cout);

    //! Measure calculation speed.
    boost::timer tm;
    int kMaxRepeatCount = 100;
    int kAnglePrecision = 100000;
    double sec1 = 0;
    double sec2 = 0;

    double tmp1 = 0;
    double omega = M_PI / 2 / (double)kAnglePrecision;

    std::stringstream dummy_output;

    tm.restart();
    for(int repeat = 0; repeat < kMaxRepeatCount; ++repeat) {
        for(int i = 0; i < kAnglePrecision; ++i) {
            tmp1 = tmp1 + sin(i * omega);
        }
    }
    sec1 = tm.elapsed();
    dummy_output << tmp1;

    tmp1 = 0;

    tm.restart();
    for(int repeat = 0; repeat < kMaxRepeatCount; ++repeat) {
        for(int i = 0; i < kAnglePrecision; ++i) {
            tmp1 = tmp1 + table.fast_sin(i * omega);
        }
    }
    sec2 = tm.elapsed();

    dummy_output << tmp1;

    std::string volatile dummy_string(dummy_output.str());

    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::cout << "elapsed sin(x):      " << std::setw(16) << std::setprecision(10) << sec1 << " sec" << std::endl;
    std::cout << "elapsed fast_sin(x): " << std::setw(16) << std::setprecision(10) << sec2 << " sec" << std::endl;
    std::cout << "average time of sin(x):      " << std::setw(10) << std::setprecision(5) << (sec1 / (kMaxRepeatCount * kAnglePrecision) * 1000000000.0) << " nanosec" << std::endl;
    std::cout << "average time of fast_sin(x): " << std::setw(10) << std::setprecision(5) << (sec2 / (kMaxRepeatCount * kAnglePrecision) * 1000000000.0) << " nanosec" << std::endl;
    std::cout << "time of fast_sin(x) / time of sin(x): " << std::setw(6) << std::setprecision(4) << (sec2 / sec1 * 100) << "%" << std::endl;
}

