#include <vector>
#include <random>
#include <array>
#include <string>

#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

#define DETERMINISTIC() true
#define RANDOMIZE_LDS_START() true

static const size_t c_numHistogramBuckets = 10;
static const size_t c_numTests = 1000; // how many times to do the tests, to get average and std dev

// what points in the input stream to dump a report
static const size_t c_reportValues[] =
{
    100,
    1000,
    10000
};
static const size_t c_numReportValues = COUNT_OF(c_reportValues);
static const size_t c_maxReportValue = c_reportValues[c_numReportValues - 1];

// ===================================== Constants =====================================

typedef std::array<float, 2> Vec2;

static const float c_goldenRatio = 1.61803398875f;
static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_root2 = 1.41421356237f;
static const float c_fractRoot2 = 0.41421356237f;

// ===================================== PDFs =====================================

namespace PDF
{
    template <typename T>
    struct Base
    {
        // Probability Function
        //
        // Returns a probability 0 to 1 for this value.
        static float PF(float x)
        {
            constexpr Vec2 MinMax = T::FuncMinMax();
            return T::Func(x) / MinMax[1];
        }

        // Cumulative Probability Function
        //
        // Returns the sum of the probabilities between x1 and x2.
        // Useful for getting the actual value a histogram bucket should have and
        // also useful for figuring out how many samples should survive rejection sampling.
        static float CumulativePF(float x1, float x2)
        {
            constexpr Vec2 MinMax = T::FuncMinMax();
            return (T::IntegratedFunc(x2) - T::IntegratedFunc(x1)) / MinMax[1];
        }

        // Inverse Probability Function
        //
        // If you want to "undo" a transformation to a pdf, you need the inverse
        // of the probability function.  This is that.
        static float InversePF(float x)
        {
            constexpr Vec2 MinMax = T::FuncMinMax();
            return MinMax[0] / T::Func(x);
        }

        // structs that derive from this class should have these functions:

        // A function that takes in an x value 0 to 1 and returns a value. Higher values are more probable.
        //static constexpr float Func(float x)

        // The minimum and maximum value that func can return between 0 and 1.
        // Useful for normalizing Func, it's inverse, and cumulative PF.
        //static constexpr Vec2 FuncMinMax()

        // indefinite integral of Func()
        //static constexpr float IntegratedFunc(float x)
    };

    // y = x+0.1
    struct Y_Eq_X_P_0_1 : public PDF::Base<Y_Eq_X_P_0_1>
    {
    private:
        friend PDF::Base<Y_Eq_X_P_0_1>;

        static constexpr float Func(float x)  
        {
            return x + 0.1f;
        }

        static constexpr Vec2 FuncMinMax()
        {
            return Vec2{ 0.1f, 1.1f };
        }

        static constexpr float IntegratedFunc(float x)
        {
            return 0.5f * x * x + 0.1f * x;
        }
    };
};

// ===================================== PMFs =====================================

// ===================================== Utils =====================================

std::mt19937 GetRNG(int seed)
{
#if DETERMINISTIC()
    std::mt19937 rng(seed);
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif
    return rng;
}

float fract(float f)
{
    return f - floor(f);
}

int FloatToItem(float f, int numItems)
{
    // if f is in [0,1], remaps to [0, numItems-1]
    return std::min(int(f * float(numItems)), numItems - 1);
}

template <size_t N>
void MakeHistogram(const std::vector<float>& stream, size_t streamCount, std::array<float, N>& histogram)
{
    std::fill(histogram.begin(), histogram.end(), 0.0f);

    for (size_t index = 0; index < streamCount; ++index)
    {
        float y = stream[index];
        histogram[FloatToItem(y, N)] += 1.0f / float(streamCount);
    }
}

float Lerp(float a, float b, float t)
{
    return a * (1.0f - t) + b * t;
}

Vec2 R2(const Vec2& input)
{
    // generalized golden ratio, from:
    // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    const float g = 1.32471795724474602596f;
    const float a1 = 1.0f / g;
    const float a2 = 1.0f / (g * g);

    return Vec2{
        fract(input[0] + a1),
        fract(input[1] + a2),
    };
}

static size_t Ruler(size_t n)
{
    size_t ret = 0;
    while (n != 0 && (n & 1) == 0)
    {
        n /= 2;
        ++ret;
    }
    return ret;
}

void Sobol(std::vector<Vec2>& values, size_t numValues)
{
    // x axis
    values.resize(numValues);
    size_t sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i)
    {
        size_t ruler = Ruler(i + 1);
        size_t direction = size_t(size_t(1) << size_t(31 - ruler));
        sampleInt = sampleInt ^ direction;
        values[i][0] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }

    // y axis
    // Code adapted from http://web.maths.unsw.edu.au/~fkuo/sobol/
    // uses numbers: new-joe-kuo-6.21201

    // Direction numbers
    std::vector<size_t> V;
    V.resize((size_t)ceil(log((double)numValues + 1) / log(2.0)));  //+1 because we are skipping index 0
    V[0] = size_t(1) << size_t(31);
    for (size_t i = 1; i < V.size(); ++i)
        V[i] = V[i - 1] ^ (V[i - 1] >> 1);

    // Samples
    sampleInt = 0;
    for (size_t i = 0; i < numValues; ++i) {
        size_t ruler = Ruler(i + 1);
        sampleInt = sampleInt ^ V[ruler];
        values[i][1] = float(sampleInt) / std::pow(2.0f, 32.0f);
    }
}

struct TestReportItem
{
    std::array<float, c_numHistogramBuckets> outputHistogram;
    std::array<float, c_numHistogramBuckets> outputHistogramAvg;
    std::array<float, c_numHistogramBuckets> outputHistogramSqdAvg;

    float survivedError;
    float survivedErrorAvg;
    float survivedErrorSqdAvg;
};

struct TestReport
{
    std::array<TestReportItem, c_numReportValues> item;
    std::string inputStreamType;
    std::string rngStreamType;
};

template <typename PDF>
void SaveTestReport(const std::vector<TestReport>& testReports, const char* baseFileName)
{
    char buffer[256];

    // show the survival error
    {
        FILE* file = nullptr;
        sprintf_s(buffer, "%s_survival.csv", baseFileName);
        fopen_s(&file, buffer, "w+t");

        // Show the error
        {
            fprintf(file, "\"Counts\"");
            for (int reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
                fprintf(file, ",\"%zu\"", c_reportValues[reportIndex]);
            fprintf(file, "\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s error\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (const TestReportItem& item : report.item)
                    fprintf(file, ",\"%f\"", abs(item.survivedError));
                fprintf(file, "\n");
            }
        }

        // Show the average error
        {
            fprintf(file, "\n\"Average Error\"\n\n\"Counts\"");
            for (int reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
                fprintf(file, ",\"%zu\"", c_reportValues[reportIndex]);
            fprintf(file, "\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s error\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (const TestReportItem& item : report.item)
                    fprintf(file, ",\"%f\"", abs(item.survivedErrorAvg));
                fprintf(file, "\n");
            }
        }

        // Show the error std dev
        {
            fprintf(file, "\n\"Error Std Dev\"\n\n\"Counts\"");
            for (int reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
                fprintf(file, ",\"%zu\"", c_reportValues[reportIndex]);
            fprintf(file, "\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s error\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (const TestReportItem& item : report.item)
                {
                    float variance = abs(item.survivedErrorSqdAvg - item.survivedErrorAvg * item.survivedErrorAvg);
                    float stddev = sqrt(variance);
                    fprintf(file, ",\"%f\"", stddev);
                }
                fprintf(file, "\n");
            }
        }

        fclose(file);
    }

    // show the histograms
    for (int reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
    {
        FILE* file = nullptr;
        sprintf_s(buffer, "%s_%i.csv", baseFileName, int(c_reportValues[reportIndex]));
        fopen_s(&file, buffer, "w+t");

        // show actual histogram
        {
            fprintf(file, "\"expected\"");
            for (size_t histogramIndex = 0; histogramIndex < c_numHistogramBuckets; ++histogramIndex)
            {
                float x1 = float(histogramIndex) / float(c_numHistogramBuckets);
                float x2 = float(histogramIndex + 1) / float(c_numHistogramBuckets);
                float expected = PDF::CumulativePF(x1, x2) / PDF::CumulativePF(0.0f, 1.0f);

                fprintf(file, ",\"%f\"", expected);
            }
            fprintf(file, "\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (float f : report.item[reportIndex].outputHistogram)
                    fprintf(file, ",\"%f\"", f);
                fprintf(file, "\n");
            }
        }

        // show histogram error
        {
            fprintf(file, "\n\"Error\"\n\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (size_t histogramIndex = 0; histogramIndex < c_numHistogramBuckets; ++histogramIndex)
                {
                    float x1 = float(histogramIndex) / float(c_numHistogramBuckets);
                    float x2 = float(histogramIndex + 1) / float(c_numHistogramBuckets);
                    float expected = PDF::CumulativePF(x1, x2) / PDF::CumulativePF(0.0f, 1.0f);

                    float f = report.item[reportIndex].outputHistogram[histogramIndex] - expected;
                    fprintf(file, ",\"%f\"", f);
                }
                fprintf(file, "\n");
            }
        }

        // show histogram average error
        {
            fprintf(file, "\n\"Average Error\"\n\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (size_t histogramIndex = 0; histogramIndex < c_numHistogramBuckets; ++histogramIndex)
                {
                    float x1 = float(histogramIndex) / float(c_numHistogramBuckets);
                    float x2 = float(histogramIndex + 1) / float(c_numHistogramBuckets);
                    float expected = PDF::CumulativePF(x1, x2) / PDF::CumulativePF(0.0f, 1.0f);

                    float f = report.item[reportIndex].outputHistogramAvg[histogramIndex] - expected;
                    fprintf(file, ",\"%f\"", f);
                }
                fprintf(file, "\n");
            }
        }

        // show histogram error std dev
        {
            fprintf(file, "\n\"Error Std Deviation\"\n\n");

            for (const TestReport& report : testReports)
            {
                fprintf(file, "\"stream %s rng %s\"", report.inputStreamType.c_str(), report.rngStreamType.c_str());
                for (size_t histogramIndex = 0; histogramIndex < c_numHistogramBuckets; ++histogramIndex)
                {
                    float x1 = float(histogramIndex) / float(c_numHistogramBuckets);
                    float x2 = float(histogramIndex + 1) / float(c_numHistogramBuckets);
                    float expected = PDF::CumulativePF(x1, x2) / PDF::CumulativePF(0.0f, 1.0f);

                    float f = abs(report.item[reportIndex].outputHistogramSqdAvg[histogramIndex] - (report.item[reportIndex].outputHistogramAvg[histogramIndex] * report.item[reportIndex].outputHistogramAvg[histogramIndex]));
                    f = sqrt(f);

                    fprintf(file, ",\"%f\"", f);
                }
                fprintf(file, "\n");
            }
        }

        fclose(file);
    }
}

template <typename PDF>
std::vector<float> RejectionSample(const std::vector<float>& inputStream, const std::vector<float>& rngStream, TestReport& report, int testIndex)
{
    // rejection sample to convert input stream into desired PDF
    std::vector<float> transformed;
    for (size_t inputIndex = 0; inputIndex < inputStream.size(); ++inputIndex)
    {
        float x = inputStream[inputIndex];
        if (rngStream[inputIndex] < PDF::PF(x))
            transformed.push_back(x);

        for (int reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            if (inputIndex + 1 == c_reportValues[reportIndex])
            {
                float alpha = 1.0f / float(testIndex + 1);

                // calculate the error of survival rate
                float survivedPercent = float(transformed.size()) / float(c_reportValues[reportIndex]);
                float survivedPercentExpected = PDF::CumulativePF(0.0f, 1.0f);
                report.item[reportIndex].survivedError = survivedPercent - survivedPercentExpected;

                // calculate the avg and squared avg of survival rate
                report.item[reportIndex].survivedErrorAvg = Lerp(report.item[reportIndex].survivedErrorAvg, report.item[reportIndex].survivedError, alpha);
                report.item[reportIndex].survivedErrorSqdAvg = Lerp(report.item[reportIndex].survivedErrorSqdAvg, report.item[reportIndex].survivedError*report.item[reportIndex].survivedError, alpha);

                // calculate histogram
                MakeHistogram(transformed, transformed.size(), report.item[reportIndex].outputHistogram);

                // calculate the avg and squared avg of histograms
                for (size_t histogramIndex = 0; histogramIndex < c_numHistogramBuckets; ++histogramIndex)
                {
                    report.item[reportIndex].outputHistogramAvg[histogramIndex] = Lerp(report.item[reportIndex].outputHistogramAvg[histogramIndex], report.item[reportIndex].outputHistogram[histogramIndex], alpha);
                    report.item[reportIndex].outputHistogramSqdAvg[histogramIndex] = Lerp(report.item[reportIndex].outputHistogramSqdAvg[histogramIndex], report.item[reportIndex].outputHistogram[histogramIndex] * report.item[reportIndex].outputHistogram[histogramIndex], alpha);
                }

                break;
            }
        }
    }
    return transformed;
}

enum ESequenceType : int
{
    WhiteWhite,
    GRWhite,
    WhiteSqrt2,
    GRSqrt2,
    SeqR2,
    SeqSobol,

    Count
};

void GenerateSequence(std::vector<float>& inputStream, std::vector<float>& rngStream, size_t numValues, ESequenceType sequenceType, std::mt19937& rng, TestReport& report, Vec2& LDSOffset)
{
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
 
    inputStream.resize(numValues);
    rngStream.resize(numValues);

    if (sequenceType == ESequenceType::SeqSobol)
    {
        report.inputStreamType = "Sobol";
        report.rngStreamType = "Sobol";

        std::vector<Vec2> inputAndRngStream(numValues);
        Sobol(inputAndRngStream, numValues);

        for (Vec2& f : inputAndRngStream)
        {
            f[0] = fract(f[0] + LDSOffset[0]);
            f[1] = fract(f[1] + LDSOffset[1]);
        }

        for (size_t index = 0; index < numValues; ++index)
        {
            inputStream[index] = inputAndRngStream[index][0];
            rngStream[index] = inputAndRngStream[index][1];
        }

    }
    else if (sequenceType == ESequenceType::SeqR2)
    {
        report.inputStreamType = "R2";
        report.rngStreamType = "R2";

        for (size_t index = 0; index < numValues; ++index)
        {
            LDSOffset = R2(LDSOffset);
            inputStream[index] = LDSOffset[0];
            rngStream[index] = LDSOffset[1];
        }
    }
    else
    {
        // make uniform input samples - either white noise or LDS
        if ((sequenceType & 1) == 0)
        {
            report.inputStreamType = "white";
            for (float& f : inputStream)
                f = dist(rng);
        }
        else
        {
            report.inputStreamType = "GR";
            for (float& f : inputStream)
            {
                LDSOffset[0] = fract(LDSOffset[0] + c_goldenRatioConjugate);
                f = LDSOffset[0];
            }
        }

        // make uniform rng samples 0 either white noise or LDS
        if (((sequenceType >> 1) & 1) == 0)
        {
            report.rngStreamType = "white";
            for (float& f : rngStream)
                f = dist(rng);
        }
        else
        {
            report.rngStreamType = "sqrt2";
            for (float& f : rngStream)
            {
                LDSOffset[1] = fract(LDSOffset[1] + c_fractRoot2);
                f = LDSOffset[1];
            }
        }
    }
}

// ===================================== Program =====================================

int main(int argc, char** argv)
{
    // Test 1 "simple transformation" - uniform to y=x+0.1
    {
        std::vector<TestReport> testReports;

        int flatTestIndex = -1;
        int percent = 0;
        printf("Test 1: uniform to y=x+0.1\n");
        for (int sequenceType = 0; sequenceType < ESequenceType::Count; ++sequenceType)
        {
            // get a test report to store data to report later
            testReports.emplace_back();
            TestReport& report = testReports.back();

            for (int testIndex = 0; testIndex < c_numTests; ++testIndex)
            {
                flatTestIndex++;

                // report percentage done
                int newPercent = int(100.0f * float(flatTestIndex) / float(c_numTests * ESequenceType::Count));
                if (percent != newPercent)
                    printf("\r%i%%", newPercent);
                percent = newPercent;

                // get white noise random number generator
                std::mt19937 rng = GetRNG(flatTestIndex);
                std::uniform_real_distribution<float> dist(0.0f, 1.0f);

                // Generate the input stream and rng stream
                Vec2 LDSOffset = Vec2{ 0.0f, 0.0f };
                #if RANDOMIZE_LDS_START()
                LDSOffset = Vec2{ dist(rng), dist(rng) };
                #endif
                std::vector<float> inputStream, rngStream;
                GenerateSequence(inputStream, rngStream, c_maxReportValue, (ESequenceType)sequenceType, rng, report, LDSOffset);

                // Do rejection sampling to convert from uniform to y=x+0.1
                std::vector<float> transformed = RejectionSample<PDF::Y_Eq_X_P_0_1>(inputStream, rngStream, report, testIndex);
            }
        }
        printf("\r100%%\n\n");

        // write out report CSVs
        SaveTestReport<PDF::Y_Eq_X_P_0_1>(testReports, "out/test1");
    }

    // TODO: test 2 is to do a multi step PDF transformation
    // TODO: test 3 is to combine the multiple steps into 1 to show survival rate etc

    system("pause");
    return 0;
}

/*

----- REJECTION SAMPLING BLOG:

! make a new repository for rejection sampling when this is done

- transmutting one pdf or pmf to another by throwing samples away.
 * you use a probability function, not pdf or pmf. difference is in normalization. PF has at most 1.0, but can be smaller.
  * want to minimize the count you throw away though.
 * could talk about how to make PDF's and PMFs from formulas
 * Full deterministic LDS has bad average error, because others have non deterministic error so will sometimes be higher or lower and over time will average to zero. FULL LDS is stuck at one sample.
 * maybe show the randomized LDS results, and also mention how you'd probably want a different LDS driving the LDS starting values instead :P
 * could talk about how the area of the square that is empty shows what % of the samples that will be lost
 * could show how putting non uniform PDF in doesn't give the right result back (intro to inverting)
 * inverting is doing 1/PF and making sure the max is < 1.
  * can't invert a probability of 0 unforuntately.
  * PDF = divide by zero then.
  * PMF = same, but could say "if you see one of these values, always take it". But, the values won't be there.
  * could "remap" values instead of throwing them away but (like inverting CDF. link to post) but that isn't rejection sampling.
  ! also relates to low probabilities in the source distribution.
   * if you have a large probability for A and a very small proabbility for B, invertint it will flop that.
   * since you then need a lot of B's you are going to have to skip a lot of A's to have the right amount
   * inverting a zero chance would mean throwing all of the list away. you still wouldn't get that value though hehe.
 * show a multi step thing like uniform to y=x back to uniform back to something else?
 * show how you can combine steps by multiplying PF together.
 * do this with PMFs and PDFs both.
 * using all 4 permutations of LDS vs white noise
  * contextualize the 4 cases.
  * white stream / white rng = naive setup?
  * white stream / LDS rng = you have some data you don't control, you want to transform.
  * LDS stream / white rng = ??
  * LDS stream / LDS rng = trying to generate a good sequence?

- can use rejection sampling with uniform PDF to choose ~N items from a long list. N is variable but scaling the PDF for the PF lets you choose how much you expect to get out.
 * could also random roll some indices. that migh thave repeats
 * could also use format preserving encryption to do a storageless shuffle and pick the first N.
 * only talking about this because it's part of the next post


----- old experiments (still valid?)

* can make a stream of values from another stream of values.
* goal: show how using LDS as RNG makes a better result.

! experiment:
 * use rejection sampling to make a desired pdf.
 * compare white noise vs LDS.
  * LDS as both input and as the rng?
 * try a couple different input and output PDFs
 ! apples to apples w/ reservoir sampling. see which does a better job!

Test 1 - uniform pdf -> rejected sampled to -> y=2x pdf
00 -> white noise as input, white noise rng
01 -> white noise as input, LDS rng
10 -> LDS as input, white noise as rng
11 -> LDS as input, LDS as rng.
! do this multiple times, show error and variance of histogram. at different sample counts too.

Test 2 & 3 - different input and output pdfs

*/