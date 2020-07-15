#include <random>
#include <vector>

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

static const float c_goldenRatio = 1.61803398875f;
static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_root2 = 1.41421356237f;
static const float c_fractRoot2 = 0.41421356237f;

// ===================================== Utils =====================================

std::mt19937 GetRNG()
{
#if DETERMINISTIC()
    static int seed = 1336;
    seed++;
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

float Lerp(float a, float b, float t)
{
    return a * (1.0f - t) + b * t;
}

// ===================================== PDFs =====================================

namespace PDF
{
    // PDF: y = 1
    struct UniformWhite
    {
        UniformWhite()
        {
            m_rng = GetRNG();
        }

        float Generate()
        {
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            return dist(m_rng);
        }

        std::mt19937 m_rng;
    };

    // PDF: y = 1
    struct UniformLDS_GR
    {
        UniformLDS_GR()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        float Generate()
        {
            m_value = fract(m_value + c_goldenRatioConjugate);
            return m_value;
        }

        float m_value;
    };

    // PDF: y = 1
    struct UniformLDS_Root2
    {
        UniformLDS_Root2()
        {
            #if RANDOMIZE_LDS_START()
            std::mt19937 rng = GetRNG();
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            m_value = dist(rng);
            #else
            m_value = 0.0f;
            #endif
        }

        float Generate()
        {
            m_value = fract(m_value + c_fractRoot2);
            return m_value;
        }

        float m_value;
    };

    // PDF:  y = (2x + 3) / 4
    // CDF:  y = (x^2 + 3x) / 4 
    // PF:   y = (2x + 3) / 5
    template <typename Generator, typename Validator>
    struct Linear
    {
        static float CDF(float x)
        {
            return (x*x + 3.0f*x) / 4.0f;
        }

        static float PF(float x)
        {
            return (2.0f * x + 3.0f) / 5.0f;
        }

        float Generate()
        {
            while (true)
            {
                float x = m_generator.Generate();
                float probability = PF(x);
                if (m_validator.Generate() <= probability)
                    return x;
            }
        }

        Generator m_generator;
        Validator m_validator;
    };

    // PDF:  y = (x^3 -10x^2 + 5x + 11) / 10.417
    // CDF:  y = 0.0959969  * (11x + 2.5x^2 - 3.33333x^3 + 0.25x^4)
    // PF:   y = (x^3 -10x^2 + 5x + 11) / 12
    template <typename Generator, typename Validator>
    struct Cubic
    {
        static float CDF(float x)
        {
            return 0.0959969f * (11.0f * x + 2.5f*x*x - 3.33333f*x*x*x + 0.25f * x*x*x*x);
        }

        static float PF(float x)
        {
            return (x*x*x - 10.0f * x*x + 5.0f * x + 11.0f) / 12.0f;
        }

        float Generate()
        {
            while (true)
            {
                float x = m_generator.Generate();
                float probability = PF(x);
                if (m_validator.Generate() <= probability)
                    return x;
            }
        }

        Generator m_generator;
        Validator m_validator;
    };
};

// ===================================== Code =====================================

template <typename TPDF>
void GenerateSequence(std::vector<float>& sequence, size_t count)
{
    TPDF pdf;
    sequence.resize(count);
    for (float& f : sequence)
        f = pdf.Generate();
}

void CalculateHistogram(const std::vector<float>& samples, size_t sampleCount, std::vector<float>& histogram)
{
    histogram.resize(c_numHistogramBuckets);
    std::fill(histogram.begin(), histogram.end(), 0.0f);
    for (size_t index = 0; index < sampleCount; ++index)
        histogram[FloatToItem(samples[index], c_numHistogramBuckets)] += 1.0f / float(sampleCount);
}

void WriteHistogram(FILE* file, const char* label, const std::vector<float>& histogram)
{
    fprintf(file, "\"%s\"", label);
    for (float f : histogram)
        fprintf(file, ",\"%f\"", f);
    fprintf(file, "\n");
}

void WriteHistogramError(FILE* file, const char* label, const std::vector<float>& histogram, const std::vector<float>& expectedHistogram)
{
    fprintf(file, "\"%s\"", label);
    for (size_t i = 0; i < c_numHistogramBuckets; ++i)
        fprintf(file, ",\"%f\"", histogram[i] - expectedHistogram[i]);
    fprintf(file, "\n");
}

void WriteHistogramStdDev(FILE* file, const char* label, const std::vector<float>& histogramAvg, const std::vector<float>& histogramSqAvg)
{
    fprintf(file, "\"%s\"", label);
    for (size_t i = 0; i < c_numHistogramBuckets; ++i)
        fprintf(file, ",\"%f\"", sqrt(abs(histogramSqAvg[i] - histogramAvg[i] * histogramAvg[i])));
    fprintf(file, "\n");
}

void HistogramCombine(const std::vector<float>& histogram, std::vector<float>& histogramAvg, std::vector<float>& histogramSqAvg, size_t sampleIndex)
{
    if (sampleIndex == 0)
    {
        histogramAvg.resize(c_numHistogramBuckets, 0.0f);
        histogramSqAvg.resize(c_numHistogramBuckets, 0.0f);
    }

    for (size_t index = 0; index < c_numHistogramBuckets; ++index)
    {
        float value = histogram[index];
        histogramAvg[index] = Lerp(histogramAvg[index], value, 1.0f / float(sampleIndex + 1));
        histogramSqAvg[index] = Lerp(histogramSqAvg[index], value*value, 1.0f / float(sampleIndex + 1));
    }
}

int main(int argc, char ** argv)
{
    // uniform to linear
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][4];
        std::vector<float> histogramAvg[c_numReportValues][4];
        std::vector<float> histogramSqAvg[c_numReportValues][4];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[4];
            GenerateSequence<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>>(samples[0], c_maxReportValue);
            GenerateSequence<PDF::Linear<PDF::UniformWhite, PDF::UniformLDS_GR>>(samples[1], c_maxReportValue);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformWhite>>(samples[2], c_maxReportValue);
            GenerateSequence<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformLDS_GR>>(samples[3], c_maxReportValue);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
            }
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/uni_lin_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);

            fclose(file);
        }
    }

    // [uniform to linear] to cubic
    {
        // calculate the expected histogram
        std::vector<float> expectedHistogram;
        for (size_t bucketIndex = 0; bucketIndex < c_numHistogramBuckets; ++bucketIndex)
        {
            float p0 = float(bucketIndex) / float(c_numHistogramBuckets);
            float p1 = float(bucketIndex + 1) / float(c_numHistogramBuckets);

            float cdf0 = PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>::CDF(p0);
            float cdf1 = PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>::CDF(p1);

            expectedHistogram.push_back(cdf1 - cdf0);
        }

        std::vector<float> histogram[c_numReportValues][4];
        std::vector<float> histogramAvg[c_numReportValues][4];
        std::vector<float> histogramSqAvg[c_numReportValues][4];

        for (size_t testIndex = 0; testIndex < c_numTests; ++testIndex)
        {
            // generate the samples
            std::vector<float> samples[4];
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformWhite>, PDF::UniformWhite>>(samples[0], c_maxReportValue);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformWhite, PDF::UniformLDS_GR>, PDF::UniformLDS_GR>>(samples[1], c_maxReportValue);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformWhite>, PDF::UniformWhite>>(samples[2], c_maxReportValue);
            GenerateSequence<PDF::Cubic<PDF::Linear<PDF::UniformLDS_Root2, PDF::UniformLDS_GR>, PDF::UniformLDS_GR>>(samples[3], c_maxReportValue);

            // calculate data for each report
            for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
            {
                size_t sampleCount = c_reportValues[reportIndex];

                // calculate the histograms
                CalculateHistogram(samples[0], sampleCount, histogram[reportIndex][0]);
                CalculateHistogram(samples[1], sampleCount, histogram[reportIndex][1]);
                CalculateHistogram(samples[2], sampleCount, histogram[reportIndex][2]);
                CalculateHistogram(samples[3], sampleCount, histogram[reportIndex][3]);

                // combine the histograms
                HistogramCombine(histogram[reportIndex][0], histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0], testIndex);
                HistogramCombine(histogram[reportIndex][1], histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1], testIndex);
                HistogramCombine(histogram[reportIndex][2], histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2], testIndex);
                HistogramCombine(histogram[reportIndex][3], histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3], testIndex);
            }
        }

        // make the reports
        for (size_t reportIndex = 0; reportIndex < c_numReportValues; ++reportIndex)
        {
            size_t sampleCount = c_reportValues[reportIndex];

            // open the file
            char buffer[256];
            sprintf_s(buffer, "out/lin_cub_%zu.csv", sampleCount);
            FILE* file = nullptr;
            fopen_s(&file, buffer, "w+t");

            // write the expected histogram
            fprintf(file, "\"Expected\"");
            for (float f : expectedHistogram)
                fprintf(file, ",\"%f\"", f);
            fprintf(file, "\n");

            // write histograms
            WriteHistogram(file, "white/white", histogram[reportIndex][0]);
            WriteHistogram(file, "white/LDS", histogram[reportIndex][1]);
            WriteHistogram(file, "LDS/white", histogram[reportIndex][2]);
            WriteHistogram(file, "LDS/LDS", histogram[reportIndex][3]);

            // write error
            fprintf(file, "\n\"Error:\"\n");
            WriteHistogramError(file, "white/white", histogram[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogram[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogram[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogram[reportIndex][3], expectedHistogram);

            // write average error
            fprintf(file, "\n\"Avg Error:\"\n");
            WriteHistogramError(file, "white/white", histogramAvg[reportIndex][0], expectedHistogram);
            WriteHistogramError(file, "white/LDS", histogramAvg[reportIndex][1], expectedHistogram);
            WriteHistogramError(file, "LDS/white", histogramAvg[reportIndex][2], expectedHistogram);
            WriteHistogramError(file, "LDS/LDS", histogramAvg[reportIndex][3], expectedHistogram);

            // write std dev
            fprintf(file, "\n\"Std Dev:\"\n");
            WriteHistogramStdDev(file, "white/white", histogramAvg[reportIndex][0], histogramSqAvg[reportIndex][0]);
            WriteHistogramStdDev(file, "white/LDS", histogramAvg[reportIndex][1], histogramSqAvg[reportIndex][1]);
            WriteHistogramStdDev(file, "LDS/white", histogramAvg[reportIndex][2], histogramSqAvg[reportIndex][2]);
            WriteHistogramStdDev(file, "LDS/LDS", histogramAvg[reportIndex][3], histogramSqAvg[reportIndex][3]);

            fclose(file);
        }
    }

    // TODO: need to account for source PDF when rejection sampling! make an interface

    return 0;
}

/*

TODO:

* work with both PDFs and PMFs both to show it working.
* do from uniform and also not from uniform.
* four combinations: white/white  white/lds  lds/white  lds/lds
* show histogram for a couple different counts: 100, 10k, 1m
* maybe do white noise tests multiple times to get average and std dev? could also show a single run.
? randomize the start of the LDS to make the multiple tests mneaningful?


Post:
* need to explain difference between PDF and PF?

PDFs:
* uniform             : done
* linear              : y = (2x + 3) / 4
* quadratic or cubic  : y = (x^3 -10x^2 + 5x + 11) / 10.417

PMFs:
???

*/