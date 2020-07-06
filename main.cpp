#include <vector>
#include <random>
#include <array>
#include <string>

#define DETERMINISTIC() true 

static const size_t c_numHistogramBuckets = 10;

typedef std::array<float, 2> Vec2;
typedef std::vector<std::vector<std::string>> CSV;

float fract(float f)
{
    return f - floor(f);
}

int FloatToItem(float f, int numItems)
{
    // if f is in [0,1], remaps to [0, numItems-1]
    return std::min(int(f * float(numItems)), numItems - 1);
}

enum SeedContext : int
{
    RejectionSampling = 1,
    ReservoirSampling = 2
};

void RejectionSamplingMain(void);

// ===================================== PDFs =====================================

/*

We are assuming x is in [0,1]

PDF - integrates to 1.0, >= 0 in the domain

PF - PDF scaled so that all y axis values are <= 1

CDF - indefinite integral of PDF

*/

namespace PDF
{
    struct Uniform
    {
        static const char* Label()
        {
            return "Uniform";
        }

        static const char* FileLabel()
        {
            return "Uniform";
        }

        static float PDF(float x)
        {
            // PDF is y=1
            return 1.0;
        }

        static float PF(float x)
        {
            // probability function can be y = 1 too
            return 1.0f;
        }

        static float CDF(float x)
        {
            // PDF integrates to y = x
            return x;
        }
    };

    struct Y_Equals_X
    {
        static const char* Label()
        {
            return "y=x";
        }

        static const char* FileLabel()
        {
            return "yeqx";
        }

        static float PDF(float x)
        {
            // y=x normalizes to y=2x
            return 2.0f * x;
        }

        static float PF(float x)
        {
            // probability function can be y=x
            return x;
        }

        static float CDF(float x)
        {
            // PDF integrates to y=x^2
            return x * x;
        }
    };

    struct Y_Equals_1_Minus_X
    {
        static const char* Label()
        {
            return "y=1-x";
        }

        static const char* FileLabel()
        {
            return "yeq1minx";
        }

        static float PDF(float x)
        {
            // y=1-x normalizes to y=2(1-x)
            return 2.0f * (1.0f - x);
        }

        static float PF(float x)
        {
            // probability function can be y=1-x
            return (1-x);
        }

        static float CDF(float x)
        {
            // PDF integrates to y=2x-x^2
            return 2.0f * x - x * x;
        }
    };

    struct Y_Equals_X_Minus_X_squared
    {
        static const char* Label()
        {
            return "y=x-x^2";
        }

        static const char* FileLabel()
        {
            return "yeqxminxsq";
        }

        static float PDF(float x)
        {
            // y=x normalizes to y=6*(x-x^2)
            return 6.0f * (x - x * x);
        }

        static float PF(float x)
        {
            // probability function can be y=4*(x-x^2)
            return 4.0f * (x - x * x);
        }

        static float CDF(float x)
        {
            // PDF integrates to y=3x^2-2x^3
            return 3.0f * x * x - 2.0f * x * x * x;
        }
    };
};

// ===================================== SEQUENCE GENERATORS =====================================

namespace Seq
{
    struct WhiteNoise2D
    {
        WhiteNoise2D(int seed)
        {
            #if DETERMINISTIC()
                std::mt19937 rng(seed);
            #else
                std::random_device rd;
                std::mt19937 rng(rd());
            #endif
            m_rng = rng;
        }

        Vec2 GetFloats01()
        {
            std::uniform_real_distribution<float> dist(0.0f, 1.0f);
            return Vec2{ dist(m_rng), dist(m_rng) };
        }

        std::mt19937 m_rng;
    };

    struct R2
    {
        R2(int seed)
        {

        }

        Vec2 GetFloats01()
        {
            // generalized golden ratio, from:
            // http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
            const float g = 1.32471795724474602596f;
            const float a1 = 1.0f / g;
            const float a2 = 1.0f / (g * g);

            m_state[0] = fract(m_state[0] + a1);
            m_state[1] = fract(m_state[1] + a2);

            return m_state;
        }

        Vec2 m_state = Vec2{ 0.0f, 0.0f };
    };
};

// ===================================== Resevoir Sampling =====================================

#pragma optimize("", off)

// Algorithm A from M. T. Chao
// https://en.wikipedia.org/wiki/Reservoir_sampling#Algorithm_A-Chao
template <typename SEQ, typename PDF>
std::vector<float> ReservoirSample(const std::vector<float>& stream, size_t desiredSamples)
{
    // initialize the returned samples with the first samples from the stream
    std::vector<float> ret(desiredSamples, 0.0f);
    float weightSum = 0.0f;
    for (size_t i = 0; i < desiredSamples; ++i)
    {
        float streamSample = stream[i];
        weightSum += PDF::PDF(streamSample);
        ret[i] = streamSample;
    }

    int accept = 0;
    int reject = 0;

    // now accept samples based on weight and probabilities, and randomly replace existing samples
    SEQ seq(SeedContext::ReservoirSampling);
    for (size_t i = desiredSamples; i < stream.size(); ++i)
    //weightSum = 0.0f;
    //for (size_t i = 0; i < stream.size(); ++i)  // TODO: temp. if leaving it, remove the code above?
    {
        float streamSample = stream[i];
        float weight = PDF::PDF(streamSample);
        weightSum += weight;

#if 1
        for (size_t j = 0; j < desiredSamples; ++j)
        {
            Vec2 p = seq.GetFloats01(); // TODO: only need 1d...
            if (p[0] < weight / weightSum)
            {
                accept++;
                ret[j] = streamSample;
            }
            else
            {
                reject++;
            }
        }
#else

        Vec2 p = seq.GetFloats01();
        if (p[0] < weight / weightSum)
        {
            accept++;
            ret[FloatToItem(p[1], int(desiredSamples))] = streamSample;
        }
        else
        {
            reject++;
        }
#endif
    }

    printf("accept = %i reject = %i  (%f%% accepted)\n", accept, reject, 100.0f * float(accept) / float(accept + reject));

    return ret;
}

#pragma optimize("", on)

// ===================================== Rejection Sampling =====================================

template <typename SEQ, typename PDF>
std::vector<float> RejectionSample(size_t desiredSamples)
{
    SEQ seq(SeedContext::RejectionSampling);

    std::vector<float> ret;
    while (ret.size() < desiredSamples)
    {
        Vec2 p = seq.GetFloats01();

        float pdf = PDF::PF(p[0]);
        if (p[1] <= pdf)
            ret.push_back(p[0]);
    }
    return ret;
}

// ===================================== CSV =====================================

void SaveCSV(const CSV& csv, const char* fileName)
{
    FILE* file = nullptr;
    fopen_s(&file, fileName, "w+t");

    for (auto& row : csv)
    {
        bool firstCell = true;
        for (auto& cell : row)
        {
            fprintf(file, "%s\"%s\"", firstCell ? "" : ",", cell.c_str());
            firstCell = false;
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

template <typename PDF>
void AddHistogramPDF(CSV& csv)
{
    size_t rowIndex = csv.size();
    csv.resize(rowIndex + 1);
    std::vector<std::string>& row = csv[rowIndex];

    row.push_back("PDF");

    char buffer[1024];
    for (int i = 0; i < c_numHistogramBuckets; ++i)
    {
        float p0 = float(i) / float(c_numHistogramBuckets);
        float p1 = float(i + 1) / float(c_numHistogramBuckets);

        sprintf_s(buffer, "%f", PDF::CDF(p1) - PDF::CDF(p0));
        row.push_back(buffer);
    }
}

void AddHistogram(CSV& csv, const char* label, const std::vector<float>& sequence, int count)
{
    std::vector<float> histogram;
    histogram.resize(c_numHistogramBuckets, 0.0f);

    for (int i = 0; i < count; ++i)
        histogram[FloatToItem(sequence[i], c_numHistogramBuckets)] += 1.0f;

    for (float& f : histogram)
        f /= float(count);

    size_t rowIndex = csv.size();
    csv.resize(rowIndex + 1);
    std::vector<std::string>& row = csv[rowIndex];

    row.push_back(label);

    char buffer[1024];
    for (int i = 0; i < c_numHistogramBuckets; ++i)
    {
        sprintf_s(buffer, "%f", histogram[i]);
        row.push_back(buffer);
    }
}

template <typename PDF>
void AddHistogramError(CSV& csv, const char* label, const std::vector<float>& sequence, int count)
{
    std::vector<float> histogram;
    histogram.resize(c_numHistogramBuckets, 0.0f);

    for (int i = 0; i < count; ++i)
        histogram[FloatToItem(sequence[i], c_numHistogramBuckets)] += 1.0f;

    for (float& f : histogram)
        f /= float(count);

    size_t rowIndex = csv.size();
    csv.resize(rowIndex + 1);
    std::vector<std::string>& row = csv[rowIndex];

    row.push_back(label);

    char buffer[1024];
    for (int i = 0; i < c_numHistogramBuckets; ++i)
    {
        float p0 = float(i) / float(c_numHistogramBuckets);
        float p1 = float(i + 1) / float(c_numHistogramBuckets);

        sprintf_s(buffer, "%f", histogram[i] - (PDF::CDF(p1) - PDF::CDF(p0)));
        row.push_back(buffer);
    }
}

// ===================================== TESTING AND MAIN =====================================

template <typename PDF>
void TestRejectionSampling()
{
    std::vector<float> samplesWhite = RejectionSample<Seq::WhiteNoise2D, PDF>(1000000);
    std::vector<float> samplesR2 = RejectionSample<Seq::R2, PDF>(1000000);

    char fileName[256];

    {
        CSV csv;
        AddHistogramPDF<PDF>(csv);
        AddHistogram(csv, "White", samplesWhite, 100);
        AddHistogram(csv, "R2", samplesR2, 100);

        csv.resize(csv.size() + 1);
        csv.resize(csv.size() + 1);
        csv.back().push_back("Error:");
        AddHistogramError<PDF>(csv, "White", samplesWhite, 100);
        AddHistogramError<PDF>(csv, "R2", samplesR2, 100);

        sprintf_s(fileName, "out/histogram_RS_%s_%i.csv", PDF::FileLabel(), 100);
        SaveCSV(csv, fileName);
    }

    {
        CSV csv;
        AddHistogramPDF<PDF>(csv);
        AddHistogram(csv, "White", samplesWhite, 10000);
        AddHistogram(csv, "R2", samplesR2, 10000);

        csv.resize(csv.size() + 1);
        csv.resize(csv.size() + 1);
        csv.back().push_back("Error:");
        AddHistogramError<PDF>(csv, "White", samplesWhite, 10000);
        AddHistogramError<PDF>(csv, "R2", samplesR2, 10000);

        sprintf_s(fileName, "out/histogram_RS_%s_%i.csv", PDF::FileLabel(), 10000);
        SaveCSV(csv, fileName);
    }

    {
        CSV csv;
        AddHistogramPDF<PDF>(csv);
        AddHistogram(csv, "White", samplesWhite, 1000000);
        AddHistogram(csv, "R2", samplesR2, 1000000);

        csv.resize(csv.size() + 1);
        csv.resize(csv.size() + 1);
        csv.back().push_back("Error:");
        AddHistogramError<PDF>(csv, "White", samplesWhite, 1000000);
        AddHistogramError<PDF>(csv, "R2", samplesR2, 1000000);

        sprintf_s(fileName, "out/histogram_RS_%s_%i.csv", PDF::FileLabel(), 1000000);
        SaveCSV(csv, fileName);
    }
}

void TestRejectionSamplingFull()
{
    /*
        Test 1 - uniform pdf -> rejected sampled to -> y=2x pdf
         00 -> white noise as input, white noise rng 
         01 -> white noise as input, LDS rng
         10 -> LDS as input, white noise as rng
         11 -> LDS as input, LDS as rng.
        ! do this multiple times, show error and variance of histogram. at different sample counts too.

        Test 2 & 3 - different input and output pdfs
    */

    // Test 1 - uniform PDF to y=2x pdf
    {
        
    }



    // TODO: maybe have it like... make a white noise stream here. make a function to rejection sample it to some other PDF (using white vs LDS). report histograms. do other tests.

    //TestRejectionSampling<PDF::Y_Equals_X>();
    //TestRejectionSampling<PDF::Uniform>();
}

int main(int argc, char **argv)
{
    RejectionSamplingMain();

    TestRejectionSamplingFull();
    return 0;

    // TODO: continue this later

    {
        std::vector<float> samplesWhite = RejectionSample<Seq::WhiteNoise2D, PDF::Uniform>(10000);

        static const size_t c_numSamplesKept = 1000;

        std::vector<float> newSamplesWhite = ReservoirSample<Seq::WhiteNoise2D, PDF::Y_Equals_X>(samplesWhite, c_numSamplesKept);
        std::vector<float> newSamplesR2 = ReservoirSample<Seq::R2, PDF::Y_Equals_X>(samplesWhite, c_numSamplesKept);

        // TODO: should we start with samplesWhite or sampleR2? maybe try/show both? at least once?

        CSV csv;
        AddHistogramPDF<PDF::Y_Equals_X>(csv);
        AddHistogram(csv, "White", newSamplesWhite, c_numSamplesKept);
        AddHistogram(csv, "R2", newSamplesR2, c_numSamplesKept);

        csv.resize(csv.size() + 1);
        csv.resize(csv.size() + 1);
        csv.back().push_back("Error:");
        AddHistogramError<PDF::Y_Equals_X>(csv, "White", newSamplesWhite, c_numSamplesKept);
        AddHistogramError<PDF::Y_Equals_X>(csv, "R2", newSamplesR2, c_numSamplesKept);

        SaveCSV(csv, "out/blah.csv");
    }

    system("pause");

    return 0;
}


/*

================= REFOCUSING - WHAT ARE WE DOING? ==================



----- unweighted single reservoir sampling?

! experiment:
 * do this N times. look at a histogram of those N.  The item at index N is the value N, so can make a histogram of the samples to show that they are evenly spaced.
 * compare white noise vs LDS. LDS should be more flat.
 * do that comparison at various times during the N.
 * also vary how many source sample count? imperfect memory? or no?

----- unweighted reservoir sampling.  take N samples with even chance from a larger stream

* the output should have the same PDF as the input i think.
* use LDS to make sure this is true
* could do the tests multiple times and take average / std dev of it.
* different techniques like imperfect memory

! experiment:
 * each time you do this you get a histogram.
 * do this N times for N histograms, show error and variance
 * compare white noise vs LDS. LDS should be less error and less variance.
 * vary source sample count and imperfect memory.

----- weighted reservoir sampling. take N samples with weighted chance from a larger stream. also: convert a stream of a PDF to another (maybe not the best algorithm though)

? what is this used for. i could see it like "take 10 lights from this list of 100000 lights"

* vanilla tests...
 * vanilla can only generate a small number of samples from a source.
 * run the test N times and look at error and variance of error.
 * use LDS as a rng and do it again
 * use LDS as input stream and do it again
 * use LDS as both rng and input stream and do it again
 * try a couple different PDFs as input and output?

* imperfect memory to get higher acceptance rate
 * your constant weight sum idea?
 * marcos idea
 * a couple different PDFs as input and output
 * LDS involved

? how does this idea compare conceptually to rejection sampling or russian roulette?
? do we really want N items selected? or do we want to spit out a stream of output for a given input?
! this may not be the best way to launder PDFs

! experiment:
 * each time you do this you get a histogram.
 * do this N times for N histograms, show error and variance
 * compare white noise vs LDS. LDS should be less error and less variance.
 * vary source sample count and imperfect memory.
 * try with various input and output pdfs

====================================================================





! marcos technique. include it in the analysis


* Average and variance of error from histogram, for low sample counts unmodified algorithm. Compare vs LDS,
* Lds as input stream? Along with as rng. both ways

* unweighted reservoir sampling should have a result that is the same pdf

? can we do reservoir sampling to output a color (correlation / anti correlation) as part of it's logic?

* there's a thing about input pdf vs output pdf overlap (integrating the product of the 2)
 * this controls output rate, even before screwing w/ the weighting i think.

* should we use other LDSs? which ones?  Sobol? GR/Sqrt2?
? what other PDFs would be interesting to show?

* make a function to generate N samples from a PDF. use regular old rejection sampling
 * also compare vs LDS later!

* Show how you can take M samples from 1 PDF and make N samples of another PDF
 * do various pdf's as input and output (could chain them to show it works? i dunno)

* is PDF() ever used?

* compare vs LDS.
 * what is the win? Do you get a more accurate histogram vs the PDF?

? what happens when not enough survive?
 * is there a link to rejection sampling

? does this work with 2d and higher as well?
? Future: could we get noise color into this?
 * maybe something about closeness to existing neighbors?
 * maybe could just try re-sorting values for desired frequency characteristics

* Tune how many histogram buckets to use?


------------------------- MARCOS CODE

    // generate outputs which follow the PDF - Method 3
    weightSum = 0.0f;
    accept = 0;
    reject = 0;
    const int memLength = 16;
    std::vector<float> outputs3(c_numOutputs);
    std::vector<int> accepts3;
    for (float input : inputs)
    {
        float weight = PDF(input);
        weightSum += weight;
        for (size_t index = 0; index < c_numOutputs; ++index)
        {
            if (dist_float_0_1(rng) <= weight / weightSum)
            {
                outputs3[index] = input;
                accept++;
            }
            else
            {
                reject++;
            }
            if (accept + reject > memLength)
            {
                weightSum *= (float)(memLength) / (float)(1 + memLength);
            }
        }
        accepts3.push_back(accept);
    }

ah yep, interesting. I need to experiment but my thought right now is that reaching a constant weight is the only way you will get a consistent histogram if you take a subsection of a stream of output.
9:11
that constant weight controls quality vs output sample rate.  The "right value" depends on how much quality vs output samples you want.  But, the output rate is also affected by how much the input PDF and output PDF overlap (the integration of them multiplied together).  Smaller overlap = lower output before tuning that constant weight value.  larger overlap = higher output before tuning.
9:13
The weight sum starting at zero and summing up to a maximum weight (in my case) may go away in the direction i'm going to be looking too, because if "consistent quality" is the desire, you'd want the start of the stream to be no different than any other part of the stream.
9:13
unknowns to figure out for sure, and things to verify...
9:14
i was also thinking it might be interesting if you didn't know / had no control over the input PDF, if maybe there could be a "temporally filtered" approximation of the input PDF, and/or just the "constant weight sum" itself, to be able to dynamically adapt to whatever the input stream was doing
9:14
anyways... some interesting things here i think & yeah, i'll share results. (and i will test the method you say and let you know how it does & compares!)
)& with a constant weight sum, i want to think about how this relates to other things like perhaps rejection sampling or Russian roulette, which are very similar in formulation)


--------------------------------------------------------------------

Note:
* unsure if I'm using LDS in ReservoirSample() correctly. we use the second value only conditionally. could be a second 1d stream maybe? not real sure.

BLOG:
* weighted reservoir sampling with replacement.  https://blog.plover.com/prog/weighted-reservoir-sampling.html
 * https://en.wikipedia.org/wiki/Reservoir_sampling#Algorithm_A-Chao
* also just rejection sampling, as a lesser thing?
* twitter thread on transmuting a PDF: https://twitter.com/CasualEffects/status/1277748266709508102?s=20
* this video and the other which takes it to multiple samples: https://www.youtube.com/watch?v=A1iwzSew5QY

* Mention how fpe can randomly select N indices out of M total, without replacement.
 * You could repeat items to get rational weighting but mapping result to index would be hard

! for color of noise, just sort for energy function
 * is there other ways? like part of the weighting could be the distance from existing samples.


*/
