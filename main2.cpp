#include <vector>
#include <random>
#include <string>

#define DETERMINISTIC() false 

static const size_t c_numHistogramBuckets = 5;

static const size_t c_numInputs = 10000;
static const size_t c_numOutputs = 5;

static int FloatToItem(float f, int numItems)
{
    // if f is in [0,1], remaps to [0, numItems-1]
    return std::min(int(f * float(numItems)), numItems - 1);
}

int main(int argc, char **argv)
{
    // PDF is for x in [0,1]
    auto PDF = [](float x)
    {
        return 2.0f * x;
    };

    auto CDF = [](float x)
    {
        return x * x;
    };

#if DETERMINISTIC()
    std::mt19937 rng;
#else
    std::random_device rd;
    std::mt19937 rng(rd());
#endif

    std::uniform_real_distribution<float> dist_float_0_1(0.0f, 1.0f);
    std::uniform_int_distribution<int> dist_int_0_nminus1(0, c_numOutputs - 1);

    // generate the uniform input sequence
    std::vector<float> inputs(c_numInputs);
    for (float& f : inputs)
        f = dist_float_0_1(rng);

    // Generate outputs which follow the PDF - Method 1
    // Algorithm A from M. T. Chao
    // https://en.wikipedia.org/wiki/Reservoir_sampling#Algorithm_A-Chao
    float weightSum = 0.0f;
    int accept = 0;
    int reject = 0;
    std::vector<float> outputs1(c_numOutputs);
    std::vector<int> accepts1;
    for (float input : inputs)
    {
        float weight = PDF(input);
        weightSum += weight;
        if (dist_float_0_1(rng) <= weight / weightSum)
        {
            outputs1[dist_int_0_nminus1(rng)] = input;
            accept++;
        }
        else
        {
            reject++;
        }
        accepts1.push_back(accept);
    }
    printf("Method 1 accepted %i, rejected %i. Accepted %0.2f%%.\n", accept, reject, 100.0f * float(accept) / float(accept+reject));

    // generate outputs which follow the PDF - Method 2
    weightSum = 0.0f;
    accept = 0;
    reject = 0;
    std::vector<float> outputs2(c_numOutputs);
    std::vector<int> accepts2;
    for (float input : inputs)
    {
        float weight = PDF(input);
        weightSum += weight;
        for (size_t index = 0; index < c_numOutputs; ++index)
        {
            if (dist_float_0_1(rng) <= weight / weightSum)
            {
                outputs2[index] = input;
                accept++;
            }
            else
            {
                reject++;
            }
        }
        accepts2.push_back(accept);
    }
    printf("Method 2 accepted %i, rejected %i. Accepted %0.2f%%.\n", accept, reject, 100.0f * float(accept) / float(accept + reject));

    // write a csv of resulting histograms so we can graph results
    typedef std::vector<std::vector<std::string>> CSV;
    CSV csv(c_numInputs + 1);
    {
        char buffer[256];

        // write the actual value, by using the CDF
        {
            csv[0].push_back("PDF");
            for (int i = 0; i < c_numHistogramBuckets; ++i)
            {
                float p0 = float(i) / float(c_numHistogramBuckets);
                float p1 = float(i + 1) / float(c_numHistogramBuckets);
                sprintf_s(buffer, "%f", CDF(p1) - CDF(p0));
                csv[i + 1].push_back(buffer);
            }
        }

        // write the histogram of method 1
        {
            std::vector<float> histogram(c_numHistogramBuckets, 0.0f);
            for (float f : outputs1)
                histogram[FloatToItem(f, c_numHistogramBuckets)] += 1.0f / float(c_numOutputs);
            csv[0].push_back("Method 1");
            for (int i = 0; i < c_numHistogramBuckets; ++i)
            {
                sprintf_s(buffer, "%f", histogram[i]);
                csv[i + 1].push_back(buffer);
            }
        }

        // write the histogram of method 2
        {
            std::vector<float> histogram(c_numHistogramBuckets, 0.0f);
            for (float f : outputs2)
                histogram[FloatToItem(f, c_numHistogramBuckets)] += 1.0f / float(c_numOutputs);
            csv[0].push_back("Method 2");
            for (int i = 0; i < c_numHistogramBuckets; ++i)
            {
                sprintf_s(buffer, "%f", histogram[i]);
                csv[i + 1].push_back(buffer);
            }
        }

        // write the accept counts over time
        {
            csv[0].push_back("");
            csv[0].push_back("Accepts Method 1");
            csv[0].push_back("Accepts Method 2");
            for (size_t i = 0; i < c_numInputs; ++i)
            {
                while (csv[i + 1].size() < 4)
                    csv[i + 1].push_back("");

                sprintf_s(buffer, "%i", accepts1[i]);
                csv[i + 1].push_back(buffer);
                sprintf_s(buffer, "%i", accepts2[i]);
                csv[i + 1].push_back(buffer);
            }
        }

        // write the csv
        FILE* file = nullptr;
        fopen_s(&file, "out/out.csv", "w+t");
        for (const auto& row : csv)
        {
            bool first = true;
            for (const std::string& col : row)
            {
                fprintf(file, "%s\"%s\"", first ? "" : ",", col.c_str());
                first = false;
            }
            fprintf(file, "\n");
        }
        fclose(file);
    }

    system("pause");

    return 0;
}
