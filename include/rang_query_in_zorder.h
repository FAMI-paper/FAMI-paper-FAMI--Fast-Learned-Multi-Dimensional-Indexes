#include <libmorton/morton.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
using namespace std;
bool isRelevant2D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy, uint_fast64_t testad)
{
    uint_fast32_t testx, testy;
    libmorton::morton2D_64_decode(testad, testx, testy);
    return testx >= startx && testx <= endx && testy >= starty && testy <= endy;
}
uint_fast64_t nextJumpIn2D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy, uint_fast64_t fromad)
{
    uint_fast64_t toad = -1;
    uint_fast64_t startad = libmorton::morton2D_64_encode(startx, starty);
    uint_fast64_t endad = libmorton::morton2D_64_encode(endx, endy);
    for (toad = fromad + 1; toad <= endad; toad++)
    {
        if (isRelevant2D(startx, endx, starty, endy, toad))
        {
            break;
        }
    }
    return toad;
}
bool isRelevant2D(uint_fast64_t startad, uint_fast64_t endad, uint_fast64_t testad)
{
    uint_fast32_t startx, starty, endx, endy, testx, testy;
    libmorton::morton2D_64_decode(startad, startx, starty);
    libmorton::morton2D_64_decode(endad, endx, endy);
    libmorton::morton2D_64_decode(testad, testx, testy);
    return testx >= startx && testx <= endx && testy >= starty && testy <= endy;
}
uint_fast64_t nextJumpIn2D(uint_fast64_t startad, uint_fast64_t endad, uint_fast64_t fromad)
{
    uint_fast64_t toad = -1;
    for (toad = fromad + 1; toad <= endad; toad++)
    {
        if (isRelevant2D(startad, endad, toad))
        {
            break;
        }
    }
    return toad;
}
// uint_fast64_t nextJumpIn2D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy, uint_fast64_t fromad)
// {
//     uint_fast32_t fromx, fromy;
//     uint_fast64_t toad = -1;
//     fromad -= 1;
//     libmorton::morton2D_64_decode(fromad, fromx, fromy);
//     fromx += 1;
//     if (fromx > endx)
//     {
//         fromx = startx;
//         fromy += 1;
//     }
//     toad = libmorton::morton2D_64_encode(fromx, fromy);
//     return toad;
// }
vector<pair<uint_fast64_t, uint_fast64_t>> divideRange2D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy)
{
    vector<pair<uint_fast64_t, uint_fast64_t>> subrange;
    uint_fast64_t startad = libmorton::morton2D_64_encode(startx, starty);
    uint_fast64_t endad = libmorton::morton2D_64_encode(endx, endy);
    if (startad == endad)
    {
        subrange.push_back(make_pair(startad, endad));
        return subrange;
    }
    uint_fast64_t substart = startad;
    for (uint_fast64_t ad = startad + 1; ad <= endad;)
    {
        if (!isRelevant2D(startx, endx, starty, endy, ad))
        {
            subrange.push_back(make_pair(substart, ad - 1));
            substart = nextJumpIn2D(startx, endx, starty, endy, ad);
            // substart = nextJumpIn2D(startx, endx, starty, endy, ad);
            ad = substart + 1;
        }
        else
        {
            ad++;
        }
    }
    subrange.push_back(make_pair(substart, endad));
    return subrange;
}
vector<pair<uint_fast64_t, uint_fast64_t>> divideRegion2D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy, int regionExp)
{
    vector<pair<uint_fast64_t, uint_fast64_t>> subranges;
    startx = startx >> regionExp << regionExp;
    endx = endx >> regionExp << regionExp;
    starty = starty >> regionExp << regionExp;
    endy = endy >> regionExp << regionExp;
    int regionLength = 1 << regionExp;
    for (uint_fast32_t y = starty; y <= endy; y += regionLength)
    {
        for (uint_fast32_t x = startx; x <= endx; x += regionLength)
        {
            vector<pair<uint_fast64_t, uint_fast64_t>> subrange = divideRange2D(x, x + regionLength - 1, y, y + regionLength - 1);
            subranges.insert(subranges.end(), subrange.begin(), subrange.end());
        }
    }
    return subranges;
}
bool isRelevant3D(uint_fast64_t startad, uint_fast64_t endad, uint_fast64_t testad)
{
    uint_fast32_t startx, starty, startz, endx, endy, endz, testx, testy, testz;
    libmorton::morton3D_64_decode(startad, startx, starty, startz);
    libmorton::morton3D_64_decode(endad, endx, endy, endz);
    libmorton::morton3D_64_decode(testad, testx, testy, testz);
    return testx >= startx && testx <= endx && testy >= starty && testy <= endy && testz >= startz && testz <= endz;
}
uint_fast64_t nextJumpIn3D(uint_fast64_t startad, uint_fast64_t endad, uint_fast64_t fromad)
{
    uint_fast64_t toad = -1;
    for (toad = fromad + 1; toad <= endad; toad++)
    {
        if (isRelevant3D(startad, endad, toad))
        {
            break;
        }
    }
    return toad;
}
vector<pair<uint_fast64_t, uint_fast64_t>> divideRange3D(uint_fast32_t startx, uint_fast32_t endx, uint_fast32_t starty, uint_fast32_t endy, uint_fast32_t startz, uint_fast32_t endz)
{
    vector<pair<uint_fast64_t, uint_fast64_t>> subrange;
    uint_fast64_t startad = libmorton::morton3D_64_encode(startx, starty, startz);
    uint_fast64_t endad = libmorton::morton3D_64_encode(endx, endy, endz);
    if (startad == endad)
    {
        subrange.push_back(make_pair(startad, endad));
        return subrange;
    }
    uint_fast64_t substart = startad;
    for (uint_fast64_t ad = startad + 1; ad <= endad;)
    {
        if (!isRelevant3D(startad, endad, ad))
        {
            subrange.push_back(make_pair(substart, ad - 1));
            substart = nextJumpIn3D(startad, endad, ad);
            ad = substart + 1;
        }
        else
        {
            ad++;
        }
    }
    subrange.push_back(make_pair(substart, endad));
    return subrange;
}
vector<pair<uint_fast32_t, uint_fast32_t>> genRange(int size)
{
    vector<pair<uint_fast32_t, uint_fast32_t>> ranges(size);
    for (int i = 0; i < size; i++)
    {
        ranges[i].first = (uint_fast32_t)rand() % 50;
        ranges[i].second = ranges[i].first + (uint_fast32_t)rand() % 50;
    }
    return ranges;
}
int main(int argc, char *argv[])
{
    int size = 1e6;
    if (argc < 9)
    {
        cout << "error arg <2D/3D> <startx> <endx> <starty> <endy> <startz> <endz>" << endl;
        return 0;
    }
    uint_fast32_t startx, endx, starty, endy, startz, endz;
    int regionExp, regionSize;
    string dim = argv[1];
    startx = atoi(argv[2]);
    endx = atoi(argv[3]);
    starty = atoi(argv[4]);
    endy = atoi(argv[5]);
    startz = atoi(argv[6]);
    endz = atoi(argv[7]);
    regionExp = atoi(argv[8]);
    srand(921224);
    vector<pair<uint_fast32_t, uint_fast32_t>> xrange = genRange(size);
    vector<pair<uint_fast32_t, uint_fast32_t>> yrange = genRange(size);
    vector<pair<uint_fast64_t, uint_fast64_t>> subrange;
    auto begin = chrono::high_resolution_clock::now();
    if (dim == "2D")
    {
        for (int i = 0; i < size; i++)
        {
            // subrange = divideRange2D(xrange[i].first, xrange[i].second, yrange[i].first, yrange[i].second);
            subrange = divideRegion2D(xrange[i].first, xrange[i].second, yrange[i].first, yrange[i].second, regionExp);
        }
        // subrange = divideRange2D(startx, endx, starty, endy);
        // subrange = divideRegion2D(startx, endx, starty, endy, regionExp);
    }
    else
    {
        subrange = divideRange3D(startx, endx, starty, endy, startz, endz);
    }
    auto end = chrono::high_resolution_clock::now();
    uint64_t time_ns = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
    cout << "total_time[s]: " << (time_ns / 1000 / 1000) / 1000.0 << endl;
    uint_fast32_t x, y, z;
    // if (dim == "2D")
    // {
    //     for (int i = 0; i < subrange.size(); i++)
    //     {
    //         cout << subrange[i].first << " " << subrange[i].second << endl;
    //         libmorton::morton2D_64_decode(subrange[i].first, x, y);
    //         cout << x << " " << y << "-";
    //         libmorton::morton2D_64_decode(subrange[i].second, x, y);
    //         cout << x << " " << y << endl;
    //     }
    // }
    // else
    // {
    //     for (int i = 0; i < subrange.size(); i++)
    //     {
    //         cout << subrange[i].first << " " << subrange[i].second << endl;
    //         libmorton::morton3D_64_decode(subrange[i].first, x, y, z);
    //         cout << x << " " << y << " " << z << "-";
    //         libmorton::morton3D_64_decode(subrange[i].second, x, y, z);
    //         cout << x << " " << y << " " << z << endl;
    //     }
    // }
}