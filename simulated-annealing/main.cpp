#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;

double getDistance(double **arrayOfCities, int a, int b);

double getDistanceOfPath(vector<int> myPath, double **arrayOfDistances)
{
    int i = 0;
    double distance = 0;
    int size = myPath.size();

    for (i = 0; i < size - 1; i++)
    {
        distance += arrayOfDistances[myPath[i]][myPath[i + 1]];
    }
    distance += arrayOfDistances[0][myPath.back()];
    return distance;
}

void twoOptSwap(vector<int> *oldPath, vector<int> *newPath, int x, int y)
{
    int size = (*oldPath).size();
    int i;

    for (i = 0; i <= x - 1; i++)
    {
        newPath->at(i) = oldPath->at(i);
    }

    int back = 0;

    for (i = x; i <= y; i++)
    {
        newPath->at(i) = oldPath->at(y - back);
        back++;
    }

    for (i = y + 1; i < size; i++)
    {
        newPath->at(i) = oldPath->at(i);
    }
}

double getDistance(double **arrayOfCities, int a, int b)
{
    double x1 = arrayOfCities[a][0];
    double y1 = arrayOfCities[a][1];
    double x2 = arrayOfCities[b][0];
    double y2 = arrayOfCities[b][1];
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int preprocessingDistance(double *arrayOfDistances, int size, vector<bool> *visited)
{
    int min = -1;
    int i = 0;
    double minValue = -1;

    for (i = 0; i < size; i++)
    {
        if (arrayOfDistances[i] != 0)
        {
            if (!visited->at(i))
            {
                minValue = arrayOfDistances[i];
                min = i;
                break;
            }
        }
    }

    if (minValue == -1)
    {
        return -1;
    }

    for (; i < size; ++i)
    {
        if (arrayOfDistances[i] < minValue && arrayOfDistances[i] != 0)
        {
            if (!visited->at(i))
            {
                minValue = arrayOfDistances[i];
                min = i;
            }
        }
    }
    return min;
}

vector<int> createInitialSolution(double **arrayOfCities, double **arrayOfDistances, int n, int runTime)
{
    vector<int> initialPath;
    int following = 0;

    initialPath.push_back(0);

    vector<bool> visitedCities(n, false);

    while (true)
    {

        visitedCities.at(following) = true;
        following = preprocessingDistance(arrayOfDistances[following], n, &visitedCities);

        if (following == -1)
        {
            following = 0;
            initialPath.push_back(following);
            break;
        }
        else
        {
            initialPath.push_back(following);
        }
    }
    return initialPath;
}

void showPath(vector<int> *finalPath)
{
    for (int i = 0; i < (*finalPath).size(); i++)
    {
        cerr << (*finalPath)[i] + 1 << endl;
    }
}

void simulatedAnnealing(vector<int> *finalPath, double **arrayOfCities, double **arrayOfDistances, int n, int runTime, time_t programInitTime)
{
    double currentTemperature = 10;
    double redFactor = 0.007;

    /*distance after preprocessing*/
    double finalLength = getDistanceOfPath(*finalPath, arrayOfDistances);

    /*best route via SA*/
    vector<int> simulatedAnnealingPath = (*finalPath);
    double simulatedAnnealingLength = finalLength;

    /*temporary route*/
    vector<int> currentRoute;
    double currentLength;

    int firstIndex;
    int secondIndex;
    double upper;
    double probability;
    double r;
    int tempCity;
    long long int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator(seed);

    uniform_real_distribution<double> rProbability(0.0, 1.0);
    uniform_int_distribution<int> indexProbability(1, n - 1);

    while (currentTemperature > 0.00001)
    {
        if (runTime < difftime(time(NULL), programInitTime))
        {
            break;
        }

        firstIndex = indexProbability(generator);
        secondIndex = indexProbability(generator);

        if (firstIndex == secondIndex)
        {
            continue;
        }
        else if (firstIndex > secondIndex)
        {
            /*swap cities*/
            tempCity = firstIndex;
            firstIndex = secondIndex;
            secondIndex = tempCity;
        }

        currentRoute = simulatedAnnealingPath;

        simulatedAnnealingLength = getDistanceOfPath(*finalPath, arrayOfDistances);

        twoOptSwap(&simulatedAnnealingPath, &currentRoute, firstIndex, secondIndex);
        currentLength = getDistanceOfPath(currentRoute, arrayOfDistances);

        if (currentLength < simulatedAnnealingLength)
        {
            simulatedAnnealingPath = currentRoute;
            simulatedAnnealingLength = currentLength;

            if (simulatedAnnealingLength < finalLength)
            {
                *finalPath = simulatedAnnealingPath;
                finalLength = simulatedAnnealingLength;
            }
        }
        else
        {
            upper = (currentLength - simulatedAnnealingLength) / currentTemperature;
            probability = exp(-upper);
            r = rProbability(generator);

            if (r < probability)
            {
                simulatedAnnealingPath = currentRoute;
                simulatedAnnealingLength = currentLength;
            }
        }
        currentTemperature = currentTemperature * (1 - redFactor);
    }
}

int main()
{
    time_t programInitTime = time(NULL);

    int n;
    cin >> n;

    double **arrayOfCities = new double *[n];

    for (int i = 0; i < n; i++)
        arrayOfCities[i] = new double[2];

    int numberOfCity;
    double xCoord;
    double yCoord;

    for (int i = 0; i < n; i++)
    {

        cin >> numberOfCity;
        cin >> xCoord;
        cin >> yCoord;

        arrayOfCities[i][0] = xCoord - 1;
        arrayOfCities[i][1] = yCoord - 1;
    }

    int runTime;
    cin >> runTime;

    double **arrayOfDistances = new double *[n];

    for (int i = 0; i < n; i++)
        arrayOfDistances[i] = new double[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            if (i == j)
                arrayOfDistances[i][j] = 0;
            else
            {
                arrayOfDistances[i][j] = getDistance(arrayOfCities, i, j);
                arrayOfDistances[j][i] = arrayOfDistances[i][j];
            }
        }
    }

    vector<int> startingPath;

    startingPath = createInitialSolution(arrayOfCities, arrayOfDistances, n, runTime);

    double finalDistance;

    if (runTime > difftime(time(NULL), programInitTime))
    {
        simulatedAnnealing(&startingPath, arrayOfCities, arrayOfDistances, n, runTime, programInitTime);
    }

    finalDistance = getDistanceOfPath(startingPath, arrayOfDistances);
    cout << finalDistance << endl;
    showPath(&startingPath);

    return 0;
}