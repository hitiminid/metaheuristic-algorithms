#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

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

void twoOpt(vector<int> *myPath, double **arrayOfDistances, int numberOfRepeats, double runTime)
{
    int size = (*myPath).size();
    double newDistance;

    vector<int> newPath = (*myPath);
    int improve = 0;
    time_t start;
    time_t end;
    start = time(NULL);
    double elapsed;

    while (improve < numberOfRepeats)
    {

        double bestResult = getDistanceOfPath((*myPath), arrayOfDistances);

        for (int i = 1; i < size - 2; i++)
        {
            for (int k = i + 1; k < size - 1; k++)
            {
                end = time(NULL);
                elapsed = difftime(end, start);
                if (elapsed >= runTime)
                {
                    return;
                }
                twoOptSwap(myPath, &newPath, i, k);
                newDistance = getDistanceOfPath(newPath, arrayOfDistances);

                if (newDistance < bestResult)
                {
                    improve = 0;
                    (*myPath) = newPath;
                    bestResult = newDistance;
                }
            }
        }
        improve++;
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

    for (; i < size; i++)
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

vector<int> createInitialSolution(double **arrayOfCities, double **arrayOfDistances, int n)
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

double getNeighbourhoodRange(double initialPathLength, int numberOfCities)
{
    double neighbourhoodRange = initialPathLength / numberOfCities;
    return neighbourhoodRange;
}

double tabuSearch(vector<int> *initialPath, double **arrayOfCities, double **arrayOfDistances, double neighbourhoodRange, int n, double currentBestPath, double runTime)
{
    int x;
    int y;
    vector<pair<int, int>> tabuList;
    vector<int> tempPath = (*initialPath);

    time_t start;
    time_t end;
    start = time(NULL);

    double elapsed;
    double distanceBetweenPoints;
    double pathLength;

    int iterations = 0;
    bool timeElapsed = false;
    int cadention = 0;

    while (iterations < 10000 && timeElapsed == false)
    {

        cadention++;
        iterations++;

        if (cadention % 25 == 0)
        {
            tabuList.clear();
        }

        x = rand() % (n - 1);

        for (y = 1; y < n - 1; y++)
        {
            distanceBetweenPoints = getDistance(arrayOfCities, x, y);
            end = time(NULL);
            elapsed = difftime(end, start);
            if (elapsed >= runTime)
            {
                timeElapsed = true;
                break;
            }

            if (x == y || x == 0)
                break;

            pair<int, int> item(x, y);

            if (!(find(tabuList.begin(), tabuList.end(), item) != tabuList.end()))
            {

                if (distanceBetweenPoints > neighbourhoodRange * 50)
                {
                    tabuList.push_back(make_pair(x, y));
                }
                else
                {

                    twoOptSwap(initialPath, &tempPath, x, y);
                    pathLength = getDistanceOfPath(tempPath, arrayOfDistances);

                    if (pathLength >= currentBestPath)
                    {
                        tabuList.push_back(make_pair(x, y));
                    }
                    else
                    {
                        currentBestPath = pathLength;
                        cout << pathLength << endl;
                        tabuList.push_back(make_pair(x, y));
                        (*initialPath) = tempPath;
                    }
                }
            }
            else
            {
                continue;
            }
        }
    }
}

int main()
{
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
    startingPath = createInitialSolution(arrayOfCities, arrayOfDistances, n);

    double abc = getDistanceOfPath(startingPath, arrayOfDistances);

    double neighbourhoodRange = getNeighbourhoodRange(abc, n);

    int numberOfRepeats = 30;
    int runTime = 300;
    twoOpt(&startingPath, arrayOfDistances, numberOfRepeats, runTime);

    abc = getDistanceOfPath(startingPath, arrayOfDistances);

    vector<int> route;

    double finalDistance;
    finalDistance = getDistanceOfPath(startingPath, arrayOfDistances);

    tabuSearch(&startingPath, arrayOfCities, arrayOfDistances, neighbourhoodRange, n, abc, runTime);

    cout << finalDistance << endl;

    for (int i = 0; i < startingPath.size(); i++)
    {
        cerr << startingPath[i] + 1 << " " << endl;
    }
    return 0;
}