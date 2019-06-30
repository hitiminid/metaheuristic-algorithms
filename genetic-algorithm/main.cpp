#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <time.h>
#include <chrono>
#include <random>
#include <algorithm>

using namespace std;

long long int seed = chrono::system_clock::now().time_since_epoch().count();
mt19937 generator(seed);

double getDistance(double **arrayOfCities, int a, int b);

bool pairComparator(pair<double, vector<int>> firstPair, pair<double, vector<int>> secondPair)
{
    return firstPair.first < secondPair.first;
}

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

vector<int> createRandomPath(int numberOfCities)
{
    int i = 0;
    vector<int> newPath;
    newPath.push_back(0);

    for (i = 1; i < numberOfCities; i++)
        newPath.push_back(i);

    newPath.push_back(0);
    shuffle(newPath.begin() + 1, newPath.end() - 1, generator);

    return newPath;
}

vector<pair<double, vector<int>>> createFirstGeneration(int numberOfCities, double **arrayOfDistances)
{
    int populationSize = (int)floor((numberOfCities * (log(numberOfCities) - 1)));
    if (populationSize > 500)
    {
        populationSize = 500;
    }

    vector<pair<double, vector<int>>> population;
    double lenghtOfNewPath;
    vector<int> newPath;

    for (int i = 0; i < populationSize; i++)
    {
        newPath = createRandomPath(numberOfCities);
        lenghtOfNewPath = getDistanceOfPath(newPath, arrayOfDistances);
        pair<double, vector<int>> newPair = make_pair(lenghtOfNewPath, newPath);
        population.insert(upper_bound(population.begin(), population.end(), newPair, pairComparator), newPair);
    }

    return population;
}

void showData(pair<double, vector<int>> *individual)
{
    double length = (*individual).first;
    vector<int> route = (*individual).second;
    cout << length << endl;
    for (int i = 0; i < route.size(); i++)
    {
        cerr << route[i] + 1 << endl;
    }
}

void crossOver(pair<double, vector<int>> *firstParent, pair<double, vector<int>> *secondParent, int numberOfCities, double **arrayOfDistances)
{
    vector<int> firstPath = (*firstParent).second;
    vector<int> secondPath = (*secondParent).second;
    int startingPosition;
    int finishPosition;
    vector<int> child = firstPath;
    uniform_real_distribution<double> parentProbability(0.0, 1.0);
    uniform_int_distribution<int> indexProbability(1, numberOfCities - 3);

    int x = indexProbability(generator);
    int y = indexProbability(generator);

    if (x > y)
    {
        startingPosition = y;
        finishPosition = x;
    }
    else if (x < y)
    {
        startingPosition = x;
        finishPosition = y;
    }
    else
    {
        startingPosition = x;
        finishPosition = y + 1;
    }

    vector<bool> used(numberOfCities, false);
    vector<int> temporary;

    for (int i = startingPosition; i <= finishPosition; i++)
    {
        used.at(firstPath.at(i)) = true;
    }

    used.at(0) = true;

    for (int i = 1; i < secondPath.size() - 1; i++)
    {
        if (!used.at(secondPath.at(i)))
        {
            temporary.push_back(secondPath.at(i));
        }
    }

    for (int i = 1; i < startingPosition; i++)
    {
        child.at(i) = temporary.at(0);
        temporary.erase(temporary.begin());
    }

    for (int i = finishPosition + 1; i < numberOfCities; i++)
    {
        child.at(i) = temporary.at(0);
        temporary.erase(temporary.begin());
    }

    double chooseParent = parentProbability(generator);
    double childLength = getDistanceOfPath(child, arrayOfDistances);

    if (chooseParent <= 0.7)
    {
        (*secondParent).second = child;
        (*secondParent).first = childLength;
    }
    else
    {
        (*firstParent).second = child;
        (*firstParent).first = childLength;
    }
}

void mutate(pair<double, vector<int>> *individual, int numberOfCities, double **arrayOfDistances)
{
    vector<int> route = (*individual).second;
    vector<int> mutatedRoute = (*individual).second;

    int firstCity;
    int secondCity;
    uniform_int_distribution<int> chooseNewIndex(2, numberOfCities - 2);

    firstCity = chooseNewIndex(generator);
    secondCity = chooseNewIndex(generator);

    if (firstCity == secondCity)
    {
        firstCity--;
    }
    else if (firstCity > secondCity)
    {
        int temp;
        temp = firstCity;
        firstCity = secondCity;
        secondCity = temp;
    }

    twoOptSwap(&route, &mutatedRoute, firstCity, secondCity);

    (*individual).second = mutatedRoute;
    (*individual).first = getDistanceOfPath(mutatedRoute, arrayOfDistances);
}

void evolvePopulation(vector<pair<double, vector<int>>> *population, double **arrayOfDistances, int numberOfCities)
{
    double mutationProbability = 0.7;
    uniform_real_distribution<double> mutationProbabilityGenerator(0.0, 1.0);
    double mutationChance;
    int populationSize = (*population).size();
    int numberOfCrosses = populationSize;

    if (populationSize % 2 != 0)
    {
        numberOfCrosses -= 2;
    }

    for (int i = 0; i < numberOfCrosses; i += 2)
    {
        crossOver(&population->at(i), &population->at(i + 1), numberOfCities, arrayOfDistances);
    }

    uniform_int_distribution<int> numberOfExtraCrossOvers(5, 15);
    uniform_int_distribution<int> chooseIndividual(0, populationSize - 1);

    int numberOfIntercourses = numberOfExtraCrossOvers(generator);
    int firstIndividual;
    int secondIndividual;

    for (int i = 0; i < numberOfIntercourses; i++)
    {

        firstIndividual = chooseIndividual(generator);
        secondIndividual = chooseIndividual(generator);

        if (firstIndividual == secondIndividual)
        {
            continue;
        }

        crossOver(&population->at(firstIndividual), &population->at(secondIndividual), numberOfCities, arrayOfDistances);
    }

    for (pair<double, vector<int>> individual : (*population))
    {

        mutationChance = mutationProbabilityGenerator(generator);

        if (mutationChance >= mutationProbability)
        {
            double length = individual.first;
            vector<int> route = individual.second;
            mutate(&individual, numberOfCities, arrayOfDistances);
        }
    }
    sort((*population).begin(), (*population).end(), pairComparator);
}

void geneticAlgorithm(vector<pair<double, vector<int>>> *population, pair<double, vector<int>> *preprocessedSolution, double **arrayOfDistances, int numberOfGenerations, int numberOfCities, int runTime, time_t programInitTime)
{
    int i = 0;
    pair<double, vector<int>> final = *preprocessedSolution;
    pair<double, vector<int>> temporary;
    pair<double, vector<int>> bestGAIndividual = final;

    while (i < numberOfGenerations)
    {

        if (runTime < difftime(time(NULL), programInitTime))
        {
            break;
        }

        evolvePopulation(population, arrayOfDistances, numberOfCities);

        temporary.first = (*population)[0].first;
        temporary.second = (*population)[0].second;

        if (temporary.first < bestGAIndividual.first)
        {

            bestGAIndividual.first = temporary.first;
            bestGAIndividual.second = temporary.second;

            if (bestGAIndividual.first < final.first)
            {
                cout << final.first << endl;
                final.first = bestGAIndividual.first;
                final.second = bestGAIndividual.second;
                (*preprocessedSolution).first = bestGAIndividual.first;
                (*preprocessedSolution).second = bestGAIndividual.second;
            }
        }
        i++;
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

    pair<double, vector<int>> result;
    result.second = createInitialSolution(arrayOfCities, arrayOfDistances, n, runTime);
    result.first = getDistanceOfPath(result.second, arrayOfDistances);

    vector<pair<double, vector<int>>> population;
    population = createFirstGeneration(n, arrayOfDistances);

    if (runTime > difftime(time(NULL), programInitTime))
    {
        geneticAlgorithm(&population, &result, arrayOfDistances, 10, n, runTime, programInitTime);
    }

    showData(&result);
    return 0;
}