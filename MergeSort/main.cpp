#include <iostream>
#include <chrono>
#include <cstdlib>

using namespace std;
using namespace std::chrono;
void merge(int leftArr[], int rightArr[], int arr[], int leftLen, int rightLen)
{
    int i = 0, l = 0, r = 0;

    // Merge elements from leftArr and rightArr
    while (l < leftLen && r < rightLen)
    {
        if (leftArr[l] < rightArr[r])
        {
            arr[i++] = leftArr[l++];
        }
        else
        {
            arr[i++] = rightArr[r++];
        }
    }

    // Add remaining elements from leftArr
    while (l < leftLen)
    {
        arr[i++] = leftArr[l++];
    }

    // Add remaining elements from rightArr/
    while (r < rightLen)
    {
        arr[i++] = rightArr[r++];
    }
}

void mergeSort(int arr[], int len)
{
    if (len <= 1)
    {
        return;
    }

    int middle = len / 2;

    // Create left and right subarrays
    int *leftArr = new int[middle];
    int *rightArr = new int[len - middle];

    // Fill the left and right subarrays
    for (int i = 0; i < middle; i++)
    {
        leftArr[i] = arr[i];
    }
    for (int i = middle; i < len; i++)
    {
        rightArr[i - middle] = arr[i];
    }

    // Recursive calls
    mergeSort(leftArr, middle);
    mergeSort(rightArr, len - middle);

    // Merge sorted subarrays
    merge(leftArr, rightArr, arr, middle, len - middle);

    // Free memory
    delete[] leftArr;
    delete[] rightArr;
}

int main(int argc, char* argv[])
{

    if (argc != 2 ){
        cerr << "Usage: " << argv[0] << "<array_size" << endl;
        return 1;

    }
    int size = atoi(argv[1]);
    int* arr = new int[size];
  
    for(int i= 0; i < size; i ++) arr[i] = rand() % 100;
    cout << "Unsorted Array: ";
    for (int i = 0; i < size; i++)
    {
        cout << arr[i] << " ";
    }

    cout << endl;\
    auto start = high_resolution_clock::now();
    mergeSort(arr, size);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Sorted Array: ";
    for (int i = 0; i < size; i++)
    {
        cout << arr[i] << " ";
    }
    cout << endl;

    cout << "Sort Time: " << duration.count() << " ms" << endl;
    cerr << "Sort Time: " << duration.count() << " ms" << endl;

    return 0;
}
