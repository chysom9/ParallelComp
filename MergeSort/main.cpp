#include <iostream>
using namespace std;

void merge(int leftArr[], int rightArr[], int arr[], int leftLen, int rightLen) {
    int i = 0, l = 0, r = 0;

    // Merge elements from leftArr and rightArr
    while (l < leftLen && r < rightLen) {
        if (leftArr[l] < rightArr[r]) {
            arr[i++] = leftArr[l++];
        } else {
            arr[i++] = rightArr[r++];
        }
    }

    // Add remaining elements from leftArr
    while (l < leftLen) {
        arr[i++] = leftArr[l++];
    }

    // Add remaining elements from rightArr
    while (r < rightLen) {
        arr[i++] = rightArr[r++];
    }
}

void mergeSort(int arr[], int len) {
    if (len <= 1) {
        return;
    }

    int middle = len / 2;

    // Create left and right subarrays
    int *leftArr = new int[middle];
    int *rightArr = new int[len - middle];

    // Fill the left and right subarrays
    for (int i = 0; i < middle; i++) {
        leftArr[i] = arr[i];
    }
    for (int i = middle; i < len; i++) {
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

int main() {
    int arr[] = {1,30,33,51,66,6, 3, 2, 4, 5, 6,100, 7, 8, 9, 10};
    int len = sizeof(arr) / sizeof(arr[0]);

    mergeSort(arr, len);

    for (int i = 0; i < len; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;

    return 0;
}
