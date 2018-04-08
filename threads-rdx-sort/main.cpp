#include <time.h>
#include <omp.h>

#include <vector>
#include <thread>
#include <iostream>
using namespace std;

#define BYTE(data, shift) (((data) >> shift) & 255)

void TestSorting(int *data1, int *data2, int n) {
	bool flag = true;

	for (int i = 0; i < n - 1; i++)
	{
		if (data1[i] > data1[i + 1]) {
			flag = false; break;
		}
		if (data1[i] != data2[i]) {
			flag = false; break;
		}
	}

	if (data1[n] != data2[n]) flag = false;

	if (flag == true)
		cout << "correct" << endl;
	else cout << "error" << endl;
}

void omp_radixSort(int* data, int n, int threads) {
	int *temp = new int[n];
	int shift = 0;

	int *counter = new int[256];
	int *counters = new int[threads * 256];

	for (int j = 0; j < 4; j++, shift += 8) {
		memset(counter, 0, 256 * sizeof(int));
		memset(counters, 0, 256 * sizeof(int) * threads);
		int *threadCounter;

#pragma omp parallel firstprivate(threadCounter)
		{
			int tid = omp_get_thread_num();
			threadCounter = counters + 256 * tid;

#pragma omp for
			for (int i = 0; i < n; i++) {
				threadCounter[BYTE(data[i], shift)]++;
			}
		}

		for (int k = 0; k < threads; k++)
		{
			threadCounter = counters + 256 * k;
			for (int i = 0; i < 256; i++) {
				counter[i] += threadCounter[i];
			}
		}

		for (int i = 1; i < 256; i++) {
			counter[i] += counter[i - 1];
		}	

		for (int i = threads - 1; i >= 0; i--) {
			threadCounter = counters + 256 * i;
			for (int k = 0; k < 256; k++) {
				counter[k] -= threadCounter[k];
				threadCounter[k] = counter[k];
			}
		}

#pragma omp parallel firstprivate(threadCounter)
		{
			int tid = omp_get_thread_num();
			threadCounter = counters + 256 * tid;

#pragma omp for 
			for (int i = 0; i < n; i++) {
				temp[threadCounter[BYTE(data[i], shift)]++] = data[i];
			}
		}

		swap(data, temp);	
	}

	delete[] temp;
	delete[] counter;
	delete[] counters;
}

void CreateCounter(int *data, int *counter, int n, int shift) {
	for (int i = 0; i < n; i++)
		counter[BYTE(data[i], shift)]++;
}

void AddCounter(int *counter, int *threadCounter) {
	for (int i = 0; i < 256; i++) {
		counter[i] += threadCounter[i];
	}
}

void Arrangement(int *data, int *temp, int *threadCounter, int n, int shift) {
	for (int i = 0; i < n; i++) {
		temp[threadCounter[BYTE(data[i], shift)]++] = data[i];
	}
}

void threads_radixSort(int* data, int n, int threads) {
	int *temp = new int[n];
	int m = n / threads; 
	std::vector<std::thread> thr;
	int shift = 0;

	int *counter = new int[256];
	int *counters = new int[threads * 256];

	for (int j = 0; j < 4; j++, shift += 8) {
		memset(counter, 0, 256 * sizeof(int));
		memset(counters, 0, 256 * sizeof(int) * threads);
		int *threadCounter;

		for (int i = 0; i < threads - 1; i++) {
			threadCounter = counters + 256 * i;
			thr.push_back(std::thread(CreateCounter, data + i*m, threadCounter, m, shift));
		}

		threadCounter = counters + 256 * (threads - 1);
		CreateCounter(data + (threads - 1)*m, threadCounter, m, shift);
		for (auto& t : thr) t.join();
		thr.clear(); 

		for (int i = 0; i < threads; i++) {
			threadCounter = counters + 256 * i;
			AddCounter(counter, threadCounter);
		}

		for (int i = 1; i < 256; i++) {
			counter[i] += counter[i - 1];
		} 

		for (int i = threads - 1; i >= 0; i--) {
			threadCounter = counters + 256 * i; 
			for (int k = 0; k < 256; k++) {
				counter[k] -= threadCounter[k];
				threadCounter[k] = counter[k];
			}
		}

		for (int i = 0; i < threads - 1; i++) {
			threadCounter = counters + 256 * i;
			thr.push_back(std::thread(Arrangement, data + i*m, temp, threadCounter, m, shift));
		}

		threadCounter = counters + 256 * (threads - 1);
		Arrangement(data + (threads - 1)*m, temp, threadCounter, m, shift);
		for (auto& t : thr) t.join();
		thr.clear();

		swap(data, temp);	
	}
	
	delete[] temp;
	delete[] counters;
	delete[] counter;
}

int main(int argc, char* argv[])
{
	double t1, t2;
	int n = atoi(argv[1]);
	int threads = atoi(argv[2]);
	int *data = new int[n];
	int *data2 = new int[n];

	srand((unsigned)time(0));
	for (int i = 0; i < n; i++) {
		data[i] = rand() | (rand() << 16);
		data2[i] = data[i];
	}

	omp_set_num_threads(threads);
	t1 = omp_get_wtime();
	omp_radixSort(data, n, threads);
	t2 = omp_get_wtime();
	cout << "omp time = " << t2 - t1 << endl;

	t1 = omp_get_wtime();
	threads_radixSort(data2, n, threads);
	t2 = omp_get_wtime();
	cout << "threads time = " << t2 - t1 << endl;
	TestSorting(data, data2, n);

	delete[] data;
	delete[] data2;

	return 0;
}