#include <time.h>
#include <omp.h>

#include <cstring>
#include <vector>
#include <thread>
#include <iostream>
using namespace std;

// проверка сортировки двух одинаковых массивов
bool TestSorting(int *data1, int *data2, int n) {
	bool flag = true;

	for (int i = 0; i < n - 1; i++) {
		if (data1[i] > data1[i + 1]) {
			flag = false; break;
		}
		if (data1[i] != data2[i]) {
			flag = false; break;
		}
	}

	if (data1[n - 1] != data2[n - 1]) flag = false;

	return flag;
}

// подсчет значений каждого байта в массиве
void CreateCounters(int *data, int *counters, int n) {
	memset(counters, 0, 256 * sizeof(int) * sizeof(int));

	unsigned char *byteP = (unsigned char*)data;
	unsigned char *dataEnd = (unsigned char*)(data + n);

	int i;
	while (byteP != dataEnd) {
		for (i = 0; i < sizeof(int); i++)
			counters[256 * i + *byteP++]++;
	}
}

// сортировка data по указанному байту
void ByteSort(int *data, int *temp, int *counter, int byte, int n) {
	unsigned char *byteP;
	int sum, c;

	sum = 0;
	for (int i = 0; i < 256; i++) {
		c = counter[i];
		counter[i] = sum;
		sum += c;
	}

	// при втором проходе по исходному массиву выполняется копирование элемента
	// во вспомогательный массив по соответствующему индексу в массиве подсчётов
	byteP = (unsigned char *)data + byte;
	for (int i = 0; i < n; i++, byteP += sizeof(int))
		temp[counter[*byteP]++] = data[i];
}

// поразрядная сортировка массива data размера n.
void radixSort(int* data, int n) {
	int *temp = new int[n];
	int *counters = new int[sizeof(int) * 256];
	int *counter;

	CreateCounters(data, counters, n);

	for (int i = 0; i < sizeof(int); i++) {
		counter = counters + 256 * i;
		ByteSort(data, temp, counter, i, n);
		swap(data, temp);
	}

	delete[] temp;
	delete[] counters;
}

// поразрядная сортировка с использованием OpenMP
void radixSort_omp(int* data, int n, int nThreads) {
	int *temp = new int[n];
	unsigned char *byteP = (unsigned char*)data;
	int byte = 0;

	int *counter = new int[256];
	int *counters = new int[nThreads * 256];

	for (int j = 0; j < sizeof(int); j++, byte++) {
		memset(counter, 0, 256 * sizeof(int));
		memset(counters, 0, 256 * sizeof(int) * nThreads);
		int *threadCounter;

#pragma omp parallel firstprivate(threadCounter, byteP)
		{
			int tid = omp_get_thread_num();
			threadCounter = counters + 256 * tid;

// каждый поток получает на обработку часть массива
// и выполняет подсчет элементов в свой массив подсчетов threadCounter.
#pragma omp for
			for (int i = 0; i < n; i++) {
				byteP = (unsigned char*)data + i*sizeof(int) + byte;
				threadCounter[*byteP]++;
			}
		}

		// c помощью массивов подсчетов со всех потоков выполняется вычисление смещений,
		// по которым будут располагаться элементы при втором проходе по массиву.
		for (int k = 0; k < nThreads; k++) {
			threadCounter = counters + 256 * k;
			for (int i = 0; i < 256; i++) {
				counter[i] += threadCounter[i];
			}
		}

		for (int i = 1; i < 256; i++) {
			counter[i] += counter[i - 1];
		}

		for (int i = nThreads - 1; i >= 0; i--) {
			threadCounter = counters + 256 * i;
			for (int k = 0; k < 256; k++) {
				counter[k] -= threadCounter[k];
				threadCounter[k] = counter[k];
			}
		}

#pragma omp parallel firstprivate(threadCounter, byteP)
		{
			int tid = omp_get_thread_num();
			threadCounter = counters + 256 * tid;

// каждый поток получает на обработку ту же часть массива, что и ранее,
// и выполняет копирование элемента во вспомогательный массив
// по соответствующему индексу в массиве смещений.
#pragma omp for
			for (int i = 0; i < n; i++) {
				byteP = (unsigned char*)data + i * sizeof(int) + byte;
				temp[threadCounter[*byteP]++] = data[i];
			}
		}

		swap(data, temp);
	}

	delete[] temp;
	delete[] counter;
	delete[] counters;
}

void CreateCounter(int *data, int *counter, int n, int byte) {
	unsigned char *byteP = (unsigned char*)data + byte;
	for (int i = 0; i < n; i++, byteP += sizeof(int))
		counter[*byteP]++;
}

void AddCounter(int *counter, int *threadCounter) {
	for (int i = 0; i < 256; i++) {
		counter[i] += threadCounter[i];
	}
}

void Arrangement(int *data, int *temp, int *threadCounter, int n, int byte) {
	unsigned char *byteP = (unsigned char*)data + byte;
	for (int i = 0; i < n; i++, byteP += sizeof(int)) {
		temp[threadCounter[*byteP]++] = data[i];
	}
}

// реализация параллельной поразрядной сортировки
// с использованием потоков C++ 11
// (повторяет реализацию с использованием OpenMP)
void radixSort_threads(int* data, int n, int nThreads) {
	int *temp = new int[n];
	int m = n / nThreads;
	std::vector<std::thread> thr;
	int byte = 0;

	int *counter = new int[256];
	int *counters = new int[nThreads * 256];

	for (int j = 0; j < sizeof(int); j++, byte++) {
		memset(counter, 0, 256 * sizeof(int));
		memset(counters, 0, 256 * sizeof(int) * nThreads);
		int *threadCounter;

		for (int i = 0; i < nThreads - 1; i++) {
			threadCounter = counters + 256 * i;
			thr.push_back(std::thread(CreateCounter, data + i*m, threadCounter, m, byte));
		}

		threadCounter = counters + 256 * (nThreads - 1);
		CreateCounter(data + (nThreads - 1)*m, threadCounter, m, byte);
		for (auto& t : thr) t.join();
		thr.clear();

		for (int i = 0; i < nThreads; i++) {
			threadCounter = counters + 256 * i;
			AddCounter(counter, threadCounter);
		}

		for (int i = 1; i < 256; i++) {
			counter[i] += counter[i - 1];
		}

		for (int i = nThreads - 1; i >= 0; i--) {
			threadCounter = counters + 256 * i;
			for (int k = 0; k < 256; k++) {
				counter[k] -= threadCounter[k];
				threadCounter[k] = counter[k];
			}
		}

		for (int i = 0; i < nThreads - 1; i++) {
			threadCounter = counters + 256 * i;
			thr.push_back(std::thread(Arrangement, data + i*m, temp, threadCounter, m, byte));
		}

		threadCounter = counters + 256 * (nThreads - 1);
		Arrangement(data + (nThreads - 1)*m, temp, threadCounter, m, byte);
		for (auto& t : thr) t.join();
		thr.clear();

		swap(data, temp);
	}

	delete[] temp;
	delete[] counters;
	delete[] counter;
}


int main(int argc, char* argv[]) {
	int n = atoi(argv[1]);
	int nThreads = atoi(argv[2]);
	int *data1 = new int[n];
	int *data2 = new int[n];

	double t1, t2;

	srand((unsigned)time(0));
	for (int i = 0; i < n; i++) {
		data1[i] = rand() | (rand() << 16);
		data2[i] = data1[i];
	}

	omp_set_num_threads(nThreads);
	t1 = omp_get_wtime();
	radixSort_omp(data1, n, nThreads);
	t2 = omp_get_wtime();
	cout << "wall clock time (omp) = " << t2 - t1 << endl;

	t1 = omp_get_wtime();
	radixSort_threads(data2, n, nThreads);
	t2 = omp_get_wtime();
	cout << "wall clock time (threads) = " << t2 - t1 << endl;

	// проверка
	if (TestSorting(data1, data2, n)) printf("\ncorrect: data size = %d, # of threads = %d\n", n, omp_get_max_threads());
    else printf("\nerror: data size = %d, # of threads = %d\n", n, omp_get_max_threads());

	delete[] data1;
	delete[] data2;

	return 0;
}
