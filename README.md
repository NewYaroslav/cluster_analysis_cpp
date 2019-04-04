# cluster_analysis_cpp
Алгоритм на C++ для кластеризации

### Описание

Данная *header-only* C++ библиотека содержит класс для работы с k-mean++ (метод к-средних++).

### Как пользоваться

```C++
#include "cluster_analysis.hpp"

std::vector<std::vector<double>> test_data;
//...
ClusterAnalysis<double> iClusterAnalysis;

for(size_t i = 0; i < test_data.size(); ++i) {
	// добавляем пример 
	std::cout << i << " " << test_data[i][0] << " " << test_data[i][1] << std::endl;
	// 
	const int USER_LABEL = 0; // пользовательская метка, вдруг пригодится, но сейчас нет
	// не будем предварительно нормализовывать данные, не нада! (ставим false в конце)
	iClusterAnalysis.add_sample(test_data[i], USER_LABEL,  false);
}

const int MAX_CLUSTER = 3; // будем искать 3 класса объектов в test_data. Почему три? Патамушто!
const int USER_SEED = 123; // зерно для функции рандома внутри метода

iClusterAnalysis.update_k_mean_pp(MAX_CLUSTER, iClusterAnalysis.EUCLIDEAN_DISTANCE, USER_SEED);

// выводим кластеры наружу в большой мир
std::cout << std::endl;
std::cout << "clusters " << std::endl;
std::cout << std::endl;
for(int k = 0; k < MAX_CLUSTER; ++k) {
	std::cout << "cluster: " << k << std::endl;
	std::vector<std::vector<double>> data;
	std::vector<int> label;
	std::vector<int> indx;
	iClusterAnalysis.get_object(k, data, label, indx);
	for(size_t i = 0; i < data.size(); ++i) {
		std::cout << indx[i] << " " << data[i][0] << " " << data[i][1] << std::endl;
	}
	std::cout << std::endl;
}

```

### Полезные ссылки

* Кластеризация: [http://proginfo.ru/clustering/](http://proginfo.ru/clustering/)
* Метод к-средних++ [https://msdn.microsoft.com/ru-ru/magazine/mt185575.aspx](https://msdn.microsoft.com/ru-ru/magazine/mt185575.aspx)
