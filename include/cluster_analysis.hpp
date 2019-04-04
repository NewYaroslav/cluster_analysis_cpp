#ifndef CLUSTER_ANALYSIS_HPP_INCLUDED
#define CLUSTER_ANALYSIS_HPP_INCLUDED

#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#define RAND_NTAB 32
#define RAND_NWUP 8
#define RAND_EPS 1.2e-7
#define RAND_RNMX (1.0 - RAND_EPS)
#define RAND_IM1 2147483563
#define RAND_IM2 2147483399
#define RAND_AM (1./RAND_IM1)
#define RAND_IMM1 (RAND_IM1-1)
#define RAND_IA1 40014
#define RAND_IA2 40692
#define RAND_IQ1 53668
#define RAND_IQ2 52774
#define RAND_IR1 12211
#define RAND_IR2 3791
#define RAND_NDIV (1 + RAND_IMM1 / RAND_NTAB)

/** \brief Класс для кластерного анализа
 */
template <typename T>
class ClusterAnalysis {
private:
    // для функции rand
    long rand_dummy; /**< Зерно */
    long dummy2 = 123456789;
    long iy = 0;
    long iv[RAND_NTAB];

    inline void rand_seed(long dum) {
        rand_dummy = dum;
    }

    float custom_rand() {
        int j;
        long k;
        //static long dummy2 = 123456789;
        //static long iy = 0;
        //static long iv[RAND_NTAB];
        float temp;
        if(rand_dummy <= 0 || !iy) {
            if(rand_dummy < 0) {
                rand_dummy = -rand_dummy;
            } else if(rand_dummy == 0) {
                rand_dummy = 1;
            }
            dummy2 = rand_dummy;
            for(j = RAND_NTAB + RAND_NWUP - 1; j >= 0; j--) {
                k = rand_dummy / RAND_IQ1;
                if((rand_dummy = RAND_IA1 * (rand_dummy - k * RAND_IQ1) - RAND_IR1*k) < 0) {
                    rand_dummy += RAND_IM1;
                }
                if(j < RAND_NTAB) {
                    iv[j] = rand_dummy;
                }
            }
            iy = iv[0];
        }
        k = rand_dummy / RAND_IQ1;
        if((rand_dummy = RAND_IA1 * (rand_dummy - k * RAND_IQ1) - RAND_IR1 * k) < 0) {
            rand_dummy += RAND_IM1;
        }
        k = dummy2 / RAND_IQ2;
        if((dummy2 = RAND_IA2 * (dummy2 - k * RAND_IQ2) - RAND_IR2 * k) < 0) {
            dummy2 += RAND_IM2;
        }
        iy = iv[j = iy / RAND_NDIV] - dummy2;
        iv[j] = rand_dummy;
        if(iy < 1) {
            iy += RAND_IMM1;
        }
        if((temp = RAND_AM * iy) > RAND_RNMX) {
            return RAND_RNMX;
        }
        return temp;
    }

    int random_number(int m, int n) {
        return (custom_rand() * (n - m + 1)) + m;
    }
public:
    enum ErrorsType {
        OK = 0,
        INVALID_PARAMETER = 1,
    };

    enum SimilarityMeasureType {
        EUCLIDEAN_DISTANCE = 0,
        PEARSON_CORRELATION,
        SPEARMAN_RANK_CORRELATION,
    };

    enum NormalizationType {
        MINMAX_0_1 = 0,
        MINMAX_1_1 = 1,
    };
private:
    /** \brief MinMax нормализация данных
     * \param in входные данные для нормализации
     * \param out нормализованный вектор
     * \param type тип нормализации (0 - нормализация данных к промежутку от 0 до 1, иначе от -1 до 1)
     * \return вернет 0 в случае успеха
     */
    int calc_min_max(std::vector<T> &in, std::vector<T> &out, int type)
    {
        if(in.size() == 0)
                return INVALID_PARAMETER;
        T max_data = *std::max_element(in.begin(), in.end());
        T min_data = *std::min_element(in.begin(), in.end());
        T ampl = max_data - min_data;
        out.resize(in.size());
        if(ampl != 0) {
            for(size_t i = 0; i < in.size(); i++) {
                out[i] = type == 0 ? (double)(in[i] - min_data) / ampl : 2.0 * ((double)(in[i] - min_data) / ampl) - 1.0;
            }
        } else {
            for(size_t i = 0; i < in.size(); i++) {
                out[i] = 0.0;
            }
        }
        return OK;
    }

    /**@brief Quick sort function (Быстрая сортировка Хоара)
     * Данная функция нужна для коэффициента корреляции Спирмена
     */
    template<class T1, class T2>
    void quick_sort(std::vector<T1> &a, std::vector<T2> &ra, long l, long r)
    {
        long i = l, j = r;
        T1 temp, p;
        p = a[ l + (r - l)/2 ];
        do {
            while(a[i] < p) i++;
            while(a[j] > p) j--;
            if (i <= j) {
                temp = a[i];
                a[i] = a[j];
                a[j] = temp;
                ra[i] = j;
                ra[j] = i;
                i++; j--;
            }
        } while(i <= j);
        if(i < r)
            quick_sort(a, ra, i, r);
        if(l < j)
            quick_sort(a, ra, l , j);
    };

    /**@brief Quick sort function (Быстрая сортировка Хоара)
     * Данная функция нужна для коэффициента корреляции Спирмена
     */
    template<class T1, class T2>
    void quick_sort(std::vector<T1> &a, std::vector<T2> &ra)
    {
            std::vector<T1> _a = a;
            quick_sort(_a, ra, 0, _a.size() - 1);
    }

    /** \brief Коэффициент корреляции Спирмена
     * Коэффициент корреляции Спирмена - мера линейной связи между случайными величинами.
     * Корреляция Спирмена является ранговой, то есть для оценки силы связи используются не численные значения, а соответствующие им ранги.
     * Коэффициент инвариантен по отношению к любому монотонному преобразованию шкалы измерения.
     * \param x первая выборка данных
     * \param y вторая выборка данных
     * \param p коэффициент корреляции Спирмена (от -1 до +1)
     * \return вернет 0 в случае успеха
     */
    int calc_spearman_rank_correlation_coefficient(std::vector<T>& x, std::vector<T> &y, double &p)
    {
            if(x.size() != y.size() || x.size() == 0) {
                return INVALID_PARAMETER;
            }
            // найдем ранги элементов
            std::vector<int> rx(x.size());
            std::vector<int> ry(y.size());
            quick_sort(x, rx);
            quick_sort(y, ry);
            T sum = 0;
            for(size_t i = 0; i < x.size(); ++i) {
                T diff = rx[i] - ry[i];
                sum += diff * diff;
            }

            double n = x.size();
            p = 1.0 - (6.0 /(n * (n - 1) * (n + 1))) * sum;
            return OK;
    }

    /** \brief Коэффициент корреляции Пирсона
     * Коэффициент корреляции Пирсона характеризует существование линейной зависимости между двумя величинами
     * \param x первая выборка данных
     * \param y вторая выборка данных
     * \param rxy коэффициент корреляции Пирсона (от -1 до +1)
     * \return вернет 0 в случае успеха
     */
    int calc_pearson_correlation_coefficient(std::vector<T>& x, std::vector<T> &y, double &rxy)
    {
        if(x.size() != y.size() || x.size() == 0) {
            return INVALID_PARAMETER;
        }
        T xm = std::accumulate(x.begin(), x.end(), T(0));
        T ym = std::accumulate(y.begin(), y.end(), T(0));
        xm /= (T)x.size();
        ym /= (T)y.size();
        T sum = 0, sumx2 = 0, sumy2 = 0;
        for(size_t i = 0; i < x.size(); ++i) {
            T dx = x[i] - xm;
            T dy = y[i] - ym;
            sum += dx * dy;
            sumx2 += dx * dx;
            sumy2 += dy * dy;
        }
        if(sumx2 == 0 || sumy2 == 0) {
                return INVALID_PARAMETER;
        }
        rxy = sum / std::sqrt(sumx2 * sumy2);
        return OK;
    }

    int calc_euclidean_distance(std::vector<T>& x, std::vector<T> &y, double &d, double max_data) {
        double sum = 0;
        if(x.size() != y.size() || x.size() == 0) {
            std::cout << "error: x.size() != y.size() || x.size() == 0; x.size() = " << x.size() << " y.size() = " << y.size() << std::endl;
            return INVALID_PARAMETER;
        }
        for(size_t i = 0; i < x.size(); ++i) {
            double diff = x[i] - y[i];
            sum = diff * diff;
        }
        d = std::sqrt(sum)/max_data;
        return OK;
    }

    /** \brief Найти расстояние между двумя образцами
     * \param x первая выборка данных
     * \param y вторая выборка данных
     * \return вернет 0 в случае успеха
     */
    int calc_distance(std::vector<T> &x, std::vector<T> &y, double &dist, int type) {
        int err = OK;
        if(type == EUCLIDEAN_DISTANCE) {
            err = calc_euclidean_distance(x, y, dist, max_euclidean_distance);
            if(err != OK)   return err;
        } else
        if(type == PEARSON_CORRELATION) {
            err = calc_pearson_correlation_coefficient(x, y, dist);
            dist /= 2.0;
            dist = 0.5 - dist;
            if(err != OK)   return err;
        } else
        if(type == SPEARMAN_RANK_CORRELATION) {
            err = calc_spearman_rank_correlation_coefficient(x, y, dist);
            dist /= 2.0;
            dist = 0.5 - dist;
            if(err != OK)   return err;
        }
        return OK;
    }

    /** \brief Класс для хранения объектов
     */
    class DataObject {
        public:
        int class_id = 0;       /**< Класс объекта */
        int label = 0;          /**< Метка объекта */
        std::vector<T> data;    /**< Данные объекта */

        DataObject();

        DataObject(int _label, std::vector<T> &_data) {
            label = _label;
            data = _data;
        }
    };

    double max_euclidean_distance = 1;  /**< Максимальное Евкливдово расстояние */
    std::vector<DataObject> objects;    /**< Массив объектов */


    int normalized(std::vector<DataObject> &objects, std::vector<std::vector<T>> &data) {
        // делаем поиск минимума и максимума
        std::vector<T> vec_max_data;
        std::vector<T> vec_min_data;
        for(size_t i = 0; i < objects.size(); ++i) {
            T max_data = *std::max_element(objects[i].data.begin(), objects[i].data.end());
            T min_data = *std::min_element(objects[i].data.begin(), objects[i].data.end());
            vec_max_data.push_back(max_data);
            vec_min_data.push_back(min_data);
        }
        T max_data = *std::max_element(vec_max_data.begin(), vec_max_data.end());
        T min_data = *std::min_element(vec_min_data.begin(), vec_min_data.end());
        T ampl = max_data - min_data;
        data.resize(objects.size());
        for(size_t i = 0; i < objects.size(); ++i) {
            data[i].resize(objects[i].data.size());
            for(size_t j = 0; j < objects[i].data.size(); ++j) {
                data[i][j] = 2.0 * ((T)(objects[i].data[j] - min_data) / ampl) - 1.0;
                std::cout << data[i][j] << " " << std::endl;
            }
            std::cout << std::endl;
        }
        return OK;
    }

    /** \brief Назначить кластеры
     * Метод Перебирает в цикле каждый элемент данных, назначая его кластеру, который сопоставлен с ближайшим текущим средним/центроидом
     * \param data данные
     * \param clustering кластеры
     * \param means средине кластеров
     * \param type тип функции для нахождения меры соответствия
     * \param result вернет true если было переназначение
     * \return вернет 0 в случае успеха
     *
     */
    int update_clustering(std::vector<std::vector<T>> &data, std::vector<int> &clustering, std::vector<std::vector<T>> &means, int type, bool &result) {
        result = false;
        std::vector<int> _clustering(data.size()); // не будем менять входной вектор, прежде чем не пройдем проверку...
        for(size_t i = 0; i < data.size(); ++i) {
            double min_distance;
            int err = calc_distance(data[i], means[0], min_distance, type);
            if(err != OK)   return err;
            _clustering[i] = 0;
            for(size_t j = 1; j < means.size(); ++j) {
                double new_min_distance;
                int err = calc_distance(data[i], means[j], new_min_distance, type);
                if(err != OK)   return err;
                if(new_min_distance < min_distance) {
                    _clustering[i] = j;
                    min_distance = new_min_distance;
                    result = true;
                }
            }
        }

        std::vector<int> num_clustering_obj(means.size());
        for(size_t i = 0; i < _clustering.size(); ++i) {
            num_clustering_obj[_clustering[i]]++;
        }
        for(size_t i = 0; i < num_clustering_obj.size(); ++i) {
            if(num_clustering_obj[i] == 0) {
                result = false;
                return OK;
            }
        }
        clustering = _clustering;
        return OK;
    }

    /** \brief
     * Метод перебирает в цикле данные, назначенные каждому кластеру, и вычисляет новое среднее/центроид для каждого кластера.
     * Этот метод установит result в false, если одно или более средних нельзя вычислить из-за того, что в кластере нет никаких элементов данных.
     * \param
     * \param
     * \return
     *
     */
    int update_means(std::vector<std::vector<T>> &data, std::vector<int> &clustering, std::vector<std::vector<T>> &means, int type, bool &result) {
        result = false;
        // переберем все кластеры
        for(size_t k = 0; k < means.size(); ++k) {
            // переберем все данные в поисках нового центра
            double min_sum = 0;
            // находим минимальное расстояние от старого центра
            for(size_t i = 0; i < data.size(); ++i) {
                if(clustering[i] != k) continue; // не наш кластер
                double min_distance;
                int err = calc_distance(data[i], means[k], min_distance, type);
                if(err != OK)   return err;
                min_sum += min_distance;
                result = true;
            }
            if(!result) return OK;
            // пробуем обновить центр
            for(size_t i = 0; i < data.size(); ++i) {
                if(clustering[i] != k) continue; // не наш кластер
                // иначе находим среднее расстояние
                double sum = 0;
                for(size_t j = 0; j < data.size(); ++j) {
                    if(clustering[j] != k || j == i) continue; // не наш кластер или тот же элемент
                    double min_distance;
                    int err = calc_distance(data[i], data[j], min_distance, type);
                    if(err != OK)   return err;
                    sum += min_distance;
                }
                if(sum < min_sum && result) {
                    min_sum = sum;
                    means[k] = data[i];
                }
            }
        }
        return OK;
    }

public:

    ClusterAnalysis() {};

    /** \brief Добавить сэмпл
     * \param data данные
     * \param label метка пользователя
     * \param нормализовать данные?
     * \return вернет 0 в случае успеха
     */
    int add_sample(std::vector<T> &data, int label, bool is_norm) {
        if(is_norm) {
            std::vector<T> norm_data;
            int err = calc_min_max(data, norm_data, MINMAX_1_1);
            if(err != OK) return err;
            objects.push_back(DataObject(label, norm_data));
        } else {
            objects.push_back(DataObject(label, data));
        }
        return OK;
    }

    /** \brief Пороговый алгоритм кластеризации
     * \param thresold пороговый уровень
     * \param type тип функции для нахождения меры соответствия
     * \return вернет 0 в случае успеха
     */
    int update_thresold_algorithm(double thresold, int type) {
        // массив центра классов
        std::vector<T> class_center(1);
        class_center[0] = objects[0].data;
        int number_of_classes = 1; // количество классов
        //
        //std::vector<int> class_of_objects(objects.size());
        //class_of_objects[0] = 0;
        objects[0].class_id = 0;
        for(int i = 1; i < objects.size(); i++) {
            double min_distance = 0;
            int err = calc_distance(class_center[0], objects[i].data, min_distance, type);
            if(err != OK)   return err;
            objects[i].class_id = 0;
            for(int j = 0; j < class_center.size(); j++) {
                double min_distance_2 = 0;
                int err = calc_distance(class_center[j], objects[i].data, min_distance_2, type);
                if(err != OK)   return err;
                if(min_distance_2 < min_distance) {
                    min_distance = min_distance_2;
                    objects[i].class_id = j;
                }
            }
            if(min_distance > thresold) {
                objects[i].class_id = class_center.size();
                class_center.resize(class_center.size() + 1);
                class_center[class_center.size() - 1] = objects[i].data;
            }
        }
        return OK;
    }

private:

    /** \brief Инициализировать первоначальные средние
     * Метод InitMeans реализует механизм инициализации k-средних++ и возвращает набор средних, которые находятся далеко друг от друга в терминах евклидова расстояния.
     * \param numClusters
     * \param data
     * \param seed
     * \return
     */
    int init_means(
            int numClusters,
            std::vector<std::vector<T>> &data,
            std::vector<std::vector<T>> &out,
            int type,
            int seed) {
        std::vector<std::vector<T>> means(numClusters); // хранит то, что вернул метод
        for(size_t i = 0; i < means.size(); ++i) {
            means[i].resize(data[0].size());
        }
        /* used содержит индексы элементов данных, назначенные в качестве начальных средних,
            чтобы можно предотвратить дублирование начальных средних*/
        std::vector<int> used;
        rand_seed(seed);
        int idx = random_number(0, data.size() - 1);
        means[0] = data[idx];
        used.push_back(idx);
        for(int k = 1; k < numClusters; ++k) {
            std::vector<double> dSquared(data.size());
            int newMean = -1;
            for(int i = 0; i < data.size(); ++i) {
                if(std::binary_search(used.begin(),used.end(),i)) continue;
                std::vector<double> distances(k);
                for (int j = 0; j < k; ++j) {
                    int err = calc_distance(data[i], means[0], distances[j], type);
                    if(err != OK)   {
                        std::cout << "error: calc distance" << std::endl;
                        return err;
                    }
                }
                double min_distances = *std::min_element(distances.begin(), distances.end());
                dSquared[i] = min_distances * min_distances;
            }
            double p = custom_rand();
            double sum = 0.0;
            for (int i = 0; i < dSquared.size(); ++i) {
                sum += dSquared[i];
            }
            double cumulative = 0.0;
            int ii = 0;
            int sanity = 0;
            while(sanity < data.size() * 2) {
                cumulative += dSquared[ii] / sum;
                if (cumulative >= p && std::binary_search(used.begin(),used.end(),ii) == false) {
                    newMean = ii; // выбранный индекс
                    used.push_back(newMean);
                    std::sort(used.begin(), used.end());
                    break;
                }
                ++ii; // следующий кандидат
                if (ii >= dSquared.size()) ii = 0; // мимо конца
                ++sanity;
            }
            means[k] = data[newMean];
        } // k, каждое остающееся среднее
        out = means;
        return OK;
    } // InitMean

    /** \brief Алгоритм k-средних ++
     *
     * \param
     * \param
     * \return
     *
     */
    int cluster(
            std::vector<DataObject> &raw_data,
            int numClusters,
            int type,
            int seed) {
        std::vector<std::vector<T>> data;
        int err = OK;
        err = normalized(raw_data, data);
        if(err != OK)   {
            std::cout << "error: normalized" << std::endl;
            return err;
        }
        bool changed = true;
        bool success = true;

        std::vector<std::vector<T>> means;
        err = init_means(numClusters, data, means, type, seed);
        if(err != OK)   {
            std::cout << "error: init means" << std::endl;
            return err;
        }

        std::vector<int> clustering(data.size());
        int maxCount = data.size() * 10;
        int ct = 0;
        while(changed == true &&
              success == true &&
              ct < maxCount) {

            err = update_clustering(data, clustering, means, type, changed);
            if(err != OK)   {
                std::cout << "error: update clustering" << std::endl;
                return err;
            }
            err = update_means(data, clustering, means, type, success);
            if(err != OK)   {
                std::cout << "error: update means" << std::endl;
                return err;
            }
            ++ct;
        }
        for(size_t i = 0; i < clustering.size(); ++i) {
            raw_data[i].class_id = clustering[i];
        }
        return OK;
    }

public:

    int update_k_mean_pp(int num_clusters, int type, int seed) {
        return cluster(objects, num_clusters, type, seed);
    }

    void get_object(int num_clusters, std::vector<std::vector<T>> &data, std::vector<int> &label, std::vector<int> &indx) {
        for(size_t i = 0; i < objects.size(); ++i) {
            if(objects[i].class_id == num_clusters) {
                data.push_back(objects[i].data);
                label.push_back(objects[i].label);
                indx.push_back(i);
            }
        }
    }
};

#endif // CLUSTER_ANALYSIS_HPP_INCLUDED
