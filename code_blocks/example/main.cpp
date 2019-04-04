#include <iostream>
#include "cluster_analysis.hpp"



int main()
{
    std::cout << "Hello world!" << std::endl;
    std::vector<std::vector<double>> test_data(20);

    test_data[0].push_back(65); test_data[0].push_back(220);
    test_data[1].push_back(73); test_data[1].push_back(160);
    test_data[2].push_back(59); test_data[2].push_back(110);
    test_data[3].push_back(61); test_data[3].push_back(120);
    test_data[4].push_back(75); test_data[4].push_back(150);
    test_data[5].push_back(67); test_data[5].push_back(240);
    test_data[6].push_back(68); test_data[6].push_back(230);
    test_data[7].push_back(70); test_data[7].push_back(220);
    test_data[8].push_back(62); test_data[8].push_back(130);
    test_data[9].push_back(66); test_data[9].push_back(210);
    test_data[10].push_back(77); test_data[10].push_back(190);
    test_data[11].push_back(75); test_data[11].push_back(180);
    test_data[12].push_back(74); test_data[12].push_back(170);
    test_data[13].push_back(70); test_data[13].push_back(210);
    test_data[14].push_back(61); test_data[14].push_back(110);
    test_data[15].push_back(58); test_data[15].push_back(100);
    test_data[16].push_back(66); test_data[16].push_back(230);
    test_data[17].push_back(59); test_data[17].push_back(120);
    test_data[18].push_back(68); test_data[18].push_back(210);
    test_data[19].push_back(61); test_data[19].push_back(130);

    ClusterAnalysis<double> iClusterAnalysis;

    for(size_t i = 0; i < test_data.size(); ++i) {
        // добавляем пример
        std::cout << i << " " << test_data[i][0] << " " << test_data[i][1] << std::endl;
        iClusterAnalysis.add_sample(test_data[i], 0,  false);
    }

    const int MAX_CLUSTER = 3;
    iClusterAnalysis.update_k_mean_pp(MAX_CLUSTER, iClusterAnalysis.EUCLIDEAN_DISTANCE,123);
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
    return 0;
}
