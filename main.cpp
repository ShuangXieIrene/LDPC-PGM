#include <iostream>
#include <vector>
#include <memory>

#include "Utilitis.h"

#define sigma std::sqrt(0.5)

using namespace std;


int main( int argc, char** argv ) {
    //num is the number of codewords
    int num = 1000;
    int cycle = 10;
    int flag;
    vector<vector<int>> input_codeword;
    vector<vector<int>> transmitted_codeword;
    vector<vector<double>> observed_codeword;
    vector<vector<vector<double>>> Prob_x;
    vector<double> prob_final_tmp(7);
    vector<vector<double>> prob_final;

    cout << "Please choose the algorithm you want to use(sum-product algorithm, "
         << endl
         << "Please enter 0 || Max product algorithm, please enter 1) : ";

    cin >> flag;

    auto ldpc = std::make_unique<Utilitis<int, double>>();


    std::vector<std::vector<int>> source;
    std::vector<std::vector<std::vector<double>>> prob_xz_all;
    for (auto i = 0; i < num; i++) {
        auto source_code = ldpc->generate_source_code();
//        cout << "source_code : " << endl;
//        ldpc->print(source_code);

        auto encoded_code = ldpc->encoder(source_code);
//        cout << "encoded_code : " << endl;
//        ldpc->print(encoded_code);


        auto observed_data = ldpc->transmit(encoded_code, 0, sigma);
//        cout << "observed_code : " << endl;
//        ldpc->print(observed_data);

        auto prob_xz = ldpc->normal_distribution(observed_data, sigma);
//        cout << "prob_xz : " << endl;
//        ldpc->print(prob_xz);

        source.push_back(encoded_code);
        prob_xz_all.push_back(prob_xz);
    }

//    std::vector<std::vector<std::vector<int>>> predict;
    for (auto i = 0; i < cycle; i++) {
        std::vector<std::vector<int>> predict;
        for (auto j = 0; j < num; j++) {
            auto prob_x = prob_xz_all[j];
            prob_x = ldpc->predict_prob(prob_x, prob_xz_all[j], 0);
//            cout << "prob_x in " << i << " iter : " << endl;
//            ldpc->print(prob_x);

            auto predict_code = ldpc->predict_value(prob_x);
//            cout << "predict_value : " << endl;
//            ldpc->print(predict_code);

            predict.push_back(predict_code);
        }
        auto error_rate = ldpc->error_rate(source, predict);
        cout << "error_rate in " << i << " iter : " << error_rate << endl;
    }
//        auto predict_code = ldpc->predict_value(prob_x);
//        cout << "predict_value : " << endl;
//        ldpc->print(predict_code);


}

//    //initialize the input codewords&transmitted codewords, input codewords are stored in a 2-dimensional array; each codeword stored in tmp, and initialize to 0
//    vector<int> tmp(7);
//    for (int i = 0; i < num; i++) {
//        input_codeword.push_back(tmp);
//    }
//    for (int j = 0; j < 7; j++) {
//        tmp.at(j) = tmp.at(j) + 1;
//    }
//
//    //encoder
//    for (int i = 0; i < num; i++) {
//        vector<double> tmp_obsvered(7);
//        transmitted_codeword.push_back(tmp);
//        //add Guassian noise to the transmitted codeword
//        for (int j = 0; j < 7; j++) {
//            tmp_obsvered.at(j) =
//                    transmitted_codeword.at(i).at(j) + Utilitis::gaussianNoise(0.0, sigma);
//        }
//        observed_codeword.push_back(tmp_obsvered);
//    }
//
//    //decoder process, calculate the probability P(X^/Z)
//    vector<vector<double>> Prob_0;
//    vector<vector<double>> Prob_1;
//    vector<double> tmp_0(7);
//    vector<double> tmp_1(7);
//    vector<double> prob_x01(2);
//    vector<vector<double>> prob_x7(7, vector<double>(2));
//    for (int i = 0; i < num; i++) {
//        for (int j = 0; j < 7; j++) {
//            double prob_0 = Utilitis::normalDistribution(0, sigma,
//                                                         observed_codeword.at(i).at(j));
//            double prob_1 = Utilitis::normalDistribution(1, sigma,
//                                                         observed_codeword.at(i).at(j));
//            //Normalize prob_0 & prob_1
//            double prob_0_norm = prob_0 / (prob_0 + prob_1);
//            double prob_1_norm = prob_1 / (prob_0 + prob_1);
//            prob_x01[0] = prob_0_norm;
//            prob_x01[1] = prob_1_norm;
//            tmp_0.at(j) = prob_0_norm;
//            tmp_1.at(j) = prob_1_norm;
//            prob_x7[j] = prob_x01;
//        }
//        Prob_0.push_back(tmp_0);
//        Prob_1.push_back(tmp_1);
//        Prob_x.push_back(prob_x7);
//    }
//
//
//    //apply sum-product and Max-product
////    vector<vector<double>> Prob_x;
//    vector<int> othervalue(3);
//    for (int m = 0; m < cycle; m++) {
//        for (int i = 0; i < num; i++) {
//            for (int j = 0; j < 7; j++) {
//                if (j == 0) {
//                    othervalue = {3, 5, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 1) {
//                    othervalue = {3, 6, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 2) {
//                    othervalue = {1, 5, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                    othervalue = {2, 6, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 3) {
//                    othervalue = {5, 6, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 4) {
//                    othervalue = {1, 3, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                    othervalue = {4, 6, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 5) {
//                    othervalue = {4, 5, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                    othervalue = {2, 3, 7};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                } else if (j == 6) {
//                    othervalue = {1, 3, 5};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                    othervalue = {4, 5, 6};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                    othervalue = {2, 3, 6};
//                    if (flag == 0) {
//                        Prob_x[i][j] = Utilitis::sumProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    } else {
//                        Prob_x[i][j] = Utilitis::maxProduct(i, j, Prob_x, othervalue, Prob_0[i][j]);
//                    }
//                }
//                 (Prob_x[i][j][0] > Prob_x[i][j][1]) ? prob_final_tmp[j] = 0 : prob_final_tmp[j] = 1;
//            }
//            prob_final.push_back(prob_final_tmp);
//            cout << "END!" <<endl;
//        }
//    }
//}











