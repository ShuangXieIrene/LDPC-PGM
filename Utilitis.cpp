//
// Created by shuangxie on 9/27/17.
//

#include "Utilitis.h"
#include <iostream>
#include <chrono>
#include <algorithm>
const auto inverse_sqrt_2pi = 0.3989422804014327;


/**
 * function generate_source_code() is used to generate 4 bit code (x) by truely random
 * @return output 4 bit code
 */
template <class T, class U>
std::vector<T> Utilitis<T, U>::generate_source_code() {
    std::vector<T> source_code;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937_64 generator (seed);
    std::uniform_int_distribution<int> distribution(0, 1);
    for (auto i = 0; i < 4; i++)
        source_code.push_back(std::move(distribution(generator)));

    return source_code;
}

template std::vector<int> Utilitis<int, double>::generate_source_code();

/**
 * function encoder is used to generate encoded code (y) from the source code by Hamming coding
 * @param _raw_code input 4 bit code
 * @return output 7 bit code
 */
template <class T, class U>
std::vector<T> Utilitis<T, U>::encoder(const std::vector<T> & _raw_code) {
    if(_raw_code.size() != 4){
        std::cout << "Error: the size of input code (x) is not 4!" << std::endl;
    }
    std::vector<int> output(_raw_code.begin(), _raw_code.end());

    T x_5 = abs(output[1] + output[2] - output[3]) % 2;
    T x_6 = abs(output[0] + output[2] - output[3]) % 2;
    T x_7 = (output[3] + output[4] + output[5]) % 2;

    output.push_back(std::move(x_5));
    output.push_back(std::move(x_6));
    output.push_back(std::move(x_7));

    if(output.size() != 7){
        std::cout << "Error: the size of transmit code (y) is not 7!" << std::endl;
    }

    // transmit it into -1 & 1
    for(auto & code : output)
        code = 1 - 2 * code;

    return output;
}

template std::vector<int> Utilitis<int, double>::encoder(const std::vector<int> &);


/**
 * Function transmit() is used to generate the transmitted code which is observed value of z
 * @tparam T
 * @tparam U
 * @param _encoded_code
 * @param _mean
 * @param _sigma
 * @return
 */
template <class T, class U>
std::vector<U> Utilitis<T, U>::transmit(const std::vector<T> & _encoded_code,
                                        const U & _mean,
                                        const U & _sigma){
    std::vector<U> output;
    for(auto & code : _encoded_code)
        output.push_back(code + gaussian_noise(_mean, _sigma));
    return output;
}

template std::vector<double> Utilitis<int, double>::transmit(const std::vector<int> & , const double & , const double & );

/**
 *
 * @tparam T
 * @tparam U
 * @param _mean
 * @param _sigma
 * @return
 */
template <class T, class U>
U Utilitis<T, U>::gaussian_noise(const U & _mean, const U & _sigma) {
    std::random_device rd;
    std::mt19937 e2(rd());
    //generate the normal distribution
    std::normal_distribution<U> dist(_mean,_sigma);
    return dist(e2);
}

template double Utilitis<int, double>::gaussian_noise(const double & _mean, const double & _sigma);

/**
 * Function normal_distribution() is used to calculate the probablity of x under condition of z
 * @tparam T
 * @tparam U
 * @param _observed_data
 * @param _sigma
 * @return
 */
template <class T, class U>
std::vector<std::vector<U>> Utilitis<T,U>::normal_distribution(const std::vector<U> & _observed_data, const U & _sigma) {
    std::vector<std::vector<U>> prob_xz;

    for (auto & prob : _observed_data){
        std::vector<U> tmp(2);
        tmp[0] = (prob - 1) / _sigma;
        tmp[0] = inverse_sqrt_2pi / _sigma * std::exp(-0.5f * tmp[0] * tmp[0]);
        tmp[1] = (prob + 1) / _sigma;
        tmp[1] = inverse_sqrt_2pi / _sigma * std::exp(-0.5f * tmp[1] * tmp[1]);
        auto sum = tmp[0] + tmp[1];
//        tmp[0] = tmp[0]/sum;
//        tmp[1] = tmp[1]/sum;
        tmp[0] = log(tmp[0]/sum);
        tmp[1] = log(tmp[1]/sum);
        prob_xz.push_back(std::move(tmp));
    }

    return prob_xz;
}
template std::vector<std::vector<double>> Utilitis<int, double>::normal_distribution(const std::vector<double> & _observed_data, const double & _sigma);

std::vector<std::vector<int>> hamming_code = {
        {0, 0, 0, 1, 1, 1, 1},
        {0, 1, 1, 0, 0, 1, 1},
        {1, 0, 1, 0, 1, 0, 1}
};

template <class T, class U>
std::vector<std::vector<U>> Utilitis<T,U>::predict_prob(const std::vector<std::vector<U>> & _prob_x,
                                                        const std::vector<std::vector<U>> & _prob_xz,
                                                        const T & _flag) {
    auto prob_x = _prob_x;

    std::vector<std::vector<std::vector<U>>> factor(hamming_code.size(),
                                                    std::vector<std::vector<U>>(_prob_x.size(),
                                                                std::vector<U>(_prob_x[0].size(), 1)));

    for(auto h = 0; h < hamming_code.size(); h++) {
        std::vector<int> index;
        for(auto i = 0; i < hamming_code[h].size(); i++){
            if(hamming_code[h][i] == 1)
                index.push_back(i);
        }

        std::vector<std::vector<std::vector<std::vector<U>>>> product;
        for(auto h = 0; h < _prob_x[index[0]].size(); h++){
            std::vector<std::vector<std::vector<U>>> product_1;
            for(auto i = 0; i < _prob_x[index[1]].size(); i++){
                std::vector<std::vector<U>> product_2;
                for(auto j = 0; j < _prob_x[index[2]].size(); j++){
                    std::vector<U> product_3;
                    for(auto k = 0; k < _prob_x[index[3]].size(); k++){
                        if((h + i + j + k) % 2 == 0){
//                            auto p = _prob_x[index[0]][h] * _prob_x[index[1]][i]
//                                     * _prob_x[index[2]][j] * _prob_x[index[3]][k];
                            auto p = _prob_x[index[0]][h] + _prob_x[index[1]][i]
                                     + _prob_x[index[2]][j] + _prob_x[index[3]][k];
                            product_3.push_back(std::move(p));
                        }
                        else{
                            product_3.push_back(log(0));
                        }
                    }
                    product_2.push_back(std::move(product_3));
                }
                product_1.push_back(std::move(product_2));
            }
            product.push_back(product_1);
        }

        std::vector<std::vector<U>> temp;
        if(_flag == 0)
            temp = std::move(sum_product(product));
        else if (_flag == 1)
            temp = std::move(max_product(product));
        else
            std::cout << "input flag error" << std::endl;


        for(auto i = 0; i < index.size(); i++){
            factor[h][index[i]] = temp[i];
        }
    }

    for(auto i = 0; i < _prob_x.size(); i++){
        for(auto j = 0; j < _prob_x[j].size(); j++){
//            prob_x[i][j] *= _prob_xz[i][j];
            prob_x[i][j] = _prob_xz[i][j];
            for(auto fac : factor){
//                prob_x[i][j] *= fac[i][j];
                prob_x[i][j] += fac[i][j];
            }
        }

//        auto sum = 0.0;
        auto sum = log(0.0);
        for(auto j = 0; j < _prob_x[j].size(); j++){
//            sum += prob_x[i][j];
            if(sum != prob_x[i][j])
                sum = LOGSUM(sum, prob_x[i][j]);
        }

        for(auto j = 0; j < _prob_x[j].size(); j++){
//            prob_x[i][j] = prob_x[i][j]/sum;
            prob_x[i][j] = prob_x[i][j] - sum;
        }
    }

    return prob_x;
};

template std::vector<std::vector<double>> Utilitis<int,double>::predict_prob(const std::vector<std::vector<double>> & _prob_x,
                                                        const std::vector<std::vector<double>> & _prob_xz,
                                                        const int & _flag);

template <class T, class U>
std::vector<std::vector<U>> Utilitis<T,U>::sum_product(const std::vector<std::vector<std::vector<std::vector<U>>>> & _product){
    std::vector<std::vector<U>> output(4);

    for(auto h = 0; h < _product.size(); h++){
//        auto sum = 0.0;
        auto sum = log(0.0);
        for(auto i = 0; i < _product[h].size(); i++){
            for(auto j = 0; j < _product[h][i].size(); j++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
                    if(sum != _product[h][i][j][k])
                        sum = LOGSUM(sum, _product[h][i][j][k]);
//                  sum += _product[h][i][j][k];
                }
            }
        }
        output[0].push_back(std::move(sum));
    }

    for(auto i = 0; i < _product[0].size(); i++){
//        auto sum = 0.0;
        auto sum = log(0.0);
        for(auto h = 0; h < _product.size(); h++){
            for(auto j = 0; j < _product[h][i].size(); j++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
//                    sum += _product[h][i][j][k];
                    if(sum != _product[h][i][j][k])
                        sum = LOGSUM(sum, _product[h][i][j][k]);
                }
            }
        }
        output[1].push_back(std::move(sum));
    }

    for(auto j = 0; j < _product[0][0].size(); j++){
//        auto sum = 0.0;
        auto sum = log(0.0);
        for(auto h = 0; h < _product.size(); h++){
            for(auto i = 0; i < _product[h].size(); i++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
//                    sum += _product[h][i][j][k];
                    if(sum != _product[h][i][j][k])
                        sum = LOGSUM(sum, _product[h][i][j][k]);
                }
            }
        }
        output[2].push_back(std::move(sum));
    }

    for(auto k = 0; k < _product[0][0][0].size(); k++){
//        auto sum = 0.0;
        auto sum = log(0.0);
        for(auto h = 0; h < _product.size(); h++){
            for(auto i = 0; i < _product[h].size(); i++){
                for(auto j = 0; j < _product[h][i].size(); j++){
//                    sum += _product[h][i][j][k];
                    if(sum != _product[h][i][j][k])
                        sum = LOGSUM(sum, _product[h][i][j][k]);
                }
            }
        }
        output[3].push_back(std::move(sum));
    }

    return output;
};

template std::vector<std::vector<double>> Utilitis<int,double>::sum_product(const std::vector<std::vector<std::vector<std::vector<double>>>> & _product);

template <class T, class U>
std::vector<std::vector<U>> Utilitis<T,U>::max_product(const std::vector<std::vector<std::vector<std::vector<U>>>> & _product){
    std::vector<std::vector<U>> output(4);

    for(auto h = 0; h < _product.size(); h++){
        auto max = 0.0;
        for(auto i = 0; i < _product[h].size(); i++){
            for(auto j = 0; j < _product[h][i].size(); j++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
                    max = std::max<U>(max, _product[h][i][j][k]);
                }
            }
        }
        output[0].push_back(std::move(max));
    }

    for(auto i = 0; i < _product[0].size(); i++){
        auto max = 0.0;
        for(auto h = 0; h < _product.size(); h++){
            for(auto j = 0; j < _product[h][i].size(); j++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
                    max = std::max<U>(max, _product[h][i][j][k]);
                }
            }
        }
        output[1].push_back(std::move(max));
    }

    for(auto j = 0; j < _product[0][0].size(); j++){
        auto max = 0.0;
        for(auto h = 0; h < _product.size(); h++){
            for(auto i = 0; i < _product[h].size(); i++){
                for(auto k = 0; k < _product[h][i][j].size(); k++){
                    max = std::max<U>(max, _product[h][i][j][k]);
                }
            }
        }
        output[2].push_back(std::move(max));
    }

    for(auto k = 0; k < _product[0][0][0].size(); k++){
        auto max = 0.0;
        for(auto h = 0; h < _product.size(); h++){
            for(auto i = 0; i < _product[h].size(); i++){
                for(auto j = 0; j < _product[h][i].size(); j++){
                    max = std::max<U>(max, _product[h][i][j][k]);
                }
            }
        }
        output[3].push_back(std::move(max));
    }

    return output;
};

template std::vector<std::vector<double>> Utilitis<int,double>::max_product(const std::vector<std::vector<std::vector<std::vector<double>>>> & _product);

template <class T, class U>
std::vector<T> Utilitis<T,U>::predict_value(const std::vector<std::vector<U>> & _prob){
    std::vector<T> output(_prob.size());
    for(auto i = 0; i < _prob.size(); i++){
        output[i] = (_prob[i][0] > _prob[i][1]) ? 0 : 1;
    }
    return output;
};

template std::vector<int> Utilitis<int,double>::predict_value(const std::vector<std::vector<double>> & _prob);

template <class T, class U>
U  Utilitis<T,U>::LOGSUM(const U & loga, const U & logb) {
    U diff = loga - logb;
    U logaplusb;
    if(diff>23)
        logaplusb = loga;
    else if (diff < -23)
        logaplusb = logb;
    else
        logaplusb = logb + log(exp(diff) + 1);
    return logaplusb;
}
template double Utilitis<int, double>::LOGSUM(const double & loga, const double & logb);

template <class T, class U>
U Utilitis<T,U>::error_rate(const std::vector<std::vector<T>> & _source_code, const std::vector<std::vector<T>> & _predict_code){
    auto output = 0.0;
    for(auto i = 0; i < _source_code.size(); i++){
        for(int j = 0 ; j < _source_code[0].size(); j++){
            if((_source_code[i][j] == 1 && _predict_code[i][j] == 1) ||
               (_source_code[i][j] == -1 && _predict_code[i][j] == 0)){
                output++;
            }
        }
    }
    output = output / (_source_code.size() * _source_code[0].size());
    return output;
}
template double Utilitis<int,double>::error_rate(const std::vector<std::vector<int>> & _source_code, const std::vector<std::vector<int>> & _predict_code);

//
//std::vector<double> Utilitis::sumProduct(int i, int j, std::vector<std::vector<std::vector<double>>> _prob_f,
//                                         std::vector<int> _otherx, double _probx_z) {
//
//
//    _prob_f[i][j][0] = _prob_f[i][j][0] * (_prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][0] +
//                                         _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][1] +
//                                         _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][1] +
//                                         _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][0] ) * _probx_z;
//
//    _prob_f[i][j][1] = _prob_f[i][j][1] * (_prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][1] +
//                                         _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][0] +
//                                         _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][0] +
//                                         _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][1] ) * _probx_z;
//
//    return _prob_f[i][j];
//
//}
//
//std::vector<double> Utilitis::maxProduct(int i, int j, std::vector<std::vector<std::vector<double>>> _prob_f,
//                                         std::vector<int> _otherx, double _probx_z) {
////    std::vector<std::vector<std::vector<double>>> prob_x;
//    auto Max = [](double _a, double _b, double _c, double _d){
//        double max = (_a > _b) ? _a : _b;
//        max = (max > _c) ? max : _c;
//        max = (max > _d) ? max : _d;
//        return max;
//    };
//
//        double prob10 = _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][0];
//        double prob20 = _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][1];
//        double prob30 = _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][1];
//        double prob40 = _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][0];
//        _prob_f[i][j][0] = _prob_f[i][j][0] * Max(prob10, prob20, prob30, prob40) * _probx_z;
//
//        double prob11 = _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][1];
//        double prob21 = _prob_f[i][_otherx[0]-1][0] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][0];
//        double prob31 = _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][0] * _prob_f[i][_otherx[2]-1][0];
//        double prob41 = _prob_f[i][_otherx[0]-1][1] * _prob_f[i][_otherx[1]-1][1] * _prob_f[i][_otherx[2]-1][1];
//        _prob_f[i][j][1] = _prob_f[i][j][1] * Max(prob11, prob21, prob31, prob41) * _probx_z;
//
//    return _prob_f[i][j];
//}
//

