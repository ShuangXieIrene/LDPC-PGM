//
// Created by shuangxie on 9/27/17.
//

#ifndef DECODER_UTILITIS_H
#define DECODER_UTILITIS_H

#include <random>
#include <ctime>
#include <cmath>
#include <chrono>
#include <iostream>

template <class T, class U>
class Utilitis {
public:
    std::vector<T> generate_source_code();
    std::vector<T> encoder(const std::vector<T> &);
    std::vector<U> transmit(const std::vector<T> &, const U &, const U &);
    std::vector<std::vector<U>> normal_distribution(const std::vector<U> &_observed_data, const U &_sigma);
    std::vector<std::vector<U>> predict_prob(const std::vector<std::vector<U>> &, const std::vector<std::vector<U>> & _prob_xz, const T &);
    std::vector<T> predict_value(const std::vector<std::vector<U>> &);
    U error_rate(const std::vector<std::vector<T>> & source_code, const std::vector<std::vector<T>> & predict_code);

    template<class A>
    inline void print(const std::vector<A> & _code){
        for(auto & code : _code) {
            std::cout << code << " ";
        }
        std::cout << std::endl;
        return;
    }

    template<class A>
    inline void print(const std::vector<std::vector<A>> & _code){
        for(auto & code : _code) {
           print(code);
        }
        std::cout << std::endl;
        return;
    }

protected:
    U gaussian_noise(const U & _mean, const U & _s);
    std::vector<std::vector<U>> sum_product(const std::vector<std::vector<std::vector<std::vector<U>>>> & );
    std::vector<std::vector<U>> max_product(const std::vector<std::vector<std::vector<std::vector<U>>>> & );
    U LOGSUM (const U & loga,const U & logb);

//    static double normalDistribution(int _flag, double _variable, double _z);
//    static std::vector<double> sumProduct(int i, int j, std::vector<std::vector<std::vector<double>>> _prob_f,
//                                          std::vector<int> _otherx, double _probx_z);
//    static std::vector<double> maxProduct(int i, int j, std::vector<std::vector<std::vector<double>>> _prob_f,
//                                          std::vector<int> _otherx, double _probx_z);

};
#endif //DECODER_UTILITIS_H
