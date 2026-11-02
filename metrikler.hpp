// Ödev 7 (görev 9): iki resmin piksel değerleri arasındaki benzerliği/farkı ölçen metrikler (sadece başlık dosyası).
// Hepsi tek kanallı 8-bit resimler için yazılmıştır.
//
//   mseImg   : Ortalama karesel hata                 mean((a - b)^2)                                   küçük = benzer
//   maeImg   : Ortalama mutlak hata (SAD / piksel)   mean(|a - b|)                                     küçük = benzer
//   psnrImg  : Tepe sinyal/gürültü oranı (dB)        10 log10(255^2 / MSE)                            büyük = benzer
//   nccImg   : Normalize çapraz korelasyon           sum(a b) / sqrt(sum(a^2) sum(b^2))               1 = benzer
//   znccImg  : Sıfır ortalamalı normalize çapraz kor. sum((a-ma)(b-mb)) / sqrt(sum((a-ma)^2) sum((b-mb)^2))  [-1, 1]
//   ssimImg  : Yapısal benzerlik (SSIM, 7x7 düzgün pencere, ortalaması)                               1 = benzer

#pragma once

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <cmath>
#include <limits>

namespace sinyal {

inline double mseImg(const cv::Mat& a, const cv::Mat& b) {
    double s = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            const double d = static_cast<double>(a.at<uchar>(y, x)) - b.at<uchar>(y, x);
            s += d * d;
        }
    return s / (static_cast<double>(a.rows) * a.cols);
}

inline double maeImg(const cv::Mat& a, const cv::Mat& b) {
    double s = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) s += std::fabs(static_cast<double>(a.at<uchar>(y, x)) - b.at<uchar>(y, x));
    return s / (static_cast<double>(a.rows) * a.cols);
}

inline double psnrImg(const cv::Mat& a, const cv::Mat& b) {
    const double mse = mseImg(a, b);
    return mse > 0.0 ? 10.0 * std::log10(255.0 * 255.0 / mse) : 999.0;
}

inline double nccImg(const cv::Mat& a, const cv::Mat& b) {
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            const double p = a.at<uchar>(y, x), q = b.at<uchar>(y, x);
            sab += p * q;
            saa += p * p;
            sbb += q * q;
        }
    const double d = std::sqrt(saa * sbb);
    return d > 0.0 ? sab / d : std::numeric_limits<double>::quiet_NaN();
}

// Resimlerden biri sabitse (varyansı 0) payda 0 olur ve ZNCC tanımsızdır (NaN döner).
inline double znccImg(const cv::Mat& a, const cv::Mat& b) {
    const double n = static_cast<double>(a.rows) * a.cols;
    double ma = 0.0, mb = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            ma += a.at<uchar>(y, x);
            mb += b.at<uchar>(y, x);
        }
    ma /= n;
    mb /= n;
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            const double p = a.at<uchar>(y, x) - ma, q = b.at<uchar>(y, x) - mb;
            sab += p * q;
            saa += p * p;
            sbb += q * q;
        }
    const double d = std::sqrt(saa * sbb);
    return d > 1e-12 ? sab / d : std::numeric_limits<double>::quiet_NaN();
}

// SSIM: her pikselin 7x7 komşuluğunda ortalama (parlaklık), varyans (kontrast) ve kovaryans (yapı) karşılaştırılır
//   SSIM = (2 mu1 mu2 + C1)(2 s12 + C2) / ((mu1^2 + mu2^2 + C1)(s1^2 + s2^2 + C2)),  C1 = (0.01*255)^2, C2 = (0.03*255)^2
inline double ssimImg(const cv::Mat& a8, const cv::Mat& b8) {
    cv::Mat a, b;
    a8.convertTo(a, CV_64F);
    b8.convertTo(b, CV_64F);
    const cv::Size win(7, 7);
    cv::Mat mu1, mu2, aa, bb, ab;
    cv::blur(a, mu1, win);
    cv::blur(b, mu2, win);
    cv::blur(a.mul(a), aa, win);
    cv::blur(b.mul(b), bb, win);
    cv::blur(a.mul(b), ab, win);
    const double C1 = (0.01 * 255.0) * (0.01 * 255.0), C2 = (0.03 * 255.0) * (0.03 * 255.0);
    double sum = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            const double m1 = mu1.at<double>(y, x), m2 = mu2.at<double>(y, x);
            const double s1 = aa.at<double>(y, x) - m1 * m1, s2 = bb.at<double>(y, x) - m2 * m2, s12 = ab.at<double>(y, x) - m1 * m2;
            sum += ((2.0 * m1 * m2 + C1) * (2.0 * s12 + C2)) / ((m1 * m1 + m2 * m2 + C1) * (s1 + s2 + C2));
        }
    return sum / (static_cast<double>(a.rows) * a.cols);
}

}  // namespace sinyal
