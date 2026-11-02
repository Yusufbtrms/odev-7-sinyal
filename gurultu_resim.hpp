// Ödev 7 (görev 8): resme gürültü ekleyen fonksiyonlar (sadece başlık dosyası).

#pragma once

#include <opencv2/core.hpp>

#include <random>
#include <stdexcept>

namespace sinyal {

// Gauss gürültüsü: her piksele ortalaması 0, standart sapması sigma olan rastgele bir sayı eklenir,
// sonuç 0..255 aralığına sınırlanır.
//   img   : tek kanallı (gri) 8-bit resim
//   sigma : gürültü kuvveti (standart sapma, parlaklık birimi; örn. 10 hafif, 40 güçlü)
inline cv::Mat addGaussianNoise(const cv::Mat& img, double sigma, unsigned seed = 1) {
    if (img.type() != CV_8UC1) throw std::invalid_argument("tek kanalli 8-bit resim gerekli");
    if (sigma < 0.0) throw std::invalid_argument("sigma negatif olamaz");
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, sigma);
    cv::Mat out(img.size(), CV_8UC1);
    for (int y = 0; y < img.rows; ++y)
        for (int x = 0; x < img.cols; ++x) out.at<uchar>(y, x) = cv::saturate_cast<uchar>(img.at<uchar>(y, x) + noise(rng));
    return out;
}

// Tuz-biber (salt and pepper) gürültüsü: piksellerin "ratio" oranı rastgele seçilir; seçilenlerin yarısı
// siyah (0, biber), yarısı beyaz (255, tuz) yapılır. Diğer pikseller değişmez.
//   ratio : bozulan piksel oranı, [0, 1] (örn. 0.05 = piksellerin %5'i)
inline cv::Mat addSaltAndPepperNoise(const cv::Mat& img, double ratio, unsigned seed = 2) {
    if (img.type() != CV_8UC1) throw std::invalid_argument("tek kanalli 8-bit resim gerekli");
    if (ratio < 0.0 || ratio > 1.0) throw std::invalid_argument("oran [0, 1] araliginda olmali");
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    cv::Mat out = img.clone();
    for (int y = 0; y < img.rows; ++y)
        for (int x = 0; x < img.cols; ++x)
            if (u(rng) < ratio) out.at<uchar>(y, x) = (u(rng) < 0.5) ? 0 : 255;
    return out;
}

}  // namespace sinyal
