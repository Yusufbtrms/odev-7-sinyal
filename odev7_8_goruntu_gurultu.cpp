// Ödev 7 - Görev 8: Resme Gauss ve tuz-biber (salt and pepper) gürültüsü ekleyip, bu resimlerin harmonik
// matrisini (2 boyutlu DFT) orijinal resmin harmonik matrisiyle karşılaştırma.
//
// Çıktılar:
//   - Gürültü kuvveti (RMS), PSNR
//   - Harmonik matrisinde yarıçap bantlarına göre ortalama genlik (orijinal / Gauss / tuz-biber)
//   - Düşük frekans enerji oranı, en güçlü harmonikler
//   - Gürültünün kendi spektrumunun ortalaması ve beyaz gürültü kuramıyla karşılaştırma: sigma * sqrt(pi / (H*W))
//   - Gürültü kuvvetine göre tarama
//   - Görüntü: ilk satır resimler, ikinci satır log spektrumlar (aynı ölçek), üçüncü satır gürültünün kendi spektrumu
//
// Kullanım: ./goruntu_gurultu [cikti_dizini] [resim_yolu=images/Cat.jpg] [sigma=20] [tuz_biber_orani=0.05]
// Derleme:  g++ -O2 odev7_8_goruntu_gurultu.cpp -o goruntu_gurultu.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp, dft2d.hpp ve gurultu_resim.hpp bu dosyayla aynı klasörde olmalıdır)

#include "dft2d.hpp"
#include "gurultu_resim.hpp"

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace fsys = std::filesystem;

static double rmsDiff(const cv::Mat& a, const cv::Mat& b) {
    double s = 0.0;
    for (int y = 0; y < a.rows; ++y)
        for (int x = 0; x < a.cols; ++x) {
            const double d = static_cast<double>(a.at<uchar>(y, x)) - b.at<uchar>(y, x);
            s += d * d;
        }
    return std::sqrt(s / (static_cast<double>(a.rows) * a.cols));
}

static double psnrDb(const cv::Mat& a, const cv::Mat& b) {
    const double rms = rmsDiff(a, b);
    return rms > 0.0 ? 20.0 * std::log10(255.0 / rms) : 999.0;
}

// Yariçap bantlarina gore (DC haric) ortalama genlik; r = sqrt(fu^2 + fv^2)
static std::vector<double> bandMeans(const cv::Mat& amp, const std::vector<std::pair<double, double>>& bands) {
    const int H = amp.rows, W = amp.cols;
    std::vector<double> sum(bands.size(), 0.0), cnt(bands.size(), 0.0);
    for (int v = 0; v < H; ++v)
        for (int u = 0; u < W; ++u) {
            if (u == 0 && v == 0) continue;
            const int fu = (u <= W / 2) ? u : u - W, fv = (v <= H / 2) ? v : v - H;
            const double r = std::sqrt(static_cast<double>(fu) * fu + static_cast<double>(fv) * fv);
            for (size_t b = 0; b < bands.size(); ++b)
                if (r >= bands[b].first && r < bands[b].second) {
                    sum[b] += amp.at<double>(v, u);
                    cnt[b] += 1.0;
                    break;
                }
        }
    std::vector<double> mean(bands.size(), 0.0);
    for (size_t b = 0; b < bands.size(); ++b) mean[b] = cnt[b] > 0 ? sum[b] / cnt[b] : 0.0;
    return mean;
}

// DC haric enerjinin |fu|,|fv| <= 16 bolgesindeki orani (%)
static double lowFreqEnergyPercent(const sinyal::Spectrum2D& S) {
    double eLow = 0.0, eAll = 0.0;
    for (int v = 0; v < S.H; ++v)
        for (int u = 0; u < S.W; ++u) {
            if (u == 0 && v == 0) continue;
            const double e = std::norm(S.X[static_cast<size_t>(v) * S.W + u]);
            eAll += e;
            const int fu = (u <= S.W / 2) ? u : u - S.W, fv = (v <= S.H / 2) ? v : v - S.H;
            if (std::abs(fu) <= 16 && std::abs(fv) <= 16) eLow += e;
        }
    return 100.0 * eLow / eAll;
}

static double meanNonDc(const cv::Mat& amp) {
    double s = 0.0;
    for (int v = 0; v < amp.rows; ++v)
        for (int u = 0; u < amp.cols; ++u)
            if (!(u == 0 && v == 0)) s += amp.at<double>(v, u);
    return s / (static_cast<double>(amp.rows) * amp.cols - 1.0);
}

// Log spektrum goruntusu, ortak olcekle: log(1 + c*genlik) / log(1 + c*ref); sifir frekans ortada
static cv::Mat logImage(const cv::Mat& amp, double c, double ref) {
    const int H = amp.rows, W = amp.cols;
    cv::Mat img(H, W, CV_8UC1);
    const double denom = std::log(1.0 + c * ref);
    for (int v = 0; v < H; ++v)
        for (int u = 0; u < W; ++u) img.at<uchar>((v + H / 2) % H, (u + W / 2) % W) = cv::saturate_cast<uchar>(255.0 * std::log(1.0 + c * amp.at<double>(v, u)) / denom);
    return img;
}

int main(int argc, char** argv) {
    const fsys::path outDir = (argc >= 2) ? fsys::path(argv[1]) : fsys::path("cikti_goruntu_gurultu");
    const std::string path = (argc >= 3) ? argv[2] : "images/Cat.jpg";
    const double sigma = (argc >= 4) ? std::stod(argv[3]) : 20.0;
    const double ratio = (argc >= 5) ? std::stod(argv[4]) : 0.05;
    fsys::create_directories(outDir);
    const int N = 256;

    cv::Mat src = cv::imread(path, cv::IMREAD_GRAYSCALE);
    if (src.empty()) {
        std::cerr << "Hata: resim okunamadi: " << path << "\n";
        return 1;
    }
    cv::Mat orig;
    cv::resize(src, orig, cv::Size(N, N), 0, 0, cv::INTER_AREA);

    const cv::Mat gauss = sinyal::addGaussianNoise(orig, sigma, 1);
    const cv::Mat sp = sinyal::addSaltAndPepperNoise(orig, ratio, 2);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "RESIM: " << path << " (" << N << "x" << N << "'e kucultuldu)\n";
    std::cout << "PARAMETRELER: Gauss sigma = " << std::setprecision(1) << sigma << ", tuz-biber orani = " << std::setprecision(3) << ratio << " (bozulan piksel sayisi ~ " << std::lround(ratio * N * N) << ")\n\n";

    // spektrumlar
    const sinyal::Spectrum2D S0 = sinyal::dft2dComplex(orig), Sg = sinyal::dft2dComplex(gauss), Ss = sinyal::dft2dComplex(sp);
    const cv::Mat amp0 = sinyal::harmonicsFromSpectrum(S0), ampG = sinyal::harmonicsFromSpectrum(Sg), ampS = sinyal::harmonicsFromSpectrum(Ss);

    // gurultunun kendi spektrumu: (gurultulu - orijinal)
    cv::Mat dG, dS;
    cv::Mat o64, g64, s64;
    orig.convertTo(o64, CV_64F);
    gauss.convertTo(g64, CV_64F);
    sp.convertTo(s64, CV_64F);
    dG = g64 - o64;
    dS = s64 - o64;
    const sinyal::Spectrum2D SnG = sinyal::dft2dComplex(dG), SnS = sinyal::dft2dComplex(dS);
    const cv::Mat ampNG = sinyal::harmonicsFromSpectrum(SnG), ampNS = sinyal::harmonicsFromSpectrum(SnS);

    // 1) gurultu kuvveti
    std::cout << "GURULTU KUVVETI (orijinale gore, kirpma sonrasi gercek degerler):\n";
    std::cout << "  " << std::left << std::setw(20) << "Resim" << std::right << std::setw(14) << "gurultu RMS" << std::setw(12) << "PSNR (dB)" << std::setw(18) << "DC (ort. parlak.)" << "\n";
    std::cout << "  " << std::left << std::setw(20) << "orijinal" << std::right << std::setw(14) << "-" << std::setw(12) << "-" << std::setprecision(3) << std::setw(18) << amp0.at<double>(0, 0) << "\n";
    const double rmsG = rmsDiff(gauss, orig), rmsS = rmsDiff(sp, orig);
    std::cout << "  " << std::left << std::setw(20) << "Gauss gurultulu" << std::right << std::setw(14) << rmsG << std::setprecision(2) << std::setw(12) << psnrDb(gauss, orig) << std::setprecision(3) << std::setw(18) << ampG.at<double>(0, 0) << "\n";
    std::cout << "  " << std::left << std::setw(20) << "tuz-biber gurultulu" << std::right << std::setw(14) << rmsS << std::setprecision(2) << std::setw(12) << psnrDb(sp, orig) << std::setprecision(3) << std::setw(18) << ampS.at<double>(0, 0) << "\n\n";

    // 2) yaricap bantlarina gore ortalama genlik
    const std::vector<std::pair<double, double>> bands = {{1, 8}, {8, 16}, {16, 32}, {32, 64}, {64, 128}, {128, 182}};
    const std::vector<double> m0 = bandMeans(amp0, bands), mg = bandMeans(ampG, bands), ms = bandMeans(ampS, bands);
    std::cout << "HARMONIK MATRISI: yaricap bantlarina gore ortalama genlik (r = sqrt(fu^2 + fv^2); DC haric)\n";
    std::cout << "  " << std::left << std::setw(14) << "Bant (r)" << std::right << std::setw(12) << "orijinal" << std::setw(12) << "Gauss" << std::setw(14) << "tuz-biber" << std::setw(14) << "Gauss/orij." << std::setw(16) << "tuz-biber/orij." << "\n";
    std::cout << std::setprecision(4);
    for (size_t b = 0; b < bands.size(); ++b)
        std::cout << "  " << std::left << std::setw(14) << (std::to_string(static_cast<int>(bands[b].first)) + " - " + std::to_string(static_cast<int>(bands[b].second))) << std::right << std::setw(12) << m0[b]
                  << std::setw(12) << mg[b] << std::setw(14) << ms[b] << std::setprecision(2) << std::setw(14) << mg[b] / m0[b] << std::setw(16) << ms[b] / m0[b] << std::setprecision(4) << "\n";

    // 3) dusuk frekans enerji orani ve en guclu harmonikler
    std::cout << "\n  DC haric enerjinin |fu|,|fv| <= 16 dusuk frekans bolgesindeki orani: orijinal = " << std::setprecision(2) << lowFreqEnergyPercent(S0) << " %, Gauss = " << lowFreqEnergyPercent(Sg)
              << " %, tuz-biber = " << lowFreqEnergyPercent(Ss) << " %\n" << std::setprecision(3);
    auto printTop = [&](const std::string& name, const cv::Mat& amp) {
        std::cout << "  En guclu 3 harmonik, " << std::left << std::setw(20) << name << std::right;
        for (const sinyal::SpecPeak& p : sinyal::topPeaks(amp, 3)) std::cout << " (fu=" << p.fu << ", fv=" << p.fv << ", " << std::setprecision(2) << p.amp << ")" << std::setprecision(3);
        std::cout << "\n";
    };
    printTop("orijinal", amp0);
    printTop("Gauss gurultulu", ampG);
    printTop("tuz-biber gurultulu", ampS);

    // 4) gurultunun kendi spektrumu ve beyaz gurultu kurami
    std::cout << "\nGURULTUNUN KENDI SPEKTRUMU (gurultulu - orijinal): ortalama genlik (DC haric)\n";
    std::cout << "  Gauss      : olculen = " << meanNonDc(ampNG) << ", beyaz gurultu kurami RMS*sqrt(pi/(H*W)) = " << rmsG * std::sqrt(sinyal::PI / (static_cast<double>(N) * N)) << "\n";
    std::cout << "  tuz-biber  : olculen = " << meanNonDc(ampNS) << ", beyaz gurultu kurami RMS*sqrt(pi/(H*W)) = " << rmsS * std::sqrt(sinyal::PI / (static_cast<double>(N) * N)) << "\n";
    {
        const std::vector<double> nb = bandMeans(ampNG, bands), nb2 = bandMeans(ampNS, bands);
        std::cout << "  Gurultu spektrumunun bant ortalamalari (dusuk bant -> yuksek bant):\n    Gauss    :";
        for (double v : nb) std::cout << " " << v;
        std::cout << "\n    tuz-biber:";
        for (double v : nb2) std::cout << " " << v;
        std::cout << "\n";
    }

    // 5) tarama
    std::cout << "\nTARAMA (ayni resim; yuksek frekans = r >= 64 bandinin ortalama genligi, orijinalde " << std::setprecision(4) << bandMeans(amp0, {{64, 182}})[0] << "):\n";
    std::cout << "  " << std::left << std::setw(26) << "Gurultu" << std::right << std::setw(12) << "PSNR (dB)" << std::setw(20) << "yuksek frek. genlik" << std::setw(22) << "dusuk frek. enerji %" << "\n";
    auto sweepRow = [&](const std::string& name, const cv::Mat& img) {
        const sinyal::Spectrum2D S = sinyal::dft2dComplex(img);
        const cv::Mat amp = sinyal::harmonicsFromSpectrum(S);
        std::cout << "  " << std::left << std::setw(26) << name << std::right << std::setprecision(2) << std::setw(12) << psnrDb(img, orig) << std::setprecision(4) << std::setw(20) << bandMeans(amp, {{64, 182}})[0]
                  << std::setprecision(2) << std::setw(22) << lowFreqEnergyPercent(S) << "\n";
    };
    for (double sg : {10.0, 20.0, 40.0}) sweepRow("Gauss sigma = " + std::to_string(static_cast<int>(sg)), sinyal::addGaussianNoise(orig, sg, 1));
    for (double r : {0.01, 0.05, 0.10}) sweepRow("tuz-biber orani = " + sinyal::fmt(r, 2), sinyal::addSaltAndPepperNoise(orig, r, 2));

    // 6) goruntuler
    double refAmp = 0.0;
    for (int v = 0; v < N; ++v)
        for (int u = 0; u < N; ++u)
            if (!(u == 0 && v == 0)) refAmp = std::max(refAmp, amp0.at<double>(v, u));
    const double c = 50.0;
    auto panelGray = [&](const cv::Mat& g, const std::string& label) {
        cv::Mat big, bgr;
        cv::resize(g, big, cv::Size(), 2.0, 2.0, cv::INTER_NEAREST);
        cv::cvtColor(big, bgr, cv::COLOR_GRAY2BGR);
        cv::putText(bgr, label, cv::Point(8, 22), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
        return bgr;
    };
    cv::Mat blank = cv::Mat::zeros(2 * N, 2 * N, CV_8UC3);
    cv::putText(blank, "gurultunun kendi spektrumu", cv::Point(20, 240), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    cv::putText(blank, "(gurultulu - orijinal)", cv::Point(20, 270), cv::FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
    cv::Mat r1, r2, r3, all;
    cv::hconcat(std::vector<cv::Mat>{panelGray(orig, "orijinal"), panelGray(gauss, "Gauss sigma=" + sinyal::fmt(sigma, 0)), panelGray(sp, "tuz-biber oran=" + sinyal::fmt(ratio, 2))}, r1);
    cv::hconcat(std::vector<cv::Mat>{panelGray(logImage(amp0, c, refAmp), "orijinal spektrum"), panelGray(logImage(ampG, c, refAmp), "Gauss spektrum"), panelGray(logImage(ampS, c, refAmp), "tuz-biber spektrum")}, r2);
    cv::hconcat(std::vector<cv::Mat>{blank, panelGray(logImage(ampNG, c, refAmp), "Gauss gurultu spektrumu"), panelGray(logImage(ampNS, c, refAmp), "tuz-biber gurultu spektrumu")}, r3);
    cv::vconcat(std::vector<cv::Mat>{r1, r2, r3}, all);
    cv::imwrite((outDir / "8_karsilastirma.png").string(), all);
    cv::imwrite((outDir / "8_gauss_gurultulu.png").string(), gauss);
    cv::imwrite((outDir / "8_tuz_biber_gurultulu.png").string(), sp);

    std::cout << "\nKaydedilenler (" << outDir.string() << "):\n"
              << "  8_karsilastirma.png         : ust: orijinal / Gauss / tuz-biber resimler; orta: log spektrumlar (ayni olcek); alt: gurultunun kendi spektrumu\n"
              << "  8_gauss_gurultulu.png, 8_tuz_biber_gurultulu.png : gurultulu resimler\n";
    return 0;
}
