// Ödev 7 - Görev 9: Gürültülü resimleri ters DFT ile orijinale döndürmeye çalışma ve benzerlik metrikleri.
//
// Bölüm A: DFT -> ters DFT tek başına resmi aynen geri verir (gürültü kalır). Bunu doğrularız.
// Bölüm B: Frekans uzayında süzme (eşik, ideal alçak geçiren, Gauss alçak geçiren) + ters DFT ile orijinale dönmeye
//          çalışırız. Her yöntem için parametre taraması yapılır, orijinalle fark ZNCC (ve diğer metrikler) ile ölçülür.
// Bölüm C: Metriklerin davranışını göstermek için, orijinalden bilinen bozulmalarla türetilmiş resimlerde
//          MSE, MAE, PSNR, NCC, ZNCC ve SSIM değerleri yan yana verilir (parlaklık/kontrast değişimi, negatif,
//          bulanıklık, kayma, döndürme, gürültü, sabit resim).
//
// Kullanım: ./geri_donus [cikti_dizini] [resim_yolu=images/Cat.jpg] [sigma=20] [tuz_biber_orani=0.05]
// Derleme:  g++ -O2 odev7_9_geri_donus_metrik.cpp -o geri_donus.exe -std=c++17 $(pkg-config --cflags --libs opencv5)
// (sinyal.hpp, dft2d.hpp, dft2d_filtre.hpp, gurultu_resim.hpp ve metrikler.hpp bu dosyayla aynı klasörde olmalıdır)

#include "dft2d_filtre.hpp"
#include "gurultu_resim.hpp"
#include "metrikler.hpp"

#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fsys = std::filesystem;

struct Row {
    std::string name;
    long long kept = 0;
    double zncc = 0, ssim = 0, psnr = 0, mse = 0, mae = 0;
    cv::Mat img;
};

static std::string num(double v, int prec) {
    if (std::isnan(v)) return "tanimsiz";
    return sinyal::fmt(v, prec);
}

static Row evaluate(const std::string& name, const cv::Mat& rec, const cv::Mat& orig, long long kept) {
    Row r;
    r.name = name;
    r.kept = kept;
    r.img = rec;
    r.zncc = sinyal::znccImg(rec, orig);
    r.ssim = sinyal::ssimImg(rec, orig);
    r.psnr = sinyal::psnrImg(rec, orig);
    r.mse = sinyal::mseImg(rec, orig);
    r.mae = sinyal::maeImg(rec, orig);
    return r;
}

static void printHeader() {
    std::cout << "  " << std::left << std::setw(26) << "Yontem / parametre" << std::right << std::setw(12) << "korunan k" << std::setw(9) << "ZNCC" << std::setw(9) << "SSIM" << std::setw(10) << "PSNR dB"
              << std::setw(10) << "MSE" << std::setw(8) << "MAE" << "\n";
}

static void printRow(const Row& r, bool best) {
    std::cout << "  " << std::left << std::setw(26) << (r.name + (best ? " *" : "")) << std::right << std::setw(12) << (r.kept > 0 ? std::to_string(r.kept) : "-") << std::setw(9) << num(r.zncc, 4)
              << std::setw(9) << num(r.ssim, 4) << std::setw(10) << num(r.psnr, 2) << std::setw(10) << num(r.mse, 1) << std::setw(8) << num(r.mae, 2) << "\n";
}

// Bir gurultulu resim icin uc sureci tara; en iyi (ZNCC'ye gore) sonucu dondur
static std::vector<Row> runNoiseCase(const std::string& title, const cv::Mat& orig, const cv::Mat& noisy) {
    std::cout << "---- " << title << " ----\n";
    const sinyal::Spectrum2D S = sinyal::dft2dComplex(noisy);

    // Bolum A: dogrulama
    {
        const cv::Mat back8 = sinyal::toImage8(sinyal::idft2dReal(S));
        double maxDiff = 0.0;
        for (int y = 0; y < noisy.rows; ++y)
            for (int x = 0; x < noisy.cols; ++x) maxDiff = std::max(maxDiff, std::fabs(static_cast<double>(back8.at<uchar>(y, x)) - noisy.at<uchar>(y, x)));
        std::cout << "  A) DFT -> ters DFT (suzme yok): gurultulu resimle en buyuk fark = " << maxDiff << " (gurultu aynen geri geldi); orijinalle ZNCC = " << num(sinyal::znccImg(back8, orig), 4)
                  << " (gurultulu resmin kendisi: " << num(sinyal::znccImg(noisy, orig), 4) << ")\n\n";
    }

    std::cout << "  B) Frekans uzayinda suzme + ters DFT ('*' = o yontemde ZNCC'ye gore en iyi parametre)\n";
    printHeader();
    std::vector<Row> best;
    {
        Row r = evaluate("gurultulu resim (suzmesiz)", noisy, orig, 0);
        printRow(r, false);
    }
    auto sweep = [&](const std::string& label, const std::vector<double>& params, const std::function<std::pair<sinyal::Spectrum2D, long long>(double)>& make) {
        std::vector<Row> rows;
        for (double p : params) {
            const auto filtered = make(p);
            rows.push_back(evaluate(label + " " + sinyal::fmt(p, 1), sinyal::toImage8(sinyal::idft2dReal(filtered.first)), orig, filtered.second));
        }
        size_t b = 0;
        for (size_t i = 1; i < rows.size(); ++i)
            if (rows[i].zncc > rows[b].zncc) b = i;
        for (size_t i = 0; i < rows.size(); ++i) printRow(rows[i], i == b);
        best.push_back(rows[b]);
    };
    sweep("esik T =", {0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 4.0}, [&](double t) {
        long long k = 0;
        sinyal::Spectrum2D y = sinyal::thresholdFilter2D(S, t, &k);
        return std::make_pair(y, k);
    });
    sweep("alcak geciren R =", {4, 8, 12, 16, 24, 32, 48, 64}, [&](double r) {
        long long k = 0;
        sinyal::Spectrum2D y = sinyal::lowPassFilter2D(S, r, &k);
        return std::make_pair(y, k);
    });
    sweep("Gauss alcak gecir. s =", {4, 8, 12, 16, 24, 32}, [&](double s) { return std::make_pair(sinyal::gaussianLowPass2D(S, s), static_cast<long long>(0)); });
    std::cout << "\n";
    return best;
}

int main(int argc, char** argv) {
    const fsys::path outDir = (argc >= 2) ? fsys::path(argv[1]) : fsys::path("cikti_geri_donus");
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

    std::cout << "RESIM: " << path << " (" << N << "x" << N << "), Gauss sigma = " << sigma << ", tuz-biber orani = " << ratio << "\n\n";

    const std::vector<Row> bestG = runNoiseCase("GAUSS GURULTULU RESIM (sigma = " + sinyal::fmt(sigma, 0) + ")", orig, gauss);
    const std::vector<Row> bestS = runNoiseCase("TUZ-BIBER GURULTULU RESIM (oran = " + sinyal::fmt(ratio, 2) + ")", orig, sp);

    // Gorsel: [orijinal | gurultulu | en iyi esik | en iyi alcak geciren | en iyi Gauss]
    auto panel = [&](const cv::Mat& g, const std::string& label, double zncc) {
        cv::Mat bgr;
        cv::cvtColor(g, bgr, cv::COLOR_GRAY2BGR);
        cv::putText(bgr, label, cv::Point(6, 18), cv::FONT_HERSHEY_SIMPLEX, 0.45, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
        if (!std::isnan(zncc)) cv::putText(bgr, "ZNCC " + sinyal::fmt(zncc, 4), cv::Point(6, 40), cv::FONT_HERSHEY_SIMPLEX, 0.45, cv::Scalar(0, 255, 255), 1, cv::LINE_AA);
        return bgr;
    };
    auto montage = [&](const cv::Mat& noisy, const std::vector<Row>& best, const std::string& noisyLabel) {
        cv::Mat row;
        cv::hconcat(std::vector<cv::Mat>{panel(orig, "orijinal", 1.0), panel(noisy, noisyLabel, sinyal::znccImg(noisy, orig)), panel(best[0].img, "esik: " + best[0].name.substr(best[0].name.find('=') + 2), best[0].zncc),
                                         panel(best[1].img, "alcak gecir.: " + best[1].name.substr(best[1].name.find('=') + 2), best[1].zncc),
                                         panel(best[2].img, "Gauss: " + best[2].name.substr(best[2].name.find('=') + 2), best[2].zncc)},
                    row);
        return row;
    };
    cv::Mat all;
    cv::vconcat(std::vector<cv::Mat>{montage(gauss, bestG, "Gauss gurultulu"), montage(sp, bestS, "tuz-biber gurultulu")}, all);
    cv::imwrite((outDir / "9_geri_donus.png").string(), all);

    // ---------------------------------------------------------------------------------------------
    // Bolum C: metriklerin davranisi
    // ---------------------------------------------------------------------------------------------
    std::cout << "---- C) METRIKLERIN DAVRANISI (referans: orijinal resim) ----\n";
    struct Variant {
        std::string name;
        cv::Mat img;
    };
    std::vector<Variant> vars;
    vars.push_back({"orijinalin kendisi", orig});
    {
        cv::Mat v;
        orig.convertTo(v, CV_8U, 0.7, 60.0);  // kazanc 0.7, kayma +60 (kirpma yok)
        vars.push_back({"kontrast 0.7x + parlaklik +60", v});
    }
    {
        cv::Mat v = 255 - orig;
        vars.push_back({"negatif (255 - resim)", v});
    }
    {
        cv::Mat v;
        cv::GaussianBlur(orig, v, cv::Size(0, 0), 1.5);
        vars.push_back({"bulaniklik (sigma 1.5)", v});
    }
    {
        cv::Mat v(N, N, CV_8UC1);
        for (int y = 0; y < N; ++y)
            for (int x = 0; x < N; ++x) v.at<uchar>(y, x) = orig.at<uchar>(y, std::max(0, x - 4));
        vars.push_back({"4 piksel saga kayma", v});
    }
    vars.push_back({"transpoz (90 derece dondurme)", orig.t()});
    vars.push_back({"Gauss gurultulu", gauss});
    vars.push_back({"tuz-biber gurultulu", sp});
    {
        const int mean = static_cast<int>(std::lround(cv::mean(orig)[0]));
        vars.push_back({"sabit resim (ortalama deger)", cv::Mat(N, N, CV_8UC1, cv::Scalar(mean))});
    }

    std::cout << "  " << std::left << std::setw(34) << "Resim" << std::right << std::setw(10) << "MSE" << std::setw(8) << "MAE" << std::setw(10) << "PSNR dB" << std::setw(10) << "NCC" << std::setw(11) << "ZNCC"
              << std::setw(9) << "SSIM" << "\n";
    for (const Variant& v : vars) {
        std::cout << "  " << std::left << std::setw(34) << v.name << std::right << std::setw(10) << num(sinyal::mseImg(v.img, orig), 1) << std::setw(8) << num(sinyal::maeImg(v.img, orig), 2) << std::setw(10)
                  << num(sinyal::psnrImg(v.img, orig), 2) << std::setw(10) << num(sinyal::nccImg(v.img, orig), 4) << std::setw(11) << num(sinyal::znccImg(v.img, orig), 4) << std::setw(9)
                  << num(sinyal::ssimImg(v.img, orig), 4) << "\n";
    }

    std::cout << "\nKaydedilenler (" << outDir.string() << "):\n"
              << "  9_geri_donus.png : ust satir Gauss, alt satir tuz-biber gurultulu; soldan saga orijinal, gurultulu, en iyi esik, en iyi alcak geciren, en iyi Gauss alcak gecirenle geri donus (ZNCC yaziyor)\n";
    return 0;
}
