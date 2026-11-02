# Ödev 7: Sinüs, Gürültü, Filtreleme ve DFT Raporu

Kodlar `sinyal.hpp` (fonksiyonlar) ve `odev7_sinyal.cpp` (ana program) dosyalarındadır. Çıktılar `cikti_sinyal/` klasöründe, ekran çıktısı `sinyal_sonuclari.txt` dosyasındadır. Görev 5 ve 6 (iki sinüsün toplamı ve frekans filtresi) için `frekans.hpp`, `odev7_5_6_frekans.cpp`, `cikti_frekans/` klasörü ve `frekans_sonuclari.txt` kullanıldı (bölüm 4). Görev 7 (2 boyutlu DFT) için `dft2d.hpp`, `odev7_7_dft2d.cpp`, `cikti_dft2d/` klasörü ve `dft2d_sonuclari.txt` kullanıldı (bölüm 5). Görev 8 (resme gürültü) için `gurultu_resim.hpp`, `odev7_8_goruntu_gurultu.cpp`, `cikti_goruntu_gurultu/` klasörü ve `goruntu_gurultu_sonuclari.txt` kullanıldı (bölüm 6). Görev 9 (ters DFT ile geri dönüş ve benzerlik metrikleri) için `dft2d_filtre.hpp`, `metrikler.hpp`, `odev7_9_geri_donus_metrik.cpp`, `cikti_geri_donus/` klasörü ve `geri_donus_sonuclari.txt` kullanıldı (bölüm 7). Görev 10 (frekans filtreleri) için `odev7_10_filtreler.cpp`, `cikti_filtreler/` klasörü ve `filtreler_sonuclari.txt` kullanıldı (bölüm 8).

## 1. Yöntem

- **Madde 1, `fillSine`:** Genlik, faz, frekans (Hz) ve örnekleme frekansı alır, vektörü `genlik * sin(2*pi*frekans*n/örnekleme + faz)` ile doldurur; vektör boyutu 2'nin kuvveti değilse hata verir (periyotla vermek için `fillSineByPeriod` da vardır). Kullanılan sinyal: N = 256, örnekleme 256 Hz (vektör 1 saniye), 8 Hz, genlik 2, faz 30°. Vektör, vektör elemanlarını nokta, sürekli sinüs dalgasını çizgi olarak gösteren grafikle çizildi.
- **Madde 2, `addNoise`:** Vektörün genliğini (en büyük mutlak değer) A bulur, her elemana `oran * A * u` ekler; u, [-1, 1] aralığında düzgün dağılımlı rastgele sayıdır. Böylece en büyük gürültü oran × A kadardır. Oran [0, 1] dışındaysa hata verir. Kullanılan oran: 0.3.
- **Madde 3:** Tek boyutlu Gauss filtresi (ağırlıklar exp(-d²/(2σ²)), toplamı 1) ve medyan filtresi; kenarlarda yansıtma. Pencere 9, Gauss σ = 2.
- **Madde 4, `dftHarmonics`:** DFT'yi doğrudan tanımından hesaplar, `X[k] = sum x[n] exp(-j 2 pi k n / N)` (sin/cos tablosu bir kere hesaplanır), k = 0..N/2 için harmonik genliklerini döndürür (DC ve N/2 harmoniği |X|/N, diğerleri 2|X|/N; böylece genlik-2 sinüs 2 olarak görünür). Faz tahmini: arg X[k] + 90° (sin için). Doğrulama olarak ters DFT yazıldı; geri elde edilen vektör ile orijinal arasındaki en büyük fark 3.1e-15 çıktı.

## 2. Sonuçlar

### Madde 2: gürültü
En büyük gürültü 0.5981 (sınır 0.6), RMS 0.3376 (kuram 0.3464), gürültülü sinyalin SNR'i 12.44 dB.

### Madde 3: filtreleme (oran 0.3, pencere 9)

| Sinyal | MSE | SNR (dB) | 8 Hz bileşen genliği (gerçek 2.0000) |
|---|---|---|---|
| Gürültülü | 0.1140 | 12.44 | 1.9882 |
| Gauss | 0.0324 | 17.90 | 1.8624 |
| Medyan | 0.0650 | 14.88 | 1.8596 |

Gürültü oranına göre MSE (pencere 9):

| Oran | Gürültülü | Gauss | Medyan |
|---|---|---|---|
| 0.1 | 0.0127 | 0.0135 | 0.0221 |
| 0.3 | 0.1140 | 0.0324 | 0.0650 |
| 0.6 | 0.4558 | 0.0929 | 0.1385 |
| 1.0 | 1.2661 | 0.2337 | 0.3160 |

Pencere boyutuna göre (oran 0.3, Gauss σ = pencere/4):

| Pencere | Gauss MSE | Gauss genlik | Medyan MSE | Medyan genlik |
|---|---|---|---|---|
| 3 | 0.0481 | 1.9707 | 0.0643 | 1.9910 |
| 5 | 0.0325 | 1.9420 | 0.0551 | 1.9543 |
| 9 | 0.0345 | 1.8445 | 0.0650 | 1.8596 |
| 15 | 0.0980 | 1.6140 | 0.1375 | 1.6643 |
| 25 | 0.4408 | 1.0998 | 0.7529 | 0.8794 |
| 41 | 1.4077 | 0.3512 | 3.3655 | 0.5069 |

Ek deney, dürtü gürültüsü (elemanların bir kısmı ±1.5·genliğe ayarlandı, pencere 9):

| Bozuk eleman oranı | Gürültülü MSE | Gauss | Medyan |
|---|---|---|---|
| 0.02 | 0.2633 | 0.0511 | 0.0201 |
| 0.05 | 0.5258 | 0.1163 | 0.0290 |
| 0.10 | 0.8165 | 0.1679 | 0.0352 |

### Madde 4: DFT harmonik vektörü (k. harmonik = k Hz)

| Sinyal | Tepe k | Tepe genlik | Faz (derece) | Taban ortalaması | Taban en büyük | Tepe / taban ortalaması |
|---|---|---|---|---|---|---|
| Temiz | 8 | 2.0000 | 30.00 | yaklaşık 2e-16 | yaklaşık 1e-5 altı | yaklaşık 8.7e15 |
| Gürültülü | 8 | 1.9882 | 29.48 | 0.03833 | 0.07965 | 51.9 |
| Gauss sonrası | 8 | 1.8624 | 29.39 | 0.01324 | 0.08083 | 140.7 |
| Medyan sonrası | 8 | 1.8596 | 29.57 | 0.02487 | 0.08017 | 74.8 |

Gürültülü vektörün en güçlü 5 harmoniği: k = 8 (1.9882), k = 2 (0.0796), k = 98 (0.0787), k = 72 (0.0773), k = 92 (0.0760).

## 3. Yorum

**Gürültü ve filtreleme (madde 3).**

- Düzgün (sınırlı) dağılımlı gürültüde Gauss filtresi medyandan yaklaşık iki kat daha düşük hata verdi (MSE 0.0324 ve 0.0650). Gauss doğrusal bir alçak geçiren filtredir; ortalama alarak rastgele gürültüyü bastırır. Medyan doğrusal olmayan bir filtredir ve bu tür gürültüde ortalamanın sağladığı kazancı sağlamaz.
- Filtreler gürültüyle birlikte sinyalin kendisini de zayıflatır: pencere 9'da 8 Hz bileşeni iki filtrede de yaklaşık %7 azaldı (2.0'dan 1.86'ya). Gauss için kendi hesabımız (kesik çekirdekle) yaklaşık 1.87 veriyor, ölçümle uyumlu.
- Gürültü azken filtre fayda yerine zarar verebilir: oran 0.1'de gürültülü sinyalin MSE'si (0.0127) Gauss (0.0135) ve medyan (0.0221) sonrasından daha düşüktü; filtrenin sinyali bozması, söndürdüğü gürültüden fazlaydı. Gürültü büyüdükçe kazanç artar: oran 1.0'da Gauss 5.4 kat, medyan 4 kat iyileştirdi.
- Pencere boyutu sinyalin periyodundan (32 örnek) çok küçük seçilmelidir. Bu denemede en iyi pencere iki filtre için de yaklaşık 5 çıktı. Pencere büyüdükçe sinyal zayıflar ve bozulur: pencere 41'de Gauss genliği 0.35'e düştü, medyanın MSE'si (3.37) sinyalin gücünden (2.0) bile büyük oldu; yani filtre sonucu, sıfır çıkışından da kötüydü.
- Dürtü gürültüsünde sonuç tersine döner: medyan, Gauss'tan 2.5 ile 4.8 kat düşük hata verdi (örneğin %5 bozuk elemanda 0.0290'a karşı 0.1163). Gauss her bozuk değeri komşularına yayar, medyan ise aykırı değeri sıralamada uca iterek atar. Medyan hatası, dürtü oranı 5 kat artınca (0.02'den 0.10'a) sadece 0.0201'den 0.0352'ye çıktı.

**DFT (madde 4).**

- Temiz sinyalin harmonik vektöründe tek tepe vardır: k = 8, genlik 2.0000, faz 30.00°; diğer bütün harmonikler sayısal hata düzeyindedir (yaklaşık 1e-16). Bunun sebebi 8 tam periyodun vektöre sığması, dolayısıyla spektral sızıntının (leakage) olmamasıdır.
- Gürültülü vektörde tepe aynı yerde ve neredeyse aynı büyüklüktedir (genlik 1.9882, yani %0.6 sapma; faz 29.48°). Gürültü ise bütün harmoniklere dağılır ve bir gürültü tabanı oluşturur: ortalama 0.0383. Beyaz gürültü için beklenen ortalama σ·√(π/N) = 0.0374'tür, ölçümle uyumlu. En güçlü gürültü harmonikleri (k = 2, 98, 72, 92) rastgele konumlardadır ve en büyüğü 0.08 civarıdır, tepenin yaklaşık %4'ü.
- Tepe/taban oranı (51.9, yaklaşık 34 dB) zaman uzayındaki SNR'den (12.44 dB) çok yüksektir. Sinyal enerjisi tek harmonikte toplanırken gürültü enerjisi 128 harmonike yayılır; bu yüzden frekans uzayında sinyali gürültüden ayırmak çok daha kolaydır.
- Filtre sonrası spektrumda Gauss gürültü tabanını 0.0383'ten 0.0132'ye (yaklaşık üçte bire), medyan 0.0249'a indirdi; ancak 8 Hz tepesi de %7 küçüldü (1.8624 ve 1.8596). Tepe/taban oranı Gauss'ta 140.7'ye, medyanda 74.8'e yükseldi. Tabanın en büyük değeri üç durumda da yaklaşık 0.08 kaldı; olası açıklama, en güçlü gürültü harmoniğinin düşük frekansta (k = 2) olması ve Gauss'un alçak geçiren yapısı yüzünden geçmesidir. Bu ayrıca kontrol edilmedi.

## 4. Görev 5 ve 6: iki sinüsün toplamı ve frekans filtresi

**Kurulum.** N = 256, örnekleme 256 Hz (her harmonik 1 Hz). Sinyal 1: 50 Hz, genlik 12; sinyal 2: 30 Hz, genlik 5; faz 0°. Her birine bağımsız düzgün dağılımlı gürültü eklendi (oran 0.3, gürültü her sinyalin kendi genliğiyle orantılı). Toplam sinyal, iki gürültülü sinyalin toplamıdır.

### 4.1 Görev 5: ayrı ve toplam DFT analizi

| Sinyal | Tepe k | Tepe genlik | Faz | Taban ortalaması | Taban en büyük |
|---|---|---|---|---|---|
| Sinyal 1 (50 Hz), gürültülü | 50 | 12.0172 | -1.97° | 0.2282 | 0.4789 |
| Sinyal 2 (30 Hz), gürültülü | 30 | 5.1670 | 0.49° | 0.0956 | 0.2210 |
| Toplam | 50 ve 30 | 12.0108 ve 5.0388 | -2.05° ve 1.75° | 0.2409 | 0.5684 |

- Gürültü tabanı sinyalin genliğiyle orantılıdır (0.2282 ve 0.0956; oran 2.39, genlik oranı 2.4).
- Toplam sinyalde gürültü RMS'i 2.1478'dir; beklenen taban ortalaması 2.1478 × sqrt(pi/256) = 0.2379 olup ölçülenle (0.2409) uyumludur. Tabanın dışındaki en büyük harmonik 0.5684'tür (k = 72).
- DFT doğrusaldır: |DFT(toplam) - DFT(sinyal 1) - DFT(sinyal 2)| en büyük değeri 1.8e-12.
- Zaman uzayında toplam sinyalin SNR'i 12.63 dB'dir; frekans uzayında en zayıf bileşen (30 Hz) tabanın 20.9 katıdır. İki bileşen ayrı harmoniklerde (k = 30 ve k = 50) net biçimde ayrışır; toplam sinyalin zaman uzayı grafiğinden bu iki frekansı gözle ayırmak zordur.

### 4.2 Görev 6: frekans filtresi

`frequencyFilter` fonksiyonu sinyalin DFT'sini alır, harmonik genliği eşiğin altında kalan harmonikleri (ve gerçek sinyalin DFT'sinin eşlenik simetrisi gereği N-k karşılıklarını) sıfırlar ve kalan harmoniklerle ters DFT yaparak sinyali yeniden üretir. Eşik 0 iken çıkışın girdiyle aynı olduğu doğrulandı (en büyük fark 3.2e-14).

Eşik 1.5 için sonuçlar (temiz sinyale göre hata):

| Sinyal | Gürültülü MSE | Filtreli MSE | Gürültülü SNR | Filtreli SNR | Korunan harmonikler |
|---|---|---|---|---|---|
| Sinyal 1 (50 Hz) | 4.1199 | 0.0857 | 12.42 dB | 29.25 dB | k = 50 |
| Sinyal 2 (30 Hz) | 0.7395 | 0.0149 | 12.28 dB | 29.24 dB | k = 30 |
| Toplam sinyal | 4.6132 | 0.1050 | 12.63 dB | 29.06 dB | k = 30, 50 |

Eşik taraması (toplam sinyal, gürültü oranı 0.3):

| Eşik | Korunan harmonik sayısı | MSE | SNR (dB) |
|---|---|---|---|
| 0.0 | 129 | 4.6132 | 12.63 |
| 0.1 | 113 | 4.5673 | 12.67 |
| 0.3 | 38 | 2.8689 | 14.69 |
| 0.6 | 2 | 0.1050 | 29.06 |
| 1.0 | 2 | 0.1050 | 29.06 |
| 1.5 | 2 | 0.1050 | 29.06 |
| 3.0 | 2 | 0.1050 | 29.06 |
| 5.5 | 1 | 12.5925 | 8.27 |
| 13.0 | 0 | 84.5000 | 0.00 |

Gürültü oranına göre (eşik 1.5, toplam sinyal):

| Oran | Gürültülü MSE | Filtreli MSE | Korunan harmonik sayısı |
|---|---|---|---|
| 0.1 | 0.5126 | 0.0117 | 2 |
| 0.3 | 4.6132 | 0.1050 | 2 |
| 0.6 | 18.4528 | 0.4199 | 2 |
| 1.0 | 51.2578 | 7.0698 | 6 |

Zaman uzayı filtreleriyle karşılaştırma (toplam sinyal, MSE): gürültülü 4.6132, Gauss (pencere 9, sigma 2) 72.4820, medyan (pencere 9) 127.7567, frekans filtresi (eşik 1.5) 0.1050.

**Yorum.**

- Eşik 1.5'te yalnızca sinyalin gerçek harmonikleri korundu ve toplam sinyalde hata 4.6132'den 0.1050'ye, SNR 12.63 dB'den 29.06 dB'ye (yaklaşık 16 dB) iyileşti; iki ayrı sinyalde de benzer iyileşme (yaklaşık 17 dB) görüldü.
- Kalan hata (0.1050), filtrenin bastıramadığı kısımdır: korunan iki harmonikteki gürültü. Hesabımıza göre 50 Hz'de -2.05° faz hatası yaklaşık 0.43 genlik birimlik bir vektör hatasına (güç 0.092), 30 Hz'de genlik ve faz hatası birlikte yaklaşık 0.16'ya (güç 0.012) karşılık gelir; toplam 0.104, ölçülen 0.1050'ye çok yakındır. Yani sinyalin kendi frekansına düşen gürültü bu yöntemle ayrılamaz.
- Doğru eşik aralığı, en büyük gürültü harmoniği (0.5684) ile en zayıf sinyal harmoniği (5.0388) arasındadır; bu aralıktaki her eşik (0.6, 1.0, 1.5, 3.0) aynı sonucu verir. Eşik düşükse gürültü harmonikleri de korunur (eşik 0.3'te 38 harmonik, hata hâlâ yüksek); eşik 5.5'te zayıf 30 Hz bileşeni de atılır ve hata (12.5925) bu bileşenin gücüne (5²/2 = 12.5) yaklaşır; eşik 13'te hiçbir harmonik kalmaz ve hata sinyalin gücüne ((12² + 5²)/2 = 84.5) eşit olur.
- Eşik 0'da filtre hiçbir şey atmaz ve sonuç gürültülü sinyalle aynıdır (MSE 4.6132).
- Gürültü oranı 0.1, 0.3 ve 0.6 iken hata gürültünün karesiyle orantılı artar (0.0117, 0.1050, 0.4199). Oran 1.0'da dört gürültü harmoniği eşiği (1.5) aştığı için 6 harmonik korundu ve hata 7.0698'e çıktı; yine de gürültülü sinyalden (51.2578) yaklaşık 7 kat iyidir. Gürültü büyüdükçe eşiğin de yükseltilmesi gerekir.
- Gauss ve medyan, bu sinyalde gürültülü sinyalden (4.6132) çok daha kötü sonuç verdi (72.4820 ve 127.7567). Sebep, 50 Hz ve 30 Hz bileşenlerinin hızlı salınımlar olmasıdır (periyotları yaklaşık 5.1 ve 8.5 örnek): 9 örneklik pencere bunları da siler ya da bozar. Gauss (sigma 2) için yaklaşık hesabımıza göre 50 Hz bileşeninin genliği 12'den 0.6 civarına düşer. Medyanın hatası sinyalin gücünden (84.5) bile büyüktür. Frekans filtresi ise sinyalin dar bantlı (iki harmonik) ve gürültünün tüm harmoniklere yayılmış olmasından yararlanır, bu yüzden aynı gürültüde çok daha başarılıdır.

## 5. Görev 7: 2 boyutlu DFT

**Fonksiyon.** `dft2dHarmonics(resim)`, gri bir resmi alır ve harmonik matrisini (H x W) döndürür. Her hücre, o frekanstaki sinüsün genliğidir (tek boyutlu `dftHarmonics` ile aynı anlam; eşlenik hücreler aynı değeri taşır); (0, 0) hücresi resmin ortalama parlaklığıdır. İstenirse faz matrisi de döndürülür. 2 boyutlu DFT ayrılabilir olduğu için iki aşamada hesaplanır: önce her satırın, sonra çıkan sonucun her sütununun tek boyutlu DFT'si (doğrudan tanım, sin/cos tabloları bir kere hesaplanır): `X[v][u] = sum f[y][x] exp(-j 2 pi (u x / W + v y / H))`. Resmin eni ve boyu 2'nin kuvveti olmalıdır. Gösterim için spektrum, sıfır frekans ortada olacak şekilde kaydırılıp logaritmik ölçekle çizildi.

**Test.** Kendi harmoniklerini bildiğimiz sentetik desenlerle ve gerçek bir resimle (256x256):

| Test | Ters DFT hatası | DC genliği (ortalama parlaklık) | En güçlü harmonikler (fu, fv, genlik) | Beklenen | DC dışı enerjinin düşük frekanstaki oranı* |
|---|---|---|---|---|---|
| A) dikey çizgiler | 4.8e-13 | 128.0000 (128.0000) | (10, 0, 99.976) | (10, 0, 100) | %100.00 |
| B) çapraz çizgiler | 6.7e-13 | 128.0000 (128.0000) | (8, 16, 99.972) | (8, 16, 100) | %100.00 |
| C) iki frekans | 9.4e-13 | 128.0000 (128.0000) | (10, 0, 59.981), (0, 24, 30.005) | (10, 0, 60), (0, 24, 30) | %79.98 |
| D) beyaz dikdörtgen 64x32 | 8.7e-13 | 7.9688 (7.9688) | (0, 1, 15.531), (0, 2, 14.350), (1, 0, 14.349) | - | %92.45 |
| E) gerçek resim (Cat.jpg) | 4.7e-12 | 153.7456 (153.7456) | (1, 0, 49.623), (0, 1, 23.436), (1, 1, 16.615) | - | %95.13 |

\* |fu| ve |fv| değerleri en fazla 16 olan bölge. Her testte Parseval teoremi doğrulandı (uzaydaki enerji, spektrumdaki enerjinin H x W'ye bölümüne eşit). Her test yaklaşık 47-52 ms sürdü.

**Yorum.**

- Sentetik desenlerde harmonik matrisi beklenen yerlerde ve beklenen genliktedir (en büyük sapma 0.03); geri kalan 0.1-0.2 genlikli küçük tepeler, piksel değerlerinin 8-bit tam sayıya yuvarlanmasından kaynaklanıyor olabilir (ayrıca doğrulanmadı).
- C'deki %79.98 oranı kuramla uyumludur: iki sinüsün enerjisi genliklerinin karesiyle orantılıdır, 60² / (60² + 30²) = %80, ve (0, 24) bileşeni (30 genlikli) |fv| ≤ 16 bölgesinin dışındadır.
- D'de DC genliği 255 x (64 x 32) / (256 x 256) = 7.96875 ile birebir uyuşur. (0, 1) ve (1, 0) harmoniklerinin genlikleri, dikdörtgenin sinc türü spektrumundan hesaplanan değerlerle (15.53 ve 14.349) örtüşür. Keskin kenarlı bir şeklin enerjisi düşük frekansta yoğundur (%92.45) ama kenarlardan kaynaklanan yüksek frekans bileşenleri de vardır.
- E'de, DC dışındaki enerjinin %95.13'ü spektrumun yalnızca yaklaşık %1.7'sini (33x33 hücre, 256x256 içinde) oluşturan düşük frekans bölgesindedir. En güçlü harmonikler (1, 0), (0, 1) ve (1, 1) gibi en düşük frekanslardadır; yani resmin kaba yapısı enerjinin büyük kısmını taşır, ince ayrıntılar (yüksek frekanslar) çok az enerji taşır. Doğal resimlerde bu beklenen bir özelliktir ve düşük geçiren filtreler ve sıkıştırma yöntemlerinin dayanağıdır.
- Doğrudan tanımdan hesaplanan 2 boyutlu DFT'nin maliyeti yaklaşık 2 x N^3'tür (256x256 için yaklaşık 34 milyon işlem). Bu yüzden gerçek resim 256x256'ya küçültüldü; 1600x1598'lik orijinal için hesap yaklaşık 250 kat daha uzun sürerdi ve boyutlar 2'nin kuvveti de değil. Hızlı Fourier dönüşümü (FFT) bu maliyeti çok düşürür; denenmedi.

## 6. Görev 8: resme Gauss ve tuz-biber gürültüsü, harmonik matrisi karşılaştırması

**Fonksiyonlar.** `addGaussianNoise(resim, sigma)`: her piksele ortalaması 0, standart sapması `sigma` olan rastgele bir sayı ekler, sonucu 0-255 arasına sınırlar (gürültü kuvveti: sigma). `addSaltAndPepperNoise(resim, oran)`: piksellerin `oran` kadarını rastgele seçer, yarısını siyah (0), yarısını beyaz (255) yapar. **Deney:** Cat.jpg, 256x256 gri; Gauss sigma = 20; tuz-biber oranı 0.05 (yaklaşık 3277 piksel). Üç resmin de 2 boyutlu DFT harmonik matrisi hesaplandı (Görev 7'deki fonksiyonla).

| Resim | Gürültü RMS | PSNR (dB) | DC genliği | DC dışı enerjinin düşük frekanstaki oranı* | (1, 0) harmoniği genliği |
|---|---|---|---|---|---|
| Orijinal | - | - | 153.746 | %95.13 | 49.62 |
| Gauss gürültülü | 19.682 | 22.25 | 153.420 | %84.78 | 49.43 |
| Tuz-biber gürültülü | 32.199 | 17.97 | 152.397 | %70.24 | 47.05 |

\* |fu| ve |fv| değerleri en fazla 16 olan bölge.

Yarıçap bantlarına göre ortalama genlik (r = sqrt(fu² + fv²), DC hariç):

| Bant (r) | Orijinal | Gauss | Tuz-biber | Gauss / orijinal | Tuz-biber / orijinal |
|---|---|---|---|---|---|
| 1 - 8 | 4.5300 | 4.5090 | 4.3061 | 1.00 | 0.95 |
| 8 - 16 | 1.0197 | 1.0244 | 0.9914 | 1.00 | 0.97 |
| 16 - 32 | 0.3517 | 0.3774 | 0.3993 | 1.07 | 1.14 |
| 32 - 64 | 0.1071 | 0.1774 | 0.2494 | 1.66 | 2.33 |
| 64 - 128 | 0.0333 | 0.1416 | 0.2257 | 4.25 | 6.78 |
| 128 - 182 | 0.0177 | 0.1377 | 0.2211 | 7.79 | 12.51 |

Gürültünün kendi spektrumu (gürültülü - orijinal): ortalama genlik (DC hariç) Gauss için ölçülen 0.136 (beyaz gürültü kuramı, gürültü RMS x sqrt(pi/(H x W)): 0.136), tuz-biber için ölçülen 0.222 (kuram 0.223). Bant ortalamaları (düşük banttan yükseğe): Gauss 0.143, 0.131, 0.136, 0.136, 0.136, 0.137; tuz-biber 0.328, 0.220, 0.219, 0.223, 0.223, 0.221.

Gürültü kuvvetine göre tarama (yüksek frekans = r >= 64 bandı; orijinalde 0.0291):

| Gürültü | PSNR (dB) | Yüksek frekans genliği | Düşük frekans enerji oranı (%) |
|---|---|---|---|
| Gauss sigma = 10 | 28.14 | 0.0765 | 92.25 |
| Gauss sigma = 20 | 22.25 | 0.1406 | 84.78 |
| Gauss sigma = 40 | 16.61 | 0.2626 | 64.55 |
| Tuz-biber oranı = 0.01 | 24.57 | 0.1093 | 88.62 |
| Tuz-biber oranı = 0.05 | 17.97 | 0.2245 | 70.24 |
| Tuz-biber oranı = 0.10 | 15.07 | 0.3095 | 54.09 |

**Yorum.**

- İki gürültü de düşük frekanslı bileşenleri (resmin kaba yapısı) büyük ölçüde korur, yüksek frekans bantlarındaki genliği ise katlar (64-128 bandında Gauss'ta 4.25, tuz-biberde 6.78 kat; en yüksek bantta 7.79 ve 12.51 kat). Düşük frekans enerji oranı gürültü arttıkça düşer (tarama tablosunda tek yönlü).
- Her iki gürültünün spektrumu da yaklaşık düzdür (beyaz); düzeyi gürültünün gücüne (RMS) bağlıdır, dağılım türüne değil. Ölçülen ortalama genlikler beyaz gürültü kuramıyla uyuşur (0.136, 0.222) ve tarama tablosunda yüksek frekans genliğinin gürültü RMS'ine oranı iki gürültüde de yaklaşık 0.0070'dir (kuram: sqrt(pi/65536) = 0.0069). Bu yüzden spektrum görüntüsüne bakarak Gauss gürültüsünü tuz-biber gürültüsünden ayırmak zordur; fark piksel uzayındadır (tuz-biber seyrek ve aşırı değerlidir) ve gürültü giderme yönteminin seçimini belirler (bkz. bölüm 3'teki dürtü gürültüsü deneyi).
- %5 tuz-biber, bu denemede sigma = 20 Gauss'tan daha yıkıcıdır (RMS 32.2'ye karşı 19.7, PSNR 17.97'ye karşı 22.25 dB; düşük frekans enerji oranı %70.24'e karşı %84.78).
- Tuz-biber, düşük frekanslı orijinal bileşenleri de zayıflatır: en güçlü harmonik (1, 0) 49.62'den 47.05'e (yaklaşık %5) ve 1-8 bandı 0.95 katına düştü. Bozulan piksellerin orijinal değerini kaybetmesi nedeniyle orijinal bileşenin yaklaşık (1 - oran) katsayısıyla çarpılması beklenir; ölçümle uyumludur. Gauss gürültüsü sıfır ortalamalı ve orijinalden bağımsız olduğu için düşük frekansları değiştirmedi (49.43; bant oranı 1.00).
- Aynı nedenle tuz-biber gürültüsünün (gürültülü - orijinal) spektrumunda en düşük bant yükselir (0.328; diğer bantlar yaklaşık 0.22): bu fark, orijinalin kendisinden bir bileşen (-oran x orijinal, yaklaşık 0.05 x 4.53 = 0.23) de taşır, kareköklü toplam 0.32 ölçülenle uyumludur. Gauss gürültüsünde böyle bir yükselme yoktur.
- DC (ortalama parlaklık) tuz-biberde 153.746'dan 152.397'ye düştü; beklenen değişim 0.05 x 0.5 x ((0 - 153.7) + (255 - 153.7)) = yaklaşık -1.3'tür.

## 7. Görev 9: ters DFT ile orijinale dönüş ve benzerlik metrikleri

### 7.1 Yöntem

Gürültülü resimler, Görev 8'deki ikisidir (Gauss sigma = 20 ve tuz-biber oranı 0.05; Cat.jpg, 256x256). Gürültülü resmin DFT'sini alıp doğrudan ters DFT yapmak resmi aynen geri verir (en büyük fark 0, ZNCC değişmedi: Gauss 0.9429, tuz-biber 0.8563); yani gürültü bu yolla kalkmaz. Orijinale yaklaşmak için DFT ile ters DFT arasında frekans uzayında süzme uygulandı (Görev 6'daki frekans filtresinin resim hali), üç yöntemle:

- **Eşik:** harmonik genliği eşiğin altında kalan harmonikler atılır (8 eşik: 0.1 ... 4.0).
- **İdeal alçak geçiren:** sıfır frekanstan uzaklığı R'den büyük harmonikler atılır (R = 4 ... 64).
- **Gauss alçak geçiren:** her harmonik exp(-r² / (2 sigma²)) ile çarpılır (sigma = 4 ... 32).

Ters DFT sonucu 0-255'e sınırlanıp yuvarlanarak 8-bit resim yapıldı ve orijinalle karşılaştırıldı. Her yöntemde parametre, orijinalle ZNCC'si en yüksek olacak şekilde seçildi (gerçek uygulamada orijinal bilinmediği için bu seçim yapılamaz).

### 7.2 Geri dönüş sonuçları

Gauss gürültülü resim:

| Yöntem (en iyi parametre) | Korunan harmonik | ZNCC | SSIM | PSNR (dB) | MSE | MAE |
|---|---|---|---|---|---|---|
| Gürültülü resim (süzmesiz) | tümü | 0.9429 | 0.3882 | 22.25 | 387.4 | 15.75 |
| Eşik T = 0.5 | 1345 | 0.9827 | 0.7137 | 27.80 | 108.0 | 7.43 |
| İdeal alçak geçiren R = 48 | 7213 | 0.9876 | 0.7529 | 29.22 | 77.7 | 6.56 |
| Gauss alçak geçiren sigma = 32 | - | 0.9901 | 0.8298 | 30.14 | 63.0 | 5.40 |

Tuz-biber gürültülü resim:

| Yöntem (en iyi parametre) | Korunan harmonik | ZNCC | SSIM | PSNR (dB) | MSE | MAE |
|---|---|---|---|---|---|---|
| Gürültülü resim (süzmesiz) | tümü | 0.8563 | 0.3498 | 17.97 | 1036.8 | 6.59 |
| Eşik T = 0.8 | 631 | 0.9718 | 0.6546 | 25.52 | 182.6 | 9.76 |
| İdeal alçak geçiren R = 32 | 3209 | 0.9810 | 0.7121 | 27.20 | 123.9 | 8.04 |
| Gauss alçak geçiren sigma = 24 | - | 0.9843 | 0.7706 | 27.63 | 112.2 | 7.32 |

**Yorum.**

- Süzme olmadan geri dönüş yoktur; süzme ile Gauss gürültüde ZNCC 0.9429'dan 0.9901'e, SSIM 0.3882'den 0.8298'e, PSNR 22.25 dB'den 30.14 dB'ye çıktı. Tuz-biberde ZNCC 0.8563'ten 0.9843'e çıktı.
- En iyi eşik, gürültünün spektrum tabanının en büyük harmoniğine yakındır. Görev 8'de gürültü spektrumunun ortalamasını Gauss için 0.136, tuz-biber için 0.222 bulmuştuk; Rayleigh dağılımlı genliklerin 32768 hücredeki beklenen en büyüğü ortalamanın yaklaşık 3.6 katıdır (yaklaşık 0.5 ve 0.8) ve en iyi eşikler de 0.5 ve 0.8 çıktı. Eşik çok düşük olursa gürültü harmonikleri kalır (Gauss gürültüde T = 0.1'de ZNCC 0.9457), çok yüksek olursa resmin zayıf harmonikleri de atılır (T = 4.0'da 0.9069).
- İki gürültüde de sıralama aynıdır: Gauss alçak geçiren > ideal alçak geçiren > eşik (ZNCC, SSIM ve PSNR'ye göre), eşik yöntemi en az harmoniği (1345 ve 631) korumasına rağmen. Yumuşak kesimin keskin kesimden iyi çıkması beklenen bir sonuçtur; keskin kesim halka biçimli artefakt (ringing) üretebilir (ürettiği çıktı görüntülerinde kontrol edilmedi).
- Tuz-biber gürültüsü daha zordur: ZNCC 0.9843'e ulaşsa da MSE (112.2) Gauss gürültüdekinin (63.0) yaklaşık iki katıdır. Olası sebep, tuz-biber darbelerinin geniş bantlı ve yüksek genlikli olması; alçak geçiren filtre bunları silmek yerine çevreye yayar (görüntü kontrol edilmedi).
- Tuz-biberde MAE süzmeyle azalmaz, artar (6.59'dan 7.32-9.76'ya); MSE ise 1036.8'den 112.2-182.6'ya düşer. Süzülmemiş gürültülü resim piksellerin %95'inde tam doğrudur, süzülmüş resim ise her yerde küçük hata taşır.

### 7.3 Piksel değerlerinin yakınlığını ölçen metrikler

Bu çalışmada altı metrik uygulandı: MSE, MAE, PSNR, NCC, ZNCC ve SSIM. Davranışlarını görmek için, orijinalden bilinen bozulmalarla türetilmiş resimler orijinalle karşılaştırıldı:

| Resim | MSE | MAE | PSNR (dB) | NCC | ZNCC | SSIM |
|---|---|---|---|---|---|---|
| Orijinalin kendisi | 0.0 | 0.00 | - | 1.0000 | 1.0000 | 1.0000 |
| Kontrast 0.7x + parlaklık +60 | 475.0 | 17.53 | 21.36 | 0.9929 | 1.0000 | 0.9469 |
| Negatif (255 - resim) | 15303.9 | 105.25 | 6.28 | 0.6565 | -1.0000 | -0.2476 |
| Bulanıklık (sigma 1.5) | 38.9 | 3.67 | 32.23 | 0.9993 | 0.9939 | 0.9076 |
| 4 piksel sağa kayma | 447.7 | 12.99 | 21.62 | 0.9916 | 0.9286 | 0.5782 |
| Transpoz (90 derece döndürme) | 6096.4 | 60.33 | 10.28 | 0.8862 | 0.0284 | 0.2805 |
| Gauss gürültülü | 387.4 | 15.75 | 22.25 | 0.9928 | 0.9429 | 0.3882 |
| Tuz-biber gürültülü | 1036.8 | 6.59 | 17.97 | 0.9808 | 0.8563 | 0.3498 |
| Sabit resim (ortalama değer) | 3137.2 | 50.15 | 13.17 | 0.9396 | tanımsız | 0.4666 |

**Metriklerin tanımı, avantaj ve dezavantajları:**

| Metrik | Ne ölçer | Avantajları | Dezavantajları |
|---|---|---|---|
| **MSE** (SSD'nin ortalaması) | Piksel farklarının karesinin ortalaması; 0 en iyi, sınırsız | Basit; türevlenebilir (optimizasyonda yaygın); büyük hatalara çok duyarlı | Parlaklık ve kontrast kaymasını çok cezalandırır (ZNCC'si 1 olan resimde 475.0); algıya uygun değil; yapıyı görmez; ölçeği sınırsız, yorumlaması zor |
| **MAE** (SAD'nin ortalaması) | Mutlak farkların ortalaması; gri seviye cinsinden | Gri seviye biriminde doğrudan yorumlanır; aykırı değerlerden MSE kadar etkilenmez | Hatanın dağılımını ayırt etmez; MSE ile ters sıralama verebilir (tuz-biber MAE 6.59 < Gauss 15.75, MSE 1036.8 > 387.4); MSE ile aynı diğer zayıflıklar |
| **PSNR** | MSE'nin dB ölçeği, 10 log10(255² / MSE); büyük iyi | Standart ve yaygın; karşılaştırması kolay | MSE ile her zaman aynı sıralamayı verir, aynı zayıflıklara sahiptir; özdeş resimde sonsuz |
| **NCC** | Ortalama çıkarılmadan normalize çapraz korelasyon; [0, 1] (pozitif piksellerde) | Kazanç (çarpan) değişiminden etkilenmez; hesaplaması basit | Ofset (parlaklık kayması) ve DC baskındır: ilişkisiz resimlerde bile yüksek çıkar (transpoz 0.8862, sabit resim 0.9396, gürültülü resimler 0.98-0.99); ayırt edici gücü düşük |
| **ZNCC** | Ortalaması çıkarılmış normalize çapraz korelasyon (Pearson korelasyonu); [-1, 1] | Parlaklık ofsetine ve kontrast kazancına değişmez (kontrast 0.7x + 60: 1.0000), aydınlanma değişimine dayanıklıdır; ilişkisiz resimde sıfıra yakın (0.0284); negatif ilişkiyi işaretle gösterir (negatif: -1.0000); şablon eşlemede standarttır | Sabit (varyansı sıfır) resimde tanımsızdır; yalnızca doğrusal ilişkiyi ölçer, parlaklık ve kontrast değişimini hiç görmez; bulanıklığa karşı duyarsızdır (0.9939); tüm resim için tek değer verdiğinden yerel yapıyı görmez; konumsal kaymada düşer ama gürültüye SSIM kadar duyarlı değildir |
| **SSIM** | Yerel parlaklık, kontrast ve yapı benzerliği; 1 en iyi, negatif olabilir | İnsan algısına daha uygundur; gürültü, bulanıklık ve kaymaya duyarlıdır (gürültülü resimlerde 0.39 ve 0.35, ZNCC 0.94 ve 0.86, PSNR 22 ve 18 dB); yerel bir harita verir | Hesabı daha ağırdır; pencere boyutu ve sabitler gibi seçimler gerektirir (burada 7x7 düzgün pencere, standart Gauss ağırlıklı 11x11'den farklıdır); yorumlaması zordur; yapısal metrik de olsa belirli bozulmalara (örneğin ofset) az duyarlıdır |

Bu çalışmada hesaplanmayan, piksel değerlerinin yakınlığını ölçmekte kullanılabilecek diğer metrikler: RMSE ve normalize RMSE (MSE'nin ölçeklenmiş hâli), SSD ve SAD (toplam hata), karşılıklı bilgi (mutual information; farklı görüntüleme kiplerindeki resimler için), histogram tabanlı uzaklıklar (konum bilgisini yok sayar), gradyan tabanlı metrikler (örneğin GMSD) ve öğrenilmiş algısal metrikler (örneğin LPIPS; eğitilmiş model gerektirir).

**Yorum (tabloya dayanarak).**

- ZNCC, parlaklık ve kontrast değişimine karşı değişmezdir (1.0000), MSE ve PSNR ise aynı resmi kötü gösterir (475.0 ve 21.36 dB). Bu, aydınlanma değişiminin önemsiz olduğu durumlarda (örneğin şablon eşleme) bir avantaj, değişimin kendisinin önemli olduğu durumlarda (örneğin gürültü gidermenin parlaklığı koruyup korumadığı) bir dezavantajdır.
- NCC pozitif piksel değerlerinde DC bileşeni yüzünden her resim için yüksek çıkar; ZNCC aynı transpoz resim için 0.0284 verir. Bu yüzden resim benzerliğinde ortalaması çıkarılmamış NCC güvenilir bir ölçüt değildir.
- ZNCC ve PSNR bulanıklığa karşı duyarsızdır (0.9939 ve 32.23 dB), SSIM daha duyarlıdır (0.9076); gürültüde ise SSIM çok daha duyarlıdır (0.3882, ZNCC 0.9429).
- Piksel piksel karşılaştıran bütün metrikler konum kaymasına duyarlıdır (4 piksel kaymada MSE 447.7, ZNCC 0.9286, SSIM 0.5782); yani karşılaştırılan resimlerin hizalı olması gerekir. 4 piksel kaymanın PSNR'si (21.62 dB) Gauss gürültününkine (22.25 dB) benzerken SSIM'leri çok farklıdır (0.5782 ve 0.3882); aynı PSNR farklı bozulma türlerini gizleyebilir.
- Sabit resimde ZNCC tanımsızdır (varyans 0); MSE (3137.2) o durumda orijinalin varyansına eşittir. Pratikte paydaya küçük bir sayı eklemek gerekir.
- Tek bir metrik bütün bozulmaları doğru yansıtmaz; hangi metriğin seçileceği, neyin önemli olduğuna bağlıdır: yapı ve algı için SSIM, aydınlanmadan bağımsız eşleme için ZNCC, piksel doğruluğu için MSE/PSNR/MAE.

## 8. Görev 10: low-pass, high-pass, band-pass ve band-reject filtreleri

### 8.1 Yöntem

Resmin 2 boyutlu DFT'sindeki her harmonik, sıfır frekanstan uzaklığı r = sqrt(fu² + fv²) ("döngü/resim"; Nyquist sınırı 128) ile belirlenen bir kazançla (maske) çarpılır, sonra ters DFT yapılır. Her filtre iki türlü kuruldu:

| Filtre | İdeal (keskin kesim) | Gauss (yumuşak geçiş) |
|---|---|---|
| Low-pass (alçak geçiren) | r <= R ise 1, değilse 0 | exp(-r² / (2 s²)), s = R / 1.1774 (r = R'de kazanç 0.5) |
| High-pass (yüksek geçiren) | r >= R ise 1, değilse 0 | 1 - low-pass |
| Band-pass (bant geçiren) | R1 <= r <= R2 ise 1, değilse 0 | exp(-((r² - r0²) / (r Wb))²), r0 = (R1 + R2) / 2, Wb = R2 - R1 |
| Band-reject (bant durduran) | 1 - band-pass | 1 - band-pass |

Kesimler: low-pass ve high-pass R = 24; band-pass ve band-reject bandı 20-50. Sıfır frekans (DC) resmin ortalama parlaklığıdır: low-pass ve band-reject bunu korur, high-pass ve band-pass sıfırlar (bu iki filtrenin çıktısı sıfır ortalamalıdır; görüntülerde gösterim için 128 ofset eklendi). Kullanılan resim: halka desenli bir zone plate (`zoneplate.png.jpeg`, 426x418, 256x256'ya küçültüldü); ortalama parlaklığı 120.52, standart sapması 93.22. Gürültülü halleri: Gauss (sigma = 20) ve tuz-biber (%5). Filtreler hem temiz resme hem iki gürültülü resme uygulandı.

### 8.2 Temiz resimde filtrelerin etkisi

| Filtre | Ortalama parlaklık (kırpmasız) | Min | Max | Kenar gücü / orijinal | Enerji oranı (çıkış varyansı / giriş) | ZNCC (orijinalle)* |
|---|---|---|---|---|---|---|
| Orijinal | 120.52 | - | - | 1.00 (kenar gücü 46.855) | 1.000 | - |
| İdeal low-pass | 120.52 | -46.0 | 311.7 | 0.38 | 0.672 | 0.824 |
| İdeal high-pass | 0.00 | -184.7 | 139.8 | 0.94 | 0.332 | - |
| İdeal band-pass | 0.00 | -132.6 | 133.8 | 0.97 | 0.382 | - |
| İdeal band-reject | 120.52 | -53.8 | 333.6 | 0.40 | 0.618 | 0.794 |
| Gauss low-pass | 120.52 | -0.0 | 255.0 | 0.36 | 0.560 | 0.898 |
| Gauss high-pass | 0.00 | -156.4 | 104.9 | 0.74 | 0.216 | - |
| Gauss band-pass | 0.00 | -107.6 | 108.1 | 0.75 | 0.241 | - |
| Gauss band-reject | 120.52 | -46.9 | 284.2 | 0.41 | 0.654 | 0.878 |

\* ZNCC yalnızca ortalaması korunan filtreler (low-pass, band-reject) için verilmiştir; kırpılmış 8-bit çıktı üzerinde hesaplanmıştır. Kenar gücü, komşu piksel farklarının mutlak değerlerinin ortalamasıdır (taşma ve halka salınımlarını da içerir).

### 8.3 Filtrelerin gürültüye etkisi

Geçen gürültü = filtre(gürültülü resim) - filtre(temiz resim) (filtreler doğrusal olduğu için filtre(gürültü)'ye eşittir). Çıkış SNR'ı = filtrelenmiş temiz resmin standart sapması / geçen gürültü RMS'i (dB). Giriş SNR'ı: Gauss gürültüde 14.29 dB (gürültü RMS 17.99), tuz-biberde 8.11 dB (gürültü RMS 36.63).

Gauss gürültülü resim:

| Filtre | Geçen gürültü RMS | Kuram (beyaz gürültü) | Geçen oran | Çıkış SNR (dB) | SNR değişimi (dB) |
|---|---|---|---|---|---|
| İdeal low-pass | 4.98 | 2.98 | 0.277 | 23.72 | +9.43 |
| İdeal high-pass | 17.29 | 17.74 | 0.961 | 9.85 | -4.44 |
| İdeal band-pass | 5.72 | 5.71 | 0.318 | 20.05 | +5.76 |
| İdeal band-reject | 17.05 | 17.06 | 0.948 | 12.67 | -1.62 |
| Gauss low-pass | 4.52 | 2.54 | 0.251 | 23.77 | +9.48 |
| Gauss high-pass | 16.99 | 17.44 | 0.945 | 8.12 | -6.17 |
| Gauss band-pass | 4.62 | 4.66 | 0.257 | 19.91 | +5.62 |
| Gauss band-reject | 16.76 | 16.74 | 0.932 | 13.06 | -1.23 |

Tuz-biber gürültülü resim:

| Filtre | Geçen gürültü RMS | Kuram (beyaz gürültü) | Geçen oran | Çıkış SNR (dB) | SNR değişimi (dB) |
|---|---|---|---|---|---|
| İdeal low-pass | 7.22 | 6.06 | 0.197 | 20.49 | +12.37 |
| İdeal high-pass | 35.92 | 36.13 | 0.980 | 3.50 | -4.62 |
| İdeal band-pass | 11.87 | 11.63 | 0.324 | 13.71 | +5.60 |
| İdeal band-reject | 34.66 | 34.74 | 0.946 | 6.51 | -1.61 |
| Gauss low-pass | 6.24 | 5.17 | 0.170 | 20.97 | +12.85 |
| Gauss high-pass | 35.29 | 35.52 | 0.963 | 1.77 | -6.34 |
| Gauss band-pass | 9.66 | 9.50 | 0.264 | 13.51 | +5.40 |
| Gauss band-reject | 34.03 | 34.10 | 0.929 | 6.91 | -1.20 |

### 8.4 Yorum

**Genel parlaklık.** Ortalama parlaklık, DC harmoniğinin kazancına eşittir. Low-pass ve band-reject DC'yi aynen korudu (ortalama 120.52), high-pass ve band-pass sıfırladı (0.00); yani ikisinin çıktısı sıfır ortalamalıdır ve gösterim için ofset gerekir. Gürültü eklemek ortalama parlaklığı belirgin biçimde değiştirmedi (Görev 8'de gösterildiği gibi).

**Kenar geçişleri.**

- Low-pass kenarları bulanıklaştırır: kenar gücü orijinalin 0.36-0.38'ine düştü (ince ayrıntı, yani yüksek frekanslar silindi). Band-reject de bu görüntüde kenar gücünü 0.40-0.41'e düşürdü; çünkü silinen 20-50 bandı bu resmin halka kenarlarının enerjisinin büyük kısmını taşıyor.
- High-pass ve band-pass kenar bilgisini korur ya da vurgular (ideal için kenar gücü oranı 0.94 ve 0.97; Gauss için 0.74 ve 0.75), ama düşük frekanslı geniş alanları (ortalama parlaklığı ve yumuşak geçişleri) kaldırır; sonuç, kenarların etrafında salınan, sıfır ortalamalı bir resimdir. İdeal filtrelerin kenar gücü, halka salınımlarının eklediği gradyanları da içerdiği için şişmiş olabilir.
- Keskin kesimli ideal filtreler 0-255 dışına taşar: ideal low-pass -46.0 ile 311.7, ideal band-reject -53.8 ile 333.6 arasına. Bu, kenarlar çevresinde halka biçimli salınımlar (ringing, Gibbs etkisi) demektir. Gauss low-pass ise tam [-0.0, 255.0] aralığında kaldı: Gauss low-pass, uzay uzayında pozitif ağırlıklı bir ortalamaya (evrişim) denk geldiği için çıktısı girdinin değer aralığını aşamaz; ideal low-pass'in uzaydaki karşılığı ise negatif loblu (sinc benzeri) olduğu için taşar. Gauss band-reject taşmayı azalttı ama sıfırlamadı (-46.9 ile 284.2). Taşmanın görsel etkisi bu raporda görüntüler üzerinden doğrulanmadı.

**Gürültü.**

- Bir filtrenin SNR'ı artırıp artırmadığı, sinyal enerjisinin ve gürültü enerjisinin hangi oranda korunduğuna bağlıdır: SNR değişimi = 10 log10(korunan sinyal enerjisi oranı / korunan gürültü enerjisi oranı). Ölçümler bu bağıntıyla ±0.02 dB içinde uyuşuyor: örneğin ideal low-pass için sinyal enerjisi oranı 0.672, geçen gürültü oranı 0.277 olduğundan 10 log10(0.672 / 0.277²) = 9.42 dB (ölçülen 9.43); high-pass için 10 log10(0.332 / 0.961²) = -4.44 dB (ölçülen -4.44); band-pass için 5.77 dB (5.76); band-reject için -1.63 dB (-1.62). Gauss filtreler için hesap 9.49, -6.16, 5.62 ve -1.23 dB, ölçülen 9.48, -6.17, 5.62 ve -1.23 dB.
- Low-pass gürültüyü en çok bastıran filtredir ve SNR'ı artırır (Gauss gürültüde +9.4 dB; tuz-biberde +12.4 dB); bedeli, ince ayrıntının ve kenarların silinmesidir. Gauss low-pass biraz daha iyi bastırır (geçen oran 0.251, ideal 0.277) ve taşma yapmaz.
- High-pass gürültünün neredeyse hepsini geçirir (%94-98) ve sinyalin yalnızca bir kısmını (enerjinin %22-33'ünü) korur, bu yüzden SNR düşer (-4.4 ile -6.3 dB).
- Band-pass yalnızca bandın içindeki gürültüyü geçirir (geçen oran 0.25-0.32); bu resimde sinyal enerjisinin %24-38'i bu banttadır, bu yüzden SNR artar (+5.4 ile +5.8 dB). Parlaklık sıfırlanır.
- Band-reject, gürültünün %93-95'ini ve sinyal enerjisinin %62-65'ini korur; SNR az azalır (-1.2 ile -1.6 dB).
- Tuz-biber gürültüsünde low-pass daha çok kazanç sağladı (+12.4 dB; Gauss gürültüde +9.4 dB), çünkü bu gürültünün enerjisi (beyaza daha yakın olduğu için) geniş bir banda yayılıyor ve düşük frekanslı küçük bir kısmı geçiyor (geçen oran 0.197, Gauss gürültüde 0.277). Bu hesap enerjiye dayanır; filtre çıktısındaki yayılmış tuz-biber noktaları (bulanık lekeler) bu rakama yansımaz ve görüntüler üzerinde kontrol edilmedi.
- Ölçülen geçen gürültü, low-pass hariç tüm filtrelerde beyaz gürültü kuramıyla uyuşur (örneğin ideal band-pass 5.72, kuram 5.71). Low-pass'te ölçülen kuramdan fazladır (4.98'e karşı 2.98). Olası sebep, gürültülü resmin 0 ve 255'te kırpılması (resimde siyah ve beyaz bölgeler var): kırpma, gürültüyü beyaz olmaktan çıkarıp düşük frekansa yaklaştırır. Bu doğrulanmadı.

**Filtre seçimi resmin spektrumuna bağlıdır.** İlk denemede (resim dosyası bulunamadığı için üretilen sentetik zone plate, enerjisi daha geniş bir banda yayılan bir desen) low-pass SNR'ı +6.51 dB (Gauss gürültüde, ideal), band-pass ise -0.19 dB değiştirmişti; sinyal enerji oranları low-pass için 0.350, band-pass için 0.093'tü. Gerçek resimde bu oranlar 0.672 ve 0.382 olduğu için kazançlar +9.43 ve +5.76 dB çıktı. Yani aynı filtre, farklı resimlerde farklı sonuç verir; sinyalin enerjisi hangi frekans bölgesinde toplanıyorsa o bölgeyi tutan filtre SNR'ı artırır.

**Özet.**

| Filtre | Gürültü | Kenar ve ayrıntı | Genel parlaklık |
|---|---|---|---|
| Low-pass | En çok azaltır | Bulanıklaştırır (ideal olanda halka/taşma) | Korur |
| High-pass | Neredeyse hiç azaltmaz (göreceli olarak artırır) | Vurgular, geniş alanları siler | Sıfırlar |
| Band-pass | Banttaki kısmı geçirir | Banttaki yapıyı ayırır | Sıfırlar |
| Band-reject | Neredeyse hiç azaltmaz | Banttaki yapıyı siler (burada kenarları da zayıflattı) | Korur |

## 9. Sınırlılıklar

- Her deney tek bir gürültü gerçekleşmesiyle (sabit tohum) yapıldı; tekrarlanan denemelerin ortalaması alınmadı.
- Gürültü düzgün dağılımlıdır (Gauss dağılımlı gürültü denenmedi); dürtü gürültüsü deneyi ödevin istediği bir madde değildir, yorumu zenginleştirmek için eklenmiştir.
- Tek frekans (8 Hz) ve tek örnekleme (256 Hz, 256 örnek) denendi. Frekans bin sınırlarına denk gelmeseydi (örneğin 8.5 Hz) spektral sızıntı olur ve temiz sinyalin spektrumu da tek tepe olmazdı; bu denenmedi.
- DFT doğrudan tanımdan (O(N²)) hesaplandı, hızlı Fourier dönüşümü (FFT) kullanılmadı.
- Grafikler (cikti_sinyal/*.png ve cikti_frekans/*.png) rapora yorumlanmadan eklenmiştir; yorumlar sayısal çıktıya dayanır.
- Görev 5 ve 6'da eşik elle seçildi (çıktıdaki en büyük gürültü harmoniğine ve en zayıf sinyal harmoniğine bakılarak). Gerçek bir uygulamada gürültü tabanının ayrıca tahmin edilmesi gerekir; gürültü ile zayıf sinyal bileşeninin genlikleri birbirine yakınsa genlik eşiği ikisini ayıramaz.
- 30 ve 50 Hz, örnekleme ve vektör uzunluğu yüzünden tam harmoniklere (k = 30 ve k = 50) denk geldi, bu yüzden spektral sızıntı yoktur. Bin sınırına denk gelmeyen frekanslar için bu sonuçlar değişir; denenmedi.
- Sinyaller sabit genlikli ve durağandır; zamanla değişen frekanslar için tek bir DFT ve tek eşik yeterli olmaz. Bu denenmedi.
- Görev 7'de gerçek resim olarak tek bir fotoğraf kullanıldı ve DFT için 256x256'ya küçültüldü (küçültme ince ayrıntıyı kaybettirir). Sentetik desenlerin frekansları tam harmoniklere denk geldi, yani spektral sızıntı yoktur. Faz matrisi fonksiyondan döndürülüyor ama bu testlerde incelenmedi.
- 2 boyutlu DFT doğrudan tanımdan hesaplandı (FFT yok); resimlerin spektrum görüntüleri (cikti_dft2d/*.png) yorumlanmadan eklenmiştir.
- Görev 8'de tek bir resim, tek bir tohum ve her gürültüden üç kuvvet denendi; gürültü fonksiyonları yalnızca tek kanallı 8-bit resimleri destekler. Gauss gürültüsünde 0-255 sınırlaması (kırpma) gürültü RMS'ini ve DC'yi çok az etkiler. Spektrum görüntüsü (cikti_goruntu_gurultu/8_karsilastirma.png) yorumlanmadan eklenmiştir; yorumlar sayısal çıktıya dayanır.
- Görev 9'da süzme parametreleri orijinal resimle ZNCC'yi en büyükleyecek şekilde seçildi (orijinal bilinmediği için gerçek bir uygulamada bu yapılamaz). Gauss gürültüde Gauss alçak geçirenin en iyi sigma'sı (32) taranan aralığın en büyük değeridir (sigma = 24'te ZNCC 0.9877, 32'de 0.9901); daha büyük sigma denenmedi, gerçek optimum daha yüksek olabilir.
- Metrikler tek bir resim ve tek bir gürültü gerçekleşmesi üzerinde, bütün resim için tek değer olarak hesaplandı (ZNCC ve NCC yerel değil küresel; SSIM 7x7 düzgün pencereyle). Bozulmalar (kontrast, negatif, bulanıklık, kayma, döndürme) tek bir şiddette denendi; farklı şiddetlerde sıralama değişebilir.
- Görev 9'daki geri dönüş görüntüsü (cikti_geri_donus/9_geri_donus.png) yorumlanmadan eklenmiştir; halka biçimli artefakt gibi görsel etkiler bu raporda doğrulanmadı.
- Görev 10'da tek bir resim, tek bir kesim değeri (R = 24; bant 20-50) ve her gürültüden tek kuvvet kullanıldı; filtre parametreleri taranmadı. Resim 426x418'den 256x256'ya küçültüldü (en-boy oranı yaklaşık %2 bozuldu); resimde başlık ya da kenar boşluğu kalıp kalmadığı doğrulanmadı, kalmışsa ortalama parlaklık ve kenar gücü etkilenmiş olabilir. Filtre çıktı görüntüleri (cikti_filtreler/*.png) rapora yorumlanmadan eklenmiştir; yorumlar sayısal çıktıya dayanır.
- Gürültünün beyaz olduğu varsayımı low-pass'te tutmadı (ölçülen geçen gürültü kuramdan büyük); kırpma etkisi doğrulanmadı. İdeal filtrelerin kenar gücü ölçümü halka salınımlarını da içerir.
