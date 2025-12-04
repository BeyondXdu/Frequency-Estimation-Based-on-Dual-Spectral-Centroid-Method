## **README — Time–Frequency Peak Detection (MATLAB)**

### **1. Software Purpose and Main Functions**

This MATLAB package implements a robust time–frequency peak detection algorithm based on a two-stage strategy:

1. iterative peak removal to obtain an unbiased noise estimate, and
2. a global Neyman–Pearson (NP) threshold for final detection.
   It includes a synthetic example consisting of two crossing chirps (100 Hz sampling, 14 s), consistent with commonly used benchmark signals.

Main functions:

* Generate crossing-chirp test signal
* Compute STFT (700 × 1400 configuration supported)
* Iterative peak removal with neighborhood dilation
* Global noise variance estimation (mean or median-based)
* Global NP-threshold peak detection
* Visualization of spectrograms, thresholds, and detected TF points

---

### **2. How to Run the Code**

**Requirements**

* MATLAB R2018a or later
* Signal Processing Toolbox
* (Optional) Image Processing Toolbox for `padarray`
* RAM ≥ 8 GB recommended

**Steps**

1. Place the code in your MATLAB working directory.
2. Run the main script:

```matlab
two_stage_global_thresholding_demo
```

Figures of the spectrogram, removed peaks, thresholds, and detection results will be displayed automatically.

---

### **3. Basic Usage Instructions**

* Adjust SNR, STFT parameters, median window length, number of iterations, and NP false-alarm rates (`alpha_init`, `alpha_final`) at the top of the script.
* To use your own signal, replace the signal-generation block with:

```matlab
x = your_signal(:);
```

Ensure the sampling rate matches `fs`.

---

## **4.Core principle (brief)**
Local noise estimates are biased near strong TF concentrations (e.g., chirp crossings).
The algorithm iteratively removes peaks, estimates noise from the remaining TF points, then applies a global NP threshold:
[
T = \hat{\sigma}^2 , (-\ln \alpha).
]
This yields stable detection across the entire TF plane.

---

### **5. Frequently Asked Questions**

**Q1. Missing `padarray`?**
Install Image Processing Toolbox or use `utils/padarray_fallback.m`.

**Q2. STFT size does not match 700×1400?**
Adjust `NFFT`, `winLen`, `hop`, or verify the signal length.

**Q3. Too many or too few detections?**
Tune `alpha_final`, `L_med`, `expand_rad`, or run more iterations.

**Q4. Crossing area still missed?**
Increase `expand_rad` or `maxIter`, or switch to median-based noise estimation.

---
