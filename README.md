# Verified MATLAB Code for [Global stability of Wright-type equations with negative Schwarzian]

This repository contains MATLAB scripts developed to support the results presented in our mathematical paper. The computations use **rigorous numerics** 
to confirm theoretical findings and are implemented using **[INTLAB](https://www.ti3.tu-harburg.de/rump/intlab/)**, a toolbox for verified numerical computations in MATLAB.

## 🧮 Purpose

The provided code rigorously verifies the conclusions of various lemmas and theorems using interval arithmetic and validated numerics. 
These computations are essential for confirming key analytical results in our paper.

## 📁 File Organization

- All files were written in MATLAB and make use of the INTLAB toolbox.
- File names correspond either to:
  - The **notation of specific mathematical functions** used in the paper, or
  - The **lemmas and theorems** whose conclusions are verified by the code.
- This naming convention allows readers to easily match each file to the relevant section of the paper.

## 🛠 Requirements

- MATLAB (tested with version XX.X or later)
- [INTLAB toolbox](https://www.ti3.tu-harburg.de/rump/intlab/) (make sure it's properly installed and initialized)

## 📜 License
This code is made publicly available for scientific and educational use.
You may reuse or modify it in line with the accompanying license (if specified).

## 📖 Citation
If you use this code in your research, please cite our paper:

To start using INTLAB in your MATLAB session:
```matlab
startintlab


▶️ How to Run
Open MATLAB.

Ensure INTLAB is installed and initialized.

Navigate to the directory containing the files.

Run the desired .m script that corresponds to the lemma, theorem, or function you wish to verify.
