# VASP Elastic Tensor Workflow

## 📌 Description

This repository provides scripts for performing **elastic tensor calculations using VASP**. It includes:

- Generation of strained structures for elastic property evaluation
- A VASP job submission script
- A parser for extracting elastic constants from output files

The scripts are simple and modular, intended for academic use in materials science projects focused on mechanical properties.

## 📂 Files

- `create_task.py` – Generate input files for different strain configurations
- `vasp_qsub_elast` – Submission script for running VASP jobs on a computing cluster
- `Read_OUT.py` – Parse VASP output files to extract the elastic tensor

## 🚀 Usage

1. Run `create_task.py` to generate strain-deformed structures.
2. Submit jobs using `vasp_qsub_elast`.
3. After VASP calculations, run `Read_OUT.py` to extract the elastic constants.

## 📄 License

For academic research use. Provided without warranty.
