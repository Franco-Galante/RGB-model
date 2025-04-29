# Information Retrieval in the Age of Generative AI: The RGB Model

This repository provides the code for Garetto, M., Cornacchia A., Galante, F., Leonardi, E., Nordio, A., & Tarable, A. (2025). Information Retrieval in the Age of Generative AI: The RGB Model. - Accepted at SIGIR 2025.

Our work presents a new quantitative model that explains the intricate interactions among different entities involved in creating, indexing, and distributing information. It also considers the role of users' reliance on generative AI tools on the spread of inaccurate information across digital platforms.


## Table of Contents

- [Abstract](#abstract)
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Usage](#usage)


## Abstract

The advent of Large Language Models (LLMs) and generative AI is fundamentally transforming information retrieval and processing on the Internet, bringing both great potential and significant concerns regarding content authenticity and reliability. This paper presents a novel quantitative approach to shed light on the complex information dynamics arising from the growing use of generative AI tools. Despite their significant impact on the digital ecosystem, these dynamics remain largely uncharted and poorly understood.
We propose a stochastic model to characterize the generation, indexing, and dissemination of information in response to new topics. This scenario particularly challenges current LLMs, which often rely on real-time Retrieval-Augmented Generation (RAG) techniques to overcome their static knowledge limitations. Our findings suggest that the rapid pace of generative AI adoption, combined with increasing user reliance, can outpace human verification, escalating the risk of inaccurate information proliferation across digital resources.
An in-depth analysis of Stack Exchange data confirms that high-quality answers inevitably require substantial time and human effort to emerge. This underscores the considerable risks associated with generating persuasive text in response to new questions and highlights the critical need for responsible development and deployment of future generative AI tools.


## Overview

The repository consists of two folders: 
- `model_simulation` which contains the C code to run a discrete-events simulator that describes the behavior of our RGB model. This allows to obtain Figures 4-7 in the paper (the plotting script are not provided but can be immediately obtained with `gnuplot`).
- `stack_exchange_analysis` which contains the scalable `Dask`-based data analysis of (some) Stack Exchange datasets. In particular, this allows to obtain Figures 8 and 9 of the paper. The analyzed domains are that of computer science (Stack Overflow) and mathematics (Math Stack Exchange).


## Prerequisites

For the `model_simulation` part, make sure you have the following:

- C compiler: You need a C compiler, in particular we will use `gcc`.
- GNU Scientific Library (GSL): Ensure GSL is installed on your system.

For the `stack_exchange_analysis`, which relies on Jupyter Lab to provide an easy interface to the `Dask` library, make sure you have the following:

- Python: You need Python 3.10 to perform the analysis.
- Stack Exchange dataset: Which is publicly available at https://archive.org/details/stackexchange thanks to the `Internet Archive` project (https://archive.org/).
- (Optional) HPC cluster with SLURM scheduler: We provide instructions to run the code locally on a Windows machine, which is suitable for most Stack Exchange datasets. Additionally, we provide guidance for executing the code on an HPC cluster (with hundreds of computational nodes), leveraging the parallelism of the `Dask` library, recommended for large communities like Stack Overflow.


## Usage

For detailed installation instructions, refer to the README files within each subfolder.


## Acknowledgments

Computational resources for this project have been provided by hpc@polito, which is a project of Academic Computing within the Department of Control and Computer Engineering at the Politecnico di Torino (http://hpc.polito.it).