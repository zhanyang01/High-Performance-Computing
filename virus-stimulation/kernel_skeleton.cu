#include <cuda_runtime.h>

#include <iostream>
#include <vector>

#include "common.h"
#include "kseq/kseq.h"

// Comments written for easier understanding and debugging purposes and revision
// purposes, the main idea is to parallelise the matching of the sample and
// signature sequences using cuda as they are independent tasks

// The cuda kernel to find the match between the sample and the signature
__global__ void matchKernel(char* device_sampleSeq, char* device_signatureSeq,
                            char* device_sampleQual, int* device_samplesSize,
                            int* device_signaturesSize,
                            double* device_matchValue, int sampleNum,
                            int signatureNum, int* device_sampleIdx,
                            int* device_signatureIdx) {
    // Initialising the thread index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int maxThreads = sampleNum * signatureNum;
    // Comparing bounds
    if (maxThreads <= idx) {
        return;
    } else {
        // Initialising the sample and signature index
        int sampleIdx = idx / signatureNum;
        int signatureIdx = idx % signatureNum;

        // Initialising the match value
        double matchValue = 0;
        int firstCharToMatch = -1;

        // Initialising the sample and signature size
        int sampleSize = device_samplesSize[sampleIdx];
        int signatureSize = device_signaturesSize[signatureIdx];

        // Initialising the sample and signature sequence
        int sampleStartingPosition = device_sampleIdx[sampleIdx];
        char* sampleSeq = &device_sampleSeq[sampleStartingPosition];

        int signatureStartingPosition = device_signatureIdx[signatureIdx];
        char* signatureSeq = &device_signatureSeq[signatureStartingPosition];

        // Initialising the sample quality
        char* sampleQual = &device_sampleQual[sampleStartingPosition];

        // Finding the match(sliding technique)
        // Iterating through the sample
        for (int i = 0; i <= sampleSize - signatureSize; i++) {
            bool isMatch = true;
            // Iterating through the signature
            for (int j = 0; j < signatureSize; j++) {
                // If the character does not match or no wildcards, break the
                // loop
                char sampleChar = sampleSeq[i + j];
                char signatureChar = signatureSeq[j];
                if (sampleChar != signatureChar && sampleChar != 'N' &&
                    signatureChar != 'N') {
                    isMatch = false;
                    break;
                }
            }
            // If match found, store the first character that matches, which is
            // the current index in the sample string
            if (isMatch) {
                firstCharToMatch = i;
                break;
            }
        }
        // If match found, calculate the match value, -1 means didnt match so
        // skip
        if (firstCharToMatch != -1) {
            for (int i = 0; i < signatureSize; i++) {
                char currentQualityAsciiChar = sampleQual[firstCharToMatch + i];
                matchValue += static_cast<double>(currentQualityAsciiChar) - 33;
            }
            matchValue = matchValue / signatureSize;
        }
        // Storing the match value
        device_matchValue[sampleIdx * signatureNum + signatureIdx] = matchValue;
    }
}

// The function to run the cuda matcher
void runMatcher(const std::vector<klibpp::KSeq>& samples,
                const std::vector<klibpp::KSeq>& signatures,
                std::vector<MatchResult>& matches) {
    // 1. Flatten the array of structs in sample and signatures by
    // allocation host memory and copying the data from the structs to the
    // host memory

    // Initialising length of array and maximum size of sample and signature
    int sampleNum = samples.size();
    int signatureNum = signatures.size();
    int numberOfPairs = sampleNum * signatureNum;
    int totalSampleSize = 0;
    int totalSignatureSize = 0;

    // Finding the max size of the sample and signature
    for (int i = 0; i < sampleNum; i++) {
        totalSampleSize += samples[i].seq.size();
    }
    for (int i = 0; i < signatureNum; i++) {
        totalSignatureSize += signatures[i].seq.size();
    }

    // Allocate host memory
    char* host_sampleSeq = (char*)malloc(totalSampleSize * sizeof(char));
    char* host_signatureSeq = (char*)malloc(totalSignatureSize * sizeof(char));
    char* host_sampleQual = (char*)malloc(totalSampleSize * sizeof(char));
    int* host_sampleSize = (int*)malloc(sampleNum * sizeof(int));
    int* host_signatureSize = (int*)malloc(signatureNum * sizeof(int));
    int* host_sampleIdx = (int*)malloc(sampleNum * sizeof(int));
    int* host_signatureIdx = (int*)malloc(signatureNum * sizeof(int));
    double* host_matchValue =
        (double*)malloc(sampleNum * signatureNum * sizeof(double));

    // Copy the data from the structs to the host memory
    int sampleStartIndexes = 0;
    for (int i = 0; i < sampleNum; i++) {
        const char* sampleSeq = samples[i].seq.c_str();
        const char* sampleQual = samples[i].qual.c_str();
        int sampleSize = samples[i].seq.size();
        host_sampleSize[i] = sampleSize;
        host_sampleIdx[i] = sampleStartIndexes;
        for (int j = 0; j < sampleSize; j++) {
            int currentSampIndex = sampleStartIndexes + j;
            host_sampleSeq[currentSampIndex] = sampleSeq[j];
            host_sampleQual[currentSampIndex] = sampleQual[j];
        }
        sampleStartIndexes += sampleSize;
    }

    int signatureStartIndexes = 0;
    for (int i = 0; i < signatureNum; i++) {
        const char* signatureSeq = signatures[i].seq.c_str();
        int signatureSize = signatures[i].seq.size();
        host_signatureSize[i] = signatureSize;
        host_signatureIdx[i] = signatureStartIndexes;
        for (int j = 0; j < signatureSize; j++) {
            int currentSigIndex = signatureStartIndexes + j;
            host_signatureSeq[currentSigIndex] = signatureSeq[j];
        }
        signatureStartIndexes += signatureSize;
    }

    // 2. Then allocate device memory and copy the data from the host memory
    // to the device memory

    // Initialising device memory variable
    char* device_sampleSeq;
    char* device_signatureSeq;
    char* device_sampleQual;
    int* device_samplesSize;
    int* device_signaturesSize;
    int* device_sampleIdx;
    int* device_signatureIdx;
    double* device_matchValue;

    // Allocate device memory using cudaMalloc
    cudaMalloc(&device_sampleSeq, totalSampleSize * sizeof(char));
    cudaMalloc(&device_signatureSeq, totalSignatureSize * sizeof(char));
    cudaMalloc(&device_sampleQual, totalSampleSize * sizeof(char));
    cudaMalloc(&device_samplesSize, sampleNum * sizeof(int));
    cudaMalloc(&device_signaturesSize, signatureNum * sizeof(int));
    cudaMalloc(&device_sampleIdx, sampleNum * sizeof(int));
    cudaMalloc(&device_signatureIdx, signatureNum * sizeof(int));
    cudaMalloc(&device_matchValue, numberOfPairs * sizeof(double));

    // Copy the data from the host memory to the device memory using
    // cudaMemcpy
    cudaMemcpy(device_sampleSeq, host_sampleSeq, totalSampleSize * sizeof(char),
               cudaMemcpyHostToDevice);
    cudaMemcpy(device_signatureSeq, host_signatureSeq,
               totalSignatureSize * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(device_sampleQual, host_sampleQual,
               totalSampleSize * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(device_samplesSize, host_sampleSize, sampleNum * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(device_signaturesSize, host_signatureSize,
               signatureNum * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(device_sampleIdx, host_sampleIdx, sampleNum * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(device_signatureIdx, host_signatureIdx,
               signatureNum * sizeof(int), cudaMemcpyHostToDevice);

    // 3. Then run the kernel to find the match
    // Initialising the block size and grid size
    int blockSize = 256;
    int gridSize = ((sampleNum * signatureNum) + blockSize - 1) / blockSize;
    matchKernel<<<gridSize, blockSize>>>(
        device_sampleSeq, device_signatureSeq, device_sampleQual,
        device_samplesSize, device_signaturesSize, device_matchValue, sampleNum,
        signatureNum, device_sampleIdx, device_signatureIdx);

    // 4. Then copy the data from the device memory to the host memory
    cudaMemcpy(host_matchValue, device_matchValue,
               sampleNum * signatureNum * sizeof(double),
               cudaMemcpyDeviceToHost);
    // printing out all the match values
    // for (int i = 0; i < sampleNum; i++) {
    //     for (int j = 0; j < signatureNum; j++) {
    //         printf("Match value: %f\n", host_matchValue[i * signatureNum
    //         + j]);
    // 5. Then copy the data from the host memory to the vector of
    // MatchResult
    for (int i = 0; i < sampleNum; i++) {
        for (int j = 0; j < signatureNum; j++) {
            if (host_matchValue[i * signatureNum + j] != 0) {
                std::string currentSampleName = samples[i].name;
                std::string currentSignatureName = signatures[j].name;
                MatchResult matchResult = {
                    currentSampleName, currentSignatureName,
                    host_matchValue[i * signatureNum + j]};
                matches.push_back(matchResult);
            }
        }
    }
    // 6. Finally, free the memory
    free(host_sampleSeq);
    free(host_signatureSeq);
    free(host_sampleQual);
    free(host_sampleSize);
    free(host_signatureSize);
    free(host_sampleIdx);
    free(host_signatureIdx);
    free(host_matchValue);
    cudaFree(device_sampleSeq);
    cudaFree(device_signatureSeq);
    cudaFree(device_sampleQual);
    cudaFree(device_samplesSize);
    cudaFree(device_signaturesSize);
    cudaFree(device_sampleIdx);
    cudaFree(device_signatureIdx);
    cudaFree(device_matchValue);
}
