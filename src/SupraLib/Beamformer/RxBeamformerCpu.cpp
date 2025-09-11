// // ================================================================================================
// // 
// // If not explicitly stated: Copyright (C) 2016, all rights reserved,
// //      Rüdiger Göbl 
// //		Email r.goebl@tum.de
// //      Chair for Computer Aided Medical Procedures
// //      Technische Universität München
// //      Boltzmannstr. 3, 85748 Garching b. München, Germany
// // 
// // ================================================================================================
// #include "RxBeamformerCuda.h"
// #include "USImage.h"
// #include "USRawData.h"
// #include "RxSampleBeamformerDelayAndSum.h"
// #include "RxSampleBeamformerDelayAndStdDev.h"
// #include "RxSampleBeamformerCoherenceFactorDelayAndSum.h"
// #include "RxSampleBeamformerTestSignal.h"
// #include "RxBeamformerCommon.h"
// #include "utilities/cudaUtility.h"
// #include <thread>
// #include <vector>
// #include <functional>

// //TODO ALL ELEMENT/SCANLINE Y positons are actually Z! Change all variable names accordingly
// namespace supra
// {
// 	RxBeamformerCuda::RxBeamformerCuda(const RxBeamformerParameters & parameters)
// 		: m_windowFunction(nullptr)
// 	{
// 		m_lastSeenDt = 0;
// 		m_numRxScanlines = parameters.getNumRxScanlines();
// 		m_rxScanlineLayout = parameters.getRxScanlineLayout();

// 		m_is3D = (m_rxScanlineLayout.x > 1 && m_rxScanlineLayout.y > 1);
// 		m_speedOfSoundMMperS = parameters.getSpeedOfSoundMMperS();
// 		m_rxNumDepths = parameters.getRxNumDepths();

// 		// create and fill new buffers - using CPU memory instead of GPU
// 		m_pRxDepths = std::unique_ptr<Container<LocationType> >(
// 			new Container<LocationType>(LocationHost, 0, parameters.getRxDepths()));

// 		m_pRxScanlines = std::unique_ptr<Container<ScanlineRxParameters3D> >(
// 			new Container<ScanlineRxParameters3D>(LocationHost, 0, parameters.getRxScanlines()));

// 		m_pRxElementXs = std::unique_ptr<Container<LocationType> >(
// 			new Container<LocationType>(LocationHost, 0, parameters.getRxElementXs()));
// 		m_pRxElementYs = std::unique_ptr<Container<LocationType> >(
// 			new Container<LocationType>(LocationHost, 0, parameters.getRxElementYs()));

// 		if (parameters.getNonlinearElementToChannelMapping())
// 		{
// 			m_elementToChannelMap = std::unique_ptr<Container<int32_t> >(
// 				new Container<int32_t>(LocationHost, 0, parameters.getElementToChannelMap()));
// 		}
// 	}

// 	RxBeamformerCuda::~RxBeamformerCuda()
// 	{
// 	}

// 	void RxBeamformerCuda::convertToDtSpace(double dt, double speedOfSoundMMperS, size_t numTransducerElements) const
// 	{
// 		if (m_lastSeenDt != dt || m_speedOfSoundMMperS != speedOfSoundMMperS)
// 		{
// 			double oldFactor = 1;
// 			double oldFactorTime = 1;
// 			if (m_lastSeenDt != 0 && m_speedOfSoundMMperS != 0)
// 			{
// 				oldFactor = 1 / (m_speedOfSoundMMperS * m_lastSeenDt);
// 				oldFactorTime = 1 / m_lastSeenDt;
// 			}

// 			double factor = 1 / oldFactor / (speedOfSoundMMperS * dt);
// 			double factorTime = 1 / oldFactorTime / dt;

// 			// Already on host, no need for memory transfers
// 			for (size_t i = 0; i < m_numRxScanlines; i++)
// 			{
// 				ScanlineRxParameters3D p = m_pRxScanlines->get()[i];
// 				p.position = p.position*factor;
// 				for (size_t k = 0; k < std::extent<decltype(p.txWeights)>::value; k++)
// 				{
// 					p.txParameters[k].initialDelay *= factorTime;
// 				}
// 				p.maxElementDistance = p.maxElementDistance*factor;
// 				m_pRxScanlines->get()[i] = p;
// 			}

// 			for (size_t i = 0; i < m_rxNumDepths; i++)
// 			{
// 				m_pRxDepths->get()[i] = static_cast<LocationType>(m_pRxDepths->get()[i] * factor);
// 			}

// 			for (size_t i = 0; i < numTransducerElements; i++)
// 			{
// 				m_pRxElementXs->get()[i] = static_cast<LocationType>(m_pRxElementXs->get()[i] * factor);
// 				m_pRxElementYs->get()[i] = static_cast<LocationType>(m_pRxElementYs->get()[i] * factor);
// 			}

// 			m_lastSeenDt = dt;
// 			m_speedOfSoundMMperS = speedOfSoundMMperS;
// 		}
// 	}

// 	// CPU implementation of the 3D beamforming kernel
// 	template <class SampleBeamformer, 
// 		bool interpolateRFlines, 
// 		bool interpolateBetweenTransmits, 
// 		bool nonlinearElementToChannelMapping,
// 		typename RFType, 
// 		typename ResultType, 
// 		typename LocationType>
// 	void rxBeamformingDTSPACE3DCpu(
// 		uint32_t numTransducerElements,
// 		vec2T<uint32_t> elementLayout,
// 		uint32_t numReceivedChannels,
// 		uint32_t numTimesteps,
// 		const RFType* RF,
// 		uint32_t numTxScanlines,
// 		uint32_t numRxScanlines,
// 		const ScanlineRxParameters3D* scanlinesDT,
// 		uint32_t numDs,
// 		const LocationType* dsDT,
// 		const LocationType* x_elemsDT,
// 		const LocationType* z_elemsDT,
// 		LocationType speedOfSound,
// 		LocationType dt,
// 		uint32_t additionalOffset,
// 		LocationType F,
// 		const WindowFunction* windowFunction,
// 		ResultType* s,
// 		const int32_t* elementToChannelMap)
// 	{
// 		// Parallelize over depth and scanlines using OpenMP-style threading
// 		const size_t numThreads = std::thread::hardware_concurrency();
// 		std::vector<std::thread> threads;
		
// 		auto processChunk = [&](size_t startIdx, size_t endIdx) {
// 			for (size_t idx = startIdx; idx < endIdx; ++idx) {
// 				int r = idx / numRxScanlines;
// 				int scanlineIdx = idx % numRxScanlines;
				
// 				if (r >= numDs || scanlineIdx >= numRxScanlines) continue;

// 				LocationType d = dsDT[r];
// 				//TODO should this also depend on the angle?
// 				LocationType aDT = squ(computeAperture_D(F, d*dt*speedOfSound) / speedOfSound / dt);
// 				ScanlineRxParameters3D scanline = scanlinesDT[scanlineIdx];

// 				LocationType scanline_x = scanline.position.x;
// 				LocationType scanline_z = scanline.position.z;
// 				LocationType dirX = scanline.direction.x;
// 				LocationType dirY = scanline.direction.y;
// 				LocationType dirZ = scanline.direction.z;
// 				vec2f maxElementDistance = static_cast<vec2f>(scanline.maxElementDistance);
// 				vec2f invMaxElementDistance = vec2f{ 1.0f, 1.0f } / min(vec2f{ sqrt(aDT), sqrt(aDT) }, maxElementDistance);

// 				float sInterp = 0.0f;

// 				int highestWeightIndex;
// 				if (!interpolateBetweenTransmits)
// 				{
// 					highestWeightIndex = 0;
// 					float highestWeight = scanline.txWeights[0];
// 					for (int k = 1; k < std::extent<decltype(scanline.txWeights)>::value; k++)
// 					{
// 						if (scanline.txWeights[k] > highestWeight)
// 						{
// 							highestWeight = scanline.txWeights[k];
// 							highestWeightIndex = k;
// 						}
// 					}
// 				}

// 				// now iterate over all four txScanlines to interpolate beamformed scanlines from those transmits
// 				for (int k = (interpolateBetweenTransmits ? 0 : highestWeightIndex);
// 					(interpolateBetweenTransmits && k < std::extent<decltype(scanline.txWeights)>::value) ||
// 					(!interpolateBetweenTransmits && k == highestWeightIndex);
// 					k++)
// 				{
// 					if (scanline.txWeights[k] > 0.0)
// 					{
// 						ScanlineRxParameters3D::TransmitParameters txParams = scanline.txParameters[k];
// 						uint32_t txScanlineIdx = txParams.txScanlineIdx;
// 						if (txScanlineIdx >= numTxScanlines)
// 						{
// 							//ERROR!
// 							continue;
// 						}
// 						float sLocal = 0.0f;
						
// 						sLocal = SampleBeamformer::template sampleBeamform3D<interpolateRFlines, nonlinearElementToChannelMapping, RFType, float, LocationType>(
// 							txParams, RF, elementLayout, numReceivedChannels, numTimesteps,
// 							x_elemsDT, z_elemsDT, scanline_x, scanline_z, dirX, dirY, dirZ,
// 							aDT, d, invMaxElementDistance, speedOfSound, dt, additionalOffset, windowFunction, nullptr, elementToChannelMap);

// 						if (interpolateBetweenTransmits)
// 						{
// 							sInterp += static_cast<float>(scanline.txWeights[k])* sLocal;
// 						}
// 						else
// 						{
// 							sInterp += sLocal;
// 						}
// 					}
// 				}
// 				s[scanlineIdx + r * numRxScanlines] = clampCast<ResultType>(sInterp);
// 			}
// 		};

// 		size_t totalWork = numDs * numRxScanlines;
// 		size_t chunkSize = (totalWork + numThreads - 1) / numThreads;
		
// 		for (size_t t = 0; t < numThreads; ++t) {
// 			size_t start = t * chunkSize;
// 			size_t end = min(start + chunkSize, totalWork);
// 			if (start < end) {
// 				threads.emplace_back(processChunk, start, end);
// 			}
// 		}
		
// 		for (auto& thread : threads) {
// 			thread.join();
// 		}
// 	}

// 	// CPU implementation of the 2D beamforming kernel
// 	template <class SampleBeamformer, 
// 		bool interpolateRFlines, 
// 		bool interpolateBetweenTransmits,
// 		bool nonlinearElementToChannelMapping,
// 		typename RFType, 
// 		typename ResultType, 
// 		typename LocationType>
// 	void rxBeamformingDTSPACECpu(
// 		size_t numTransducerElements,
// 		size_t numReceivedChannels,
// 		size_t numTimesteps,
// 		const RFType* RF,
// 		size_t numTxScanlines,
// 		size_t numRxScanlines,
// 		const ScanlineRxParameters3D* scanlinesDT,
// 		size_t numDs,
// 		const LocationType* dsDT,
// 		const LocationType* x_elemsDT,
// 		LocationType speedOfSound,
// 		LocationType dt,
// 		uint32_t additionalOffset,
// 		LocationType F,
// 		const WindowFunction* windowFunction,
// 		ResultType* s,
// 		const int32_t* elementToChannelMap)
// 	{
// 		// Parallelize over depth and scanlines using threading
// 		const size_t numThreads = std::thread::hardware_concurrency();
// 		std::vector<std::thread> threads;
		
// 		auto processChunk = [&](size_t startIdx, size_t endIdx) {
// 			for (size_t idx = startIdx; idx < endIdx; ++idx) {
// 				int r = idx / numRxScanlines;
// 				int scanlineIdx = idx % numRxScanlines;
				
// 				if (r >= numDs || scanlineIdx >= numRxScanlines) continue;

// 				LocationType d = dsDT[r];
// 				//TODO should this also depend on the angle?
// 				LocationType aDT = computeAperture_D(F, d*dt*speedOfSound) / speedOfSound / dt;
// 				ScanlineRxParameters3D scanline = scanlinesDT[scanlineIdx];
// 				LocationType scanline_x = scanline.position.x;
// 				LocationType dirX = scanline.direction.x;
// 				LocationType dirY = scanline.direction.y;
// 				LocationType dirZ = scanline.direction.z;
// 				LocationType maxElementDistance = static_cast<LocationType>(scanline.maxElementDistance.x);
// 				LocationType invMaxElementDistance = 1 / min(aDT, maxElementDistance);

// 				float sInterp = 0.0f;

// 				int highestWeightIndex;
// 				if (!interpolateBetweenTransmits)
// 				{
// 					highestWeightIndex = 0;
// 					float highestWeight = scanline.txWeights[0];
// 					for (int k = 1; k < std::extent<decltype(scanline.txWeights)>::value; k++)
// 					{
// 						if (scanline.txWeights[k] > highestWeight)
// 						{
// 							highestWeight = scanline.txWeights[k];
// 							highestWeightIndex = k;
// 						}
// 					}
// 				}

// 				// now iterate over all four txScanlines to interpolate beamformed scanlines from those transmits
// 				for (int k = (interpolateBetweenTransmits ? 0 : highestWeightIndex);
// 					(interpolateBetweenTransmits && k < std::extent<decltype(scanline.txWeights)>::value) ||
// 					(!interpolateBetweenTransmits && k == highestWeightIndex);
// 					k++)
// 				{
// 					if (scanline.txWeights[k] > 0.0)
// 					{
// 						ScanlineRxParameters3D::TransmitParameters txParams = scanline.txParameters[k];
// 						uint32_t txScanlineIdx = txParams.txScanlineIdx;
// 						if (txScanlineIdx >= numTxScanlines)
// 						{
// 							//ERROR!
// 							continue;
// 						}

// 						float sLocal = 0.0f;
// 						sLocal = SampleBeamformer::template sampleBeamform2D<interpolateRFlines, nonlinearElementToChannelMapping, RFType, float, LocationType>(
// 							txParams, RF, numTransducerElements, numReceivedChannels, numTimesteps,
// 							x_elemsDT, scanline_x, dirX, dirY, dirZ,
// 							aDT, d, invMaxElementDistance, speedOfSound, dt, additionalOffset, windowFunction, elementToChannelMap);

// 						if (interpolateBetweenTransmits)
// 						{
// 							sInterp += static_cast<float>(scanline.txWeights[k])* sLocal;
// 						}
// 						else
// 						{
// 							sInterp += sLocal;
// 						}
// 					}
// 				}
// 				s[scanlineIdx + r * numRxScanlines] = clampCast<ResultType>(sInterp);
// 			}
// 		};

// 		size_t totalWork = numDs * numRxScanlines;
// 		size_t chunkSize = (totalWork + numThreads - 1) / numThreads;
		
// 		for (size_t t = 0; t < numThreads; ++t) {
// 			size_t start = t * chunkSize;
// 			size_t end = min(start + chunkSize, totalWork);
// 			if (start < end) {
// 				threads.emplace_back(processChunk, start, end);
// 			}
// 		}
		
// 		for (auto& thread : threads) {
// 			thread.join();
// 		}
// 	}

// 	template <class SampleBeamformer, typename RFType, typename ResultType, typename LocationType>
// 	void rxBeamformingDTspaceCpu3D(
// 		bool interpolateRFlines,
// 		bool interpolateBetweenTransmits,
// 		size_t numTransducerElements,
// 		vec2s elementLayout,
// 		size_t numReceivedChannels,
// 		size_t numTimesteps,
// 		const RFType* RF,
// 		size_t numTxScanlines,
// 		size_t numRxScanlines,
// 		const ScanlineRxParameters3D* scanlines,
// 		size_t numZs,
// 		const LocationType* zs,
// 		const LocationType* x_elems,
// 		const LocationType* y_elems,
// 		LocationType speedOfSound,
// 		LocationType dt,
// 		uint32_t additionalOffset,
// 		LocationType F,
// 		const WindowFunction* windowFunction,
// 		ResultType* s,
// 		const int32_t* elementToChannelMap)
// 	{
// 		if (interpolateRFlines)
// 		{
// 			if (interpolateBetweenTransmits)
// 			{
// 				if(elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, true, true, true>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, true, true, false>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 			else {
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, true, false, true>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, true, false, false>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 		}
// 		else {
// 			if (interpolateBetweenTransmits)
// 			{
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, false, true, true>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, false, true, false>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 			else {
// 				if (elementToChannelMap)
// 				{

// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, false, false, true>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACE3DCpu<SampleBeamformer, false, false, false>(
// 						(uint32_t)numTransducerElements, static_cast<vec2T<uint32_t>>(elementLayout),
// 						(uint32_t)numReceivedChannels, (uint32_t)numTimesteps, RF,
// 						(uint32_t)numTxScanlines, (uint32_t)numRxScanlines, scanlines,
// 						(uint32_t)numZs, zs, x_elems, y_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 		}
// 	}

// 	template <class SampleBeamformer, typename RFType, typename ResultType, typename LocationType>
// 	void rxBeamformingDTspaceCpu(
// 		bool interpolateRFlines,
// 		bool interpolateBetweenTransmits,
// 		size_t numTransducerElements,
// 		size_t numReceivedChannels,
// 		size_t numTimesteps,
// 		const RFType* RF,
// 		size_t numTxScanlines,
// 		size_t numRxScanlines,
// 		const ScanlineRxParameters3D* scanlines,
// 		size_t numZs,
// 		const LocationType* zs,
// 		const LocationType* x_elems,
// 		LocationType speedOfSound,
// 		LocationType dt,
// 		uint32_t additionalOffset,
// 		LocationType F,
// 		const WindowFunction* windowFunction,
// 		ResultType* s,
// 		const int32_t* elementToChannelMap)
// 	{
// 		if (interpolateRFlines)
// 		{
// 			if (interpolateBetweenTransmits)
// 			{
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, true, true, true>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, true, true, false>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 			else {
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, true, false, true>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, true, false, false>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 		}
// 		else {
// 			if (interpolateBetweenTransmits)
// 			{
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, false, true, true>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, false, true, false>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 			else {
// 				if (elementToChannelMap)
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, false, false, true>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 				else
// 				{
// 					rxBeamformingDTSPACECpu<SampleBeamformer, false, false, false>(
// 						numTransducerElements, numReceivedChannels, numTimesteps, RF,
// 						numTxScanlines, numRxScanlines, scanlines,
// 						numZs, zs, x_elems, speedOfSound, dt, additionalOffset, F, windowFunction, s, elementToChannelMap);
// 				}
// 			}
// 		}
// 	}

// 	template <typename ChannelDataType, typename ImageDataType>
// 	shared_ptr<USImage> RxBeamformerCuda::performRxBeamforming(
// 		RxBeamformerCuda::RxSampleBeamformer sampleBeamformer,
// 		shared_ptr<const USRawData> rawData,
// 		double fNumber,
// 		double speedOfSoundMMperS,
// 		WindowType windowType,
// 		WindowFunction::ElementType windowParameter,
// 		bool interpolateBetweenTransmits,
// 		int32_t additionalOffset) const
// 	{
// 		//Ensure the raw-data are on the CPU
// 		auto gRawData = rawData->getData<ChannelDataType>();
// 		if (gRawData->isGPU())
// 		{
// 			gRawData = std::make_shared<Container<ChannelDataType> >(LocationHost, *gRawData);
// 		}

// 		size_t numelOut = m_numRxScanlines*m_rxNumDepths;
// 		shared_ptr<Container<ImageDataType> > pData = std::make_shared<Container<ImageDataType> >(ContainerLocation::LocationHost, 0, numelOut);

// 		double dt = 1.0 / rawData->getSamplingFrequency();

// 		if (!m_windowFunction || m_windowFunction->getType() != windowType || m_windowFunction->getParameter() != windowParameter)
// 		{
// 			m_windowFunction = std::unique_ptr<WindowFunction>(new WindowFunction(windowType, windowParameter, m_windowFunctionNumEntries));
// 		}

// 		auto beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 		auto beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 		switch (sampleBeamformer)
// 		{
// 		case DelayAndSum:
// 			beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 			beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 			break;
// 		case DelayAndStdDev:
// 			beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerDelayAndStdDev, ChannelDataType, ImageDataType, LocationType>;
// 			beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerDelayAndStdDev, ChannelDataType, ImageDataType, LocationType>;
// 			break;
// 		case CoherenceFactorDelayAndSum:
// 			beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerCoherenceFactorDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 			beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerCoherenceFactorDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 			break;
// 		case TestSignal:
// 			beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerTestSignal, ChannelDataType, ImageDataType, LocationType>;
// 			beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerTestSignal, ChannelDataType, ImageDataType, LocationType>;
// 			break;
// 		case INVALID:
// 		default:
// 			beamformingFunction3D = &rxBeamformingDTspaceCpu3D<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 			beamformingFunction2D = &rxBeamformingDTspaceCpu<RxSampleBeamformerDelayAndSum, ChannelDataType, ImageDataType, LocationType>;
// 		}


// 		convertToDtSpace(dt, speedOfSoundMMperS, rawData->getNumElements());
// 		if (m_is3D)
// 		{
// 			beamformingFunction3D(
// 				true,
// 				interpolateBetweenTransmits,
// 				rawData->getNumElements(),
// 				rawData->getElementLayout(),
// 				rawData->getNumReceivedChannels(),
// 				rawData->getNumSamples(),
// 				gRawData->get(),
// 				rawData->getNumScanlines(), // numTxScanlines
// 				m_numRxScanlines,			// numRxScanlines
// 				m_pRxScanlines->get(),
// 				m_rxNumDepths, m_pRxDepths->get(),
// 				m_pRxElementXs->get(),
// 				m_pRxElementYs->get(),
// 				static_cast<LocationType>(m_speedOfSoundMMperS),
// 				static_cast<LocationType>(dt),
// 				additionalOffset,
// 				static_cast<LocationType>(fNumber),
// 				m_windowFunction.get(),
// 				pData->get(),
// 				(m_elementToChannelMap ? m_elementToChannelMap->get() : nullptr)
// 				);
// 		}
// 		else {
// 			beamformingFunction2D(
// 				true,
// 				interpolateBetweenTransmits,
// 				rawData->getNumElements(),
// 				rawData->getNumReceivedChannels(),
// 				rawData->getNumSamples(),
// 				gRawData->get(),
// 				rawData->getNumScanlines(), // numTxScanlines
// 				m_numRxScanlines,			// numRxScanlines
// 				m_pRxScanlines->get(),
// 				m_rxNumDepths, m_pRxDepths->get(),
// 				m_pRxElementXs->get(),
// 				static_cast<LocationType>(m_speedOfSoundMMperS),
// 				static_cast<LocationType>(dt),
// 				additionalOffset,
// 				static_cast<LocationType>(fNumber),
// 				m_windowFunction.get(),
// 				pData->get(),
// 				(m_elementToChannelMap ? m_elementToChannelMap->get() : nullptr)
// 				);
// 		}

// 		if (rawData->getImageProperties() != m_lastSeenImageProperties)
// 		{
// 			m_lastSeenImageProperties = rawData->getImageProperties();
// 			shared_ptr<USImageProperties> newProps = std::make_shared<USImageProperties>(*m_lastSeenImageProperties);
// 			newProps->setScanlineLayout(m_rxScanlineLayout);
// 			newProps->setNumSamples(m_rxNumDepths);
// 			newProps->setImageState(USImageProperties::RF);
// 			m_editedImageProperties = std::const_pointer_cast<const USImageProperties>(newProps);
// 		}

// 		auto retImage = std::make_shared<USImage>(
// 			vec2s{ m_numRxScanlines, m_rxNumDepths },
// 			pData,
// 			m_editedImageProperties,
// 			rawData->getReceiveTimestamp(),
// 			rawData->getSyncTimestamp());

// 		return retImage;
// 	}

// 	template
// 	shared_ptr<USImage> RxBeamformerCuda::performRxBeamforming<int16_t, int16_t>(
// 		RxBeamformerCuda::RxSampleBeamformer sampleBeamformer,
// 		shared_ptr<const USRawData> rawData,
// 		double fNumber,
// 		double speedOfSoundMMperS,
// 		WindowType windowType,
// 		WindowFunction::ElementType windowParameter,
// 		bool interpolateBetweenTransmits,
// 		int32_t additionalOffset) const;
// 	template
// 	shared_ptr<USImage> RxBeamformerCuda::performRxBeamforming<int16_t, float>(
// 		RxBeamformerCuda::RxSampleBeamformer sampleBeamformer,
// 		shared_ptr<const USRawData> rawData,
// 		double fNumber,
// 		double speedOfSoundMMperS,
// 		WindowType windowType,
// 		WindowFunction::ElementType windowParameter,
// 		bool interpolateBetweenTransmits,
// 		int32_t additionalOffset) const;
// 	template
// 	shared_ptr<USImage> RxBeamformerCuda::performRxBeamforming<float, int16_t>(
// 		RxBeamformerCuda::RxSampleBeamformer sampleBeamformer,
// 		shared_ptr<const USRawData> rawData,
// 		double fNumber,
// 		double speedOfSoundMMperS,
// 		WindowType windowType,
// 		WindowFunction::ElementType windowParameter,
// 		bool interpolateBetweenTransmits,
// 		int32_t additionalOffset) const;
// 	template
// 	shared_ptr<USImage> RxBeamformerCuda::performRxBeamforming<float, float>(
// 		RxBeamformerCuda::RxSampleBeamformer sampleBeamformer,
// 		shared_ptr<const USRawData> rawData,
// 		double fNumber,
// 		double speedOfSoundMMperS,
// 		WindowType windowType,
// 		WindowFunction::ElementType windowParameter,
// 		bool interpolateBetweenTransmits,
// 		int32_t additionalOffset) const;
// }