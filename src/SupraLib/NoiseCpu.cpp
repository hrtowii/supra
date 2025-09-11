// // ================================================================================================
// // 
// // CPU implementation replacing CUDA noise generation functionality
// // 
// // ================================================================================================

// #include "NoiseCuda.h"
// #include <random>
// #include <algorithm>
// #include <cmath>

// using namespace std;

// namespace supra
// {
// 	namespace NoiseCpuInternal
// 	{
// 		typedef NoiseCuda::WorkType WorkType;

// 		// CPU implementation of the noise processing
// 		template <typename InputType, typename OutputType>
// 		void processKernel(
// 			const InputType* inputImage,
// 			const WorkType* noiseAdditiveUniform,
// 			const WorkType* noiseAdditiveGauss,
// 			const WorkType* noiseMultiplicativeUniform,
// 			const WorkType* noiseMultiplicativeGauss,
// 			vec3s size,
// 			WorkType additiveUniformMin, WorkType additiveUniformRange,
// 			WorkType multiplicativeUniformMin, WorkType multiplicativeUniformRange, 
// 			OutputType* outputImage)
// 		{
// 			size_t width = size.x;
// 			size_t height = size.y;
// 			size_t depth = size.z;

// 			for (size_t z = 0; z < depth; z++)
// 			{
// 				for (size_t y = 0; y < height; y++)
// 				{
// 					for (size_t x = 0; x < width; x++)
// 					{
// 						size_t idx = x + y*width + z*width*height;
						
// 						// Get the input pixel value and cast it to working type
// 						WorkType inPixel = inputImage[idx];
// 						WorkType noiseAdditiveUniformVal = additiveUniformMin + additiveUniformRange * noiseAdditiveUniform[idx];
// 						WorkType noiseAdditiveGaussVal = noiseAdditiveGauss[idx];
// 						WorkType noiseMultiplicativeUniformVal = multiplicativeUniformMin + multiplicativeUniformRange * noiseMultiplicativeUniform[idx];
// 						WorkType noiseMultiplicativeGaussVal = noiseMultiplicativeGauss[idx];

// 						// Add the noise
// 						WorkType value = inPixel * noiseMultiplicativeUniformVal * noiseMultiplicativeGaussVal + noiseAdditiveUniformVal + noiseAdditiveGaussVal;

// 						// Store the output pixel value with clamping
// 						outputImage[idx] = clampCast<OutputType>(value);
// 					}
// 				}
// 			}
// 		}
// 	}

// 	// Simple CPU-based noise correlation using basic averaging filter
// 	shared_ptr<Container<NoiseCuda::WorkType> > NoiseCuda::makeNoiseCorrelated(const shared_ptr<const Container<WorkType>>& in,
// 		size_t width, size_t height, size_t depth)
// 	{
// 		auto out = make_shared<Container<WorkType> >(LocationHost, 0, width*height*depth);
		
// 		// Simple 5x5 averaging filter to simulate spatial correlation
// 		const int filterSize = 5;
// 		const int filterOffset = filterSize / 2;
		
// 		for (size_t z = 0; z < depth; z++)
// 		{
// 			for (size_t y = 0; y < height; y++)
// 			{
// 				for (size_t x = 0; x < width; x++)
// 				{
// 					WorkType sum = 0;
// 					int count = 0;
					
// 					for (int fy = -filterOffset; fy <= filterOffset; fy++)
// 					{
// 						for (int fx = -filterOffset; fx <= filterOffset; fx++)
// 						{
// 							int nx = static_cast<int>(x) + fx;
// 							int ny = static_cast<int>(y) + fy;
							
// 							if (nx >= 0 && nx < static_cast<int>(width) && 
// 								ny >= 0 && ny < static_cast<int>(height))
// 							{
// 								sum += in->get()[nx + ny*width + z*width*height];
// 								count++;
// 							}
// 						}
// 					}
					
// 					if (count > 0)
// 						out->get()[x + y*width + z*width*height] = sum / count;
// 					else
// 						out->get()[x + y*width + z*width*height] = in->get()[x + y*width + z*width*height];
// 				}
// 			}
// 		}

// 		return out;
// 	}

// 	template <typename InputType, typename OutputType>
// 	shared_ptr<Container<OutputType> > NoiseCuda::process(const shared_ptr<const Container<InputType>>& imageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated)
// 	{
// 		size_t width = size.x;
// 		size_t height = size.y;
// 		size_t depth = size.z;
// 		size_t numel = width * height * depth;

// 		// Ensure data is on CPU
// 		auto inImageData = imageData;
// 		if (inImageData->isGPU())
// 		{
// 			inImageData = make_shared<Container<InputType> >(LocationHost, *inImageData);
// 		}
		
// 		// Prepare output memory
// 		auto outImageData = make_shared<Container<OutputType> >(LocationHost, 0, numel);

// 		// Generate noise using standard C++ random facilities
// 		std::random_device rd;
// 		std::mt19937 gen(rd());
// 		std::uniform_real_distribution<WorkType> uniformDist(0.0f, 1.0f);
// 		std::normal_distribution<WorkType> normalDist(0.0f, 1.0f);
		
// 		// Prepare noise containers
// 		auto noiseAdditiveUniform = make_shared<Container<WorkType> >(LocationHost, 0, numel);
// 		auto noiseAdditiveGauss = make_shared<Container<WorkType> >(LocationHost, 0, numel);
// 		auto noiseMultiplicativeUniform = make_shared<Container<WorkType> >(LocationHost, 0, numel);
// 		auto noiseMultiplicativeGauss = make_shared<Container<WorkType> >(LocationHost, 0, numel);

// 		// Generate random noise
// 		for (size_t i = 0; i < numel; i++)
// 		{
// 			noiseAdditiveUniform->get()[i] = uniformDist(gen);
// 			noiseAdditiveGauss->get()[i] = additiveGaussMean + additiveGaussStd * normalDist(gen);
// 			noiseMultiplicativeUniform->get()[i] = uniformDist(gen);
// 			noiseMultiplicativeGauss->get()[i] = multiplicativeGaussMean + multiplicativeGaussStd * normalDist(gen);
// 		}

// 		// Apply spatial correlation if requested
// 		if (additiveUniformCorrelated)
// 		{
// 			noiseAdditiveUniform = makeNoiseCorrelated(noiseAdditiveUniform, width, height, depth);
// 		}
// 		if (additiveGaussCorrelated)
// 		{
// 			noiseAdditiveGauss = makeNoiseCorrelated(noiseAdditiveGauss, width, height, depth);
// 		}
// 		if (multiplicativeUniformCorrelated)
// 		{
// 			noiseMultiplicativeUniform = makeNoiseCorrelated(noiseMultiplicativeUniform, width, height, depth);
// 		}
// 		if (multiplicativeGaussCorrelated)
// 		{
// 			noiseMultiplicativeGauss = makeNoiseCorrelated(noiseMultiplicativeGauss, width, height, depth);
// 		}
		
// 		// Process the image with noise
// 		NoiseCpuInternal::processKernel<InputType, OutputType>(
// 			inImageData->get(),
// 			noiseAdditiveUniform->get(),
// 			noiseAdditiveGauss->get(),
// 			noiseMultiplicativeUniform->get(),
// 			noiseMultiplicativeGauss->get(),
// 			size,
// 			additiveUniformMin, additiveUniformMax - additiveUniformMin,
// 			multiplicativeUniformMin, multiplicativeUniformMax - multiplicativeUniformMin,
// 			outImageData->get());

// 		return outImageData;
// 	}

// 	// Template instantiations
// 	template
// 		shared_ptr<Container<uint8_t> > NoiseCuda::process<int16_t, uint8_t>(const shared_ptr<const Container<int16_t> >& inImageData, vec3s size, 
// 			WorkType additiveUniformMin, WorkType additiveUniformMax,
// 			WorkType additiveGaussMean, WorkType additiveGaussStd,
// 			WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 			WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 			bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 			bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<uint8_t> > NoiseCuda::process<float, uint8_t>(const shared_ptr<const Container<float> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<uint8_t> > NoiseCuda::process<uint8_t, uint8_t>(const shared_ptr<const Container<uint8_t> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<float> > NoiseCuda::process<int16_t, float>(const shared_ptr<const Container<int16_t> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<float> > NoiseCuda::process<float, float>(const shared_ptr<const Container<float> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<float> > NoiseCuda::process<uint8_t, float>(const shared_ptr<const Container<uint8_t> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<int16_t> > NoiseCuda::process<int16_t, int16_t>(const shared_ptr<const Container<int16_t> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<int16_t> > NoiseCuda::process<float, int16_t>(const shared_ptr<const Container<float> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// 	template
// 	shared_ptr<Container<int16_t> > NoiseCuda::process<uint8_t, int16_t>(const shared_ptr<const Container<uint8_t> >& inImageData, vec3s size,
// 		WorkType additiveUniformMin, WorkType additiveUniformMax,
// 		WorkType additiveGaussMean, WorkType additiveGaussStd,
// 		WorkType multiplicativeUniformMin, WorkType multiplicativeUniformMax,
// 		WorkType multiplicativeGaussMean, WorkType multiplicativeGaussStd,
// 		bool additiveUniformCorrelated, bool additiveGaussCorrelated,
// 		bool multiplicativeUniformCorrelated, bool multiplicativeGaussCorrelated);
// }